/*
 * Original work: Copyright (c) 2017, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file epll.cpp
 * @brief EPLL denoising functions
 *
 * @author Thibaud Ehret <ehret.thibaud@gmail.com>
 **/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <string>

#include "epll.h"
#include "LibMatrix.h"
#include "../Utilities/Utilities.h"
#include <ctime>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

void EPLLhalfQuadraticSplit(Video<float>& noiseI, Video<float>& finalI, Video<float>& origI, bool partialPSNR, float noiseSD, int patchsize, int patchsizeChannels, std::vector<float> betas, int T, int pas, std::vector<Model>& models)
{
	int N = patchsize*patchsize;
	int lambda = N;

	finalI = noiseI;

	for(int i = 0; i < betas.size(); ++i)
	{
		for(int it = 0; it < T; ++it)
		{
			Video<float> tempI(noiseI.sz, 0.);
			float sigma = noiseSD/std::sqrt(betas[i]);
			sigma /= 255.;
			aprxMAPGMM(finalI, tempI, sigma, patchsize, patchsizeChannels, pas, models);

			for(int x = 0; x < finalI.sz.width; ++x)	
			for(int y = 0; y < finalI.sz.height; ++y)	
			for(int t = 0; t < finalI.sz.frames; ++t)	
			for(int c = 0; c < finalI.sz.channels; ++c)	
			{
				finalI(x,y,t,c) = (lambda*noiseI(x,y,t,c) + betas[i]*N*tempI(x,y,t,c))/(lambda + betas[i]*N);
			}

			if(partialPSNR)
			{
				float psnr, rmse;
				VideoUtils::computePSNR(finalI, origI, psnr, rmse);
				float psnr2, rmse2;
				VideoUtils::computePSNR(tempI, origI, psnr2, rmse2);
				printf("Beta = \t%f, Iter =\t%d  PSNR =\t%f RMSE =\t%f (temporary image PSNR =\t%f RMSE =\t%f)\n", betas[i], it, psnr, 255.*rmse, psnr2, 255.*rmse2);
			}
		}
	}
}

void aprxMAPGMM(Video<float>& noiseI, Video<float>& tempI, float sigma, int ps, int psc, int step, std::vector<Model>& models)
{
	int N = ps*ps*psc;
	float sigma2 = sigma*sigma;
	Video<float> count(noiseI.sz);

	std::vector<float> tempPatch(N);
	std::vector<float> patch(N);
	std::vector<float> tempVects;

	Video<int> mask(noiseI.sz.width, noiseI.sz.height, noiseI.sz.frames, noiseI.sz.channels, 0);

	//! Compute the mask of patches that need denoising 
	int nbP = 0;
	for(int x = 0; x <= noiseI.sz.width-ps; ++x)	
	for(int y = 0; y <= noiseI.sz.height-ps; ++y)	
	for(int t = 0; t < noiseI.sz.frames; ++t)	
	for(int c = 0; c <= noiseI.sz.channels-psc; ++c)	
	{
		if((x % step == 0) || (x == noiseI.sz.width-ps))
		if((y % step == 0) || (y == noiseI.sz.height-ps))
		{
			mask(x,y,t,c) = 1;
			++nbP;
		}
	}

	std::vector<float> patches(N*nbP);
	std::vector<float> means(nbP);

	// Prep the different objects in models
	for(int m = 0; m < models.size(); ++m)
	{
		models[m].logdet = 0;
		for(int v = 0; v < models[m].eigVals.size(); ++v)
			models[m].logdet += std::log(models[m].eigVals[v] + sigma2);

		float *invSqrtCov = models[m].invSqrtCov.data();
		float *eigv = models[m].eigVects.data();
		for (unsigned k = 0, plouf = 0; k < models[m].r; ++k)
		for (unsigned i = 0; i < N; ++i, ++plouf)
			*invSqrtCov++ = (*eigv++) / std::sqrt(sigma2 + models[m].eigVals[k]);
	}

	// Compute the denoised patches
	int k = 0;
	for(int x = 0; x < noiseI.sz.width; ++x)	
	for(int y = 0; y < noiseI.sz.height; ++y)	
	for(int t = 0; t < noiseI.sz.frames; ++t)	
	for(int c = 0; c < noiseI.sz.channels; ++c)	
	{
		if(mask(x,y,t,c) == 1)
		{
			// Compute the DC component
			means[k] = 0.f;
			for(int dc = 0; dc < psc; ++dc)
			for(int dx = 0; dx < ps; ++dx)
			for(int dy = 0; dy < ps; ++dy)
				means[k] += noiseI(x+dx,y+dy,t,c+dc);

			means[k] /= (float)N;

			//load patch into vector
			for(int dc = 0, d = 0; dc < psc; ++dc)	
			for(int dx = 0;        dx < ps; ++dx)	
			for(int dy = 0;        dy < ps; ++dy, ++d)	
				patches[k + nbP*d] = (noiseI(x+dx,y+dy,t,c+dc) - means[k]);
			++k;
		}
	}

	//Compute the weights for all patches at the same time
	float csta = N * std::log(2*PI);
	std::vector<std::vector<float> > weights(models.size(), std::vector<float>(nbP));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD) \
	shared(weights, models)
#endif
	for(int m = 0; m < models.size(); ++m)
		loggausspdf(patches, N, nbP, models[m], csta, weights[m]);

	// Compute the denoised patches
	k = 0;
	for(int x = 0; x < noiseI.sz.width; ++x)	
	for(int y = 0; y < noiseI.sz.height; ++y)	
	for(int t = 0; t < noiseI.sz.frames; ++t)	
	for(int c = 0; c < noiseI.sz.channels; ++c)	
	{
		if(mask(x,y,t,c) == 1)
		{
			// Select the best Gaussian
			int best = 0;
			float bestv = weights[0][k];
			for(int m = 1; m < models.size(); ++m)
			{
				if(weights[m][k] > bestv)
				{
					best = m;
					bestv = weights[m][k];
				}
			}

			for(int d = 0; d < N; ++d)
				patch[d] = patches[k + nbP*d];

			//! Z' = X'*U
			productMatrix(tempPatch,
					patch,
					models[best].eigVects,
					1, models[best].r, N,
					false, false);

			//! U * W
			std::vector<float> eigVecs(models[best].eigVects);
			float *eigv = eigVecs.data();
			for (unsigned k = 0; k < models[best].r; ++k)
			for (unsigned i = 0; i < N; ++i)
				*eigv++ *= (models[best].eigVals[k] / (models[best].eigVals[k] + sigma2));

			//! hX' = Z'*(U*W)'
			productMatrix(patch,
					tempPatch,
					eigVecs,
					1, N, models[best].r,
					false, true);

			for(int d = 0; d < N; ++d)
				patch[d] += means[k];

			for(int dc = 0, d = 0; dc < psc; ++dc)	
			for(int dx = 0;        dx < ps; ++dx)	
			for(int dy = 0;        dy < ps; ++dy, ++d)	
			{
				tempI(x+dx,y+dy,t,c+dc) += patch[d];
				count(x+dx,y+dy,t,c+dc)++;
			}
			++k;
		}
	}

//	std::vector<std::vector<int> > assignement(models.size(), std::vector<int>());
//	// Assign each patch to its best model
//	k = 0;
//	for(int x = 0; x < noiseI.sz.width; ++x)	
//	for(int y = 0; y < noiseI.sz.height; ++y)	
//	for(int t = 0; t < noiseI.sz.frames; ++t)	
//	for(int c = 0; c < noiseI.sz.channels; ++c)	
//	{
//		if(mask(x,y,t,c) == 1)
//		{
//			// Select the best Gaussian
//			int best = 0;
//			float bestv = weights[0][k];
//			for(int m = 1; m < models.size(); ++m)
//			{
//				if(weights[m][k] > bestv)
//				{
//					best = m;
//					bestv = weights[m][k];
//				}
//			}
//			assignement[best].push_back(k);
//			++k;
//		}
//	}
//
//	std::vector<std::vector<float> > processPatches(models.size());
//	for(int m = 0; m < models.size(); ++m)
//	processPatches[m].resize(N*assignement[m].size());
//
//#ifdef _OPENMP
//#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD) \
//	shared(weights, models)
//#endif
//	for(int m = 0; m < models.size(); ++m)
//	{
//		if(assignement[m].size() == 0)
//			continue;
//		
//		// load the patches into a group
//		// denoise them
//		for(int p = 0; p < assignement[m].size(); ++p)
//		for(int d = 0; d < N; ++d)
//			processPatches[m][p + assignement[m].size()*d] = patches[assignement[m][p] + nbP*d];
//
//		std::vector<float> tempPatches(assignement.size() * models[m].r);
//
//		//! Z' = X'*U
//		//productMatrix(tempPatches,
//		//		processPatches[m],
//		//		models[m].eigVects,
//		//		assignement[m].size(), models[m].r, N,
//		//		false, false);
//
//		////! U * W
//		//std::vector<float> eigVecs(models[m].eigVects);
//		//float *eigv = eigVecs.data();
//		//for (unsigned k = 0; k < models[m].r; ++k)
//		//	for (unsigned i = 0; i < N; ++i)
//		//		*eigv++ *= (models[m].eigVals[k] / (models[m].eigVals[k] + sigma2));
//
//		////! hX' = Z'*(U*W)'
//		//productMatrix(processPatches[m],
//		//		tempPatches,
//		//		eigVecs,
//		//		assignement[m].size(), N, models[m].r,
//		//		false, true);
//	}
//
//	for(int m = 0; m < models.size(); ++m)
//	{
//		for(int p = 0; p < assignement[m].size(); ++p)
//		{
//			for(int d = 0; d < N; ++d)	
//				patches[assignement[m][p] + nbP*d] = processPatches[m][p + assignement[m].size()*d];
//		}
//	}
//
//	k = 0;
//	for(int x = 0; x < noiseI.sz.width; ++x)	
//	for(int y = 0; y < noiseI.sz.height; ++y)	
//	for(int t = 0; t < noiseI.sz.frames; ++t)	
//	for(int c = 0; c < noiseI.sz.channels; ++c)	
//	{
//		if(mask(x,y,t,c) == 1)
//		{
//			// Aggregate back all elements from patches
//			for(int dx = 0, d = 0; dx < ps; ++dx)	
//			for(int dy = 0;        dy < ps; ++dy)	
//			for(int dc = 0;        dc < psc; ++dc, ++d)	
//			{
//				tempI(x+dx,y+dy,t,c+dc) += patches[k + nbP*d] + means[k];
//				count(x+dx,y+dy,t,c+dc)++;
//			}
//			++k;
//		}
//	}

	for(int x = 0; x < noiseI.sz.width; ++x)	
	for(int y = 0; y < noiseI.sz.height; ++y)	
	for(int t = 0; t < noiseI.sz.frames; ++t)	
	for(int c = 0; c < noiseI.sz.channels; ++c)	
		tempI(x,y,t,c) /= count(x,y,t,c);

	count.saveVideo("count.tiff", 1, 1);
}

void loggausspdf(std::vector<float>& patches, int dim, int nbP, Model& model, float csta, std::vector<float>& results)
{
	for(int p = 0; p < nbP; ++p)
		results[p] = 0.;

	std::vector<float> output(patches.size());

	productMatrix(output,
			patches,
			model.invSqrtCov,
			nbP, model.r, dim,
			false, false);

	for(int p = 0; p < nbP; ++p)
	{
		for(int d = 0; d < dim; ++d)
			results[p] -= output[p + nbP*d]*output[p + nbP*d];

		results[p] -= (csta + model.logdet); 
		results[p] /= 2.;
		results[p] += model.logweight;
	}
}
