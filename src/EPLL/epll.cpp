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

	finalI = noiseI;

	for(int i = 0; i < betas.size(); ++i)
	{
		for(int it = 0; it < T; ++it)
		{
			Video<float> tempI(noiseI.sz, 0.);
			float sigma = noiseSD/std::sqrt(betas[i]);
			sigma /= 255.;

			// Compute the maximum a posteriori estimation of the set of patches in finalI for the given pre-learned GMM model
			aprxMAPGMM(finalI, tempI, sigma, patchsize, patchsizeChannels, pas, models);

			// Update the final estimate by combining the noisy image and the current estimate
			for(int x = 0; x < finalI.sz.width; ++x)	
			for(int y = 0; y < finalI.sz.height; ++y)	
			for(int t = 0; t < finalI.sz.frames; ++t)	
			for(int c = 0; c < finalI.sz.channels; ++c)	
				finalI(x,y,t,c) = (noiseI(x,y,t,c) + betas[i]*tempI(x,y,t,c))/(1. + betas[i]);

			// Print partial PSNR to track the improvement
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

	// Compute the mask of patches that need denoising 
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

	// Prepare the different objects in models:
	// 1/ Compute the log of the determinant
	// 2/ Precompute the matrix to fast compute the inverse of the covariance matrix
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
			// Compute the DC component (average patch)
			means[k] = 0.f;
			for(int dc = 0; dc < psc; ++dc)
			for(int dx = 0; dx < ps; ++dx)
			for(int dy = 0; dy < ps; ++dy)
				means[k] += noiseI(x+dx,y+dy,t,c+dc);

			means[k] /= (float)N;

			// Load patch into the specific vector while removing the DC component
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
		// If this patch requires denoising
		if(mask(x,y,t,c) == 1)
		{
			// Select the best Gaussian based on the probability estimated for each Gaussian of the GMM (weights of the GMM are taken into account)
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

			// Extract the patch from the list of patches
			for(int d = 0; d < N; ++d)
				patch[d] = patches[k + nbP*d];

			// Compute the MAP patch for the chosen Gaussian 
			productMatrix(tempPatch,
					patch,
					models[best].eigVects,
					1, models[best].r, N,
					false, false);
			std::vector<float> eigVecs(models[best].eigVects);
			float *eigv = eigVecs.data();
			for (unsigned k = 0; k < models[best].r; ++k)
			for (unsigned i = 0; i < N; ++i)
				*eigv++ *= (models[best].eigVals[k] / (models[best].eigVals[k] + sigma2));
			productMatrix(patch,
					tempPatch,
					eigVecs,
					1, N, models[best].r,
					false, true);

			// Add back the DC component
			for(int d = 0; d < N; ++d)
				patch[d] += means[k];

			// Aggregate the result on the result image
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

	// Finish the aggregation by averaging all contribution
	for(int x = 0; x < noiseI.sz.width; ++x)	
	for(int y = 0; y < noiseI.sz.height; ++y)	
	for(int t = 0; t < noiseI.sz.frames; ++t)	
	for(int c = 0; c < noiseI.sz.channels; ++c)	
		tempI(x,y,t,c) /= count(x,y,t,c);
}

void loggausspdf(std::vector<float>& patches, int dim, int nbP, Model& model, float csta, std::vector<float>& results)
{
	// Initialize the results
	for(int p = 0; p < nbP; ++p)
		results[p] = 0.;

	std::vector<float> output(patches.size());

	// For all patches provided: compute the log of the probability for the given Gaussian distribution
	// First compute x.U.sqrt(D)
	productMatrix(output,
			patches,
			model.invSqrtCov,
			nbP, model.r, dim,
			false, false);

	for(int p = 0; p < nbP; ++p)
	{
		// Compute xSx^T=(x.U.sqrt(D)) . (x.U.sqrt(D))^T
		for(int d = 0; d < dim; ++d)
			results[p] -= output[p + nbP*d]*output[p + nbP*d];

		results[p] -= (csta + model.logdet); 
		results[p] /= 2.;
		results[p] += model.logweight;
	}
}
