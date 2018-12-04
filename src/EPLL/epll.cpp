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

void EPLLhalfQuadraticSplit(std::vector<float>& noiseI, std::vector<float>& finalI, std::vector<float>& origI, ImageSize& imSize, bool partialPSNR, float noiseSD, int patchsize, int patchsizeChannels, std::vector<float> betas, int T, int pas, std::vector<Model>& models)
{
	finalI = noiseI;

	for(int i = 0; i < betas.size(); ++i)
	{
		for(int it = 0; it < T; ++it)
		{
			std::vector<float> tempI(finalI.size(), 0.);
			float sigma = noiseSD/std::sqrt(betas[i]);
			sigma /= 255.;

			// Compute the maximum a posteriori estimation of the set of patches in finalI for the given pre-learned GMM model
			aprxMAPGMM(finalI, tempI, imSize, sigma, patchsize, patchsizeChannels, pas, models);

			// Update the final estimate by combining the noisy image and the current estimate
			for(int k = 0; k < finalI.size(); ++k)	
				finalI[k] = (noiseI[k] + betas[i]*tempI[k])/(1. + betas[i]);

			// Print partial PSNR to track the improvement
			if(partialPSNR)
			{
				float psnr, rmse;
				computePsnr(finalI, origI, psnr, rmse, 1.);
				float psnr2, rmse2;
				computePsnr(tempI, origI, psnr2, rmse2, 1.);
				printf("Beta = \t%f, Iter =\t%d  PSNR =\t%f RMSE =\t%f (temporary image PSNR =\t%f RMSE =\t%f)\n", betas[i], it, psnr, 255.*rmse, psnr2, 255.*rmse2);
			}
		}
	}
}

void aprxMAPGMM(std::vector<float>& noiseI, std::vector<float>& tempI, ImageSize& imSize, float sigma, int ps, int psc, int step, std::vector<Model>& models)
{
	int N = ps*ps*psc;
	float sigma2 = sigma*sigma;
	std::vector<float> count(noiseI.size());

	std::vector<float> tempPatch(N);
	std::vector<float> patch(N);
	std::vector<float> tempVects;

	std::vector<int> mask(imSize.width*imSize.height*imSize.nChannels, 0);


	// Compute the mask of patches that need denoising 
	int nbP = 0;
	for(int c = 0; c <= imSize.nChannels-psc; ++c)	
	for(int y = 0; y <= imSize.height-ps; ++y)	
	for(int x = 0; x <= imSize.width-ps; ++x)	
	{
		if((x % step == 0) || (x == imSize.width-ps))
		if((y % step == 0) || (y == imSize.height-ps))
		{
			mask[x*imSize.nChannels + y*imSize.width*imSize.nChannels + c] = 1;
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
		for (unsigned k = 0; k < models[m].r; ++k)
		for (unsigned i = 0; i < N; ++i)
			*invSqrtCov++ = (*eigv++) / std::sqrt(sigma2 + models[m].eigVals[k]);
	}

	// Extract patches
	int k = 0;
	for(int y = 0; y < imSize.height; ++y)	
	for(int x = 0; x < imSize.width; ++x)	
	for(int c = 0; c < imSize.nChannels; ++c)	
	{
		if(mask[x*imSize.nChannels + y*imSize.width*imSize.nChannels + c] == 1)
		{
			// Compute the DC component (average patch)
			means[k] = 0.f;
			for(int dy = 0; dy < ps; ++dy)
			for(int dx = 0; dx < ps; ++dx)
			for(int dc = 0; dc < psc; ++dc)
				means[k] += noiseI[(x+dx)*imSize.nChannels + (y+dy)*imSize.width*imSize.nChannels + c+dc];

			means[k] /= (float)N;

			// Load patch into the specific vector while removing the DC component
			for(int dc = 0, d = 0; dc < psc; ++dc)	
			for(int dx = 0;        dx < ps; ++dx)	
			for(int dy = 0;        dy < ps; ++dy, ++d)	
				patches[k + nbP*d] = (noiseI[(x+dx)*imSize.nChannels + (y+dy)*imSize.width*imSize.nChannels + c+dc] - means[k]);
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
	for(int y = 0; y < imSize.height; ++y)	
	for(int x = 0; x < imSize.width; ++x)	
	for(int c = 0; c < imSize.nChannels; ++c)	
	{
		// If this patch requires denoising
		if(mask[x*imSize.nChannels + y*imSize.width*imSize.nChannels + c] == 1)
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
				tempI[(x+dx)*imSize.nChannels + (y+dy)*imSize.width*imSize.nChannels + c+dc] += patch[d];
				count[(x+dx)*imSize.nChannels + (y+dy)*imSize.width*imSize.nChannels + c+dc]++;
			}
			++k;
		}
	}

	// Finish the aggregation by averaging all contribution
	for(int k = 0; k < tempI.size(); ++k)	
		tempI[k] /= count[k];
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
