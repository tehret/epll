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

/**
 * @brief Apply the half-quadratic splitting denoising method
 *
 * @param noisyI : Image to denoise
 * @param finalI : Output image
 * @param origI  : Non-noisy image (used to compute partial PSNRs)
 * @param imSize : Image size
 * @param partialPSNR : Print partial PSNR for each iteration and beta
 * @param noiseSD : Standard deviation of the noise
 * @param patchSize : Size of the patch
 * @param patchChannels : Number of channels to use for the denoising
 *                        (defined by the model)
 * @param betas  : Sequence of beta coefficients used for the H-Q splitting
 * @param T      : Number of iterations
 * @param step   : Step of the spatial grid of processed patches
 * @param models : patch Gaussian mixture model
 **/
void EPLLhalfQuadraticSplit(
		std::vector<float>& noisyI,
		std::vector<float>& finalI,
		std::vector<float>& origI,
		ImageSize& imSize,
		bool partialPSNR,
		float noiseSD,
		int patchSize,
		int patchChannels,
		std::vector<float> betas,
		int T,
		int step,
		std::vector<Model>& models)
{
	finalI = noisyI;

	for(int i = 0; i < betas.size(); ++i)
	{
		for(int it = 0; it < T; ++it)
		{
			std::vector<float> aggMAPI(finalI.size(), 0.);
			float sigma = noiseSD/std::sqrt(betas[i]);
			sigma /= 255.;

			// aggregated MAP patch estimates according to GMM a-priori model
			// (denoted Z in the IPOL paper)
			aprxMAPGMM(finalI, aggMAPI, imSize, sigma, patchSize, patchChannels,
			           step, models);

			// Update image estimate by combining the noisy image and the patches MAP
			for(int k = 0; k < finalI.size(); ++k)	
				finalI[k] = (noisyI[k] + betas[i]*aggMAPI[k])/(1. + betas[i]);

			// Print partial PSNR
			if(partialPSNR)
			{
				float psnr, rmse, psnr2, rmse2;
				computePsnr(finalI , origI, psnr , rmse , 1.);
				computePsnr(aggMAPI, origI, psnr2, rmse2, 1.);
				printf("beta = % 5.1f  iter = % 2d - image estimate: PSNR = %5.2f RMSE = %6.3f\n", betas[i], it, psnr, 255.*rmse);
				printf("                 - aggregated patch MAPs: PSNR = %5.2f RMSE = %6.3f\n", psnr2, 255.*rmse2);
			}
		}
	}
}

/**
 * @brief Computes an image resulting from the aggregation of the approximative
 * MAP estimates given the noisy images. The a priori distribution is given by
 * the pre-trained patch GMM.
 *
 * @param noisyI : Noisy input image
 * @param aggMAPI : Output image with the aggregated patch MAP estimates
 * @param imSize : Image size
 * @param sigma : Standard deviation of the noise
 * @param ps : Patch size
 * @param pc : Number of channels to use for the denoising (defined by the model)
 * @param step : Step of the spatial grid of processed patches
 * @param models : patch Gaussian mixture model
 **/
void aprxMAPGMM(
		std::vector<float>& noisyI,
		std::vector<float>& aggMAPI,
		ImageSize& imSize,
		float sigma,
		int ps,
		int pc,
		int step,
		std::vector<Model>& models)
{
	int N = ps*ps*pc;
	float sigma2 = sigma*sigma;
	std::vector<float> count(noisyI.size());

	std::vector<float> tempPatch(N);
	std::vector<float> patch(N);
	std::vector<float> tempVects;

	std::vector<int> mask(imSize.width*imSize.height*imSize.nChannels, 0);


	// Compute the mask of patches that need denoising 
	int nbP = 0;
	for(int c = 0; c <= imSize.nChannels-pc; ++c)	
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
		for (unsigned k = 0; k < models[m].rank; ++k)
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
			for(int dc = 0; dc < pc; ++dc)
				means[k] += noisyI[(x+dx)*imSize.nChannels + (y+dy)*imSize.width*imSize.nChannels + c+dc];

			means[k] /= (float)N;

			// Load patch into the specific vector while removing the DC component
			for(int dc = 0, d = 0; dc < pc; ++dc)	
			for(int dx = 0;        dx < ps; ++dx)	
			for(int dy = 0;        dy < ps; ++dy, ++d)	
				patches[k + nbP*d] = (noisyI[(x+dx)*imSize.nChannels + (y+dy)*imSize.width*imSize.nChannels + c+dc] - means[k]);
			++k;
		}
	}

	//Compute the weights for all patches at the same time
	std::vector<std::vector<float> > weights(models.size(), std::vector<float>(nbP));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD) \
	shared(weights, models)
#endif
	for(int m = 0; m < models.size(); ++m)
		logGausspdf(patches, N, nbP, models[m], weights[m]);

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
					1, models[best].rank, N,
					false, false);
			std::vector<float> eigVecs(models[best].eigVects);
			float *eigv = eigVecs.data();
			for (unsigned k = 0; k < models[best].rank; ++k)
			for (unsigned i = 0; i < N; ++i)
				*eigv++ *= (models[best].eigVals[k] / (models[best].eigVals[k] + sigma2));
			productMatrix(patch,
					tempPatch,
					eigVecs,
					1, N, models[best].rank,
					false, true);

			// Add back the DC component
			for(int d = 0; d < N; ++d)
				patch[d] += means[k];

			// Aggregate the result on the result image
			for(int dc = 0, d = 0; dc < pc; ++dc)	
			for(int dx = 0;        dx < ps; ++dx)	
			for(int dy = 0;        dy < ps; ++dy, ++d)	
			{
				aggMAPI[(x+dx)*imSize.nChannels + (y+dy)*imSize.width*imSize.nChannels + c+dc] += patch[d];
				count[(x+dx)*imSize.nChannels + (y+dy)*imSize.width*imSize.nChannels + c+dc]++;
			}
			++k;
		}
	}

	// Finish the aggregation by averaging all contribution
	for(int k = 0; k < aggMAPI.size(); ++k)	
		aggMAPI[k] /= count[k];
}

/**
 * @brief Compute the Gaussiian log pdf for different patches, the weight of
 * the given Gaussian inside the GMM is taken into account
 *
 * @param patches : Patches for which the probability needs to be computed
 * @param dim : Dimension of a patch
 * @param nbP : Number of patches
 * @param model : Model of the Gaussian to be used
 * @param result : Log Gaussian probability for each patch
 **/
void logGausspdf(
		std::vector<float>& patches,
		int dim,
		int nbP,
		Model& model,
		std::vector<float>& results)
{

	std::vector<float> output(model.rank * nbP);

	// Compute  R^{-1/2}Q^T x_i for all patches x_i in patches
	/* NOTE: patches is a nbP x dim matrix (in column-major layout). Therefore each
	 * patch x_i is a row in the matrix patches. model.invSqrtCov contains Q R^{-1/2}.
	 * In order to avoid unnecessary traspositions, we compute patches * Q R^{-1/2}.
	 */
	productMatrix(output,
			patches,
			model.invSqrtCov,
			nbP, model.rank, dim,
			false, false);

	// Initialize the results
	float constant = dim * std::log(2*PI);
	for(int p = 0; p < nbP; ++p)
	{
		// squared L2 norm R^{-1/2}Q^T x_i
		results[p] = 0.;
		for(int d = 0; d < dim; ++d)
			results[p] -= output[p + nbP*d]*output[p + nbP*d];

		// add logdet, constants, etc
		results[p] -= (constant + model.logdet);
		results[p] /= 2.;
		results[p] += model.logweight;
	}
}
