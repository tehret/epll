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

#include <cmath>

#include "epll.h"
#include "LibMatrix.h"

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
	int pdim = ps*ps*pc;
	float sigma2 = sigma*sigma;

	int W = imSize.width, H = imSize.height, C = imSize.nChannels;

	// Compute the mask of patches that need denoising
	std::vector<int> mask(imSize.width*imSize.height*imSize.nChannels, 0);
	int nbP = 0; // number of patches
	for(int y = 0; y <= H - ps; ++y)
	for(int x = 0; x <= W - ps; ++x)
	{
		if((x % step == 0) || (x == W-ps))
		if((y % step == 0) || (y == H-ps))
		{
			for(int c = 0; c <= C - pc; ++c)
			{
				mask[y*W*C + x*C + c] = 1;
				++nbP;
			}
		}
	}

	// Prepare the different objects in models:
	for(int m = 0; m < models.size(); ++m)
	{
		// Compute the log of the determinant
		models[m].logdet = 0;
		for(int v = 0; v < models[m].eigVals.size(); ++v)
			models[m].logdet += std::log(models[m].eigVals[v] + sigma2);

		// Precompute the matrix to fast compute the inverse of the covariance matrix
		/* NOTE: in the notation of the IPOL article, this corresponds to
		 * R_m^{-1/2} Q_m^T. The matrix Q_m contains the eigenvectors as columns,
		 * and R_m is a diagonal matrix consisting of the addition of the eigenvalues
		 * plus the noise variance. The product R_m^{-1/2} Q_m^T amounts to scaling
		 * each eigenvector by corresponding coefficient in the diagonal of R_m^{-1/2}.
		 * The matrix Q_m is stored in models[m].eigVects in column-major layout,
		 * meaning that the first pdim elements are the first eigenvector, the next
		 * pdim components the second, etc.
		 */
		float *invSqrtCov = models[m].invSqrtCov.data();
		float *eigv = models[m].eigVects.data();
		for (int k = 0; k < models[m].rank; ++k) // for each eigenvector ...
		for (int i = 0; i < pdim; ++i)           // ... scale all its entries
			*invSqrtCov++ = (*eigv++) / std::sqrt(sigma2 + models[m].eigVals[k]);
	}

	// Extract patches
	std::vector<float> patchesDC(nbP);     // DC components of the patches
	std::vector<float> patches(pdim*nbP);
	int k = 0;
	for(int y = 0; y < H; ++y)
	for(int x = 0; x < W; ++x)
	for(int c = 0; c < C; ++c)
	if(mask[y*W*C + x*C + c] == 1)
	{
		// Compute the DC component (average value of patch)
		patchesDC[k] = 0.f;
		for(int dy = 0; dy < ps; ++dy)
		for(int dx = 0; dx < ps; ++dx)
		for(int dc = 0; dc < pc; ++dc)
			patchesDC[k] += noisyI[(y+dy)*W*C + (x+dx)*C + c+dc];

		patchesDC[k] /= (float)pdim;

		// Load patch into the specific vector while removing the DC component
		/* NOTE: patches is a nbP x pdim matrix stored in column-major layout */
		for(int dc = 0, d = 0; dc < pc; ++dc)
		for(int dx = 0;        dx < ps; ++dx)
		for(int dy = 0;        dy < ps; ++dy, ++d)
			patches[k + nbP*d] = (noisyI[(y+dy)*W*C + (x+dx)*C + c+dc] - patchesDC[k]);

		++k;
	}

	// Compute the conditional mixing weights for all patches
	std::vector<std::vector<float> > mixingWeights(models.size(), std::vector<float>(nbP));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(NTHREAD) \
	shared(mixingWeights, models)
#endif
	for(int m = 0; m < models.size(); ++m)
		logGausspdf(patches, pdim, nbP, models[m], mixingWeights[m]);

	// Compute MAP estimate for each patch
	std::vector<float> aggCount(noisyI.size(), 0.f);
	k = 0;
	for(int y = 0; y < H; ++y)
	for(int x = 0; x < W; ++x)
	for(int c = 0; c < C; ++c)
	if(mask[y*C*W + x*C + c] == 1) // only process some patches
	{
		// Select the best Gaussian component (highest conditional mixing weight)
		int best = 0;
		float bestv = mixingWeights[0][k];
		for(int m = 1; m < models.size(); ++m)
		{
			if(mixingWeights[m][k] > bestv)
			{
				best = m;
				bestv = mixingWeights[m][k];
			}
		}

		// Extract the patch from the list of patches
		std::vector<float> patch(pdim);
		for(int d = 0; d < pdim; ++d)
			patch[d] = patches[k + nbP*d];

		// Compute the MAP of the patch according to best Gaussian component

		// First project patch over the eigenvectors (only r of them):
		// tempPatch <-- Q_{kmax}^T * patch, where Q_{kmax} is the eigenvector matrix,
		std::vector<float> tempPatch(models[best].rank);
		productMatrix(tempPatch,                  // A^T*B
		              models[best].eigVects,      // A=eigVects, pdim x rank
		              patch,                      // B=patch,    pdim x 1
		              models[best].rank, 1, pdim, // dimensions
		              true, false);               // transpose A, don't transpose B

		// Apply shrinkage coefficients to patch components:
		// tempPatch <-- S_{kmax} * tempPatch, with S_{kmax} the diagonal shrinkage matrix
		for (unsigned k = 0; k < models[best].rank; ++k)
			tempPatch[k] *= (models[best].eigVals[k] / (models[best].eigVals[k] + sigma2));

		// Reconstruct patch from shrunk components: patch = Q_{kmax} * tempPatch
		productMatrix(patch,                       // A*B
		              models[best].eigVects,       // A=eigVects,  pdim x rank
		              tempPatch,                   // B=tempPatch, rank x 1
		              pdim, 1, models[best].rank,  // dimensions
		              false, false);               // don't transpose A nor B

		// Add back the DC component
		for(int d = 0; d < pdim; ++d)
			patch[d] += patchesDC[k];

		// Aggregate the result on the result image
		for(int dc = 0, d = 0; dc < pc; ++dc)
		for(int dx = 0;        dx < ps; ++dx)
		for(int dy = 0;        dy < ps; ++dy, ++d)
		{
			aggMAPI [(y+dy)*W*C + (x+dx)*C + c+dc] += patch[d];
			aggCount[(y+dy)*W*C + (x+dx)*C + c+dc] ++;
		}
		++k;
	}

	// Normalize the aggregated image
	for(int k = 0; k < aggMAPI.size(); ++k)
		aggMAPI[k] /= aggCount[k];
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
