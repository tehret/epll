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

#ifndef EPLL_H_INCLUDED
#define EPLL_H_INCLUDED

#define NTHREAD 16
#define PI 3.1415926535

#include <vector>
#include "../Utilities/LibImages.h"


/**
 * @brief Structure containing a Gaussian model used for the pdf estimation.
 *
 * @param logweight: Log of the weight associated to this Gaussian in the GMM
 * @param eigVects: Eigenvectors of the covariance matrix of the Gaussian
 * @param eigVals: Eigenvalues of the covariance matrix of the Gaussian
 * @param invSqrtCov: Stores the precomputation of the inverse of the square
 *                    root of the covariance matrix for the computation of the log pdf
 * @param logdet: Precomputation of the log of the determinant
 * @param rank: Max rank kept for the covariance matrix;
 **/
struct Model{
	float logweight;
	std::vector<float> eigVects;
	std::vector<float> eigVals;
	std::vector<float> invSqrtCov;
	float logdet;
	int rank;
};

/**
 * @brief Apply the half-quadratic splitting denoising method
 *
 * @param noiseI : Image to denoise
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
 * @param models : Gaussian mixture model
 **/
void EPLLhalfQuadraticSplit(
		std::vector<float>& noiseI,
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
		std::vector<Model>& models);

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
		std::vector<Model>& models);

/**
 * @brief Compute the Gaussiian log pdf for different patches, the weight of the given Gaussian inside the GMM is taken into account
 *
 * @param patches : Store the patches for which the probability needs to be computed
 * @param dim : Dimension of a patch
 * @param nbP : Number of patches
 * @param model : Model of the Gaussian to be used 
 * @param csta : Precomputed constant set to dim + log(2*PI)
 * @param result : Store the probability for each patch
 **/
void loggausspdf(std::vector<float>& patches, int dim, int nbP, Model& model, float csta, std::vector<float>& result);

#endif
