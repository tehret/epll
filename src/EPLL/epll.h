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

#define NTHREAD 32
#define PI 3.1415926535

#include <vector>
#include "../Utilities/LibImages.h"


/**
 * @brief Structure containing a Gaussian model used for the pdf estimation.
 *
 * @param logweight: Log of the weight associated to this Gaussian in the GMM
 * @param eigVects: Eigenvectors of the covariance matrix of the Gaussian
 * @param eigVals: Eigenvalues of the covariance matrix of the Gaussian
 * @param invSqrtCov: Stores the precomputation of the inverse of the square root of the covariance matrix for the computation of the log pdf
 * @param logdet: Precomputation of the log of the determinant
 * @param r: Max rank kept for the covariance matrix;
 **/
struct Model{
	float logweight;
	std::vector<float> eigVects;
	std::vector<float> eigVals;
	std::vector<float> invSqrtCov;
	float logdet;
	int r;
};

/**
 * @brief Apply the half-quadratic splitting denoising method
 *
 * @param noiseI : Image that need denoising
 * @param finalI : Stores the output image
 * @param origI : Non-noisy image to compute PSNRs
 * @param partialPSNR : Print partial PSNR to show the evolution when going through iterations and betas 
 * @param noiseSD : Standard deviation of the noise
 * @param patchsize : Spatial size of the patch
 * @param patchsizeChannels : Number of channel to use for the denoising (defined by the model)
 * @param betas : Stores the serie of beta coefficient used for the splitting
 * @param T : Number of iteration 
 * @param pas : Step of the spatial grid used to extract the patches in the image
 * @param models : Stores the Gaussian mixture model (one Gaussian per entry of the array)
 **/
void EPLLhalfQuadraticSplit(std::vector<float>& noiseI, std::vector<float>& finalI, std::vector<float>& origI, bool partialPSNR, float noiseSD, int patchsize, int patchsizeChannels, std::vector<float> betas, int T, int pas, std::vector<Model>& models);

/**
 * @brief Compute the approximative maximum a priori image reconstructed using the Gaussian mixture model represented by models
 *
 * @param noiseI : Image that need denoising
 * @param tempI : Store the results of the reconstruction
 * @param sigma : Standard deviation of the noise
 * @param ps : Spatial size of the patch
 * @param psc : Number of channel to use for the denoising (defined by the model)
 * @param step : Step of the spatial grid used to extract the patches in the image
 * @param models : Stores the Gaussian mixture model (one Gaussian per entry of the array)
 **/
void aprxMAPGMM(std::vector<float>& noiseI, std::vector<float>& tempI, float sigma, int ps, int psc, int step, std::vector<Model>& models);

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
