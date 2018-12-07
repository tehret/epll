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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>

#include <string>
#include <sstream>

#include "EPLL/LibImages.h"
#include "EPLL/epll.h"
#include "EPLL/LibMatrix.h"
#include "cmd_option.h"

using namespace std;

/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 * @author THIBAUD EHRET  <ehret.thibaud@gmail.com>
 **/

int main(int argc, char **argv)
{
	clo_usage("EPLL denoising");
	clo_help(" NOTE: Input (<) and output (>).\n");

	//! Paths to input/output images
	using std::string;
	const string input_path = clo_option("-i"    , ""              , "< input image");
	const string noisy_path = clo_option("-nsy"  , "noisy.tiff"    , "> noisy image");
	const string diff_path  = clo_option("-diff" , "diff.tiff"     , "> difference image");
	const string final_path = clo_option("-deno" , "denoised.tiff" , "> denoised image");
	const string cov_model_path = clo_option("-cmodel" , "covs_opt.txt", "< model containing the covariance of the Gaussians");
	const string w_model_path   = clo_option("-wmodel" , "w_opts.txt"  , "< model containing the weigths of the Gaussians");

	//! General parameters
	const float sigma       = clo_option("-sigma", 0, "< standard deviation of the noise");
	const bool  add_noise   = clo_option("-add", true, "< add noise of given standard deviation sigma");
	const int   patch_size  = clo_option("-ps", 8, "< patch size. It must be the same than the one of the model");
	const int   patch_size_channels = clo_option("-psc", 1, "< number of channel of the model
			                                       (1 for grayscale or 3 for color). It must be the same than the one of the model");
	const int   step        = std::min(patch_size, clo_option("-st", 1, "< step size"));
	const int   iter        = clo_option("-T", 1, "< nb iter");
	const bool  partialPSNR = clo_option("-psnr", true, "< print partial PSNR");

	const bool  changeBasis = clo_option("-yuv", false, "< use YUV instead of RGB colorspace (if color image)");
	const int   rank        = clo_option("-r", 100, "< maximum rank used for the covariance of the Gaussians (as a percentage)");


	int firstFrame = 1, lastFrame = 1, frameStep = 1;
	bool verbose = false;

    std::vector<float> betas;
    if(sigma < 30)
        betas = {1,4,8,16,32,64};
    else
        betas = {1,2,8,13,32,64};

	//! Check inputs
	if (input_path == "")
	{
		fprintf(stderr, "%s: no input images.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	//! Declarations
	std::vector<float> original, noisy, final, diff;
	ImageSize imSize;

	//! Load input videos
	loadImage(input_path.c_str(), original, imSize, verbose);


	//! Add noise
	if(add_noise && sigma)
	{
		addNoise(original, noisy, sigma, verbose);

		//! Save noisy image
		saveImage(noisy_path.c_str(), noisy, imSize);
	}
	else
		noisy = original;

	if(patch_size <= 0)
		return EXIT_FAILURE;

	if(changeBasis)
	{
		transformColorSpace(noisy, imSize, true);
		transformColorSpace(original, imSize, true);
	}

	//! Loads the covs and the weights
	FILE* wfile = fopen(w_model_path.c_str(), "r");
	FILE* covfile = fopen(cov_model_path.c_str(), "r");

	int pdim = patch_size*patch_size*patch_size_channels;
	std::vector<Model> models;
	Model tmp_model;
	float w;
	int kpt = 0;

	//! Diagonalize the covariance matrices for the following computations. If speed is necessary, 
	//! these can be pre-saved diagonalized. In practice it this computation is quite fast anyway 
	//! and is in no way a bottleneck.
	while(fscanf(wfile, "%f,", &w) != EOF)
	{
		tmp_model.logweight = std::log(w);
		tmp_model.eigVects.resize(pdim*pdim);
		tmp_model.invSqrtCov.resize(pdim*pdim);
		tmp_model.eigVals.resize(pdim);
		tmp_model.rank = std::min(std::max(pdim*rank/100,1),pdim);

		std::vector<float> covMat(pdim*pdim);
		for(int d = 0; d < pdim*pdim; ++d)
			fscanf(covfile, "%f,", &covMat[d]);

		int info = matrixEigs(covMat, pdim, tmp_model.rank, tmp_model.eigVals, tmp_model.eigVects);
		models.push_back(tmp_model);
	}
	fclose(wfile);
	fclose(covfile);


	printf("Loaded %d models\n", models.size());

	//! Denoising
	if (verbose) printf("Running EPLL on the noisy image\n");

	//! Computation is done on images ranging in [0,1] contrary to the usual [0,255]
	for(int i = 0; i < noisy.size(); ++i)	
	{
		noisy[i] /= 255;
		original[i] /= 255;
	}

	//! Run denoising algorithm
	EPLLhalfQuadraticSplit(noisy, final, original, imSize, partialPSNR, sigma, patch_size,
	                       patch_size_channels, betas, iter, step, models);

	//! Compute PSNR and RMSE
	float final_psnr = -1, final_rmse = -1;
	computePsnr(original, final, final_psnr, final_rmse, 1.);

	printf("final PSNR =\t%f\tRMSE =\t%f\n", final_psnr, 255.*final_rmse);

	//! Compute Difference
	computeDiff(original, final, diff, sigma);

	//! Go back to the usual range [0,255] before saving the images
	for(int i = 0; i < noisy.size(); ++i)	
	{
		final[i] *= 255.;
		diff[i] *= 255.;
	}

    //! Go back to RGB colorspace if necessary
	if(changeBasis)
	{
		transformColorSpace(final, imSize, false);
		transformColorSpace(diff, imSize, false);
	}

	//! Save output images
	saveImage(final_path.c_str(), final,  imSize);
	saveImage(diff_path.c_str(), diff,   imSize);

	printf("Done\n");
	return EXIT_SUCCESS;
}
