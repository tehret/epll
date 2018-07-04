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

#include <string>
#include <sstream>

#include "Utilities/Utilities.h"
#include "Utilities/LibImages.h"
#include "EPLL/epll.h"
#include "EPLL/LibMatrix.h"
#include "Utilities/cmd_option.h"

using namespace std;

/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 *
 * @author THIBAUD EHRET  <ehret.thibaud@gmail.com>
 **/

int main(int argc, char **argv)
{
	clo_usage("EPLL denoising");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input/output sequences
	using std::string;
	const string  input_path = clo_option("-i"    , ""         , "< input image");
	const string  noisy_path = clo_option("-nsy"    , "noisy.tiff"         , "> noisy image");
	const string  diff_path = clo_option("-diff"    , "diff.tiff"         , "> difference image");
	const string  final_path = clo_option("-deno" , "denoised.tiff" , "> denoised image");
	const string  cov_model_path = clo_option("-cmodel" , "covs_opt.txt", "< model containing the covariance of the Gaussians");
	const string  w_model_path = clo_option("-wmodel" , "w_opts.txt", "< model containing the weigths of the Gaussians");

	//! General parameters
	const float sigma       = clo_option("-sigma", 0, "< standard deviation of the noise");
	const bool  add_noise = clo_option("-add", true, "< add noise of given standard deviation sigma");
	const int   patch_size  = clo_option("-ps",    8, "< patch size. It must be the same than the one of the model");
	const int   patch_size_channels  = clo_option("-psc",   1, "< number of channel of the model (1 for a learning in grayscale and 3 for one in color). It must be the same than the one of the model");
	const int   step        = std::max(patch_size, clo_option("-st",    1, "< step size"));
	const int   iter        = clo_option("-T",     1, "< nb iter");
	const bool  partialPSNR = clo_option("-psnr", true, "< print partial PSNR");

	const bool  changeBasis = clo_option("-yuv", false, "< change the RGB basis to YUV if color image");
	const int   rank        = clo_option("-r", 1000000000, "< maximum rank used for the covariance of the Gaussians");


	int firstFrame = 1, lastFrame = 1, frameStep = 1;
	bool verbose = false;

	// Defaults betas suggested in the original article
	std::vector<float> betas{1,4,8,16,32}; 

	//! Check inputs
	if (input_path == "")
	{
		fprintf(stderr, "%s: no input sequence.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	//! Declarations
	std::vector<float> original, noisy, final, diff;
	ImageSize imSize;

	//! Load input videos
	loadImage(input_path, original, imSize, verbose);


	//! Add noise
	if(add_noise && sigma)
	{
		addNoise(original, noisy, sigma, verbose);

		//! Save noisy image
		saveImage(noisy_path, noisy, imSize);
	}
	else
		noisy = original;

	if(patch_size <= 0)
	       return 0;	

	if(changeBasis)
	{
		transformColorSpace(noisy, imSize, true);
		transformColorSpace(original, imSize, true);
	}

	// Loads the covs and the weights
	
	FILE* wfile = fopen(w_model_path.c_str(), "r");
	FILE* covfile = fopen(cov_model_path.c_str(), "r");

	int N = patch_size*patch_size*patch_size_channels;
	std::vector<Model> models;
	Model current;
	float w;
	int kpt = 0;

	// Diagonalize the covariance matrices for the following computations. If speed is necessary, 
	// these can be pre-saved diagonalized. In practice it this computation is quite fast anyway 
	// and is in no way a bottleneck.
	while(fscanf(wfile, "%f,", &w) != EOF)
	{
		current.logweight = std::log(w);
		current.eigVects.resize(N*N);
		current.invSqrtCov.resize(N*N);
		current.eigVals.resize(N);
		current.r = std::min(N, rank);

		std::vector<float> covMat(N*N);
		for(int d = 0; d < N*N; ++d)
			fscanf(covfile, "%f,", &covMat[d]);

		int info = matrixEigs(covMat, N, current.r, current.eigVals, current.eigVects);
		models.push_back(current);
	}
	fclose(wfile);
	fclose(covfile);


	printf("Loaded %d models\n", models.size());

	//! Denoising
	if (verbose) printf("Running EPLL on the noisy image\n");

	// Computation is done on images ranging in [0,1] contrary to the usual [0,255]
	for(int x = 0; x < noisy.sz.width; ++x)	
	for(int y = 0; y < noisy.sz.height; ++y)	
	for(int t = 0; t < noisy.sz.frames; ++t)	
	for(int c = 0; c < noisy.sz.channels; ++c)	
	{
		noisy(x,y,t,c) /= 255;
		original(x,y,t,c) /= 255;
	}
	

	//! Run denoising algorithm
	EPLLhalfQuadraticSplit(noisy, final, original, partialPSNR, sigma, patch_size, patch_size_channels, betas, iter, step, models);

	//! Compute PSNR and RMSE
	float final_psnr = -1, final_rmse = -1;
	computePsnr(original, final, final_psnr, final_rmse, 1.);

	printf("final PSNR =\t%f\tRMSE =\t%f\n", final_psnr, 255.*final_rmse);

	//! Compute Difference
	computeDiff(original, final, diff, sigma);

	// Go back to the usual range [0,255] before saving the images
	for(int x = 0; x < noisy.sz.width; ++x)	
	for(int y = 0; y < noisy.sz.height; ++y)	
	for(int t = 0; t < noisy.sz.frames; ++t)	
	for(int c = 0; c < noisy.sz.channels; ++c)	
	{
		final(x,y,t,c) *= 255.;
		diff(x,y,t,c) *= 255.;
	}


	if(changeBasis)
	{
		transformColorSpace(final, imSize, false);
		transformColorSpace(diff, imSize, false);
	}

	//! Save output sequences
	saveVideo(final, final_path, firstFrame, frameStep);
	saveVideo( diff_path, firstFrame, frameStep);

	printf("Done\n");
	return EXIT_SUCCESS;
}
