/*
 * Original work: Copyright (c) 2016, Thibaud Ehret <ehret.thibaud@gmail.com>
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
#include "Utilities/cmd_option.h"
#include "Utilities/LibImages.h"

using namespace std;

/**
 * @file   compute_results.cpp
 * @brief  Compute the PSNR and difference between two images
 *
 * @author THIBAUD EHRET <ehret.thibaud@gmail.com>
 **/

int main(int argc, char **argv)
{
	clo_usage("Compute PSNR and difference");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input images
	using std::string;
	const string  input_path   = clo_option("-i" , "" , "< input sequence");
	const string  inbsc_path   = clo_option("-r" , "" , "< input sequence 2");
	const float   sigma        = clo_option("-s" , 20., "< noise std");
	const string  outpsnr_path = clo_option("-f" , "psnr.txt" , "> where to write the PSNR informations");
	const string  out_path     = clo_option("-o1", "original.png" , "> where to write the original image (used for tiff to png transformation)");
	const string  outdiff_path = clo_option("-o2", "diff.png" , "> where to write the diff image");
	const string  comment      = clo_option("-m" , "" , "> comment before printing the result, used for formatting");


	//! Declarations
	std::vector<float> original, final, diff;
	ImageSize imSize;

	//! Load input image
	loadImage(input_path.c_str(), original, imSize);
	loadImage(inbsc_path.c_str(),    final, imSize);

	saveImage(out_path.c_str(),      final, imSize);

	float final_psnr = -1, final_rmse = -1;
	computePsnr(original, final, final_psnr, final_rmse, 255.);
	computeDiff(original, final, diff, sigma);
	saveImage(outdiff_path.c_str(), diff, imSize);

	ofstream file;
	file.open(outpsnr_path, ios::app);
	file << comment;
	file << "\tPSNR =\t" << final_psnr << ",\tRMSE =\t"<< final_rmse << endl;
	file.close();

	return EXIT_SUCCESS;
}
