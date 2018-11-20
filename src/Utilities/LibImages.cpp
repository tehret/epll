/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
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
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibImages.h"
#include "io_png.h"
#include "Utilities.h"
#include "mt19937ar.h"

#include <unistd.h> // getpid
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

extern "C" {
#include "iio.h"
}

using namespace std;

/**
 * @brief Load image, check the number of channels.
 *
 * @param p_name : name of the image to read;
 * @param o_im : vector which will contain the image : R, G and B concatenated;
 * @param o_imSize : will contain the size of the image;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_SUCCESS if the image has been loaded, EXIT_FAILURE otherwise.
 **/
int loadImage(
	const char* p_name
,	std::vector<float> &o_im
,	ImageSize &o_imSize
,	const bool p_verbose
){
	//! read input image
	if (p_verbose) cout << endl << "Read input image...";

	float *imTmp = NULL;
	//size_t w, h, c;
	//imTmp = read_png_f32(p_name, &w, &h, &c);
	int w, h, c;
	imTmp =  iio_read_image_float_vec(p_name, &w, &h, &c);

	if (!imTmp) {
		cout << "error :: " << p_name << " not found or not a correct png image" << endl;
		return EXIT_FAILURE;
	}

	if (p_verbose) cout << "done." << endl;

	//! test if image is really a color image and exclude the alpha channel
	if (c > 2) {
		unsigned k = 0;
		while (k < w * h && imTmp[k] == imTmp[w * h + k] && imTmp[k] == imTmp[2 * w * h + k])
			k++;

		c = (k == w * h ? 1 : 3);
	}

	//! Some image informations
	if (p_verbose) {
		cout << "image size :" << endl;
		cout << " - width          = " << w << endl;
		cout << " - height         = " << h << endl;
		cout << " - nb of channels = " << c << endl;
	}

	//! Initializations
	o_imSize.width      = w;
	o_imSize.height     = h;
	o_imSize.nChannels  = c;
	o_imSize.wh         = w * h;
	o_imSize.whc        = w * h * c;
	o_im.resize(w * h * c);
	for (unsigned k = 0; k < w * h * c; k++)
		o_im[k] = imTmp[k];

	free(imTmp);

	return EXIT_SUCCESS;
}

/**
 * @brief write image.
 *
 * @param p_name : path+name+extension of the image;
 * @param i_im : vector which contains the image;
 * @param p_imSize : size of the image;
 *
 * @return EXIT_SUCCESS if the image has been saved, EXIT_FAILURE otherwise
 **/
int saveImage(
    const char* p_name
,   std::vector<float> const& i_im
,   const ImageSize &p_imSize
){
    //! Allocate Memory
    float* imTmp = new float[p_imSize.whc];

    unsigned c = p_imSize.nChannels;
    unsigned w = p_imSize.width;
    unsigned h = p_imSize.height;

    //! Check for boundary problems
    for (unsigned k = 0; k < p_imSize.whc; k++)
        imTmp[k] = i_im[k];

    iio_save_image_float_vec(p_name, imTmp, w, h, c);
    //! Free Memory
    delete[] imTmp;

    return EXIT_SUCCESS;
}

/**
 * @brief add noise to img.
 *
 * @param i_im : original noise-free image;
 * @param o_imNoisy = im + noise;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return none.
 **/
void addNoise(
    std::vector<float> const& i_im
,   std::vector<float> &o_imNoisy
,   const float p_sigma
,   const bool p_verbose
){
    if (p_verbose) {
        cout << "Add noise [sigma = " << p_sigma << "] ...";
    }

	//! Initialization
    o_imNoisy = i_im;
    mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) getpid());

    //! Add noise
    for (unsigned k = 0; k < i_im.size(); k++) {
        const double a = mt_genrand_res53();
        const double b = mt_genrand_res53();

        o_imNoisy[k] += p_sigma * (float) (sqrtl(-2.0l * log(a)) * cos(2.0l * M_PI * b));
    }

    if (p_verbose) {
        cout << "done." << endl;
    }
}

/**
 * @brief Compute PSNR and RMSE between i_im1 and i_im2
 *
 * @param i_im1 : pointer to an allocated array of pixels;
 * @param i_im2 : pointer to an allocated array of pixels;
 * @param o_psnr  : will contain the PSNR;
 * @param o_rmse  : will contain the RMSE;
 * @param p_imageName: name of the image;
 * @param p_verbose: if true, print values of PSNR and RMSE.
 *
 * @return EXIT_FAILURE if both images haven't the same size.
 **/
int computePsnr(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   float &o_psnr
,   float &o_rmse
,   const float i_max = 1.
){
    if (i_im1.size() != i_im2.size()) {
        return EXIT_FAILURE;
    }
    float sum = 0.f;
    for (unsigned k = 0; k < i_im1.size(); k++)
        sum += (i_im1[k] - i_im2[k]) * (i_im1[k] - i_im2[k]);

    o_rmse = sqrtf(sum / (float) i_im1.size());
    o_psnr = 20.f * log10f(i_max / o_rmse);

    return EXIT_SUCCESS;
}

/**
 * @brief Compute a difference image between i_im1 and i_im2.
 *
 * @param i_im1: reference image;
 * @param i_im2: image to compare;
 * @param o_imDiff: will contain the difference;
 * @param p_sigma : standard deviation of the noise;
 * @param p_min, p_max : range of data (usually [0, 255]);
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if i_im1 and i_im2 don't have the same size.
 **/
int computeDiff(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   std::vector<float> &o_imDiff
,   const float p_sigma
,   const float p_max
){
    if (i_im1.size() != i_im2.size()) {
        return EXIT_FAILURE;
    }

    const unsigned size = i_im1.size();
    if (o_imDiff.size() != size) {
        o_imDiff.resize(size);
    }

    for (unsigned k = 0; k < size; k++) {
        o_imDiff[k] =  (i_im1[k] - i_im2[k] + p_sigma) * p_max / (2.f * p_sigma);
    }

    return EXIT_SUCCESS;
}

/**
 * @brief Transform the color space of an image, from RGB to YUV, or vice-versa.
 *
 * @param io_im: image on which the transform will be applied;
 * @param p_imSize: size of io_im;
 * @param p_isForward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
 *
 * @return none.
 **/
void transformColorSpace(
	std::vector<float> &io_im
,	const ImageSize p_imSize
,	const bool p_isForward
){
	//! If the image as only one channel, do nothing
	if (p_imSize.nChannels == 1) return;

	//! Initialization
	const unsigned width  = p_imSize.width;
	const unsigned height = p_imSize.height;
	const unsigned chnls  = p_imSize.nChannels;
	const unsigned wh     = width * height;
	vector<float> imTmp(wh * chnls);

	//! RGB to YUV
	if (p_isForward) {
		if (chnls == 3) {
			const unsigned red   = 0;
			const unsigned green = 1;
			const unsigned blue  = 2;
			const float a = 1.f / sqrtf(3.f);
			const float b = 1.f / sqrtf(2.f);
			const float c = 2.f * a * sqrtf(2.f);

			for (unsigned k = 0; k < wh; k++) {
				//! Y channel
				imTmp[k*chnls + red] = a * (io_im[k*chnls + red] + io_im[k*chnls + green] + io_im[k*chnls + blue]);

				//! U channel
				imTmp[k*chnls + green] = b * (io_im[k*chnls + red] - io_im[k*chnls + blue]);

				//! V channel
				imTmp[k*chnls + blue] = c * (0.25f * io_im[k*chnls + red ] - 0.5f * io_im[k*chnls + green]
				                      + 0.25f * io_im[k*chnls + blue]);
			}
		}
		else { //! chnls == 4
			const unsigned Gr = 0;
			const unsigned R  = 1;
			const unsigned B  = 2;
			const unsigned Gb = 3;
			const float a = 0.5f;
			const float b = 1.f / sqrtf(2.f);

			for (unsigned k = 0; k < wh; k++) {
				imTmp[k*chnls + Gr] = a * ( io_im[k*chnls + Gr] + io_im[k*chnls + R ] +
				                      io_im[k*chnls + B ] + io_im[k*chnls + Gb]);
				imTmp[k*chnls + R ] = b * ( io_im[k*chnls + R ] - io_im[k*chnls + B ]);
				imTmp[k*chnls + B ] = a * (-io_im[k*chnls + Gr] + io_im[k*chnls + R ] +
				                      io_im[k*chnls + B ] - io_im[k*chnls + Gb]);
				imTmp[k*chnls + Gb] = b * (-io_im[k*chnls + Gr] + io_im[k*chnls + Gb]);
			}
		}
	}
	//! YUV to RGB
	else {
		if (chnls == 3) {
			const unsigned red   = 0;
			const unsigned green = 1;
			const unsigned blue  = 2;
			const float a = 1.f / sqrtf(3.f);
			const float b = 1.f / sqrtf(2.f);
			const float c = a / b;

			for (unsigned k = 0; k < wh; k++) {
				//! R channel
				imTmp[k*chnls + red  ] = a * io_im[k*chnls + red] + b * io_im[k*chnls + green]
				                               + c * 0.5f * io_im[k*chnls + blue];
				//! G channel
				imTmp[k*chnls + green] = a * io_im[k*chnls + red] - c * io_im[k*chnls + blue];

				//! B channel
				imTmp[k*chnls + blue ] = a * io_im[k*chnls + red] - b * io_im[k*chnls + green]
				                               + c * 0.5f * io_im[k*chnls + blue];
			}
		}
		else {	//! chnls == 4
			const unsigned Gr = 0;
			const unsigned R  = 1;
			const unsigned B  = 2;
			const unsigned Gb = 3;
			const float a = 0.5f;
			const float b = 1.f / sqrtf(2.f);
			for (unsigned k = 0; k < wh; k++) {
				imTmp[k*chnls + Gr] = a * io_im[k*chnls + Gr] - a * io_im[k*chnls + B] - b * io_im[k*chnls + Gb];
				imTmp[k*chnls + R ] = a * io_im[k*chnls + Gr] + b * io_im[k*chnls + R] + a * io_im[k*chnls + B];
				imTmp[k*chnls + B ] = a * io_im[k*chnls + Gr] - b * io_im[k*chnls + R] + a * io_im[k*chnls + B];
				imTmp[k*chnls + Gb] = a * io_im[k*chnls + Gr] - a * io_im[k*chnls + B] + b * io_im[k*chnls + Gb];
			}
		}
	}

	io_im = imTmp;
}

/**
 * @brief Clip the values of im_im between val_min and val_max.
 *
 * @param io_im: image on which the cipping will be applied;
 * @param val_min: minimum value;
 * @param val_max: maximum value.
 *
 * @return none.
 **/
void clip(
    std::vector<float> &io_im
,   const float val_min
,   const float val_max
)
{
	for(int i = 0; i < io_im.size(); ++i)	
		io_im[i] = std::min(std::max(io_im[i], val_min), val_max);
}
