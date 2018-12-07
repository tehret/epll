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
#ifndef LIB_IMAGES_H_INCLUDED
#define LIB_IMAGES_H_INCLUDED

#include <vector>
#include <string>

/**
 * @brief Structure containing size informations of an image.
 *
 * @param width     : width of the image;
 * @param height    : height of the image;
 * @param nChannels : number of channels in the image;
 * @param wh        : equal to width * height. Provided for convenience;
 * @param whc       : equal to width * height * nChannels. Provided for convenience.
 **/
struct ImageSize
{
	unsigned width;
	unsigned height;
	unsigned nChannels;
	unsigned wh;
	unsigned whc;
};

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
,   std::vector<float> &o_im
,   ImageSize &o_imSize
,   const bool p_verbose = false
);

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
);

/**
 * @brief add noise to img.
 *
 * @param i_im : original noise-free image;
 * @param o_imNoisy = im + noise;
 * @param p_sigma : standard deviation of the noise.
 * @param p_verbose : if true, print some informations.
 *
 * @return none.
 **/
void addNoise(
    std::vector<float> const& i_im
,   std::vector<float> &o_imNoisy
,   const float p_sigma
,   const bool p_verbose
);

/**
 * @brief Compute PSNR and RMSE between i_im1 and i_im2
 *
 * @param i_im1 : pointer to an allocated array of pixels;
 * @param i_im2 : pointer to an allocated array of pixels;
 * @param o_psnr: will contain the PSNR;
 * @param o_rmse: will contain the RMSE;
 * @param i_max : range of the image, usually either 1 or 255;
 *
 * @return EXIT_FAILURE if both images haven't the same size.
 **/
int computePsnr(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   float &o_psnr
,   float &o_rmse
,   const float i_max
);

/**
 * @brief Compute a difference image between i_im1 and i_im2.
 *
 * @param i_im1: reference image;
 * @param i_im2: image to compare;
 * @param o_imDiff: will contain the difference;
 * @param p_sigma: standard deviation of the noise;
 * @param p_max : range of data (usually 1 or 255));
 *
 * @return EXIT_FAILURE if i_im1 and i_im2 don't have the same size.
 **/
int computeDiff(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   std::vector<float> &o_imDiff
,   const float p_sigma
,   const float p_max = 255.
);


/**
 * @brief Transform the color space of an image, from RGB to OPP, or vice-versa.
 *
 * @param io_im: image on which the transform will be applied;
 * @param p_imSize: size of io_im;
 * @param p_isForward: if true, go from RGB to OPP, otherwise go from OPP to RGB.
 *
 * @return none.
 **/
void transformColorSpace(
	std::vector<float> &io_im
,	const ImageSize p_imSize
,	const bool p_isForward
);

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
);

#endif // LIB_IMAGES_H_INCLUDED
