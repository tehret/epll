IMPLEMENTATION OF THE EPLL IMAGE DENOISING ALGORITHM
====================================================

* Author    : EHRET Thibaud <ehret.thibaud@gmail.com>
* Copyright : (C) 2018 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see gpl.txt

OVERVIEW
--------

This source code provides an implementation of EPLL developped in "D. Zoran and Y. Weiss, From learning models of natural image patches to whole image restoration, ICCV 2011".
It also extends the original work to color images, and provides a script to use
it in conjunction with the multiscaler of 

[BLA BLA BLA](link)

This code is part of an [IPOL](http://www.ipol.im/) publication. Plase cite it
if you use this code as part of your research.

COMPILATION
-----------

The code is compilable on Unix/Linux and hopefully on Mac OS (not tested!). 

Compilation: requires the cmake and make programs.

Dependencies: FFTW3, CBLAS, LAPACKE, OpenMP [can be disabled], iio (which requires libpng, libtiff and libjpeg).
 
USAGE
-----

1. Download the code package and extract it. Go to that directory. 

2. Configure and compile the source code using cmake and make. 
It is recommended that you create a folder for building:

UNIX/LINUX/MAC:
$ mkdir build; cd build
$ cmake ..
$ make

Binaries will be created in build/bin folder.

NOTE: By default, the code is compiled with OpenMP multithreaded
parallelization enabled (if your system supports it). 
The number of threads used by the code is defined in epll/epll.h.

3. Usage instruction:
Running './main_epll -i image.png -cmodel ../../models/sigma_gs_original.txt -wmodel ../../models/w_gs_original.txt -sigma 20' computes the result using the method with the original Gaussian mixture models after adding an additive white Gaussian noise of standard deviation 20. The result is available in denoised.tiff.
'./main_epll --help' list all available options.
The multiscale version is also available using "denoising_multiscale.sh". An example of usage is './denoising_multiscale.sh input.png 20 output.png "-cmodel ../../models/sigma_gs_original.txt -wmodel ../../models/w_gs_original.txt"' to denoise a degraded version of input.png with an additive white Gaussian noise o stndard deviation 20, the result is available in output.png. 

4. Three models are provided with the code in the models folder. We provide the original grayscale model "gs_homemade" and one we trained "gs_homemade". These models can directly be used with color images (in this case the denoising will be channel by channel, setting the parameter -yuv to true means that the denoising is done in the YUV colorspace instead of RGB). We also provide a color model "color_homemade" that requires setting the parameter -psc to 3 to denoise color images.

5. This project contains the following source files:
	src/main_epll.cpp
	src/denoising_multiscale.sh
	src/epll/epll.h
	src/epll/epll.cpp
    src/epll/LibMatrix.h
    src/epll/LibMatrix.cpp
	src/Utilities/iio.h
	src/Utilities/iio.c
	src/Utilities/LibImages.h
	src/Utilities/LibImages.cpp
	src/Utilities/cmd_option.h
	src/Utilities/comparators.h
	src/Utilities/mt19937ar.h
	src/Utilities/mt19937ar.c
	src/Utilities/Utilities.h
	src/Utilities/Utilities.cpp

It also contains the source code of the multiscaler tool from https://github.com/npd/multiscaler/ in multiscaler.


6. The files that have been reviewed for IPOL publication are
TODO ONCE REVIEW IS DONE


# ABOUT THIS FILE

Copyright 2018 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
