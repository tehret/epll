% EPLL: an Image Denoising Method Using a Gaussian Mixture Model Learned on a Large Set of Patches

# ABOUT

* Author    : EHRET Thibaud <ehret.thibaud@gmail.com>
* Copyright : (C) 2018 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see gpl.txt

Version 1.0, released on November 28, 2018

# OVERVIEW

This source code provides an implementation of EPLL developped in "D. Zoran and Y. Weiss, From learning models of natural image patches to whole image restoration, ICCV 2011", which has been studied in details in http://www.ipol.im

# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux and Mac OS. 

- Compilation. 
Automated compilation requires the cmake and make programs.

- Library. 
This code requires the libpng, libtiff, libjpeg, CBLAS and lapack library.

- Image format. 
Only the TIFF, JPEG and PNG formats are supported. 
 
-------------------------------------------------------------------------
Usage:
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

4. This project contains the following source files:
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


5. The files that have been reviewed for IPOL publication are
TODO ONCE REVIEW IS DONE


# ABOUT THIS FILE

Copyright 2018 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
