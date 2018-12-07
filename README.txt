% EPLL: an Image Denoising Method Using a Gaussian Mixture Model Learned on a
%       Large Set of Patches

* Author    : EHRET Thibaud <ehret.thibaud@gmail.com>
* Copyright : (C) 2018 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see gpl.txt

OVERVIEW
--------

This source code provides an implementation of EPLL developped in
"D. Zoran and Y. Weiss, From learning models of natural image patches to whole
image restoration, ICCV 2011".
It also extends the original work to color images, and provides a script to use
it in conjunction with the multiscaler of N. Pierazzo and G. Facciolo
(https://github.com/npd/multiscaler).

This code is part of an IPOL publication (http://www.ipol.im/). Plase cite it
if you use this code as part of your research.

COMPILATION
-----------

The code is compilable on Unix/Linux and hopefully on Mac OS (not tested!). 

Compilation: requires the cmake and make programs.

Dependencies: CBLAS, LAPACKE, OpenMP [can be disabled]. 
For image i/o we use Enric Meinhardt's iio (https://github.com/mnhrdt/iio),
which requires libpng, libtiff and libjpeg.
 
Configure and compile the source code using cmake and make.  It is recommended
that you create a folder for building:

UNIX/LINUX/MAC:
$ mkdir build; cd build
$ cmake ..
$ make

Binaries will be created in `build/bin folder`.

NOTE: By default, the code is compiled with OpenMP multithreaded
parallelization enabled (if your system supports it). 
The number of threads used by the code is defined in `epll/epll.h`.

USAGE
-----

The following commands have to be run from the `build/bin` folder:

List all available options:</br>
$ ./main_epll --help

There are 4 mandatory input arguments:
* -i     : the input image
* -cmodel: the GMM covariances
* -wmodel: the GMM weights
* -sigma : the noise level (additive white Gaussian noise, AGWN)

We provide three trained models in the `models/` folder:
* `gs_original:` GMM for grayscale patches proposed in the work of Zoran and Weiss
* `gs_homemade:` our GMM for grayscale patches
* `color_homemade:` our GMM for RGB patches

-----

The following command denoises an image using the GMMs trained in the original
article.

$ ./main_epll -i image.png -cmodel ../../models/sigma_gs_original.txt -wmodel ../../models/w_gs_original.txt -sigma 20

Noise of std. dev. 20 is added by the program before denoising (use option
`-add 0` if the input image already has noise). The result is
available in `denoised.tiff`. If the input image is color, the grayscale model
is applied channel by channel.  One can choose between two colorspaces: RGB
(default) and OPP by setting `-opp 1`.

-----

If the model trained for RGB patches is used, then we must specify that the model 
has three channels by passing the argument `-psc 3`:

$ ./main_epll -i image.png -cmodel ../../models/sigma_color_homemade.txt -wmodel ../../models/w_color_homemade.txt -psc 3 -sigma 20

NOTE: the `color_homemade` GMM was trained in the RGB colorspace. It should not be 
used with `-opp 1`.

-----

To use the multiscaler, we provide the BASH script `denoising_multiscale.sh`.
Examples are:

$ ./denoising_multiscale.sh input.png 20 output.png "-cmodel ../../models/sigma_gs_original.txt -wmodel ../../models/w_gs_original.txt -opp 1" 

$ ./denoising_multiscale.sh input.png 20 output.png "-cmodel ../../models/sigma_color_original.txt -wmodel ../../models/w_color_original.txt -psc 3"


FILES
-----

This project contains the following source files:

	main function:            src/main_epll.cpp
	multiscaler script:       src/denoising_multiscale.sh
	command line parsing:     src/cmd_option.h
	epll implementation:      src/EPLL/epll.h
	                          src/EPLL/epll.cpp
	basic matrix operations:  src/EPLL/LibMatrix.h
	                          src/EPLL/LibMatrix.cpp
	image i/o:                src/EPLL/iio.h
	                          src/EPLL/iio.c
	image container:          src/EPLL/LibImages.h
	                          src/EPLL/LibImages.cpp
	random number genernator: src/EPLL/mt19937ar.h
	                          src/EPLL/mt19937ar.c
	GMM models:               models/sigma_color_homemade.txt
	                          models/w_color_homemade.txt
	                          models/sigma_gs_homemade.txt
	                          models/w_gs_homemade.txt
	                          models/sigma_gs_original.txt
	                          models/w_gs_original.txt

It also contains the source code of the multiscaler
(https://github.com/npd/multiscaler/) in folder `src/multiscaler`.

ABOUT THIS FILE
---------------

Copyright 2018 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
