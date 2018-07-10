#!/bin/bash
# All rights reserved.
#
# This program is free software: you can use, modify and/or
# redistribute it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later
# version. You should have received a copy of this license along
# this program. If not, see <http://www.gnu.org/licenses/>.

# dctdenoising_multi.sh 
#        noisy.tif
#        40
#        out.tif      # also creates out.tif.non.tif
#        "-w 8 -1"    # dctdenoising params (optional)
#        3            # LEVELS of pyramid (default: 4 = auto) (optional)
#        2            # R_PYR pyramid ratio: 2 (optional)
#        0.7          # PAR_PYR recomposition ratio : 0.7 (optional)

#set -x
#set -u
#set -e

if [ $# -lt 4 ]; then
   echo "$0 noisy.tif sigma out.tif"
   echo "    \"-w 8 -1\"  # dctdenoising params (optional)"
   echo "    4            # LEVELS of pyramid (default: 4) (optional)"
   echo "    2            # R_PYR pyramid ratio: 2 (optional), 1.5 also possible"
   echo "    0.7          # PAR_PYR recomposition ratio : 0.7 (optional)"
   exit 1
fi

INPUT=$1
NOISE=$2
OUTPUT=$3
EPLL_ARGS=$4
LEVELS=4
R_PYR=2   
PAR_PYR=0.7

if [ -n "$5" ]; then
   LEVELS=$5
fi
if [ -n "$6" ]; then
   R_PYR=$6
fi
if [ -n "$7" ]; then
   PAR_PYR=$7
fi

DEC_ARGS="-r ${R_PYR}"
MS_ARGS="-c ${PAR_PYR}"

# Create the noisy version of the image
main_epll -i ${INPUT} -sigma $NOISE -ps 0

# Decompose the noisy image into the LEVELS scales requested
decompose noisy.tiff level_ ${LEVELS} .tiff ${DEC_ARGS}

# Denoise each scale independently, while taking into account the different noise level for a given scale
for ((lvl=LEVELS-1; lvl>=0; --lvl))
do
  sigma=$(bc <<< "scale=2; $NOISE / ${R_PYR}^$lvl")
  main_epll -i level_${lvl}.tiff -sigma ${sigma} -deno level_${lvl}_non.tiff ${EPLL_ARGS} -add false -psnr false
done

# Compose the final result from the independent denoised results 
recompose level_ ${LEVELS} _non.tiff $OUTPUT ${MS_ARGS}

