#!/bin/bash
#set -x
# dctdenoising_multi.sh 
#        noisy.tif
#        40
#        out.tif      # also creates out.tif.non.tif
#        "-w 8 -1"    # dctdenoising params (optional)
#        3            # LEVELS of pyramid (default: -1 = auto) (optional)
#        2            # R_PYR pyramid ratio: 2 (optional)
#        0.7          # PAR_PYR recomposition ratio : 0.7 (optional)

export PATH=$PATH:/bin:/usr/local/bin:/usr/bin:./

if [ $# -lt 4 ]; then
   echo "$0 noisy.tif sigma out.tif"
   echo "    \"-w 8 -1\"  # dctdenoising params (optional)"
   echo "    -1           # LEVELS of pyramid (default: -1 = auto) (optional)"
   echo "    2            # R_PYR pyramid ratio: 2 (optional), 1.5 also possible"
   echo "    0.7          # PAR_PYR recomposition ratio : 0.7 (optional)"
   exit 1
fi

INPUT=$1
NOISE=$2
OUTPUT=$3
EPLL_ARGS=$4
LEVELS=-1
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

# determine the levels based on the image size
if [ $LEVELS -eq -1 ]; then
   PIXELS=$(./num_pixels $INPUT)
   echo $PIXELS
   if [ ${PIXELS} -lt 500000 ]; then
     LEVELS=1
   elif [ ${PIXELS} -lt 2000000 ]; then
     LEVELS=2
   elif [ ${PIXELS} -lt 8000000 ]; then
     LEVELS=3
   else
     LEVELS=4
   fi
   echo "Scales: $LEVELS (auto)"
else
   echo "Scales: $LEVELS (requested)"
fi


DEC_ARGS="-r ${R_PYR}"
MS_ARGS="-c ${PAR_PYR}"


./main_epll -i ${INPUT} -sigma $NOISE ${EPLL_ARGS} -ps 0

./decompose noisy.tiff level_ ${LEVELS} .tiff ${DEC_ARGS}

lvl=$((LEVELS-1))
sigma=$(bc <<< "scale=2; $NOISE / ${R_PYR}^$lvl")
./main_epll -i level_${lvl}.tiff -sigma ${sigma} -deno level_${lvl}_non.png ${EPLL_ARGS} -add false -psnr false


for ((lvl=LEVELS-2; lvl>=0; --lvl))
do
  sigma=$(bc <<< "scale=2; $NOISE / ${R_PYR}^$lvl")
  ./main_epll -i level_${lvl}.tiff -sigma ${sigma} -deno level_${lvl}_non.png ${EPLL_ARGS} -add false -psnr false
done

wait
./recompose level_ ${LEVELS} _non.png $OUTPUT ${MS_ARGS}
./comp_psnr -i ${INPUT} -r ${OUTPUT}
