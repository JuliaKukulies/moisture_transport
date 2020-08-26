#! /bin/bash

## This bash script preprocesses data from a 9 km WRF downscaling over the Tibetan Plateau
 

for y in {1983..2007}
do
    for m in {05..09}
    do
	DIR="/media/DataP2/PlevTP/CF${y}/${m}/"
	# check if directory exists
	echo "${DIR}"
	if [ -d "${DIR}" ]; then
	    # calculate monthly mean from hourly files 
	    cdo -select,name=WaterFlx, precip_g ${DIR}post_wrfout_d01_${y}-${m}-*nc data/selected.nc
	    cdo timmean data/selected.nc data/wrfout_TP9km_CF_2D_${y}_${m}.nc
	    rm data/selected.nc
	fi
    done
done





