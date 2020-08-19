#! /bin/bash

## This bash script preprocesses data from a 9 km WRF downscaling over the Tibetan Plateau

for y in {1979..2002}
do
    for m in {05..09}
    do
	DIR="/media/DataP3/PlevTP/CF${y}/${m}/"
	# check if directory exists
	echo "${DIR}"
	if [ -d "${DIR}" ]; then
	    # calculate monthly mean from hourly files 
	    cdo -select,name=u_tr_P,v_tr_p,r_v_p,Z_p,T_p ${DIR}wrfout_d01_TP9km_CF_${y}-${m}-*nc selected.nc
	    cdo timmean selected.nc data/wrfout_TP9km_${y}_${m}.nc
	    rm selected.nc
	fi
    done
done





