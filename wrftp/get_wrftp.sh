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
	    cdo -ensmean -select,name=u_tr_p,v_tr_p,r_v_p,Z_p,T_p,precip_g ${DIR}wrfout_d01_TP9km_CF_${y}-${m}-*nc data/wrfout_TP9km_${y}_${m}.nc
	fi
    done
done





