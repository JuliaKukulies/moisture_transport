#! /bin/bash

## This bash script preprocesses data from a 9 km WRF downscaling over the Tibetan Plateau 

## remember to activate conda basic env before CDO can be used on server!! 


for y in {1983..2007}
do
    for m in {05..09}
    do
	for d in {01..31}
	do
		for h in {00..23}
		do
		DIR="/media/DataP2/PlevTP/CF${y}/${m}/"
		# check if directory exists
		echo "${DIR}"
		if [ -d "${DIR}" ]; then
			cdo -select,name=r_v_p,r_solid_p,r_rain_p,r_cloud_p,u_tr_p,v_tr_p,Z_p,T_p ${DIR}post_wrfout_d01_${y}-${m}-${d}_${h}.nc data/hourly/wrfout_TP9km_CF_${y}-${m}-${d}-${h}.nc
		fi
		done
	done
    done
done





