## This python script downloads hourly data from ERA5 at pressure levels, calculates deviations from its monthly mean and saves the monthly mean of the sum ofmean qV and its primes.

import numpy as np
import xarray as xr 
import glob
import cdsapi
import os
import atmotrans

# use pansat package for hourly download 
from datetime import datetime
from pansat.products.reanalysis.era5 import ERA5Product

# specific gas constant for dry air 
R = 287.058

pressure_vars = ['temperature', 'geopotential','specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content',  'specific_humidity', 'u_component_of_wind','v_component_of_wind']
srfc_vars = ['surface_pressure']
srfc_vars_monthly = ['vertical_integral_of_northward_water_vapour_flux', 'vertical_integral_of_eastward_water_vapour_flux']
domain = [10, 60, 50, 130]


# temporary directory for hourly files 
#os.mkdir('tmpdir')

# create product instance
srfc_hourly = ERA5Product('hourly','surface', srfc_vars, domain)
pressure_hourly = ERA5Product('hourly','pressure', pressure_vars,domain)
srfc_monthly = ERA5Product('monthly', 'surface', srfc_vars_monthly, domain)
pressure_monthly = ERA5Product('monthly','pressure', pressure_vars, domain)

for year in np.arange(1979,2020):
    for month in np.arange(5,10):
        if month == 6 or month == 9:
            d= 30
        else:
            d= 31

        if month < 10:
            m = '0' + str (month)
        else:
            m = str(month)
            
        t_0 = datetime(year, month, 1, 0)
        t_1 = datetime(year, month, d, 23)

        # download also monthly mean data 
        monthly = pressure_monthly.download(t_0, t_1, destination = 'tmpdir/monthly')
        srfc_monthly.download(t_0, t_1, destination = 'tmpdir/monthly')
        
        if len(glob.glob('tmpdir/reanalysis-era5-single-levels_' + str(year) + m + '*surface_pressure*.nc')) > 700:
            pressure_files = glob.glob('tmpdir/reanalysis-era5-pressure-levels_'+ str(year)+m +'???*geopotential*.nc')
            srfc_files = glob.glob('tmpdir/reanalysis-era5-single-levels_'+ str(year) + m + '*_surface_pressure*.nc')
            monthly = glob.glob('tmpdir/monthly/reanalysis-era5-pressure-levels*'+ str(year)+m +'*geopotential*.nc')


        if os.path.isfile('tmpdir/processed/qu-int' +str(year)+m+ '.nc') == False: 
            srfc_files = srfc_hourly.download(t_0, t_1, destination= 'tmpdir')
            pressure_files = pressure_hourly.download(t_0, t_1, destination = 'tmpdir')

            # hourly files of one month  
            srfc_files.sort()
            pressure_files.sort()

            assert len(srfc_files) == len(pressure_files)
            assert len(pressure_files) > 0 

            qu_integral = np.zeros((201,321))
            qv_integral = np.zeros((201,321))

            qu_prim = np.zeros((37,201,321))
            qv_prim = np.zeros((37,201,321))

            mfd = np.zeros((201,321))

            # loop through all files in one month 
            for i in np.arange(len(pressure_files)):
                # one hour file at the time 
                data= xr.open_dataset(pressure_files[i])
                srfcdata = xr.open_dataset(srfc_files[i])
                sp = srfcdata.sp.values[0] /100
                meandata= xr.open_dataset(monthly[0])

                # extract monthly means
                u_mean = meandata.u.values[0]
                v_mean = meandata.v.values[0]
                q_mean = meandata.q.values[0] + meandata.ciwc.values[0] + meandata.clwc.values[0]

                # extract variables from xarray dataset 
                t = data.t.values[0]
                u = data.u.values[0]
                v= data.v.values[0]
                z= data.z.values[0]
                # include cloud particles in total moisture 
                ciwc= data.ciwc.values[0]
                clwc = data.clwc.values[0]
                q= data.q.values[0] + ciwc + clwc 
                pressure = data.level.values

                # convert specific humidity to absolute humidity in kg/m3
                for plev in np.arange(37):
                    p_d = (pressure[plev] * 100)/(R*t[plev])
                    q[plev] *= p_d

                # components of moisture flux  
                qu = q*u
                qv = q*v

                # set geopotential to 0, where surface pressure < 1000 hpa (needed for column integration) 
                coords = np.where(sp < 1000)

                for i, ilat in enumerate(coords[0]):
                    ilon = coords[1][i]
                    sp_value = sp[ilat,ilon]
                    idx, pl = atmotrans.find_nearest_idx(pressure, sp_value)
                    if sp_value > pl:
                        idx = idx + 1     
                    # set q value below ground to 0 
                    qu[idx:36, ilat, ilon] = 0
                    qv[idx:36, ilat, ilon] = 0

                # integral of hourly values
                qu_integral += atmotrans.column_integration(np.flip(qu, axis= 0), np.flip(z, axis = 0), ax = 0)    
                qv_integral += atmotrans.column_integration(np.flip(qv, axis= 0), np.flip(z, axis = 0), ax = 0)

                # calculate hourly deviations from monthly mean
                q_dev = q - q_mean 
                u_dev = u - u_mean 
                v_dev = v - v_mean 

                qu_prim += (q_dev * u_dev)
                qv_prim += (q_dev * v_dev)

                data.close()
                srfcdata.close()

            # save hourly-based integral of qV for one month
            xr.DataArray(qu_integral/744).to_netcdf('tmpdir/processed/qu-int' + str(year) + m + '.nc')
            xr.DataArray(qv_integral/744).to_netcdf('tmpdir/processed/qv-int' +str(year) + m + '.nc')

            # save monthly average of primes (hourly deviations from monthly mean)
            xr.DataArray(qu_prim/744).to_netcdf('tmpdir/processed/qu-prim' + str(year) + m + '.nc')
            xr.DataArray(qv_prim/744).to_netcdf('tmpdir/processed/qv-prim' +str(year) + m + '.nc')


            # remove hourly files in month to save space  
            for i in np.arange(len(pressure_files)):
                os.remove(pressure_files[i])
                os.remove(srfc_files[i])

            print('hourly integral calculated for year' +str(year) + 'and month ' + m)
