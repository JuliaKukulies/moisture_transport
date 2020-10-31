## This python script downloads hourly data from ERA5 at pressure levels, calculates deviations from its monthly mean and saves the monthly mean of the sum of mean qV and its primes.

import numpy as np
import xarray as xr 
import glob
import cdsapi
import os


# use pansat package for hourly download 
from datetime import datetime
from pansat.products.reanalysis.era5 import ERA5Product

pressure_vars = ['temperature', 'geopotential','ciwc', 'clwc' 'specific_humidity', 'u_component_of_wind','v_component_of_wind']
srfc_vars = ['surface_pressure']

# temporary directory for hourly files 
#os.mkdir('tmpdir')

# create product instance
srfc_hourly = ERA5Product('hourly','surface', variables= srfc_vars, domain = ['60', '50', '10','130'])
pressure_hourly = ERA5Product('hourly','pressure' variables = pressure_vars, domain = ['60', '50', '10','130'])
                              
for year in np.arange(1979,2020):
    for month in np.arange(5,10):
        t_0 = datetime(year, month, 1, 0)
        t_1 = datetime(year, month, 31, 23)
    
        srfc_files = srfc_hourly.download(t_0, t_1, destination= 'tmpdir')
        pressure_files = pressure_hourly.download(t_0, t_1, destination = 'tmpdir')

        # hourly files 
        srfc_files.sort()
        pressure_files.sort()
        
        qu_integral = np.zeros((201,321))
        qv_integral = np.zeros((201,321))
        mfd = np.zeros((201,321))

        
        for i, f in enumerate(pressure_files): 
            # one hour file at the time 
            data= xr.open_dataset(f)
            srfcdata = xr.open_dataset(srfc_files[i])
            sp = srfcdata.sp.values[0] /100

            # extract variables from xarray dataset 
            t = data.t.values[0]
            u = data.u.values[0]
            v= data.v.values[0]
            z= data.z.values[0]
            # include cloud particles in total moisture 
            ciwc= data.ciwc.values[0]
            clwc = data.clwc.values[0]
            q= data.q.values[0] + ciwc + clwc 
            pressure = data.levels.values

            # convert specific humidity to absolute humidity in kg/m3
            for plev in np.arange(37):
                p_d = (pressure[plev] * 100)/(R*t[plev])
                q[plev] *= p_d

            # monthly mean 
            qu = q*u
            qv = q*v

            # set geopotential to 0, where surface pressure < 1000 hpa 
            coords = np.where(sp < 1000)

            for i, ilat in enumerate(coords[0]):
                ilon = coords[1][i]
                sp_value = sp[ilat,ilon]
                idx, pl = find_nearest_idx(pressure, sp_value)
                if sp_value > pl:
                    idx = idx + 1     
                # set q value below ground to 0 
                qu[idx:36, ilat, ilon] = 0
                qv[idx:36, ilat, ilon] = 0

            # integral of hourly values
            qu_integral += column_integration(np.flip(qu, axis= 0), np.flip(z, axis = 0), ax = 0)    
            qv_integral += column_integration(np.flip(qv, axis= 0), np.flip(z, axis = 0), ax = 0)    

            data.close()
            srfcdata.close()

            # save hourly-based integral for one months
            if month < 10:
                m = '0' + str (m)
            else:
                m = str(month)
                              
            xr.DataArray(qu_integral/744).to_netcdf('tmpdir/qu-int' + str(year) + m +'.nc')
            xr.DataArray(qv_integral/744).to_netcdf('tmpdir/qv-int' +str(year) + m+ '.nc')

        # remove hourly files to save space  
        for i, f in enumerate(pressure_files):
            os.remove(f)
            os.remove(srfc_files[i])

        print('hourly integral calculated for year' +str(year) + 'and month ' + m)

