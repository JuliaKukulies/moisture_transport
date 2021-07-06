## This python script downloads hourly data from ERA5 at pressure levels, calculates deviations from its monthly mean and saves the monthly mean of q and V as well as the product of deviations/primes at hourly, 6-hourly and daily time scale. These output files can be used to decompose the total moisture divergence into its main components (e.g. mean flow and transient eddies). 

import numpy as np
import xarray as xr 
import glob
import os
# use pansat package for hourly download 
from datetime import datetime
from pansat.products.reanalysis.era5 import ERA5Monthly, ERA5Hourly

pressure_vars = ['specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content',  'specific_humidity', 'u_component_of_wind','v_component_of_wind']
domain = [10, 60, 50, 130]

# create product instance
pressure_hourly = ERA5Hourly('pressure', pressure_vars, domain)
pressure_monthly = ERA5Monthly('pressure', pressure_vars,domain)
nr_levels = 37
nr_lats = 201
nr_lons = 321 


for year in np.arange(1980,2020):
    for month in np.arange(7,10):
        if month == 6 or month == 9:
            day= 30
            hours = 30*24
        else:
            day= 31
            hours = 31*24 

        if month < 10:
            m = '0' + str (month)
        else:
            m = str(month)
            
        t_0 = datetime(year, month, 1, 0)
        t_1 = datetime(year, month, day, 23)

        # check whether processed files not already exist 
        if os.path.isfile('tmpdir/processed/qu_prime' +str(year)+ m + '.nc') == False:
            # download monthly and hourly data only if files are not already there 
            if len(glob.glob('tmpdir/hourly/reanalysis-era5-pressure-levels_' + str(year) + m + '*.nc')) > 500:
                hourly_files = glob.glob('tmpdir/hourly/reanalysis-era5-pressure-levels_'+ str(year)+ m +'*.nc')
                monthly_files = glob.glob('tmpdir/monthly/reanalysis-era5-pressure-levels-monthly-means_'+ str(year) + m + '*.nc')
                print(len(hourly_files), ' files found for', str(year), m)
            else:
                hourly_files = pressure_hourly.download(t_0, t_1, destination = 'tmpdir/hourly')
                monthly_files = pressure_monthly.download(t_0, t_1, destination = 'tmpdir/monthly')
                print('all files downloaded.')
                
            # sort hourly files of one month  
            hourly_files.sort()
            assert len(hourly_files) > 500
            assert len(monthly_files) == 1 

            print('starting preprocessing of hourly files for :',str(year), m)
            # initialize arrays for calculations
            
            qu_prim = np.zeros((nr_levels, nr_lats, nr_lons))
            qv_prim = np.zeros((nr_levels, nr_lats, nr_lons))
            qu_prim6 = np.zeros((nr_levels, nr_lats, nr_lons))
            qv_prim6 = np.zeros((nr_levels, nr_lats, nr_lons))           
            qu_prim_d = np.zeros((nr_levels, nr_lats, nr_lons))
            qv_prim_d = np.zeros((nr_levels, nr_lats, nr_lons))
             
            # loop through all files in one month (on hourly basis)
            for i in np.arange(len(hourly_files)):
                # control when 6 hours have passed 
                if np.mod(i,6) == 0:
                    if i!= 0:
                        # calculate deviations from monthly mean
                        q_dev6 = q6 / 6 - q_mean 
                        u_dev6 = u6 / 6 - u_mean 
                        v_dev6 = v6 / 6 - v_mean 
                        qu_prim6 += (q_dev6 * u_dev6)
                        qv_prim6 += (q_dev6 * v_dev6)
                                    
                    # initiate 6-hourly variables
                    u6= np.zeros((nr_levels, nr_lats, nr_lons))
                    v6= np.zeros((nr_levels, nr_lats, nr_lons))
                    q6= np.zeros((nr_levels, nr_lats, nr_lons))

                # when one day (24 hours) have passed 
                if np.mod(i,24)== 0:
                    if i!= 0:
                        # calculate deviations from monthly mean
                        q_dev_d = q_d / 24 - q_mean 
                        u_dev_d = u_d / 24 - u_mean 
                        v_dev_d = v_d / 24 - v_mean 
                        qu_prim_d += (q_dev_d * u_dev_d)
                        qv_prim_d += (q_dev_d * v_dev_d)
                    
                    # initiate daily variables
                    u_d= np.zeros((nr_levels, nr_lats, nr_lons))
                    v_d= np.zeros((nr_levels, nr_lats, nr_lons))
                    q_d= np.zeros((nr_levels, nr_lats, nr_lons))

                ###########################processing hourly files ################################################ 
                data= xr.open_dataset(hourly_files[i])
                meandata= xr.open_dataset(monthly_files[0])
                
                # extract monthly means
                u_mean = meandata.u.values[0]
                v_mean = meandata.v.values[0]
                q_mean = meandata.q.values[0] + meandata.ciwc.values[0] + meandata.clwc.values[0]

                # extract variables from hourly data  
                u = data.u.values[0]
                v = data.v.values[0]
                               
                # include cloud particles and frozen particles in total atmospheric moisture budget  
                ciwc = data.ciwc.values[0]
                clwc = data.clwc.values[0]
                q = data.q.values[0] + ciwc + clwc
                pressure = data.level.values

                # add to 6-hourly and daily eddies 
                q_d += q
                u_d += u
                v_d += v 
                q6 += q
                u6 += u
                v6 += v
                
                # calculate hourly deviations from monthly mean
                q_dev = q - q_mean 
                u_dev = u - u_mean 
                v_dev = v - v_mean 

                qu_prim += (q_dev * u_dev)
                qv_prim += (q_dev * v_dev)

                data.close()

            meandata.close()
                
                
            # save monthly average of primes (hourly deviations from monthly mean)
            xr.DataArray(data = qu_prim/hours, name = 'primeproduct_qu').to_netcdf('tmpdir/processed/qu_prim' + str(year) + m + '.nc')
            xr.DataArray(data = qv_prim/hours, name = 'primeproduct_qv').to_netcdf('tmpdir/processed/qv_prim' +str(year) + m + '.nc')
            
            # save monthly average of 6-hourly eddies
            xr.DataArray(data = qu_prim6/(hours/6), name = 'primeproduct_qu').to_netcdf('tmpdir/processed/qu-prim-6hr-' + str(year) + m + '.nc')
            xr.DataArray(data = qv_prim6/(hours/6), name = 'primeproduct_qv' ).to_netcdf('tmpdir/processed/qv-prim-6hr-' +str(year) + m + '.nc')
           
            # save monthly average of daily eddies 
            xr.DataArray(data= qu_prim_d/(hours/24), name = 'primeproduct_qu').to_netcdf('tmpdir/processed/qu-prim-daily-' + str(year) + m + '.nc')
            xr.DataArray(data = qv_prim_d/(hours/24), name = 'primeproduct_qv').to_netcdf('tmpdir/processed/qv-prim-daily-' +str(year) + m + '.nc')
                                    
            # remove hourly files in month to save space           
            print('all subseasonal eddies calculated for year' +str(year) + 'and month ' + m)

            for f in hourly_files:
                os.remove(f)
