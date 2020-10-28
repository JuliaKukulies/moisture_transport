## This python script downloads hourly data from ERA5 at pressure levels, calculates deviations from its monthly mean and saves the monthly mean of the sum of mean qV and its primes.


import numpy as np
import xarray as xr 
import glob
import cdsapi


# download hourly files




files = glob.glob('cache/1979/era5-1979-??-hour??.nc')
quvec = np.zeros(q.shape)
qvvec = np.zeros(q.shape)

# get mean data 
f = 'cache/era5_1979_may.nc'
# get variables 
flds = xr.open_dataset(f)
u_mean = flds.u[0]
v_mean = flds.v[0]
# specific humidity (water vapour in g/kg)
q_mean = flds.q[0] 

# convert specific humidity to absolute humidity in kg/m3
for plev in np.arange(37):
    p_d = (pressure[plev] * 100)/(R*temp[plev])
    q_mean[plev] *= p_d



# loop through hourly files 
for f in files:
    # one hour file at the time 
    data= xr.open_dataset(f)
    u = data.u.values[0]
    v= data.v.values[0]
    q= data.q.values[0]
        
    # convert specific humidity to absolute humidity in kg/m3
    for plev in np.arange(37):
        p_d = (pressure[plev] * 100)/(R*temp[plev])
        m = q[plev] *p_d
        q[plev] = m 
        
    # calculate deviations 
    u_prime = u - u_mean.values
    v_prime = v - v_mean.values
    q_prime = q - q_mean.values

    #  calculate product of primes 
    uq = u_prime * q_prime 
    vq = v_prime * q_prime 
        
    # monthly mean 
    qu = q*u
    qv = q*v

    # save monthly mean deviation from mean 
    quvec += (qu + uq)
    qvvec += (qv + vq)

    # save as xarray 
    xr.DataArray(quvec/744).to_netcdf('cache/qu-primes_1979_05.nc')
    xr.DataArray(qvvec/744).to_netcdf('cache/qv-primes_1979_05.nc')



# rmeove hourly files 




    
