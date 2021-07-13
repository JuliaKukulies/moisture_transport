''' This python module collects functions for the calculation of atmospheric water vapor transport and moisture divergence. '''

# import libraries 
import numpy as np 
import xarray as xr
import metpy
from metpy import calc
import scipy as sp
import scipy.signal
from scipy.signal import convolve
from metpy.units import units

############################# CONSTANTS##############################


# assign unit to grid spacing
Rad = 6371*1000
# gravitational accelration 
g = 9.8 
# density of water in kg/m3 
pw= 997
# density for dry air 
pd = 1.225 
# specific gas constant for dry air 
R = 287.058
# constant for unit in mm per day 
C= -1/(g*pw)
c= -1/(g)

# gas constant for water vapour in J K-1 kg-1
Rvap= 461

# constants for Tetens formula (for saturation over water)
c1= 611.21
c2= 17.502
c3= 32.19
# freezing point
T0 = 273.16 

############################# BASIC CALCULATIONS ##########################


def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def geopotential_to_height(z):
    """ This function converts geopotential heights to geometric heights. This approximation takes into account the varying gravitational force with heights, but neglects latitudinal vairations.
    Parameters:
    ------------
    z(float) : (1D or multi-dimenstional) array with geopotential heights
    Returns:
    ----------
    geometric_heights : array of same shape containing altitudes in metres
    """
    g = 9.80665 # standard gravity 
    Re = 6.371 * 10**6  # earth radius
    geometric_heights   = (z*Re) / (g * Re - z)
    return geometric_heights 



def get_surface_humidity(temperature, spressure):
    '''
    This function calculates near-surface humidity for ERA5 
    based on the 2m dew point temperature and surafce pressure. 

    Args: 
    temperature: 2D array with 2m dew point temperature 
    spressure: 2D array with surface pressure values in hpa

    Returns:
    q_sat: near surface humidity in kg/kg 

    '''
    
    #### define constants #### 

    # gas constants for dry air and water vapour in J K-1 kg-1
    Rdry= 287
    Rvap= 461
    # constants for Tetens formula (for saturation over water)
    c1= 611.21
    c2= 17.502
    c3= 32.19
    # freezing point
    T0 = 273.16 
    
    spressure = spressure*100
    e_sat = c1* np.exp( c2 * ((temperature - T0)/ (temperature - c3)))
    
    q_sat = ((Rdry / Rvap) * e_sat ) / (spressure - (1- Rdry/Rvap) * e_sat )
    return q_sat




def colint_pressure(values,pressure_levels):
    """ This function calculates the column-integrated water vapor
    in kg/m2 from specific humidity (kg/kg) at different hpa levels.

    """
    return np.trapz(values, pressure_levels, axis = 0)* g



def column_integration_height(values, z, ax = None ):
    """This functions calculates the column-integrated value of a given atmospheric variable at different pressure levels
    Parameters:
    -----------
    values(float): 1D or multi-dimensional array with values of atmospheric variable at different pressure levels
    z(int): array with geopotential heights for values
    axis = axis along which to integrated. The default is 0.
    Returns:
    --------
    colint(float): array with column-integrated values of variable (dimension reduced by 1)
    """
    # convert geopotential to geometric heights in meters 
    geometric_heights   = geopotential_to_height(z)

    if ax == None:
        ax = 0

    # integration of column values
    colint = np.trapz(values, x= geometric_heights, axis =ax )

    return colint



def column_integration(levels, sp, var):
    """

    This function integrates over vertical pressure levels in ERA5 after setting values
    that are below surface pressure to 0 and replacing the surface pressure with 
    values above the maximum pressure level 1000 hpa with extrapolated values. 

    Args: 
    levels: 1D array that contains pressure coordinates
    sp: 2D field with surface pressures 
    war: 3D field with variable to integrate 

    Returns:
    colint: integrated 2D field 

    """
    from scipy import interpolate
    import wrf 
    coords = np.where(sp < 100000)
    pressure = np.zeros(var.shape)
    for i, ilat in enumerate(coords[0]):
        ilon = coords[1][i]
        sp_value = sp[ilat,ilon]

        pressure[:,ilat,ilon]= levels
        pressure[levels.size-1, :, :] =  sp
        idx, pl = find_nearest_idx(levels, sp_value)

        # function for extrapolation/ interpolation: 
        x_vals = levels
        y_vals= var[:,ilat,ilon]
        f = interpolate.interp1d(x_vals, y_vals, fill_value = "extrapolate", kind = 'cubic')

        # set q value below ground to 0 
        if sp_value < 1000:
            pressure[idx, ilat,ilon] = sp_value
            var[idx, ilat,ilon] = f(sp_value)
            var[idx+1:levels.size, ilat, ilon] =  0

        if sp_value > 1000:
            var[levels.size-1, ilat, ilon] = f(sp_value)

    colint = colint_pressure(var, pressure)

    return colint




def column_integration_q(data, sp, var, temp):
    """

    This function integrates specific humidity  over vertical pressure levels in ERA5 
    after interpolating to levels inbetween two pressure levels, setting levels 
    below surface pressure to 0 and replacing the surface pressure with 
    values above the maximum pressure level 1000 hpa with the calculated surface humidity. 

    Args: 
    data: xarray that contains pressure coordinates
    sp: 2D field with surface pressures 
    q: 3D field of specific humidity at pressure levels 
    temp: 2D field with 2m dew point temperature 

    Returns:
    colint: integrated 2D field 

    """
    from scipy import interpolate
    import wrf 
    coords = np.where(sp < 10000)
    pressure = np.zeros((37,201,321))
    for i, ilat in enumerate(coords[0]):
        ilon = coords[1][i]
        sp_value = sp[ilat,ilon]

        pressure[:,ilat,ilon]= data.level.values
        pressure[36, :, :] =  sp
        idx, pl = find_nearest_idx(data.level.values, sp_value)

        # function for extrapolation/ interpolation: 
        x_vals = data.level.values
        y_vals= qu[:,ilat,ilon]
        f = interpatm.olate.interp1d(x_vals, y_vals, fill_value = "extrapolate", kind = 'cubic')

        # set q value below ground to 0 
        if sp_value < 1000:
            if sp_value > pl:
                idx = idx + 1  
            pressure[idx, ilat,ilon] = sp_value
            var[idx:37, ilat, ilon] =  0

        if sp_value > 1000:
            var[36, ilat, ilon] = f(sp_value)

    colint = atm.colint_pressure(var, pressure)

    return colint





def total_integrated_moisture_flx(qu, qv):
    """
    Returns 2D field with total column-integrated water vapour flux, given: 

    qu: 2D field with column-integrated moisture flux u -component 
    qv: 2D field with column-integrated moisture flux v -component    

    """
    return np.sqrt(qu **2 + qv **2)



def weighted_mean(arr, data):
    """

    This function calculates the area-weighted mean over a 2D field. 
    
    Args: 
    arr: 2D array with variable that should be averaged
    data: xarray with latitude and longitude coordinates

    Returns: 
    weighted_mean: scalar that is the weighted mean over the area

    """
    dataset=xr.DataArray(arr,  dims= {'latitude':data.latitude[:-1].values, 'longitude':data.longitude[:-1].values})
    weights = np.cos(np.deg2rad(data.latitude[:-1]))
    weights.name = "weights"
    data_weighted = dataset.weighted(weights)
    weighted_mean = data_weighted.mean(("latitude", "longitude"), skipna= True)
    return weighted_mean.values


def weighted_mean_timeseries(arr, data):
    """
    This function calculates the area-weighted mean over a 3D field. 
    
    Args: 
    arr: 3D array with variable that should be averaged over space
    data: xarray with latitude and longitude coordinates

    Returns: 
    weighted_mean: timeseries with weighted means over the area

    """   
    dataset=xr.DataArray(arr,  dims= {'time':data.latitude[:-1].values,'latitude':data.latitude[:-1].values, 'longitude':data.longitude[:-1].values})
    weights = np.cos(np.deg2rad(data.latitude[:-1]))
    weights.name = "weights"
    data_weighted = dataset.weighted(weights)
    weighted_mean = data_weighted.mean(("latitude", "longitude"), skipna= True)
    return weighted_mean.values










############################## MOISTURE DIVERGENCE##############################################

def get_spacing(lats, lons):
    '''This functions calculates the grid spacing in meter.


    Args:
    lats(numpy array): 1D array with latitudes
    lons(numpy array): 1D array with longitudes


    Returns:
    dlat(numpy array): latitude spacings in m 
    dlon(numpy array): longitude spacings in m 

'''

    # creating 2D fields for lats and lons 
    latitudes = np.stack([lats]*np.shape(lons)[0], axis = 1)
    longitudes = np.stack([lons]*np.shape(lats)[0], axis = 0)

    # convert lats and lons to cartesian coordinates 
    x = Rad * np.cos(np.radians(latitudes)) * np.cos(np.radians(longitudes))
    y = Rad * np.cos(np.radians(latitudes)) * np.sin(np.radians(longitudes))
    z = Rad *np.sin(np.radians(latitudes))

    # stack to get 3D array 
    cartesian = np.stack([x, y, z], axis = 2)
    # pythagorean theorem to get distances in meter
    dlat = np.sqrt(np.sum((cartesian[2:, :] - cartesian[:-2,:]) ** 2, axis=-1))
    dlon = np.sqrt(np.sum((cartesian[:, 2:] - cartesian[:,:-2]) ** 2, axis=-1))

    return dlat, dlon


def get_delta(lats, Rad):
    ## grid spacings 
    dx = 2*np.pi*Rad * (0.25/360)
    # latitude dependent 
    dy = 2*np.pi*Rad *(0.25/360) * np.cos(np.nanmean(lats))*(-1)
    dx = dx * units.meters
    dy = dy * units.meters

    return dx, dy 


def derivative_u(quint,dlon):
    """

    This function calculates the derivative in v direction in the spectral space using FFT.


    Args:

    quint(np.array): 2D field of integrated water vapor flux in u direction 
    dlon(np.array): 2D field of longitude spacings (accounting for different distances dependent on latitude)

    Returns:

    2D field with first derivative of vertically integrated water vapor flux in u direction 

    """
    quint_padded = np.hstack([np.fliplr(quint[:-1, :-1]), quint[:-1, :-1], np.fliplr(quint[:-1, :-1])]) 
    f_quint = np.fft.fft(quint_padded, axis=1)

    m, n = f_quint.shape
    m2 = m // 2
    n2 = n // 2
        
    f_lon = (2.0 * np.pi * np.fft.fftfreq(n, d= dlon[:-1,[0]]/2) )
    f_lon[:, n2] = 0.0
    #f_lon = np.broadcast_to(f_lon.reshape(1, -1), (m, n)) 
    
    df_quint_dx = f_quint.copy() * 1j * f_lon 
    
    d_n = 50
    df_quint_dx[:, n2 - d_n : n2 + d_n + 1] *= 0.0
    
    real = np.fft.ifft(df_quint_dx, axis = 1).real

    return real[ :, quint[:-1, :-1].shape[1]: quint[:-1, :-1].shape[1] * 2]





def derivative_v(qvint,dlat):
    """

    This function calculates the derivative in v direction in the spectral space using FFT.


    Args:

    qvint(np.array): 2D field of integrated water vapor flux in v direction 
    dlat(np.array): 2D field of latitude spacings 

    Returns:

    2D field with first derivative of vertically integrated water vapor flux in v direction 

    """
    qvint_padded = np.vstack([np.flipud(qvint[:-1, :-1]), qvint[:-1, :-1], np.flipud(qvint[:-1, :-1])]) 
    f_qvint = np.fft.fft(qvint_padded, axis=0)

    m, n = f_qvint.shape
    m2 = m // 2
    n2 = n // 2
        
    f_lat = 2.0 * np.pi * np.fft.fftfreq(m, d= dlat[0,0]/2)
    f_lat[m2] = 0.0
    f_lat = np.broadcast_to(f_lat.reshape(-1, 1), (m, n)) 
    
    df_qvint_dy = f_qvint.copy() * -1j * f_lat
    
    d_m = 60
    df_qvint_dy[m2 - d_m : m2 + d_m + 1, :] *= 0.0
    real = np.fft.ifft(df_qvint_dy, axis = 0).real

    return real[qvint[:-1, :-1].shape[0]: qvint[:-1, :-1].shape[0] * 2, :]


def dy_dlat(y, dlat):
    '''This functions calculates the derivative along latitudes of a variable y using a finite differential method.

    Args:
    y(numpy array): atmospheric variable, e.g. u, v, q, qu, qv and so on
    dlat(np.array): grid spacing in meter, should be more or less constant'''

    k_lat = np.array([[-1], [0], [1]])
    result = convolve(y, k_lat, mode="valid") / dlat
    return result 


def dy_dlon(y, dlon):
    '''This functions calculates the derivative along longitudes of a variable y using a finite differential method.

    Args:
    y(numpy array): atmospheric variable, e.g. u, v, q, qu, qv and so on
    dlon(np.array): grid spacing in meter, should be varying dependent on latitude '''
    k_lon= np.array([[1, 0, -1]])
    result = convolve(y, k_lon, mode="valid") / dlon
    return result




def get_surface_values(field, nlat, nlon,nlev, surface_pressures,pressure_levels):
    """
    This function reduces a 3D meteorological field to two dimensions given the surface pressures. 

    Args:
    field(numpy array): 3 dimensionsal field of meteorological variable
    nlat(int): number of latitudes (2nd dimension)
    nlon(int): number of longitudes (3rd dimension)
    nlev(int): number of levels (1st dimension )
    surface_pressures(numpy array): 2D field with surface pressure values, must have same latitude and longitudes as field
    pressure_levels(numpy array): 1D array with pressure levels of field variable 


    Returns: 2D array with values at surface 
    """
    for lat in np.arange(nlat):
        for lon in np.arange(nlon):
            idx,pr = find_nearest_idx(pressure_levels,surface_pressures[lat,lon])
            field[np.arange(nlev)!=idx,lat,lon] = 0
            
    return np.nansum(field,axis = 0 )




def divergence(data,qu,qv):
    """

    This function calculates the divergence of a given flux.

    Args:
    data: xarray dataset containing coordinate references 
    qu: u-component of 2D flux field (e.g. moisture flux)
    qv: v-component of 2D flux field 

    Returns: 2D field with divergence of the flux. 
    
    """
    import wrf
    dlat, dlon = get_spacing(data.latitude.values, data.longitude.values)
    udiff= derivative_u(qu, dlon)
    vdiff= derivative_v(qv, dlat)
    conv_total = (udiff + vdiff)

    return wrf.smooth2d(-(vdiff + udiff)*86400 , passes = 3)



def correct_column_integration(data, sp, q, u, v, u10, v10):
    """

    This function performs an interpolation and extrapolation before 
    the vertical column integration over pressure coordinates. Humidity 
    values are interpolated onto the surface pressure field (extrapolated 
    when surface pressure is higher
    than the lowest pressure level in the model). 

    Args: 
    data: xarray dataset with data and coordinate references
    sp: 2D field with corresponding surafce pressures
    q: humidity field on pressure levels  
    u: u wind field on pressure levels 
    v: v wind field on pressure levels 

    u10: 2D field with surface wind u-component
    v10: 2D field with surface wind v-component

    Returns: 
    colint: total column water vapour 
    qu: vertically integrated moisture eastward flux 
    qv: vertically integrated moisture northward flux 


    """
    from scipy import interpolate
    import wrf 
    
    coords = np.where(sp < 10000)
    pressure = np.zeros((q.shape))
    # get humidity value for surface pressures 
    surface_humidity = wrf.interplevel(q, pressure, sp)

    for i, ilat in enumerate(coords[0]):
        ilon = coords[1][i]
        sp_value = sp[ilat,ilon]

        pressure[:,ilat,ilon]= data.level.values
        pressure[36] =  sp
        idx, pl = find_nearest_idx(data.level.values, sp_value)

        # function for extrapolation/ interpolation: 
        x_vals = data.level.values
        y_vals= q[:,ilat,ilon]
        f = interpolate.interp1d(x_vals, y_vals, fill_value = "extrapolate", kind = 'cubic')

        # set q value below ground to 0 
        if sp_value < 1000:
            if sp_value > pl:
                idx = idx + 1  
            q[idx, ilat,ilon] = surface_humidity[ilat,ilon]
            u[idx, ilat,ilon] = u10[ilat,ilon]
            v[idx, ilat,ilon] = v10[ilat,ilon]
            pressure[idx, ilat,ilon] = sp_value
            q[idx:37, ilat, ilon] =  0

        if sp_value > 1000:
            q[36, ilat, ilon] = f(sp_value)
            u[36, ilat, ilon ] = u10[ilat,ilon]
            v[36, ilat, ilon ] = v10[ilat,ilon]

    colint = colint_pressure(q, pressure)
    qu = colint_pressure(q*u, pressure)
    qv= colint_pressure(q*v, pressure)
    
    return colint, qu, qv 






