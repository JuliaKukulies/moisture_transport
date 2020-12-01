''' This python module collects functions for the calculation of atmospheric moisture divergence. '''


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


def column_integration(values, z, ax = None ):
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



############################## MOISTURE DIVERGENCE##############################################

def get_spacing(lats, lons):
    '''This functions calculates the grid spacing in meter.


    Args:
    lats(numpy array): 1D array with latitudes
    lons(numpy array): 1D array with longitudes'''

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


def derivative_u(quint):
    """

    This function calculates the derivative in v direction in the spectral space using FFT.


    Args:

    quint(np.array): 2D field of integrated water vapor flux in u direction 

    Returns:

    2D field with first derivative of vertically integrated water vapor flux in u direction 

    """
    quint_padded = np.hstack([np.fliplr(quint[:-1, :-1]), quint[:-1, :-1], np.fliplr(quint[:-1, :-1])]) 
    f_quint = np.fft.fft(quint_padded, axis=1)

    m, n = f_quint.shape
    m2 = m // 2
    n2 = n // 2
        
    f_lon = (2.0 * np.pi * np.fft.fftfreq(n)) 
    f_lon[n2] = 0.0
    f_lon = np.broadcast_to(f_lon.reshape(1, -1), (m, n)) 
    
    df_quint_dx = f_quint.copy() * 1j * f_lon 
    
    d_n = 50
    df_quint_dx[:, n2 - d_n : n2 + d_n + 1] *= 0.0
    
    return np.fft.ifft(df_quint_dx, axis = 1).real





def derivative_v(qvint):
    """

    This function calculates the derivative in v direction in the spectral space using FFT.


    Args:

    qvint(np.array): 2D field of integrated water vapor flux in v direction 

    Returns:

    2D field with first derivative of vertically integrated water vapor flux in v direction 

    """
    qvint_padded = np.vstack([np.flipud(qvint[:-1, :-1]), qvint[:-1, :-1], np.flipud(qvint[:-1, :-1])]) 
    f_qvint = np.fft.fft(qvint_padded, axis=0)

    m, n = f_qvint.shape
    m2 = m // 2
    n2 = n // 2
        
    f_lat = 2.0 * np.pi * np.fft.fftfreq(m)
    f_lat[m2] = 0.0
    f_lat = np.broadcast_to(f_lat.reshape(-1, 1), (m, n)) 
    
    df_qvint_dy = f_qvint.copy() * -1j * f_lat
    
    d_m = 60
    df_qvint_dy[m2 - d_m : m2 + d_m + 1, :] *= 0.0
    
    return np.fft.ifft(df_qvint_dy, axis = 0).real












def dy_dlat(y, dlat):
    '''This functions calculates the horizontal divergence of a variable y.

    Args:
    y(numpy array): atmospheric variable, e.g. u, v, q, qu, qv and so on
    dlat(np.array): grid spacing in meter, should be more or less constant'''

    k_lat = np.array([[-1], [0], [1]])
    result = convolve(y, k_lat, mode="valid") / dlat
    return result 


def dy_dlon(y, dlon):
    '''This functions calculates the horizontal divergence of a variable y.

    Args:
    y(numpy array): atmospheric variable, e.g. u, v, q, qu, qv and so on
    dlon(np.array): grid spacing in meter, should be varying dependent on latitude '''
    k_lon= np.array([[1, 0, -1]])
    result = convolve(y, k_lon, mode="valid") / dlon
    return result


