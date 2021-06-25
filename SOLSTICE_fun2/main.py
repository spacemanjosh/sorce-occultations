#import libraries
import netCDF4 as nc
import os
import numpy as np
import astropy
from astropy.time import Time



def plot_irradianceSOL(file_name):

    with nc.Dataset(file_name, 'r', format='NETCDF4') as file_grp:

        #assign names to variables required
        irradiance = file_grp.variables['irradiance']
        time = file_grp.variables['microsecondssincegpsepoch'][:]
        wavelength = file_grp.variables['wavelength']

        #convert time matrix from microseconds to seconds
        time /= 1000000

        #convert to julian dates
        time = Time(time, format='gps')
        time.format='jd'
        print(time)



#call back my function
if __name__ == '__main__':
    file_name = 'C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_fun2\\venv\\dataset\\solstice_occultation20040706T181400.000.nc'
    plot_irradianceSOL(file_name)