# import libraries
import netCDF4 as nc
import os
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.time import Time
from matplotlib.dates import DateFormatter
import pandas as pd


def plot_irradianceSOL(file_name):
    with nc.Dataset(file_name, 'r', format='NETCDF4') as file_grp:

        # assign names to variables required
        irradiance = file_grp.variables['irradiance']
        time = file_grp.variables['microsecondssincegpsepoch'][:]
        wavelength = file_grp.variables['wavelength']
        tanheight = file_grp.variables['tanht']

        # convert time matrix from microseconds to seconds
        time /= 1000000

        # find dy of wavelength
        dx = 1
        y = wavelength[:]
        dy = np.diff(y) / dx

        """" Find wavelength boundaries for each scan to pull out each individual irradiance scan """

        # preallocate a matrix for resulting values of 'dy' loop
        scanboundaries = []
        scanboundaries.append([0, 0])

        # begin loop to find where curve's slope is negative using derivative of wavelength
        for i in range(1, dy.size):
            if (dy[i] * dy[i - 1]) < 0:
                scanboundaries[-1][1] = i - 1
                scanboundaries.append([i, dy.size - 1])

        """Plot each individual irradiance scan"""

        # set a counter for this set of plots
        count = 0
        # pre allocate lists for resulting values
        irradiancelist = []
        wavelengthlist = []
        timelist = []
        tanheightlist = []

        # begin loop using only values within our scan boundaries
        for scanboundary in scanboundaries:
            waveboundary = wavelength[scanboundary[0]:scanboundary[1]]
            irrboundary = irradiance[scanboundary[0]:scanboundary[1]]
            timeboundary = time[scanboundary[0]:scanboundary[1]]
            tanheightboundary = tanheight[scanboundary[0]:scanboundary[1]]

            irradiancebin = []
            wavelengthbin = []
            timebin = []
            tanheightbin = []

            # find all irradiance values in each 0.1 wavelength sections
            for wave in np.arange(120.5, 122.0, 0.1):
                indices = np.where((waveboundary >= wave - 0.05) & (wave + 0.05 >= waveboundary))
                countindices = indices[0].size
                wavelengthbin.append(wave)

                # if no values found, add nan and append to bin
                if countindices == 0:
                    irradiancebin.append(np.nan)
                    timebin.append(np.nan)
                    tanheightbin.append(np.nan)
                # if values found then average all values within 0.1 bin and add to irradiancebin
                else:
                    irradiancebin.append(np.mean(irrboundary[indices]))
                    timebin.append(np.mean(timeboundary[indices]))
                    tanheightbin.append(np.mean(tanheightboundary[indices]))

            # make bins into arrays
            irradiancebin = np.array(irradiancebin, dtype=np.double)
            wavelengthbin = np.array(wavelengthbin, dtype=np.double)
            timebin = np.array(timebin, dtype=np.double)
            tanheightbin = np.array(tanheightbin, dtype=np.double)

            # make a list of all arrays
            irradiancelist.append(irradiancebin)
            wavelengthlist.append(wavelengthbin)
            timelist.append(timebin)
            tanheightlist.append(tanheightbin)

            plt.figure()

            plt.plot(waveboundary, irrboundary, color='red', linestyle='dashed', label='irradiance over wavelength')
            plt.plot(wavelengthbin, irradiancebin, color='green', label='average irradiance over wavelength')
            plt.legend(loc=9)
            plt.yscale("log")
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Irradiance")
            plt.savefig(
                "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_practice\\venv\\plots--solstice_occultation20040706T181400\\"
                "individual_scans\\Scan" + str(count) + ".png")
            plt.close()
            count += 1

        wavelengthlist = np.vstack(wavelengthlist)
        irradiancelist = np.vstack(irradiancelist)
        timelist = np.vstack(timelist)
        tanheightlist = np.vstack(tanheightlist)

        """ Plot time series for each individual wavelength"""

        # loop through list of irradiance and time arrays to plot each element
        for i in np.arange(wavelengthlist[0, :].size):
            plt.plot(timelist[:, i], irradiancelist[:, i])
            plt.title('Wavelength ' + "{:.2f}".format(wavelengthlist[0, i]))
            plt.xlabel("Time (GPSseconds)")
            plt.ylabel("Irradiance")
            plt.savefig(
                "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_practice\\venv\plots--solstice_occultation20040706T181400\\"
                "time_series_plot\\wavelength-" + "{:.2f}".format(wavelengthlist[0, i]) + ".png")
            plt.close()

        # loop through list of irradiance and tangent height arrays to plot each element
        for i in np.arange(wavelengthlist[0, :].size):
            plt.plot(tanheightlist[:, i], irradiancelist[:, i])
            plt.title('Wavelength ' + "{:.2f}".format(wavelengthlist[0, i]))
            plt.xlabel("Tangent Height (km)")
            plt.ylabel("Irradiance")
            plt.savefig(
                "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_practice\\venv\\heightplot\\wavelength-" + "{:.2f}".format(
                    wavelengthlist[0, i]) + ".png")
            plt.close()

            x = np.where((tanheightlist[:, i] >= 200) & (250 >= tanheightlist[:, i]))
            if x[0].size > 0:
                ave_irradiance = np.mean(irradiancelist[:, i][x])
                ########print(ave_irradiance)
                #########print(wavelengthlist[0, i])

                extinction_ratio = irradiancelist[:, i] / ave_irradiance

                plt.plot(tanheightlist[:, i], extinction_ratio)
                plt.title("Extinction Ratio at Wavelength-" + "{:.2f}".format(wavelengthlist[0, i]))
                plt.xlabel("Tangent Height (km)")
                plt.ylabel("Extinction Ratio")
                plt.savefig(
                    "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_practice\\venv\\plots--solstice_occultation20040706T181400\\"
                    "extinction_ratio_plot\\wavelength " + "{:.2f}".format(wavelengthlist[0, i]) + ".png")
                plt.close()

        # convert to julian dates
        time = Time(time, format='gps')
        time.format = 'jd'


# call back my function
if __name__ == '__main__':
    file_name = 'C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_practice\\venv\\SOLdataset' \
                '\\solstice_occultation20040706T181400.000.nc '
    plot_irradianceSOL(file_name)
