# import libraries
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from astropy.time import Time
from matplotlib.dates import DateFormatter
from scipy import interpolate


# insert link here
def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, 'same')
    return y_smooth


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

            plt.rcParams["figure.figsize"] = [11,5]

            plt.plot(waveboundary, irrboundary, color='red', linestyle='dashed', label='irradiance over wavelength')
            plt.plot(wavelengthbin, irradiancebin, color='green', label='average irradiance over wavelength')
            plt.legend(loc=1)
            plt.yscale("log")
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Irradiance")
            plt.savefig(
                "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\plots--solstice_occultation20040706T181400\\"
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
            # convert to julian dates
            x = np.where(~np.isnan(timelist[:, i]))

            time = Time(timelist[:, i][x], format='gps')
            time.format = 'jd'
            formatter = DateFormatter("%d/%m/%y %H: %M: %S:")

            fig, ax = plt.subplots()

            plt.plot_date(time.plot_date, irradiancelist[:, i][x])
            ax.xaxis.set_major_formatter(formatter)
            ax.xaxis.set_tick_params(rotation=10, labelsize=5)

            plt.title('Wavelength ' + "{:.2f}".format(wavelengthlist[0, i]))
            plt.xlabel("Time (GPSseconds)")
            plt.ylabel("Irradiance")
            plt.savefig(
                "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\plots--solstice_occultation20040706T181400\\"
                "time_series_plots\\wavelength-" + "{:.2f}".format(wavelengthlist[0, i]) + ".png")
            plt.close()

        # read O2 cross section file and assign names to variables
        # data to be used for interpolation
        cross_section = 'C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\O2_cross_section\\o2_cross_section.ncdf'
        file_handle = nc.Dataset(cross_section, mode='r')
        sigma = np.squeeze(file_handle.variables['sigma'])
        wave = np.squeeze(file_handle.variables['wave'])

        # interpolate sigma for wavelengths 120.5-121.9
        f = interpolate.interp1d(wave, sigma)
        sigmanew = f(wavelengthlist[0, :])
        plt.plot(wavelengthlist[0, :], sigmanew, '-')

        # loop through list of irradiance and tangent height arrays to plot each element
        for i in np.arange(wavelengthlist[0, :].size):
            plt.plot(tanheightlist[:, i], irradiancelist[:, i])
            plt.title('Wavelength ' + "{:.2f}".format(wavelengthlist[0, i]))
            plt.xlabel("Tangent Height (km)")
            plt.ylabel("Irradiance")
            plt.savefig(
                "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\plots--solstice_occultation20040706T181400\\"
                "tanheight_plots\\wavelength-" + "{:.2f}".format(wavelengthlist[0, i]) + ".png")
            plt.close()

            x = np.where((tanheightlist[:, i] >= 200) & (250 >= tanheightlist[:, i]))
            if x[0].size > 0:
                ave_irradiance = np.mean(irradiancelist[:, i][x])

                extinction_ratio = irradiancelist[:, i] / ave_irradiance

                plt.plot(tanheightlist[:, i], extinction_ratio)
                plt.title("Extinction Ratio at Wavelength-" + "{:.2f}".format(wavelengthlist[0, i]))
                plt.xlabel("Tangent Height (km)")
                plt.ylabel("Extinction Ratio")
                plt.savefig(
                    "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\plots--solstice_occultation20040706T181400\\"
                    "extinction_ratio_plots\\wavelength " + "{:.2f}".format(wavelengthlist[0, i]) + ".png")
                plt.close()

                # column density equation, derived from extinction ratio equation
                column_density = -np.log(extinction_ratio) / sigmanew[i]

                # plot column density against tangent height in log scale
                plt.plot(tanheightlist[:, i], column_density)
                plt.title("Column Density at Wavelength-" + "{:.2f}".format(wavelengthlist[0, i]))
                plt.yscale("log")
                plt.xlabel("Tangent Height (km)")
                plt.ylabel("Column Density (n/cm^2)")
                plt.savefig(
                    "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\plots--solstice_occultation20040706T181400\\"
                    "column_density_logscale\\wavelength " + "{:.2f}".format(wavelengthlist[0, i]) + ".png")
                plt.close()

                # plot column density against tangent height not in log scale
                plt.plot(tanheightlist[:, i], column_density)
                plt.title("Column Density at Wavelength-" + "{:.2f}".format(wavelengthlist[0, i]))
                plt.xlabel("Tangent Height (km)")
                plt.ylabel("Column Density (n/cm^2)")
                plt.savefig(
                    "C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\plots--solstice_occultation20040706T181400\\"
                    "column_density\\wavelength " + "{:.2f}".format(wavelengthlist[0, i]) + ".png")
                plt.close()


        plt.rcParams["figure.figsize"] = [15,7]
        alltanheightlist = []
        column_densitylist = []
        for i in np.arange(wavelengthlist[0, :].size):
            x = np.where((tanheightlist[:, i] >= 200) & (250 >= tanheightlist[:, i]))
            if x[0].size > 0:
                ave_irradiance = np.mean(irradiancelist[:, i][x])
                extinction_ratio = irradiancelist[:, i] / ave_irradiance

                # column density equation, derived from extinction ratio equation
                column_density = -np.log(extinction_ratio) / sigmanew[i]
                column_densitylist.extend(column_density)
                alltanheightlist.extend(tanheightlist[:, i])
                plt.plot(tanheightlist[:, i], column_density, '.', label="{:.2f}".format(wavelengthlist[0, i]) + "nm")
                plt.yscale("log")
                plt.xlabel("Tangent Height")
                plt.ylabel("Column Density")
                plt.savefig("C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\scatter_plot\\plot-20040706T181400.png")

        alltanheightlist = np.array(alltanheightlist, dtype=np.double)
        x = alltanheightlist.argsort()
        alltanheightlist = alltanheightlist[x]
        column_densitylist = np.array(column_densitylist, dtype=np.double)
        column_densitylist = column_densitylist[x]

        x = np.where(~np.isnan(column_densitylist))
        alltanheightlist = alltanheightlist[x]
        column_densitylist = column_densitylist[x]


        plt.plot(alltanheightlist, smooth(column_densitylist, 50), 'g-')
        plt.legend(loc=0)
        plt.savefig("C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\scatter_plot\\plot-20040706T181400.png")


# call back my function
if __name__ == '__main__':
    file_name = 'C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_current_code\\venv\\SOLdatasets' \
                '\\solstice_occultation20040706T181400.000.nc'
    plot_irradianceSOL(file_name)
