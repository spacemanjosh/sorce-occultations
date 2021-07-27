import netCDF4 as nc
import matplotlib.pyplot as plt



if __name__ == '__main__':
    file_name = 'C:\\Users\\Computer\\PycharmProjects\\SOLSTICE_fun\\venv\\solsticedata\\solstice_occultation20040706T181400.000.nc'
    with nc.Dataset(file_name, 'r', format='NETCDF4') as file_grp:
        variables = file_grp.variables
        #print(variables)

        # locate variables
        irradiance = file_grp.variables['irradiance']
        print(irradiance)
        time = file_grp.variables['microsecondssincegpsepoch']
        print (time)

        plt.plot(time, irradiance)
        plt.show()






