from DEM_generator import DEM_generator
import netCDF4 as nc
"""
This file will give an example of how to run DEM generator with given files both a grid and a transect file
 to use either comment the other

"""
# gridded bathymetry data run through scaleC
filelist = ['http://134.164.129.55/thredds/dodsC/FRF/survey/gridded/FRF_20160726_1121_FRF_NAVD88_LARC_GPS_UTC_v20170320_grid_latlon.nc']
# transect bathymetry data run though scaleC
# filelist = ['/home/spike/repos/scalecInterp_python/FRF_20160726_1121_FRF_NAVD88_LARC_GPS_UTC_v20170320.nc']  # files with NEW data that are in background grid


ncfile = nc.Dataset(filelist[0])  # open netCDF file
x0, y0 = 1600, 1600  # north east corner of grid
x1, y1 = 0, -200     # south west corner of grid
dict = {'x0': x0,    #gp.FRFcoord(x0, y0)['Lon'],  # -75.47218285,
        'y0': y0,    #gp.FRFcoord(x0, y0)['Lat'],  #  36.17560399,
        'x1': x1,    #gp.FRFcoord(x1, y1)['Lon'],  # -75.75004989,
        'y1': y1,    #gp.FRFcoord(x1, y1)['Lat'],  #  36.19666112,
        'lambdaX': 10,  # grid spacing in x  -  Here is where CMS would hand array of variable grid spacing
        'lambdaY': 10,  # grid spacing in y
        'msmoothx': 100,  # smoothing length scale in x
        'msmoothy': 200,  # smoothing length scale in y
        'msmootht': 1,   # smoothing length scale in Time
        'filterName': 'hanning',
        'nmseitol': 0.75,
        'grid_coord_check': 'FRF',
        'grid_filename': '',  # should be none if creating background Grid!  becomes best guess grid
        'data_coord_check': 'FRF',
        'xFRF_s': ncfile['xFRF'][:],
        'yFRF_s': ncfile['yFRF'][:],
        'Z_s': ncfile['elevation'][:]
        }

out = DEM_generator(dict)

from matplotlib import pyplot as plt
plt.figure()
plt.subplot(221)
plt.title('Zi')
plt.pcolor(out['x_out'], out['y_out'], out['Zi'] )
plt.colorbar()
plt.subplot(222)
plt.title('MSEi')
plt.pcolor(out['x_out'], out['y_out'], out['MSEi'])
plt.colorbar()
plt.subplot(223)
plt.title('NMSEi')
plt.pcolor(out['x_out'], out['y_out'], out['NMSEi'])
plt.colorbar()
plt.subplot(224)
plt.title('MSRi')
plt.pcolor(out['x_out'], out['y_out'], out['MSRi'])
plt.colorbar()
plt.tight_layout()
