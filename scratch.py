from DEM_generator import DEM_generator
import netCDF4 as nc
import numpy as np
import os

# Call dataBuilder to construct data in necessary format for interpolation
# filelist = list_files(datapath, datatype)  # creates a list of files to be interpolated
filelist = ['http://134.164.129.55/thredds/dodsC/FRF/survey/gridded/FRF_20160726_1121_FRF_NAVD88_LARC_GPS_UTC_v20170320_grid_latlon.nc'] #''/home/spike/repos/scalecInterp_python/FRF_20160726_1121_FRF_NAVD88_LARC_GPS_UTC_v20170320.nc',]  # files with NEW data that are in background grid
bathy = nc.Dataset(filelist[0])

xV = bathy['xFRF'][:]
yV = bathy['yFRF'][:]
dataZ = bathy['elevation'][:]
dataZ = dataZ[0, :]
dataX, dataY = np.meshgrid(xV, yV)

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
        'xFRF_s': dataX,
        'yFRF_s': dataY,
        'Z_s': dataZ,
        }

out = DEM_generator(dict)

from matplotlib import pyplot as plt
fig_loc = 'C:\Users\dyoung8\Desktop\David Stuff\Projects\CSHORE\Bathy Interpolation\Test Figures'
fig_name = 'test.png'
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
plt.savefig(os.path.join(fig_loc, fig_name))
plt.close()
