# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 09:33:00 2014

@author: jwlong
"""

import numpy as np
import os, time
from scipy.io import netcdf
from dataBuilder import dataBuilder, gridBuilder
# from gridBuilder import gridBuilder
from subsampleData import subsampleData
from scalecInterpolation import scalecInterpTilePerturbations
import datetime as DT



# Define inputs (will eventually come from ScienceBase user interface)
toolkitpath = '' # D:\\CDI_DEM\\geoprocessing'        # Path to the interpolation toolkit codes - should be local to repo
savepath = ''  # ''D:\\CDI_DEM\\geoprocessing'        # Path to the final output directory for saving the DEM
datapath = '' #FRF_20170227_1131_FRF_NAVD88_LARC_GPS_UTC_v20170320.nc'     # Path to the raw data files
# datatype = 'mat'                                      # Type of data to be analyzed (file extension; e.g. 'las' for lidar tile files)
                                                      #     ['las', 'laz', 'nc', 'txt', 'mat']
x0 = -75.47218285                                     # Minimum x-value of the output grid (origin)
x1 = -75.75004989                                              # Maximum x-value of the output grid
y0 = 36.17560399                                      # Minimum y-value of the output grid (origin)
y1 = 36.19666112                                      # Maximum y-value of the output grid
lambdaX = 10                                        # Grid spacing in the x-direction
lambdaY = 10                                         # Grid spacing in the y-direction
msmoothx = 100                                        # Smoothing length scale in the x-direction
msmoothy = 200                                        # Smoothing length scale in the y-direction
msmootht = 1                                          # Smoothing length scale in time
filtername = 'hanning'                                # Name of the filter type to smooth the data
                                                      #      ['hanning', 'linloess', 'quadloess', 'boxcar', si']
nmseitol = 0.75                                       # Normalized error tolerance the user will tolerate in the final grid
                                                      #      (0 - (no error) to 1 (no removal of bad points))
grid_coord_check = 'LL'                               # ['LL' or 'UTM'] - Designates if the grid supplied by the user (if one exists)
                                                      #      is in UTM or lat-lon coordinates
grid_filename = ' '                                   # Name of the grid file (if supplied)
outFname = 'TestOutput.nc'
###########################
data_coord_check = 'LL' #, 'NCSP']                               # ['LL' or 'UTM'] - Designates if the data supplied by the user
                                                      #      is in UTM or lat-lon coordinates


def DEM_generator(dict):
    """
    :param dict:
    x0                          # Minimum x-value of the output grid (origin)
    y0                          # Minimum y-value of the output grid
    x1                          # Maximum x-value of the output grid
    y1                          # Maximum y-value of the output grid
    lambdaY                     # grid spacing in the y-direction
    lambdaX                     # Grid spacing in the x-direction
    msmoothx                    # Smoothing length scale in the x-direction
    msmoothy                    # Smoothing length scale in the y-direction
    msmootht                    # Smoothing length scale in time
    filtername                  # Name of the filter type to smooth the data
                                #      ['hanning', 'linloess', 'quadloess', 'boxcar', si']
    nmseitol                    # Normalized error tolerance the user will tolerate in the final grid
                                #      (0 - (no error) to 1 (no removal of bad points))

    :return: dict with keys:
        zi, the depth estimate
        msei, the mean square interpolation error estimate (units of z)
        nmsei, the normalized mean square error
        msri, the mean square residuals
    """
    x0 = dict['x0']             # Minimum x-value of the output grid (origin)
    y0 = dict['y0']             # Minimum y-value of the output grid
    x1 = dict['x1']             # Maximum x-value of the output grid
    y1 = dict['y1']             # Maximum y-value of the output grid
    lambdaY = dict['lambdaY']   # grid spacing in the y-direction
    lambdaX = dict['lambdaX']   # Grid spacing in the x-direction
    msmoothx = dict['msmoothx'] # Smoothing length scale in the x-direction
    msmoothy = dict['msmoothy'] # Smoothing length scale in the y-direction
    msmootht = dict['msmootht'] # Smoothing length scale in time
    filtername = dict['filterName']  # Name of the filter type to smooth the data
                                     #      ['hanning', 'linloess', 'quadloess', 'boxcar', si']
    nmseitol = dict['nmseitol']  # Normalized error tolerance the user will tolerate in the final grid
                                                      #      (0 - (no error) to 1 (no removal of bad points))

    #### data checks ###########3
    filters = ['hanning', 'linloess', 'quadloess', 'boxcar', 'si']
    assert filtername in filters, 'Check filter name, not appropriate for current DEM generator function'
    ####################################################################
    # ############################### Load Data ########################
    ####################################################################
    t = DT.datetime.now()
    x, z = dataBuilder(filelist, data_coord_check)  # function loads x, y, z data and concatenates in long array
    s = np.ones((np.size(x[:,1]),1))
    print 'loading time is %s seconds' % (DT.datetime.now() - t)

    ####################################################################
    # Call grid builder to make a grid based on x,y min and max values #
    ####################################################################

    x_grid, y_grid = gridBuilder(x0, x1, y0, y1, lambdaX, lambdaY, grid_coord_check, grid_filename)
    t_grid = np.zeros_like((x_grid))  # dummy values, for time?
    xi = np.array([x_grid.flatten(), y_grid.flatten(), t_grid.flatten()]).T  # grid locations, flatten make row-major style
    xsm = msmoothx*np.ones_like(x_grid)
    ysm = msmoothy*np.ones_like(y_grid)
    tsm = msmootht*np.ones_like(t_grid)
    lx = np.array([xsm.flatten(), ysm.flatten(), tsm.flatten()]).T  # smoothing array , flatten takes row-major style

    N, M = np.shape(x_grid)
    x_out = x_grid[0,:].copy()  # utm output for x and y
    y_out = y_grid[:,0].copy()

    del x_grid, y_grid, t_grid
    #####################################################################
    # subsample the data   ##############################################
    #####################################################################
    DXsmooth = np.array([msmoothx, msmoothy, msmootht])/4
    DXsmooth[2] = 1
    t = DT.datetime.now()
    Xi, zprime, si = subsampleData(x, z, s, DXsmooth)
    # What's returned here
    print 'subsampling time is %s seconds' % (DT.datetime.now() - t)

    #####################################################################
    # Send it all into scalecinterpolation  -  Here is where the interpolation takes place
    #####################################################################
    t = DT.datetime.now()
    print 'Interpolating'
    zi, msei, nmsei, msri = scalecInterpTilePerturbations(Xi, zprime, si, xi, lx, filtername, nmseitol)
    print 'Interpolating time is %s seconds' % (DT.datetime.now() - t)

    # save the ouput
    # reshape
    zi = np.reshape(zi, (M,N))  # what are these errors ???
    msei = np.reshape(msei, (M,N))
    nmsei = np.reshape(nmsei, (M,N))
    msri = np.reshape(msri, (M,N))

       # zi, the estimate
       # msei, the mean square interpolation error estimate (units of z)
       # nmsei, the normalized mean square error
       # msri, the mean square residuals
    out = {'Zi': zi,
           'MSEi': msei,
           'NMSEi': nmsei,
           'MSRi': msri}
    return out







# ########## scratch
from sblib import geoprocess as gp

# Call dataBuilder to construct data in necessary format for interpolation
# filelist = list_files(datapath, datatype)  # creates a list of files to be interpolated
filelist = ['FRF_20170518_1134_FRF_NAVD88_LARC_GPS_UTC_v20170525.nc',]  # files with NEW data that are in background grid

x0, y0 = 1600, 1600  # north east corner of grid
x1, y1 = 0, -200     # south west corner of grid
dict = {'x0': gp.FRFcoord(x0, y0)['Lon'],  # -75.47218285,
        'y0': gp.FRFcoord(x0, y0)['Lat'],  #  36.17560399,
        'x1': gp.FRFcoord(x1, y1)['Lon'],  # -75.75004989,
        'y1': gp.FRFcoord(x1, y1)['Lat'],  #  36.19666112,
        'lambdaX': 10,  # grid spacing in x  -  Here is where CMS would hand array of variable grid spacing
        'lambdaY': 10,  # grid spacing in y
        'msmoothx': 100,  # smoothing length scale
        'msmoothy': 200,  # smoothing length scale
        'filtername': 'hanning',
        'nmseitol': 0.75,
        'grid_coord_check': 'LL',
        'grid_filename': '',  # should be none if creating background Grid!  becomes best guess grid
        'data_coord_check': 'LL',
        'outFname': 'TestOutput.nc'}



######################################################################
# open a new netCDF file for writing.
######################################################################
# save file to output
ncfile = netcdf.netcdf_file(outFname, 'w')
ncfile.history = 'Topo/Bathy DEM created using XXX'

# create the lat and lon dimensions.
ncfile.createDimension('y_utm', N)
ncfile.createDimension('x_utm', M)

# Define the coordinate variables. They will hold the coordinate
# information, that is, the latitudes and longitudes.
y_utm = ncfile.createVariable('y_utm','float32',('y_utm',))
x_utm = ncfile.createVariable('x_utm','float32',('x_utm',))

# Assign units attributes to coordinate var data. This attaches a
# text attribute to each of the coordinate variables, containing the
# units.
y_utm.units = 'meters'
x_utm.units = 'meters'

# write data to coordinate vars.
y_utm[:] = y_out
x_utm[:] = x_out

# create the pressure and temperature variables 
z = ncfile.createVariable('z','float32',('x_utm','y_utm'))
mse = ncfile.createVariable('mean_square_error','float32',('x_utm','y_utm'))
nmse = ncfile.createVariable('normalized_mean_square_error','float32',('x_utm','y_utm'))
msr = ncfile.createVariable('mean_square_residuals','float32',('x_utm','y_utm'))

# set the units attribute.
z.units =  'meters'
mse.units =  'meters'
nmse.units =  'meters'
msr.units =  'meters'

z.long_name =  'gridded elevation (datum NAVD88)'
mse.long_name =  'mean square error of the gridded surface'
nmse.long_name =  'normalized mean square error of the gridded surface'
msr.long_name =  'mean square residual of the gridded surface'

# write data to variables.
z[:] = zi
mse[:] = msei
nmse[:] = nmsei
msr[:] = msri

# having trouble getting these inputs to write correctly
# also write the inputs used
#ncfile.createDimension('smoothing_params',3)
#ncfile.createVariable('smoothing_params','i',('smoothing_params',))
#smoothing_params[:] = (msmoothx, msmoothy, msmootht)
#
#ncfile.createDimension('filtername',1)
#ncfile.createVariable('filtername','S1',('filtername',))
#filtername = filtername
#
#ncfile.createDimension('nmseitol',1)
#ncfile.createVariable('nmseitol','f4',('nmseitol',))
#nmseitol[:] = nmseitol
#nmseitol.long_name =  'normalized mean square error tolerance'

# close the file.
ncfile.close()
print '*** SUCCESS writing DEM to netcdf file!'  

