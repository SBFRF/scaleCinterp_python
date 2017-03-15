# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 09:33:00 2014

@author: jwlong
"""

import numpy as np
import os, time
from scipy.io import netcdf
from dataBuilder import dataBuilder
from list_files import list_files
from gridBuilder import gridBuilder
from subsampleData import subsampleData
from scalecInterpolation import scalecInterpTilePerturbations

WorkingDir = "%scratchworkspace%"

# inputs from ArcGIS WPS
#x0 = arcpy.GetParameterAsText(0)
#x1 = arcpy.GetParameterAsText(1)
#y0 = arcpy.GetParameterAsText(2)
#y1 = arcpy.GetParameterAsText(3)
#lambdaX = arcpy.GetParameterAsText(4)
#lambdaY = arcpy.GetParameterAsText(5)
#msmoothx = arcpy.GetParameterAsText(6)
#msmoothy = arcpy.GetParameterAsText(7)
#msmootht = arcpy.GetParameterAsText(8)
#nmseitol = arcpy.GetParameterAsText(9)


# Define inputs (will eventually come from ScienceBase user interface)
toolkitpath = 'D:\\CDI_DEM\\geoprocessing'            # Path to the interpolation toolkit codes
savepath = 'D:\\CDI_DEM\\geoprocessing'               # Path to the final output directory for saving the DEM
datapath = 'D:\\CDI_DEM\\2010_dauphin_lidar_test'     # Path to the raw data files
datatype = 'las'                                      # Type of data to be analyzed (file extension; e.g. 'las' for lidar tile files)
                                                      #     ['las', 'laz', 'nc', 'txt', 'mat']
x0 = -88.36                                           # Minimum x-value of the grid
x1 = -88                                              # Maximum x-value of the grid
y0 = 30.19                                            # Minimum y-value of the grid
y1 = 30.26                                            # Maximum y-value of the grid
lambdaX = 40                                          # Grid spacing in the x-direction
lambdaY = 100                                         # Grid spacing in the y-direction
msmoothx = 100                                        # Smoothing length scale in the x-direction
msmoothy = 200                                        # Smoothing length scale in the y-direction
msmootht = 1                                          # Smoothing length scale in time
filtername = 'hanning'                                # Name of the filter type to smooth the data
                                                      #      ['hanning', 'linloess', 'quadloess', 'boxcar', si']
nmseitol = 0.75                                       # Normalized error tolerance the user will tolerate in the final grid
                                                      #      (0 - (no error) to 1 (no removal of bad points))
grid_coord_check = 'LL'                               # ['LL' or 'UTM'] - Designates if the grid supplied by the user (if one exists)
                                                      #      is in UTM or lat-lon coordinates
data_coord_check = 'LL'                               # Name of the grid file (if supplied)
grid_filename = ' '                                   # ['LL' or 'UTM'] - Designates if the data supplied by the user
                                                      #      is in UTM or lat-lon coordinates

# Call dataBuilder to construct data in necessary format for interpolation
filelist = list_files(datapath, datatype)

os.chdir(datapath)
x, z = dataBuilder(filelist, data_coord_check)
s = np.ones((np.size(x[:,1]),1))
lfile = np.shape(filelist)

# Call grid builder to make a grid based on x,y min and max values
os.chdir(toolkitpath)
x_grid, y_grid = gridBuilder(x0, x1, y0, y1, lambdaX, lambdaY, grid_coord_check, grid_filename)
t_grid = np.zeros((x_grid.shape))
xi = np.array([x_grid.flatten(1), y_grid.flatten(1), t_grid.flatten(1)]).T
xsm = msmoothx*np.ones(x_grid.shape)
ysm = msmoothy*np.ones(y_grid.shape)
tsm = msmootht*np.ones(t_grid.shape)
lx = np.array([xsm.flatten(1), ysm.flatten(1), tsm.flatten(1)]).T

N, M = np.shape(x_grid)
x_out = x_grid[1,:].copy()
y_out = y_grid[:,1].copy()

del x_grid, y_grid, t_grid

# subsample the data
DXsmooth = np.array([msmoothx,msmoothy,msmootht])/4
DXsmooth[2] = 1
t = time.time()
Xi, zprime, si = subsampleData(x,z,s,DXsmooth)
elapsed = time.time() - t
print 'subsampling time is %d seconds' % elapsed

# Send it all into scalecinterpolation  -  Here is where the interpolation takes place
t = time.time()
print 'Interpolating'
zi, msei, nmsei, msri = scalecInterpTilePerturbations(Xi, zprime, si, xi, lx, filtername, nmseitol)
elapsed = time.time() - t
print 'Interpolating time is %d seconds' % elapsed

# save the ouput
os.chdir(savepath)
# reshape
zi = np.reshape(zi, (M,N))
msei = np.reshape(msei, (M,N))
nmsei = np.reshape(nmsei, (M,N))
msri = np.reshape(msri, (M,N))

# open a new netCDF file for writing.
ncfile = netcdf.netcdf_file('DEM_output.nc', 'w')
ncfile.history = 'Topo/Bathy DEM created using XXX'

# create the lat and lon dimensions.
ncfile.createDimension('y_utm',N)
ncfile.createDimension('x_utm',M)

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

