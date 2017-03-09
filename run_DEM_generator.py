# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 09:33:00 2014

@author: jwlong
"""

import arcpy
import numpy as np
import os

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
toolkitpath = 'D:\\CDI_DEM\\geoprocessing'
savepath = 'D:\\CDI_DEM\\geoprocessing'
datapath = 'D:\\CDI_DEM\\2010_dauphin_lidar_test'
datatype = 'las'
x0 = -88.36
x1 = -88
y0 = 30.19
y1 = 30.26
lambdaX = 40
lambdaY = 100
msmoothx = 100
msmoothy = 200
msmootht = 1
filtername = 'hanning'
nmseitol = 0.75
grid_coord_check = 'LL'
data_coord_check = 'LL'
grid_filename = ' '

# Call dataBuilder to construct data in necessary format for interpolation
from list_files import list_files
filelist = list_files(datapath, datatype)

from dataBuilder import dataBuilder
os.chdir(datapath)
x, z = dataBuilder(filelist, data_coord_check)
s = np.ones((np.size(x[:,1]),1))
lfile = np.shape(filelist)

# Call grid builder to make a grid based on x,y min and max values
os.chdir(toolkitpath)
from gridBuilder import gridBuilder
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

del x_grid
del y_grid
del t_grid

# subsample the data
from subsampleData import subsampleData
DXsmooth = np.array([msmoothx,msmoothy,msmootht])/4
DXsmooth[2] = 1
import time
t = time.time()
Xi, zprime, si = subsampleData(x,z,s,DXsmooth)
elapsed = time.time() - t
print 'subsampling time is %d seconds' % elapsed

# Send it all into scalecinterpolation
from scalecInterpolation import scalecInterpTilePerturbations
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
from scipy.io import netcdf 
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

