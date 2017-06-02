# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 09:33:00 2014

@author: jwlong
"""
import numpy as np
import os
import time
import list_files

#WorkingDir = "%scratchworkspace%"



# Define inputs (will eventually come from ScienceBase user interface)
toolkitpath = 'D:\\CDI_DEM\\geoprocessing'
savepath = 'D:\\CDI_DEM\\geoprocessing'
datapath = 'D:\\CDI_DEM\\2010_dauphin_lidar_test'
datatype = 'las'

# these are required -- don't run unless a user provides these
x0 = -88.36
x1 = -88
y0 = 30.19
y1 = 30.26
dx = 40
dy = 100

# make some defaults based on dx and dy
msmoothx = 2*dx
msmoothy = 2*dy
msmootht = 1
Lysmooth = 0
splinebc = 1
splinedxm = 2
splinelc = 1 
splinebcpert = 1
targetvar = 1

# flag for making alongshore boundaries uniform
if (Lysmooth == 0):
    isflow = 0
elif    (Lysmooth > 0):
    isflow = 1
elif (Lysmooth < 0):
    Lysmooth = abs(Lysmooth)
    isflow = 2
    print 'setting spline bc to 1,10 for periodic bc'
    splinebcpert = np.array([1,10])

# set everything for invoking a spline
if not splinebc:
    splinebc = -1
if (splinebc >= 0):
    # check other splines
    if not splinelc:
        splinelc = 1
    if not splinedxm:
        splinedxm = 2
    if not splinebcpert:
        splinebcpert = 1
    
# see if user suppied a previous grid; got a grid from ScienceBase; or wants a NOAA DEM?
priorgrid_filename = 'northern_gulf_coast_navd_88.grd'
print 'no prior DEM provided; loading NOAA DEM by default'        
isBestGuess = 1

filtername = 'hanning'
nmseitol = 0.75

# would be better to have the program check this automatically
grid_coord_check = 'LL'
data_coord_check = 'LL'

# allow users to provide a grid of their own but someow specify what it
# needsto look like
grid_filename = ' '

#################### end inputs ###########################################

# Call dataBuilder to construct data in necessary format for interpolation
os.chdir(toolkitpath)
filelist = list_files.list_files(datapath, datatype)

from dataBuilder import dataBuilder
os.chdir(datapath)
t = time.time()
x, z = dataBuilder(filelist, data_coord_check)
s = np.ones((np.size(x[:,1]),1))
lfile = np.shape(filelist)
elapsed = time.time() - t
print 'data building time is %d seconds' % elapsed

# Call grid builder to make a grid based on x,y min and max values
os.chdir(toolkitpath)
from gridBuilder import gridBuilder
x_grid, y_grid = gridBuilder(x0, x1, y0, y1, dx, dy, grid_coord_check, grid_filename)
t_grid = np.zeros((x_grid.shape))
xi = np.array([x_grid.flatten(1), y_grid.flatten(1), t_grid.flatten(1)]).T
xsm = msmoothx*np.ones(x_grid.shape)
ysm = msmoothy*np.ones(y_grid.shape)
tsm = msmootht*np.ones(t_grid.shape)
lx = np.array([xsm.flatten(1), ysm.flatten(1), tsm.flatten(1)]).T

N, M = np.shape(x_grid)
x_out = x_grid[1,:].copy()
y_out = y_grid[:,1].copy()

# load the priod DEM
# need to add switch so the user doesn't have to load NOAA DEM
os.chdir(toolkitpath)
from dataBuilder import loadNOAAdem
os.chdir(datapath)
Xprior, Yprior, Zprior = loadNOAAdem(priorgrid_filename, x0, x1, y0, y1)

# get smoothing scales based on grid
Lxg = 2 * dx
Lyg = 2 * dy

# must not be finer than data smoothing
Lxg1 = Lxg*np.ones(xsm.shape)
idd = np.where(Lxg1 < msmoothx)
if (np.size(idd, axis=0) > 0):
    Lxg1[idd] = msmoothx
    Lxg = Lxg1 #np.reshape(Lxg1, (np.shape(Lxg)))
    
Lyg1 = Lyg*np.ones(ysm.shape)
idd = np.where(Lyg1 < msmoothy)
if (np.size(idd, axis=0) > 0):
    Lyg1[idd] = msmoothy
    Lyg = Lyg1 #np.reshape(Lyg1, (np.shape(Lyg)))
del Lxg1, Lyg1
    
# smoothing scales should be smooth
#LXG = 0.5*LXG + 0.25*[LXG(:,2:end), LXG(:,end-1)] + 0.25*[LXG(:,2), LXG(:,1:end-1)];
tmp1 = np.reshape(Lxg[:,-1], (np.size(Lxg[:,-1],axis=0), 1) )
tmp2 = np.reshape(Lxg[:,1], (np.size(Lxg[:,1],axis=0), 1) )
Lxg = 0.5 * Lxg + 0.25 * np.concatenate((Lxg[:,0:-1], tmp1), axis=1) + 0.25 * np.concatenate((tmp2, Lxg[:,0:-1]), axis=1)
   
#LYG = 0.5*LYG + 0.25*[LYG(2:end,:); LYG(end-1,:)] + 0.25*[LYG(2,:); LYG(1:end-1,:)];
tmp3 = np.reshape(Lyg[-1,:], (1, np.size(Lyg[-1,:], axis=0) ))
tmp4 = np.reshape(Lyg[1,:], (1, np.size(Lyg[1,:], axis=0) ))
Lyg = 0.5 * Lyg + 0.25 * np.concatenate((Lyg[0:-1,:], tmp3), axis=0) + 0.25 * np.concatenate((tmp4, Lyg[0:-1,:]), axis=0)
del tmp1, tmp2, tmp3, tmp4
    
# need a regular grid for tile/spline    
#mindxg = min(np.ravel(Dxg))
#mindyg = min(np.ravel(Dyg))
    
# reset splinedxm
tmp5 = splinedxm
# if for spline
if (splinebc > 0):    
    splinedxm = []
    # splines (if used) can depend on smoothing scales
    splinedxm.append(max(np.array([tmp5, np.fix(0.25 * min(Lxg.flatten(1)) / dx)])))
    splinedxm.append(max(np.array([splinedxm[0], np.fix(0.25 * min(Lyg.flatten(1)) / dy)])))
    # taper for smoothing at lateral boundaries
    # tmp6 = np.reshape(Dyg[0,:], (np.size(Dyg[0,:], axis=0), 1)).T
    #Nysmooth = np.fix(Lysmooth / max(np.concatenate((tmp6, Dyg[:,:]),axis=0).flatten(1))) # returns an array of size with the value of Nysmooth inside
    Nysmooth = np.fix(Lysmooth / dy)
    # Nysmooth = np.reshape(Nysmooth, (1)) # Nysmooth comes out to be of size zero and therefore can't be parsed... but can be made to be of size 1, which can be parsed
    # Nysmooth = Nysmooth[0] # Once the size has been changed, the value can be extracted
    if (Nysmooth < 1):
        Nysmooth = 1
del tmp5 #, tmp6
 
# Ignore this...we require this as input for this phase
# check for implied simple grid
#if (all(gridStruct['Xg'][0,:] == gridStruct['Xg'][-1,:]) and all(gridStruct['Yg'][:,0] == gridStruct['Yg'][:,0])): # the second term of this conditional might have a typo... a[:,0] == a[:,0] will always be true... 
#if (all(gridStruct['Xg'][0,:] == gridStruct['Xg'][-1,:]) and all(gridStruct['Yg'][:,0] == gridStruct['Yg'][:,-1])): # This conditional term makes more sense to me... 
#    # it is simple grid        
#    Xreg = gridStruct['Xg']
#    Yreg = gridStruct['Yg']
#else:
#    # input grid is not simple rectangular
#    print 'input grid is not simple rectangular grid, generating regular grid \n'
#    xreg = np.arange(gridStruct['Xg'][0,0], (gridStruct['Xg'][-1,-1] + mindxg), mindxg) # tack on one extra to be inclusive
#    yreg = np.arange(gridStruct['Yg'][0,0], (gridStruct['Yg'][-1,-1] + mindyg), mindyg)
#    Xreg, Yreg = np.meshgrid(xreg,yreg)

Xreg = x_grid
Yreg = y_grid
del x_grid, y_grid

Ny, Nx = np.shape(Xreg)
    
# and copy the variable smoothing scales to the local grid
from scipy.interpolate import interp2d # May be deleted... maybe
from scipy.interpolate import CloughTocher2DInterpolator
interpObj1 = CloughTocher2DInterpolator(np.array([Xreg.flatten(1), Yreg.flatten(1)]).T, Lxg.flatten(1), tol = 1e-6)
msmoothx = interpObj1(Xreg, Yreg)
interpObj2 = CloughTocher2DInterpolator(np.array([Xreg.flatten(1), Yreg.flatten(1)]).T, Lyg.flatten(1), tol = 1e-6)
msmoothy = interpObj2(Xreg, Yreg)
    
idd = np.nonzero(np.isnan(msmoothx))
if (np.size(idd) > 0):
    msmoothx[idd] = max(np.ravel(msmoothx))
idd = np.nonzero(np.isnan(msmoothy))
if (np.size(idd) > 0):
    msmoothy[idd] = max(np.ravel(msmoothy))
 
# need time points for output    
#Treg = np.array([])
#msmootht = np.array([])
# catch the input variant
#if ('Tg' in gridStruct):
#    gridStruct['tg'] = gridStruct['Tg']
#    del gridStruct['Tg']
#if ('tg' in gridStruct and gridStruct['tg'] != None):
#    Treg = np.tile(gridStruct['tg'], np.shape(Xreg))
#    msmootht = np.tile(msmootht, np.shape(Xreg))
#else:
    # Declare the following variables to be VERY small as opposed to passing them as zeros or nans... which will cause problems
#Treg = 0.0001 * np.ones(np.shape(Xreg), float) 
#msmootht = 0.0001 * np.ones(np.shape(msmoothx), float)
        
Treg = np.tile(0, np.shape(Xreg))   # need this to be te time of te grid
msmootht = np.tile(msmootht, np.shape(Xreg))
    
# subsample the data
os.chdir(toolkitpath)
from subsampleData import subsampleData
DXsmooth = np.array([np.min(msmoothx),np.min(msmoothy),np.min(msmootht)])/4
DXsmooth[2] = 1
t = time.time()
Xi, zprime, si = subsampleData(x,z,s,DXsmooth)
elapsed = time.time() - t
print 'subsampling time is %d seconds' % elapsed
#del x
#del z

tmp = zprime[~np.isnan(zprime)]
if (np.size(tmp) == 0):
    print 'ERROR: no input data found'
else:
    Xi = Xi[~np.isnan(zprime).any(1),:]
    s = s[~np.isnan(zprime)]
    zprime = zprime[~np.isnan(zprime)]

zprime = zprime[:,np.newaxis]
s = s[:,np.newaxis]
    
# remove background best guess if exists    
if (isBestGuess):
    print 'best guess found, removing it before interpolating data... \n'
    # get the best guess at the input points
    interpObj = CloughTocher2DInterpolator(np.array([Xprior.flatten(1), Yprior.flatten(1)]).T, Zprior.flatten(1), tol=1e-6)     
    zbgdata = interpObj(Xi[:,0], Xi[:,1])
    zbgdata = zbgdata[:,np.newaxis]
    
#    zbgdata = np.reshape(zbgdata, (len(zbgdata),1))
    # get the detrended data        
    zprime = zprime - zbgdata
#    del zbgdata
    # get the best guess at the output points
    Zbg = interpObj(Xreg, Yreg)
else: #### CONDITIONS TO REACH THIS else HAVE NOT BEEN MET WITH CURRENT TESTING #####
    # no background, create dummy bg on output points
    Zbg = 0
        
    # if we are doing flow model, may need an alongshore uniform constraint
    if (isflow > 0):
        if ((2 * Nysmooth) > (Ny / 1)):
            print 'ERROR: grid not big enough to smooth ends'
        else:
            print 'WARNING: generating a smooth background grid'
        # get enough data for at least a trend surface on top and bottom boundaries
        # use only data near end: id = find(data(:,2)>(Yreg(1)-10*msmooth(2)) & data(:,2)<(Yreg(Nysmooth,1)+10*msmooth(2)));
        bcsmooth = np.array([max(msmoothx.flatten(1)), min(msmoothy.flatten(1)) * Lysmooth / dy])
        idd = np.arange(1, (np.size(Xi, axis=0) + 1))
        from scalecInteprolation import scalecInterp
        Zitrend1, Eitrend1, NEitrend1, MSRitrend1 = scalecInterp(Xi[idd, 0:2], zprime[idd], s[idd], np.array([Xreg[0,:].T, Yreg[0,:].T]).T, bcsmooth, filtername, 0.5)#, nargout=4) # MAY NEED REVISION
        wb = (1 - Eitrend1) / (targetvar + Eitrend1)
        
        from bsplineFunctions import bspline_pertgrid
        Zbstrend1 = bspline_pertgrid(Zitrend1, wb, splinebc, splinelc, splinedxm[0])
        Zitrend2, Eitrend2, NEitrend2, MSRitrend2 = scalecInterp(Xi[idd, 0:2], zprime[idd], s[idd], np.array([Xreg[0,:].T, Yreg[-2,:].T]), bcsmooth, filtername, 0.5)#, nargout=4) # MAY NEED REVISION
        wb = 1 - Eitrend2 / (targetvar + Eitrend2)
        Zbstrend2 = bspline_pertgrid(Zitrend2, wb, splinebc, splinelc, splinedxm[0])
        
        # check for periodic conditions            
        if (isflow == 2):
            print 'making identical lateral boundaries for periodic condition'
            Zbstrend1 = 0.5 * (Zbstrend1 + Zbstrend2)
            Zbstrend2 = Zbstrend1
            
        # make trend
        Zbg = np.zeros(np.shape(Xreg))
        Zbg[0:Nysmooth,:] = np.tile(Zbstrend1, Nysmooth, 1)
        Zbg[(-1 - Nysmooth + 1):-1,:] = np.tile(Zbstrend2,(Nysmooth,1))
        # fill in the middle
        Zbg[((Nysmooth + 1) - 1):(-1 -Nysmooth),:] = interp2d(np.tile(Xreg[0,:],(2,1))), np.tile(np.array([Yreg[(Nysmooth - 1),0], Yreg[(Ny - Nysmooth + 1 - 1),0]]),(1,Nx)), np.array([Zbstrend1,Zbstrend2]), Xreg[((Nysmooth + 1) - 1):(-1 - Nysmooth),:], Yreg[((Nysmooth + 1)-1):(-1-Nysmooth),:]
        # remove this trend
        isBestGuess = 1 # update flag that we have bg
        interpObj = interp2d(Xreg, Yreg, Zbg)
        zprime = zprime - interpObj(Xi[:,0], Xi[:,1])
        idd = np.nonzero(not np.isnan(zprime))
        zprime = zprime[idd]
        Xi = Xi[idd,:]
        s = s[idd]
        
# Send it all into scalecinterpolation
from scalecInterpolation import scalecInterpTilePerturbations
t = time.time()
print 'Interpolating'
zi, msei, nmsei, msri = scalecInterpTilePerturbations(Xi, zprime, si, xi, lx, filtername, nmseitol)
elapsed = time.time() - t
print 'Interpolating time is %d seconds' % elapsed

# check output
idd = np.nonzero(np.isnan(zi))
if (np.size(idd) > 0):
    print ' WARNING: found', np.size(idd, axis=1) ,' nans in Zi' 
    zi[idd] = 0
    nmsei[idd] = 1
Zi = np.reshape(zi[:,np.newaxis], np.shape(Xreg.T)).T
NEi = np.reshape(nmsei, np.shape(Xreg.T)).T
Ei = np.reshape(msei, np.shape(Xreg.T)).T
MSRi = np.reshape(msri, np.shape(Xreg.T)).T
    
# post process with splines?
if (splinebc > 0):
    # may need to fix flow constraints and spline
    # set up the default spline weight based on INTERPOLATION errors
    wb = np.array([1 - Ei / (targetvar + Ei)]) # based on error requirements
    
    # if this is flow, force back to the boundary conditions with spline weights and data mask        
    if (isflow > 0): ##### ALL LINES BEYOND THIS POINT HAVE NOT BEEN TESTED #####
        print ('splining flow \n')
        from supportingMethods import makeWBFlow
        wbflow = makeWBFlow(Yreg, Nysmooth, dy)
        # spline force gradient to zero
        wb = wb * wbflow
        wb = wb * np.array([1 - NEi]) # force perturbations to zero based on sample errors
        # weights damp both the data directly and the spline output
        # splinebcpert used
        print  'using splinebcpert = \n', splinebcpert
        from bsplineFunctions import bspline_pertgrid
        Zbs = bspline_pertgrid(Zi*wb, wb, splinebcpert, splinelc, splinedxm)
    else:
        print 'splining generic \n'
        # weights are simpler and we don't let the attack the data
        # splinebc used (NOT splinebcpert)
        from bsplineFunctions import bspline_pertgrid
        Zbs = bspline_pertgrid(Zi, wb, splinebc, splinelc, splinedxm)
        
    # not supposed to be nans, but it is possible        
    idd = np.nonzero(np.isnan(Zbs))
    if (np.size(idd) > 0):
        print ' WARNING: found ', np.size(idd, axis=0) ,' nans in Zbs '
        Zbs[idd]=0
        
    # put the trend back in now
    Zbs = Zbs + Zbg
    
# Put the trend back into scalecInterp product    
Zi += Zbg
    
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

