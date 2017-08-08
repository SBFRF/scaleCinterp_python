# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 09:33:00 2014

@author: jwlong, edited Spicer Bak 8/4/17
"""
import numpy as np
from dataBuilder import dataBuilder, gridBuilder
from subsampleData import subsampleData
from scalecInterpolation import scalecInterpTilePerturbations
import datetime as DT



def DEM_generator(dict):
    """
    This is a function that takes in data with the below keys, and will generate a cartesian grid if there is not
    one already built.  if there is a grid given it will grid to the same locations.  If there is not it will generate
    a new grid with the x0,y0 and x1, y1, lambdaX and lambdaY keys.  This function will import data, subsample data,
    and grid data using scaleCinterp from plant 2002.  This will not incorporate to a background grid.

    :param dict:
    x0                          # Minimum x-value of the output grid (origin)
    y0                          # Minimum y-value of the output grid
    x1                          # Maximum x-value of the output grid
    y1                          # Maximum y-value of the output grid
    grid_filename               # full filepath of the existing grid (which this will build upon), if it exists
    lambdaY                     # grid spacing in the y-direction
    lambdaX                     # Grid spacing in the x-direction
    msmoothx                    # Smoothing length scale in the x-direction
    msmoothy                    # Smoothing length scale in the y-direction
    msmootht                    # Smoothing length scale in time
    filtername                  # Name of the filter type to smooth the data
                                #      ['hanning', 'linloess', 'quadloess', 'boxcar', si']
    nmseitol                    # Normalized error tolerance the user will tolerate in the final grid
                                #      (0 - (no error) to 1 (no removal of bad points))
    xFRF_s                      # survey xFRF coordinates
    yFRF_s                      # survey yFRF coordinates
    Z_s                         # survey bottom elevations

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
    try:
        grid_filename = dict['grid_filename']
    except:
        pass
    lambdaY = dict['lambdaY']        # grid spacing in the y-direction
    lambdaX = dict['lambdaX']        # Grid spacing in the x-direction
    msmoothx = dict['msmoothx']      # Smoothing length scale in the x-direction
    msmoothy = dict['msmoothy']      # Smoothing length scale in the y-direction
    msmootht = dict['msmootht']      # Smoothing length scale in time
    filtername = dict['filterName']  # Name of the filter type to smooth the data
                                     #      ['hanning', 'linloess', 'quadloess', 'boxcar', si']
    nmseitol = dict['nmseitol']      # Normalized error tolerance the user will tolerate in the final grid
                                     #      (0 - (no error) to 1 (no removal of bad points))

    xFRF_s = dict['xFRF_s']          # survey xFRF coordinates
    yFRF_s = dict['yFRF_s']          # survey yFRF coordinates
    Z_s = dict['Z_s']                # survey bottom elevations

    # if these are 2D, convert them
    if xFRF_s.ndim == 2 and yFRF_s.ndim == 2:
        xFRF_s = np.reshape(xFRF_s, (1, np.shape(xFRF_s)[0]*np.shape(xFRF_s)[1]))
        yFRF_s = np.reshape(yFRF_s, (1, np.shape(yFRF_s)[0] * np.shape(yFRF_s)[1]))
        Z_s = np.reshape(Z_s, (np.shape(Z_s)[0] * np.shape(Z_s)[1]))
        xFRF_s = xFRF_s[0]
        yFRF_s = yFRF_s[0]
        Z_s = Z_s[0]
    else:
        pass


    #### data checks ###########3
    filters = ['hanning', 'linloess', 'quadloess', 'boxcar', 'si']
    assert filtername in filters, 'Check filter name, not appropriate for current DEM generator function'
    ####################################################################
    # ############################### Load Data ########################
    ####################################################################
    t = DT.datetime.now()
    # I use my dictionary instead of the dataBuilder function from plant's code !!!!!
    # x, z = dataBuilder(filelist, data_coord_check='FRF')
    x = np.array([xFRF_s, yFRF_s, np.zeros(xFRF_s.size)]).T
    z = Z_s[:, np.newaxis]
    s = np.ones((np.size(x[:,1]),1))     # TODO estimate measurement error from the crab and incorporate to scripts
    print 'loading time is %s seconds' % (DT.datetime.now() - t)
    assert x.shape[0] > 1, 'Data Did not Load!'
    ####################################################################
    # Call grid builder to make a grid based on x,y min and max values #
    ####################################################################
    x_grid, y_grid = gridBuilder(x0, x1, y0, y1, lambdaX, lambdaY, dict['grid_coord_check'], dict['grid_filename'])
    t_grid = np.zeros_like((x_grid))  # Interpolate in time -- Not Developed Yet, but place holder there
    xi = np.array([x_grid.flatten(), y_grid.flatten(), t_grid.flatten()]).T  # grid locations, flatten make row-major style
    # now make smoothing array same shape as  xi
    xsm = msmoothx*np.ones_like(x_grid)
    ysm = msmoothy*np.ones_like(y_grid)
    tsm = msmootht*np.ones_like(t_grid)
    lx = np.array([xsm.flatten(), ysm.flatten(), tsm.flatten()]).T  # smoothing array , flatten takes row-major style

    N, M = np.shape(x_grid)
    x_out = x_grid[0,:].copy()  # grid coordinate for output
    y_out = y_grid[:,0].copy()  # grid coordinate for output

    del x_grid, y_grid, t_grid
    #####################################################################
    # subsample the data   ##############################################
    #####################################################################
    DXsmooth = np.array([msmoothx, msmoothy, msmootht])/4
    DXsmooth[2] = 1  # this hard codes a time smoothing of 1 (units unclear?)
    t = DT.datetime.now()
    Xi, zprime, si = subsampleData(x, z, s, DXsmooth)
    # a plot to compare original data to subsampled data
    # from matplotlib import pyplot as plt
    # plt.figure()
    # plt.subplot(211)
    # plt.plot(x[:,0], x[:,1], '.', label='Raw')
    # plt.plot(Xi[:,0], Xi[:,1], '.', label='SubSampled')
    # plt.legend()
    # plt.subplot(212)
    # plt.plot(np.sqrt(x[:, 0]**2 + x[:, 1]**2), z, '.', label='raw')
    # plt.plot( np.sqrt(Xi[:,0]**2 + Xi[:,1]**2), zprime, '.', label='subsampled')
    # plt.legend()
    # plt.close()

    # What's returned here
    print 'subsampling time is %s seconds' % (DT.datetime.now() - t)

    #####################################################################
    # Send it all into scalecinterpolation  -  Here is where the interpolation takes place
    #####################################################################
    t = DT.datetime.now()
    print 'Interpolating'
    zi, msei, nmsei, msri = scalecInterpTilePerturbations(Xi, zprime, si, xi, lx, filtername, nmseitol)
    print 'Interpolating time is %s seconds' % (DT.datetime.now() - t)

    # save the ouput and reshape
    # reshape
    zi = np.reshape(zi, (M, N)).T           # zi, the estimate
    msei = np.reshape(msei, (M, N)).T       # msei, the mean square interpolation error estimate (units of z)
    nmsei = np.reshape(nmsei, (M, N)).T     # nmsei, the normalized mean square error
    msri = np.reshape(msri, (M, N)).T       # msri, the mean square residuals

    # package to dictionary
    out = {'Zi': zi,
           'MSEi': msei,
           'NMSEi': nmsei,
           'MSRi': msri,
           'x_out': x_out,
           'y_out': y_out}
    return out


