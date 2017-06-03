# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 11:25:54 2014

@author: jwlong
"""
import numpy as np  
import pyproj
import scipy.io as sio


def gridBuilder(x0, x1, y0, y1, dx, dy, grid_coord_check, grid_filename, EPSG=26918):
    """
        uses default to EPSG code of US EST coast (FRF location)
    This function could also use key word generation
    :param x0:
    :param x1:
    :param y0:
    :param y1:
    :param dx:
    :param dy:
    :param grid_coord_check:
    :param grid_filename:
    :param EPSG:
    :return:
    """
    if (grid_filename.strip() == ''):
        # build grid in UTM using dx dy, and 2 corners of grid(x0, y0)
        if (grid_coord_check.strip() == 'LL'):
            # must convert to UTM (meters)
            UTM16N=pyproj.Proj("+init=EPSG:%d" %EPSG) #:32616") # UTM coords, zone 16N, WGS84 datum
            x0, y0 = UTM16N(x0, y0)
            x1, y1 = UTM16N(x1, y1)
       
        # x0 = np.round(x0, decimals=0)  # why are these rounded ?  this moves the origin of the grid
        # x1 = np.round(x1, decimals=0)
        # y0 = np.round(y0, decimals=0)
        # y1 = np.round(y1, decimals=0)
        numGridPointsX = np.abs((x1 - x0)/dx)  # this assumes finite difference grid (points are vertex located)
        numGridPointsY = np.abs((y1 - y0)/dy)
        x_grid, y_grid = np.meshgrid(np.linspace(x0, x1, numGridPointsX), np.linspace(y0, y1, numGridPointsY))
        pass
    else:
        try:
            gridFile = sio.loadmat(grid_filename) # Currently only works with MAT file
            x_grid = gridFile['xgrid'] 
            y_grid = gridFile['ygrid']
        except IOError:
            print 'The file,', grid_filename, ', does not exist in the path. Please try again.' 
            
    return x_grid, y_grid