# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 11:25:54 2014

@author: jwlong
"""
import numpy as np  
import pyproj
    
def gridBuilder(x0, x1, y0, y1, dx, dy, grid_coord_check, grid_filename):

    if (grid_filename == ' '):
        if (grid_coord_check == 'LL'):
            # must convert to UTM (meters)
            UTM16N=pyproj.Proj("+init=EPSG:32616") # UTM coords, zone 16N, WGS84 datum
            x0, y0 = UTM16N(x0, y0)
            x1, y1 = UTM16N(x1, y1)
       
        x0 = np.round(x0, decimals=0)
        x1 = np.round(x1, decimals=0)  
        y0 = np.round(y0, decimals=0)  
        y1 = np.round(y1, decimals=0)  
        numGridPointsX = (x1 - x0)/dx
        numGridPointsY = (y1 - y0)/dy
        x_grid, y_grid = np.meshgrid(np.linspace(x0, x1, numGridPointsX),np.linspace(y0, y1, numGridPointsY))
        pass
    else:
        import scipy.io as sio
        try:
            gridFile = sio.loadmat(grid_filename) # Currently only works with MAT file
            x_grid = gridFile['xgrid'] 
            y_grid = gridFile['ygrid']
        except IOError:
            print 'The file,', grid_filename, ', does not exist in the path. Please try again.' 
            
    return x_grid, y_grid