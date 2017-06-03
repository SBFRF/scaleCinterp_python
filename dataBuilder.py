# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 15:32:36 2014

@author: szthompson
"""

import scipy.io as sio
import numpy as np
import netCDF4 as nc
import pyproj

""" BUILD #1: """ 
def readInDataSet(filename):
    dataX, dataY, dataZ = [], [], []              
    # Handle NetCDF files
    if filename.endswith('nc'):
        try:
            ncfile = nc.Dataset(filename)
            dataX = ncfile['x'][:]
            dataY = ncfile['y'][:]
            dataZ = ncfile['z'][:]
        except IndexError:  # assume FRF ncfile keys
            dataX = ncfile['lon'][:]
            dataY = ncfile['lat'][:]
            dataZ = ncfile['elevation']
        except IOError:
            print '1 - cannot open', filename, 'it may not exist.  Please check path' 
    
    # Handle LAZ files

    elif filename.endswith('.laz'):
        print('is a LAZ file, this is not compatible now')
        # import laszip 
    # Handle LAS files
    elif filename.endswith('.las'):
        from liblas import file
        try:
            f = file.File(filename, mode = 'r')
            print 'Reading LAS'
            for p in f:
                if p.classification == 2:
                    dataX.append(p.x) 
                    dataY.append(p.y)
                    dataZ.append(p.z)
        except IOError: 
            print ('2 - cannot open', filename, 'it may not exist.  Please check path')
        
    elif filename.endswith('.xyz'):
        # open file, read lines
        f = open(filename)
        lines = f.readlines()
        f.close()
        x, y, z, = [], [], []
        # parse x y z of the lines
        for line in lines:
            split = line.split()
            dataX.append(float(split[0]))
            dataY.append(float(split[1]))
            dataZ.append(float(split[2]))

    # Handle Ascii files
    elif filename.endswith('.txt'):
        # Modified from... http://stackoverflow.com/questions/16155494/python-parse-lines-of-input-file
        try:
            with open(filename, 'r') as f:
                for line in f: # Parse the columnated data on each line
                    if line.find(" "): # Each data value on each line is seperated by a space        
                        info = line.split() # Split the data into variables based on seperation criteria: the space
                        #print info[0], info[1], info[2]
                        dataX.append(float(info[0]))
                        dataY.append(float(info[1]))
                        dataZ.append(float(info[2]))
        except IOError: 
            print '3 - cannot open', filename, 'it may not exist.  Please check path'
    
    # Handle Mat files
    elif filename.endswith('.mat'):
        try:
            matFile = sio.loadmat(filename)
            #dataX = matFile['lidar']['E'][0][0]
            #dataY = matFile['lidar']['N'][0][0]
            #dataZ = matFile['lidar']['Z'][0][0]
            dataX = matFile['xi'][:,0]
            dataY = matFile['xi'][:,1]
            dataZ = matFile['xi'][:,2]
        except IOError: 
            print '4 - cannot open', filename, 'it may not exist.  Please check path'
    else:
        print 'The file extension of,', filename,', is not supported. Please try again'
        print 'Supported file extension: \n.nc\n.laz\n.las\n.txt\n.mat'
        
    # Reshape the data... for compatibility
    dataX = np.reshape(dataX, (len(dataX),))
    dataY = np.reshape(dataY, (len(dataY),))
    dataZ = np.reshape(dataZ, (len(dataZ),))

    # Need to handle any NaNs?
    
    return dataX, dataY, dataZ
  
def dataBuilder(filelist, data_coord_check, EPSG=26918):
    """
    this function reads data sets and converts to
    This function assumes UTM zone 18N at the FRF in NAD83
    :param filelist:
    :param data_coord_check:
    :param EPSG:
    :return:
    """
    tempX, tempY, tempZ = [], [], [] 
    for files in filelist:
        dataX, dataY, dataZ = readInDataSet(files)
        tempX = np.concatenate((tempX,dataX))
        tempY = np.concatenate((tempY,dataY))
        tempZ = np.concatenate((tempZ,dataZ))
    
    x = np.array([tempX, tempY, np.zeros(tempX.size)]).T
    z = tempZ[:,np.newaxis]
    
    if (data_coord_check == 'LL'):
        UTM16N=pyproj.Proj("+init=EPSG:%s" % EPSG ) # 32618") # UTM coords, zone 18N, WGS84 datum  #
        [xutm,yutm] = UTM16N(tempX, tempY) # convert to UTM coord
        x = np.array([xutm, yutm, np.zeros(xutm.size)]).T
        # import utm
        # [xutm, yutm, zoneNum, zoneLet] = utm.from_latlon(tempX, tempY)
        # x = np.array([xutm, yutm, np.zeros(xutm.size)]).T
    return x, z

# load NOAA DEM
def loadNOAAdem(filename, x0, x1, y0, y1):
    from scipy.io import netcdf 
    f = netcdf.netcdf_file(filename, 'r')
    xtmp = f.variables['x'][:]
    ytmp = f.variables['y'][:]
    ztmp = f.variables['z'][:]
    
    [xtmp,ytmp] = np.meshgrid(xtmp,ytmp)
    xtmp = xtmp.flatten(1) 
    ytmp = ytmp.flatten(1) 
    ztmp = ztmp.flatten(1) 

    Xprior = xtmp[np.where((xtmp > x0) & (xtmp < x1) & (ytmp > y0) & (ytmp < y1))]
    Yprior = ytmp[np.where((xtmp > x0) & (xtmp < x1) & (ytmp > y0) & (ytmp < y1))]
    Zprior = ztmp[np.where((xtmp > x0) & (xtmp < x1) & (ytmp > y0) & (ytmp < y1))]
    
    import pyproj
    UTM16N=pyproj.Proj("+init=EPSG:32616") # UTM coords, zone 16N, WGS84 datum
    [Xprior, Yprior] = UTM16N(Xprior,Yprior)
    
    return Xprior, Yprior, Zprior