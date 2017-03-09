# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 14:28:56 2014

@author: jwlong
"""
import os
def list_files(datapath, datatype):
    # returns a list of names (with extension, without full path) of all files 
    # in folder path
    filelist = []
    listing = os.listdir(datapath)
    for files in listing:
        if files.endswith(datatype):
           filelist.append(files)
    return filelist 