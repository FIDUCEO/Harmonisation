""" Functions to read netCDF file of a sensor pair.
Use N-D arrays from xarray package to read and use netCDF data. """

import numpy as np
import pandas as pd
import xarray as xr
import os

print 'Current directory:', os.getcwd() # print current working directory
# Folder where the data is
datadir = "D:\Projects\FIDUCEO\Data\Simulated\Data"
print 'Setting the working directory to', datadir
os.chdir(datadir) # change working directory to datadir

# read a netCDF file
ncfile = 'm02_n15.nc'

ncdata = xr.open_dataset(ncfile) # read netCDF file 
print ncdata.data_vars # print xarray variables from netCDF file
Im = ncdata.lm # matchups index matrix 

didx = [0, 5, 6, 7, 8, 9] # columns with data in H, Ur, Us
Hdata = ncdata.H[:,didx] # data variables
# add column-names to Hdata 
# add K values to Hdata xarray

Hrnd = ncdata.Ur[:,didx] # random uncertainties
# add K random uncertainty, label columns
