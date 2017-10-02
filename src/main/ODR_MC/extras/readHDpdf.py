""" Read netCDF of a sensor pair using netCDF4 package to read the file, 
and data frame from pandas to store harmonisation data. """

# import numpy as np
import netCDF4 as nc
import pandas as pd
import os

def rHDpair(folder, filename):
    print 'Setting the working directory to', folder
    os.chdir(folder) # change working directory to folder

    print 'Opening netCDF file', filename
    ncid = nc.Dataset(filename,'r')

    H = ncid.variables['H'][:,:] # harmonisation variables; empty vars included
    Ur = ncid.variables['Ur'][:,:] # random uncertainty for H vars
    Us = ncid.variables['Us'][:,:] # systematic uncertainty for H vars
    K = ncid.variables['K'][:] # K adjustment values
    Ks = ncid.variables['Ks'][:] # systematic uncertainty for K values
    Kr = ncid.variables['Kr'][:] # matchup random uncertainty
    mutime = ncid.variables['CorrIndexArray'][:]

    ncid.close()

    """ Store netCDF data in data frames:
    - Hdata for harmonisation variables including K values
    - Hrnd for random uncertainties for each variable
    - Hsys for systemtic uncertainties for each variable """
    didx = [0, 5, 6, 7, 8, 9] # columns with data in H, Ur, Us
    
    """ Create data frame of harmonisation variables from H matrix and K """
    Hdata = pd.DataFrame(data=H[:,didx],columns=['Lref','Cspace2','CICT2','CEarth2','LICT2','Torb2'])
    Hdata['Adj'] = K # add K values to the data frame
    del H, K     # delete the big arrays included in Hdata
    
    """ Create data frame of random uncertainties from Ur and Kr """
    Hrnd = pd.DataFrame(data=Ur[:,didx],columns=['Lref','Cspace2','CICT2','CEarth2','LICT2','Torb2'])
    Hrnd['Adj'] = Kr # add K random uncertainty to the data frame
    del Ur, Kr
    
    """ Create data frame of systematic uncertainties from Us and Ks """
    Hsys = pd.DataFrame(data=Us[:,didx],columns=['Lref','Cspace2','CICT2','CEarth2','LICT2','Torb2'])
    Hsys['Adj'] = Ks # add K random uncertainty to the data frame
    del Us, Ks
    
    return Hdata, Hrnd, Hsys
