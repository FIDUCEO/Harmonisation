""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 06-10-2016
    Last update: 22-01-2017
Read matchup data from a netCDF file of a sensors pair, 
extract netCDF variables to be used in harmonisation, 
calculate the other harmonisation variables from netCDF vars. 
The read function rHDpair uses netCDF4 package to read the netCDF file, 
and numpy ndarrays to store harmonisation data. """

import netCDF4 as nc
import numpy as np
import os
import csv
# import datetime

""" Read netcdf of a sensors' pair to harmonisation data arrays """
def rHDpair(folder, filename):     
    pfn = os.path.join(folder, filename) # filename with path  
    print 'Opening netCDF file', pfn
    ncid = nc.Dataset(pfn,'r')
    
    Im = ncid.variables['lm'][:] # matchup index array
    H = ncid.variables['H'][:,:] # harmonisation variables; empty vars included
    Ur = ncid.variables['Ur'][:,:] # random uncertainty for H vars
    Us = ncid.variables['Us'][:,:] # systematic uncertainty for H vars
    K = ncid.variables['K'][:] # evaluated K adjustment values
    Kr = ncid.variables['Kr'][:] # matchup random uncertainty
    Ks = ncid.variables['Ks'][:] # systematic uncertainty for K values
    corIdx = ncid.variables['CorrIndexArray'][:] # matchup time; internal format
    corLen = ncid.variables['corrData'][:] # length of averaging window

    ncid.close()   
    # Extract non-empty columns
    if Im[0,0] == -1: 
        rspair = 1 # reference-sensor pair
        didx = [0, 5, 6, 7, 8, 9] # non-empty columns in H, Ur, Us
        H = H[:,didx] 
        Ur = Ur[:,didx] 
        Us = Us[:,didx] 
    else:
        rspair = 0 # sensor-sensor pair
    
    # compile ndarrays of harmonisation data
    nor = Im[0,2]
    noc = H.shape[1] + 1 # plus one column for K data
    # harmonisation variables
    Hdata = np.zeros((nor, noc)) 
    Hdata[:,:-1] = H
    Hdata[:,noc-1] = K # adjustment values in last column
    # random uncertainties 
    Hrnd = np.zeros((nor, noc)) 
    Hrnd[:,:-1] = Ur
    Hrnd[:,noc-1] = Kr
    # systematic uncertainty
    Hsys = np.zeros((nor, noc)) 
    Hsys[:,:-1] = Us
    Hsys[:,noc-1] = Ks
    
    return rspair,Im,Hdata,Hrnd,Hsys,corIdx,corLen

""" Extract sensors from a list of netCDF files """
def sensors(nclist):
   fsl = [] # list of sensor labels in nclist
   for netcdf in nclist:
       s1 = netcdf[0:3]
       s2 = netcdf[4:7]
       
       if s1 not in fsl:
           fsl.append(s1)           
       if s2 not in fsl:
           fsl.append(s2)
   return fsl

""" Get input calibration coefficients of sensors in netCDF filelist """
def sInCoeff(csvfolder, nclist):
   fsl = sensors(nclist) # list of sensors in nc files 
   
   inCc = {} # dictionary of fsl sensors input coefficients
   cfn = os.path.join(csvfolder, 'CalCoeff.csv')
   with open(cfn, 'rb') as f:
       reader = csv.reader(f)
       reader.next() # skip header
       for row in reader:
           sl = row[0] # sensor label
           coefs = np.array(row[1:5]).astype('float')
           if sl in fsl:
               inCc[sl] = coefs
   return inCc

""" Create harmonisation data from netcdf files in the list.
Data matrix with fixed number of columns, and number of rows equal to the
number of matchups observations for the whole series. """                
def rHData(folder, nclist):
    fsl = sensors(nclist) # list of sensors; same order as coeffs array
    # initialise harmonisation data arrays 
    Im = np.empty([0,3], dtype=int) # matrix index for pairs
    Hdata = np.empty([0,12], order = 'C')
    Hrnd = np.empty([0,12], order = 'C')
    Hsys = np.empty([0,12], order = 'C')
    Hsen = np.empty([0,2], dtype=int, order = 'C')

    # loop through the netCDFs filelist 
    for ncfile in nclist:
        # read harmonisation data, variables and uncertainties
        rsp,lm,Hd,Hr,Hs,corIdx,corLen = rHDpair(folder, ncfile)
        print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
        
        nor = lm[0,2] # no of matchups for the current pair
        tH = np.zeros((nor, 12)) # temporary matrix for Hdata
        tHr = np.ones((nor, 12)) # for random err; fill empty values with 1 
        tHs = np.zeros((nor, 12)) # temporary for systematic error
        tsen = np.zeros((nor, 2), dtype=int) # array with sensor indices
        if rsp: # reference-sensor pair
            tH[:,0] = Hd[:,0] # reference radiance
            tH[:,1].fill(1.) # space counts 1st sensor
            tH[:,3].fill(1.) # Earth counts 1st sensor
            tH[:,6:12] = Hd[:,1:7] # 2nd sensor calibration data and K values
            tHr[:,0] = Hr[:,0] # reference radiance random error
            tHr[:,6:12] = Hr[:,1:7] # 2nd sensor random error data
            tHs[:,0] = Hr[:,0] # ref. radiance systematic error (if <>0)
            tHs[:,6:12] = Hr[:,1:7] # 2nd sensor systematic error data
            sl1 = 'm02'# 1st sensor label
        else: # sensor-sensor pair
            tH[:,1:12] = Hd[:,0:11] # 1st and 2nd sensor calibration data
            tHr[:,1:12] = Hr[:,0:11] # 1st and 2nd sensor random error
            tHs[:,1:12] = Hr[:,0:11] # 1st and 2nd sensor syetmatic error 
            sl1 = 'n' + str(lm[0,0]) # 1st sensor label
        sidx1 = fsl.index(sl1) # index in the list of sensors
        sl2 = 'n' + str(lm[0,1]) # 2nd sensor label
        sidx2 = fsl.index(sl2) # index in the list of sensors
        tsen[:,0] = sidx1 # 1st sensor's index in coeffs array 
        tsen[:,1] = sidx2 # 2nd sensor's index
        
        # add data of the current pair to harmonisation arrays
        Im = np.append(Im, lm, axis=0)
        Hdata = np.append(Hdata, tH, axis=0)
        Hrnd = np.append(Hrnd, tHr, axis=0)
        Hsys = np.append(Hsys, tHs, axis=0)
        Hsen = np.append(Hsen, tsen, axis=0)
        
    return Im, Hdata, Hrnd, Hsys, Hsen

""" Create harmonisation data for the whole series specific for AVHRR:
extract non-zero columns from the matrix of systematic errors and add 
as variables to the Hdata matrix. """
def ravhrrHD(folder, nclist):
    fsl = sensors(nclist) # list of sensors; same order as coeffs array
    # initialise harmonisation data arrays 
    Im = np.empty([0,3], dtype=int) # matrix index for pairs
    Hdata = np.empty([0,16], order = 'C')
    Hrnd = np.empty([0,16], order = 'C')
    Hsen = np.empty([0,2], dtype=int, order = 'C')

    # loop through the netCDFs filelist 
    for ncfile in nclist:
        # read harmonisation data, variables and uncertainties
        rsp,lm,Hd,Hr,Hs,corIdx,corLen = rHDpair(folder, ncfile)
        print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
        
        nor = lm[0,2] # no of matchups for the current pair
        tH = np.zeros((nor, 12)) # temporary matrix for Hdata
        tHr = np.ones((nor, 12)) # for random err; fill empty values with 1 
        tHs = np.zeros((nor, 12)) # temporary for systematic error
        tsen = np.zeros((nor, 2), dtype=int) # array with sensor indices
        if rsp: # reference-sensor pair
            tH[:,0] = Hd[:,0] # reference radiance
            tH[:,1].fill(1.) # space counts 1st sensor
            tH[:,3].fill(1.) # Earth counts 1st sensor
            tH[:,6:12] = Hd[:,1:7] # 2nd sensor calibration data and K values
            tHr[:,0] = Hr[:,0] # reference radiance random error
            tHr[:,6:12] = Hr[:,1:7] # 2nd sensor random error data
            tHs[:,0] = Hr[:,0] # ref. radiance systematic error (if <>0)
            tHs[:,6:12] = Hr[:,1:7] # 2nd sensor systematic error data
            sl1 = 'm02'# 1st sensor label
        else: # sensor-sensor pair
            tH[:,1:12] = Hd[:,0:11] # 1st and 2nd sensor calibration data
            tHr[:,1:12] = Hr[:,0:11] # 1st and 2nd sensor random error
            tHs[:,1:12] = Hr[:,0:11] # 1st and 2nd sensor syetmatic error 
            sl1 = 'n' + str(lm[0,0]) # 1st sensor label

        # merge systematic error data into variables matrix
        selector = [4,5,9,10] # Lict and To systematic err of both sensors
        tHd = np.concatenate((tH, tHs[:,selector]), axis=1)

        # fill array with sensors' indices 
        sidx1 = fsl.index(sl1) # index in the list of sensors
        sl2 = 'n' + str(lm[0,1]) # 2nd sensor label
        sidx2 = fsl.index(sl2) # index in the list of sensors
        tsen[:,0] = sidx1 # 1st sensor's index in coeffs array 
        tsen[:,1] = sidx2 # 2nd sensor's index
        
        # add data of the current pair to harmonisation arrays
        Im = np.append(Im, lm, axis=0)
        Hdata = np.append(Hdata, tHd, axis=0)
        Hrnd = np.append(Hrnd, tHr, axis=0)
        Hsen = np.append(Hsen, tsen, axis=0)
        
    return Im, Hdata, Hrnd, Hsen
