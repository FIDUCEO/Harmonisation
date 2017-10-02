""" FIDUCEO FCDR harmonisation 
    Author:       Arta Dilo / NPL MM, Sam Hunt / NPL ECO
    Date created: 06-10-2016
    Last update:  15-03-2017
    Version:      8.0
Read matchup data from a netCDF file of a sensors pair, 
extract netCDF variables to be used in harmonisation, 
calculate the other harmonisation variables from netCDF vars. 
The read function rHDpair uses netCDF4 package to read the netCDF file, 
and numpy ndarrays to store harmonisation data. """

from datetime import datetime as dt
from netCDF4 import Dataset
from numpy import append, array, empty, mean as npmean, ones, zeros, where
from os.path import join as pjoin
from csv import reader as csvread


""" Extract sensors from a list of netCDF files """
def satSens(nclist):
   fsl = [] # list of satellite sensor labels in nclist
   for netcdf in nclist:
       s1 = netcdf[0:3]  # label of 1st sensor 
       s2 = netcdf[4:7] # label of 2nd sensor 
       
       if s1 not in fsl: # if not yet in the list of satellites
           fsl.append(s1)  # add to the list    
       if s2 not in fsl:
           fsl.append(s2)
   return fsl

""" Get input calibration coefficients of sensors in netCDF filelist """
def sInCoeff(csvfolder, nclist):
   fsl = satSens(nclist) # list of sensors in nc files 
   
   inCc = {} # dictionary of fsl sensors input coefficients
   cfn = pjoin(csvfolder, 'CalCoeff.csv')
   with open(cfn, 'rb') as f:
       reader = csvread(f)
       reader.next() # skip header
       for row in reader:
           sl = row[0] # sensor label
           coefs = array(row[1:5]).astype('float')
           if sl in fsl:
               inCc[sl] = coefs
   return inCc


""" Read netcdf of a sensors' pair to harmonisation data arrays """
def rHDpair(folder, filename):     
    pfn = pjoin(folder, filename) # filename with path  
    print '\nOpening netCDF file', pfn
    ncid = Dataset(pfn,'r')
    
    Im = ncid.variables['lm'][:] # matchup index array
    H = ncid.variables['H'][:,:] # harmonisation variables; empty vars included
    Ur = ncid.variables['Ur'][:,:] # random uncertainty for H vars
    Us = ncid.variables['Us'][:,:] # systematic uncertainty for H vars
    K = ncid.variables['K'][:] # evaluated K adjustment values
    Kr = ncid.variables['Kr'][:] # matchup random uncertainty
    Ks = ncid.variables['Ks'][:] # SRF uncertainty for K values
    corIdx = ncid.variables['CorrIndexArray'][:] # matchup time; internal format
    corLen = ncid.variables['corrData'][:] # length of averaging window
    # arrays of Cspace and Cict random uncertainty for 51 mav scanlines per each matchup
    CsUr = ncid.variables['cal_Sp_Ur'][:,:] # array for space counts
    CictUr = ncid.variables['cal_BB_Ur'][:,:] # array for ICT counts
    
    #print '\ncorrData value for calculating pixel-to-pixel correlation', corLen[0]
    ncid.close()   

    ''' Compile ndarrays of harmonisation data '''
    nor = Im[0,2] # number of matchups in the pair
    
    if Im[0,0] == -1: # reference-sensor pair
        rspair = 1 
        
        # Extract non-empty columns in data matrices, H, Ur, Us
        # 0 [Lref], 5 [Cspace], 6 [Cict], 7 [CEarth], 8 [Lict], 9 [To]
        didx = [0, 5, 6, 7, 8, 9] # non-empty columns in H, Ur, Us
        H = H[:,didx] 
        Ur = Ur[:,didx] 
        Us = Us[:,didx] 
        
        # create data and uncertainty arrays
        noc = H.shape[1] + 1 # plus one column for K data

        # data variables
        Hdata = zeros((nor, noc)) 
        Hdata[:,:-1] = H
        Hdata[:,noc-1] = K # adjustment values in last column        
        # random uncertainties 
        Hrnd = zeros((nor, noc)) 
        Hrnd[:,:-1] = Ur
        Hrnd[:,noc-1] = Kr
        # systematic uncertainty
        Hsys = zeros((nor, noc)) 
        Hsys[:,:-1] = Us
        Hsys[:,noc-1] = Ks
        
    else: # sensor-sensor pair
        rspair = 0 
        noc = H.shape[1] + 2 # plus two columns for K and Lref
        
        Hdata = zeros((nor, noc)) 
        Hdata[:,1:11] = H
        Hdata[:,noc-1] = K # adjustment values in last column       
        # random uncertainties 
        Hrnd = zeros((nor, noc)) 
        Hrnd[:,1:11] = Ur
        Hrnd[:,noc-1] = Kr
         # systematic uncertainty
        Hsys = zeros((nor, noc)) 
        Hsys[:,1:11] = Us
        Hsys[:,noc-1] = Ks
       
    return rspair,Im,Hdata,Hrnd,Hsys,corIdx,corLen, CsUr,CictUr
    #return rspair,Im,Hdata,Hrnd,Hsys,corIdx,corLen

""" Reset the systematic error matrix to constant values such that the ODR 
problem setting corresponds to Peter's LS optimisation problem. 
- Lref, Cspace, Cict, CEarth have 0 systematic uncertainty
- Lict and To have a constant systematic uncertainty 
- K random uncertainty from SRF shifting is stored in the Hs matrix """
def resetHs(Hs, rspair):
    Hsys = zeros(Hs.shape)
    nor = Hs.shape[0]

    # set Lict and To systematic to mean of corresponding Hs column
    # columns [4] and [5] are respectively Lict and To of series sensor in 
    # ref-sensor pair, and Lict and To of the 1st sensor in sensor-sensor pair
    sLict = npmean(Hs[:,4]) # mean Lict through all matchups
    sTo = npmean(Hs[:,5])  # mean To 
    Hsys[:,4] = sLict * ones(nor)
    Hsys[:,5] = sTo  * ones(nor)
        
    if rspair: # reference sensor pair, i.e. Hs has 7 columns     
        # keep K random uncertainty of SRF shifting   
        Hsys[:,6] = Hs[:,6] 

    else: # sensor-sensor pair, Hs has 12 columns
        Hsys[:,11] = Hs[:,11] # keep K random uncertainty from SRF shifting
        
        # set 2nd sensor' Lict and To systematic to mean of Hs values
        sLict = npmean(Hs[:,9]) 
        sTo = npmean(Hs[:,10]) 
        Hsys[:,9] = sLict * ones(nor)
        Hsys[:,10] = sTo * ones(nor)
       
    return Hsys # return the new set of sytematic uncertainties


""" Create harmonisation data from netcdf files in the list.
Data matrix with fixed number of columns (12) and number of rows equal to 
number of matchup observations for the whole series. """                
def rHData(folder, nclist):
    fsl = satSens(nclist) # list of sat.sensors; same order as coeffs array
    # initialise harmonisation data arrays 
    Im = empty([0,3], dtype=int) # matrix index for pairs
    Hdata = empty([0,12], order = 'C')
    Hrnd = empty([0,12], order = 'C')
    Hsys = empty([0,12], order = 'C')
    tidx = empty([0]) # array of matchups' scanline time
    psens = empty([0,2], dtype=int) # sensor labels in a pair

    # loop through the netCDFs filelist 
    for ncfile in nclist:
        # read harmonisation data, variables and uncertainties
        rsp,lm,Hd,Hr,Hs,corIdx,corLen,CsUar,CictUar = rHDpair(folder, ncfile)
        #rsp,lm,Hd,Hr,corIdx = rHDpair(folder, ncfile)
        print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
        
        nor = lm[0,2] # no of matchups for the current pair
        tH = zeros((nor, 12))
        tHr = ones((nor, 12)) 
        tHs = zeros((nor, 12))
        tps = zeros((nor, 2), dtype=int) 
        if rsp: # reference-sensor pair
            tH[:,0] = Hd[:,0] # reference radiance
            tH[:,1].fill(1.) # space counts 1st sensor
            tH[:,3].fill(1.) # Earth counts 1st sensor
            tH[:,6:12] = Hd[:,1:7] # 2nd sensor calibration data
            tHr[:,0] = Hr[:,0] # reference radiance random uncertainty
            tHr[:,6:12] = Hr[:,1:7] # 2nd sensor random uncertainty data
            tHs[:,6:12] = Hs[:,1:7] # 2nd sensor systematic uncertainty
            sl1 = 'm02'# 1st sensor label
        else: # sensor-sensor pair
            tH = Hd # data variables 
            tHr[:,1:12] = Hr[:,1:12] # 1st and 2nd sensor random uncertainty 
            tHs = Hs # systematic uncertainty 
            if lm[0,0] < 10:
                sl1 = 'n0' + str(lm[0,0]) # 1st sensor label
            else:
                sl1 = 'n' + str(lm[0,0]) # 1st sensor label            
        sidx1 = fsl.index(sl1) # index in the list of sensors

        if lm[0,1] < 10:
            sl2 = 'n0' + str(lm[0,1]) # 2nd sensor label
        else:
            sl2 = 'n' + str(lm[0,1]) # 2nd sensor label        
        sidx2 = fsl.index(sl2) # index in the list of sensors
        
        tps[:,0] = sidx1 # 1st sensor's index in coeffs array 
        tps[:,1] = sidx2 # 2nd sensor's index
        
        # add data of the current pair to harmonisation arrays
        Im = append(Im, lm, axis=0)
        Hdata = append(Hdata, tH, axis=0)
        Hrnd = append(Hrnd, tHr, axis=0)
        Hsys = append(Hsys, tHs, axis=0)
        psens = append(psens, tps, axis=0)
        tidx = append(tidx, corIdx)
        
    return Im, Hdata, Hrnd, Hsys, psens, tidx, corLen


""" Find rows of sensor matchups in H matrices of multiple pairs. Arguments to fun
- Im: matchups index array
- slab: 2 digits label of a sensor
- cidx: column in Im, i.e. 1st (cidx=0) or 2nd (cidx=1) sensor in a pair """
def sliceHidx(Im, slab, cidx=1):
    
    bidx = where(Im[:,1]==slab) # entries in Im for sensor slab
    sidx = bidx[0][0] # first occurance of slab in Im 
    
    startMU = 0 # start index of sensor matchups
    for i in range(sidx): # loop through Im rows until 1st slab occurance
        startMU += Im[i,2]
    endMU = startMU + Im[sidx,2]
    
    return startMU, endMU # return indices of first and last matchup record 

                
""" Sam's function for time conversion:
Return unix time in seconds (from 1970) from AVHRR time in seconds (from 1975)"""
def conv2date(inTime):
    # Calculate difference from AVHRR start time and unix start time in seconds
    start_time_AVHRR = dt(1975, 1, 1)
    start_time_unix = dt(1970, 1, 1)
    time_diff = (start_time_AVHRR - start_time_unix).total_seconds()

    # Convert to time from 1975 to date
    outTime = [dt.fromtimestamp(time+time_diff) for time in inTime]

    return outTime