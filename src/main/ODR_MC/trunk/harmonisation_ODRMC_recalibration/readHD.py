""" FCDR sensor recalibration modules
    Project:        H2020 FIDUCEO
    Author:         Arta Dilo \NPL MM
    Reviewer:       Peter Harris \NPL MM, Sam Hunt \NPL ECO
    Date created:   06-10-2016
    Last update:    24-03-2017
    Version:        10.0
Read matchup data from a netCDF file of a sensors pair, 
extract netCDF variables to be used in harmonisation, 
calculate the other harmonisation variables from netCDF vars. 
The read function rHDpair uses netCDF4 package to read the netCDF file, 
and numpy ndarrays to store harmonisation data. """

from netCDF4 import Dataset
from numpy import array, mean as npmean, ones, zeros
from os.path import join as pjoin
from csv import reader as csvread


""" Extract sensors from a list of netCDF files """
def satSens(nclist):
   fsl = [] # list of sensor labels in nclist
   for netcdf in nclist:
       s1 = netcdf[0:3]  # 1st sensor' label
       s2 = netcdf[4:7] # 2nd sensor' label
       
       if s1 not in fsl: # if not yet in the list of sensors
           fsl.append(s1)  # add to the list    
       if s2 not in fsl:
           fsl.append(s2)
   return fsl

""" Get input calibration coefficients of sensors in netCDF filelist """
def sInCoeff(csvfolder, nclist, notd):
   fsl = satSens(nclist) # list of sensors in nc files 
   
   if notd: # known calib. coeffs for non-/ time dependency
       csvname = 'CalCoeff_notd.csv'
   else:
       csvname = 'CalCoeff.csv'
   cfn = pjoin(csvfolder, csvname)
   
   inCc = {} # dictionary of fsl sensors input coefficients
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
    print 'Opening netCDF file', pfn
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
