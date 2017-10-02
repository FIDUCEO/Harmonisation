""" Functions to read a netCDF file and set of netCDF files and create
harmonisation data for the ODR package,
functions to run ODR an a reference-sensor pair,
function to calculate correlation between observations. """

import os
import netCDF4 as nc
import numpy as np
import scipy.odr as odr


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

""" Create harmonisation data from netcdf files in the list """                
def rHData(folder, nclist):
    fsl = sensors(nclist) # list of sensors; same order as coeffs array
    # initialise harmonisation data arrays 
    Im = np.empty([0,3], dtype=int) # matrix index for pairs
    Hdata = np.empty([0,14], order = 'F')
    Hrnd = np.empty([0,12], order = 'F')

    # loop through the netCDFs filelist 
    for ncfile in nclist:
        # read harmonisation data, variables and uncertainties
        rsp,lm,Hd,Hr,Hs,corIdx,corLen = rHDpair(folder, ncfile)
        print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
        
        nor = lm[0,2] # no of matchups for the current pair
        tH = np.zeros((nor, 14))
        tHr = np.zeros((nor, 12))
        if rsp: # reference-sensor pair
            tH[:,0] = Hd[:,0] # reference radiance
            tH[:,1].fill(1.) # space counts 1st sensor
            tH[:,3].fill(1.) # Earth counts 1st sensor
            tH[:,6:12] = Hd[:,1:7] # 2nd sensor calibration data
            tHr[:,0] = Hr[:,0] # reference radiance random uncertainty
            tHr[:,6:12] = Hr[:,1:7] # 2nd sensor random uncertainty data
            sl1 = 'm02'# 1st sensor label
        else: # sensor-sensor pair
            tH[:,1:12] = Hd[:,0:11] # 1st and 2nd sensor calibration data
            tHr[:,1:12] = Hr[:,0:11] # 1st and 2nd sensor random uncertainty 
            sl1 = 'n' + str(lm[0,0]) # 1st sensor label
        sidx1 = fsl.index(sl1) # index in the list of sensors
        sl2 = 'n' + str(lm[0,1]) # 2nd sensor label
        sidx2 = fsl.index(sl2) # index in the list of sensors
        tH[:,12] = sidx1 # 1st sensor's index in coeffs array 
        tH[:,13] = sidx2 # 2nd sensor's index
        
        # add data of the current pair to harmonisation arrays
        Im = np.append(Im, lm, axis=0)
        Hdata = np.append(Hdata, tH, axis=0)
        Hrnd = np.append(Hrnd, tHr, axis=0)
        
    return Im, Hdata, Hrnd
    
""" AVHRR measurement model """
def fcnAVHRRp(coef, data):
    a0 = coef[0] # AVHRR model coefficients
    a1 = coef[1]
    a2 = coef[2]
    a3 = coef[3]    
    # transposed ndarrays 
    Cs = data[0,:] # space counts
    Cict = data[1,:] # ICT counts
    CE = data[2,:] # Earth counts
    Lict = data[3,:] # ICT radiance
    To = data[4,:] # orbit temperature
    
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) + a3*To
    #print 'Current iteration coefficients:', [a0,a1,a2, a3]
    return LE  
    
""" Read harmonisation data and perform ODR fit """
def pairODR(Hdata, Hrnd, beta0, wgt):
    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    Xvars = Hdata[:,1:6].transpose() # transpose data matrix
    Yvar = Hdata[:,0] + Hdata[:,6] # reference radiance + adjustment values
    # errors in X and Y variables
    VX = (Hrnd[:,1:6]**2).transpose()
    VY = Hrnd[:,0]**2 + Hrnd[:,6]**2

    # data for (weighted) ODR 
    if wgt: # weighted ODR
        odrData = odr.Data(Xvars, Yvar, wd=1./VX, we = 1./VY)
        print 'Running weighted ODR for a sensors pair'
    else: # unweighted
        odrData = odr.Data(Xvars, Yvar)
        print 'Running unweighted ODR for a sensors pair'

    # model for odr fit
    odrMod = odr.Model(fcnAVHRRp)
    odrFit = odr.ODR(odrData, odrMod, beta0)
    
    mFit = odrFit.run()    
    return mFit # return ODR output

# Function to give error correlation between scan lines
# Operated on a case by case basis (does not assume inputs are arrays)
# USAGE:
#    corr_val = return_correlation(CorrIndexArray,corrData,i,j)
# where
#    CorrIndexArray  - Data from file (this name)
#    corrData        - Auxil data from file (this name)
#    i               - central scanline of averaging
#    j               - outer scanline of interest
# 
def return_correlation(index,corr_array,cent_pos,req_pos):

    diff = abs(index[cent_pos]-index[req_pos])
    if diff > corr_array[0]:
        return 0.
    else:
        return 1.-(diff/corr_array[0])


if __name__ == '__main__':
    
    datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
    ncfile = 'm02_n15.nc' # netCDF file
    initC = [-10., -4.e-3, 1.e-5, 1.e-3] #  initialise fit coefficients
    
    # read harmonisation data: variables and uncertainties
    rsp,Im,Hd,Hr,Hs,corIdx,corLen = rHDpair(datadir, ncfile)
    print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'

    # apply odr to data using random errors
    podr = pairODR(Hd, Hr, initC, 1)
    print 'NOAA-15 calibration coefficients and covariance:'
    podr.pprint()
