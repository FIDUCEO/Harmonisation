""" Perform model fitting for reference-sensor pairs """

import readHD as rhd
from harFun import *
#import os

#print 'Current directory', os.getcwd() # print current working directory
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
#filelist = ["m02_n19.nc", "m02_n17.nc", "m02_n16.nc", "m02_n15.nc"]
filelist = ["m02_n15.nc"]
nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist

slist = rhd.sensors(filelist) # list of sensors in filelist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print 'Sensors in the file list', slist
print 'Input calibration coefficients', inCoef

initC = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
for ncfile in filelist:
    s2 = ncfile[4:7]
    initC[3] = inCoef[s2][3] # set a3 to the input value
    invCoef = (1./inCoef[s2]).tolist() # inverse coefficients to use for scaling
    print 'Input coefficients for NOAA-', s2, ':', inCoef[s2]
    print 'Inverted values to use for scaling beta:', invCoef
    
    # read harmonisation data: variables and uncertainties
    rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
    print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
    
    VX = Hr[:,1:6]**2
    VY = Hr[:,0]**2 + Hr[:,6]**2
    print 'Range of inverse Y variance:', max(1./VY)
    print 'Range of inverse X variances:'
    for cidx in range(0,5):
        print max(1./VX[:,cidx])
    
    # apply ODR to data without error' weighting 
    podr = pairODR(Hdata, Hrnd, initC, 0)
    print 'ODR output for', s2
    podr.pprint()
    
    ## apply weighted ODR using X & Y errors
    #podr = pairODR(Hd, Hr, initC, 1)
    #print 'Weighted ODR output for', s2, 'by fixing a3'
    #podr.pprint()
    
#    # apply OLS to harmonisation data
#    pols = pairOLS(Hd, Hr, initC, 0)
#    print 'OLS output for', s2
#    pols.pprint()
#
#    # apply WLS to harmonisation data
#    pwls = pairOLS(Hd, Hr, initC, 1)
#    print 'WLS output for', s2
#    pwls.pprint()

# need to run low level odr with full_output=1 for delta and epsilon errors
