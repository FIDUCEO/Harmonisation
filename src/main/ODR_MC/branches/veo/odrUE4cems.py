#!/usr/bin/env python

import numpy as np
import readHD as rhd
import harFun as har

##---- Create the data that is needed for each MC trial; to be run once ------
filelist = ["m02_n15.nc"]
datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS

beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
slist = rhd.sensors(filelist) # list of sensors in filelist
ncfile = filelist[0] # netCDF file to work with 
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 

fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb, Hs) # run odr with random & syst weights
b0odr = podr.beta # odr fit coefficients

##------------ an MC trial, each to be run in parallel at 500 CEMS nodes ------
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
sigma = np.sqrt(Hr**2 + Hs**2) # standard uncertainty of X and Y vars

# generate errors for MC trial and add to X & Y best estimates
errStr = np.random.normal(sigma) 
Xdt = X.T + errStr[:,1:6] # X variables for odr in the MC trial
Ydt = Y + errStr[:,0] + errStr[:,6] # Y var: Lref + K

# run ODR on new X & Y vals and Hr&Hs weights
mcodr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb, Hs) # ODR result
b0 = mcodr.beta # beta coeffs for the MC trial

# OUTPUT i need from each CEMS job: Calibration Coefficients 
print b0
##----------