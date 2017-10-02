#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 01-03-2017
Evaluate ODR fit coefficients uncertainty via Monte Carlo (MC) trials. 
Use the same error structure with weights in ODR, i.e. combined random and 
systematic uncertainty. Systematic effect is set constant for all the matchups 
in order to have comparable results with Peter & Sam's. """

import numpy as np
from datetime import datetime as dt
import readHD as rhd
import harFun as har

st = dt.now() # start of script execution

filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]
datadir = "/home/ad6/Data" # in eoserver
ncfile = filelist[0] # netCDF file to work with 

nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

# read data from the netCDF file
#rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(oldir, ncfile)
rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 
print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'
print 'Range of values of CEarth random uncertainties [',Hr[:,3].min(axis=0), ',', Hr[:,3].max(axis=0),']'

# perform odr fit and extract output
fixb = [1,1,1,0] # fix a3 coefficient
#podr = har.odrP(Hd, Hr, beta, fixb, Hs) # random & systematic uncert. for weights
podr = har.odrP(Hd, Hr, beta, fixb) # weights from random uncertainties
print '\nODR results on Jon data, weights from random uncertainty'
podr.pprint()

# extract ODR output and eval. sigma matrix for the MC
b0odr = podr.beta # odr fit coefficients
sd0odr = podr.sd_beta # standard error of fit coefficients 
cov0odr = podr.cov_beta # odr evaluated covariance matrix
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
sigma = np.sqrt(Hr**2)# + Hs**2) # standard uncertainty of X and Y vars
print '\nGenerate MC data with error distribution from input random uncertainty.'

# MC runs ODR on new data: best estimate + full correlation error draw
notr = 500 # number of MC trials
b0 = np.empty([notr, len(beta)]) # array to store beta vals from MC

''' Run notr MC trials '''
for i in range(notr):
    
    # generate errors for MC trial and add to odr best estimates
    errStr = np.random.normal(loc=0.0, scale=sigma) # errors in MC trial i
    Xdt = X.T + errStr[:,1:6] # X variables for odr in the MC trial
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable for odr 
    
    # run ODR on new X & Y vals with Hr (and Hs) weights 
    #mcodr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb, Hs)
    mcodr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb) # weights from random uncert.
    # store fit coefficients
    b0[i] = mcodr.beta

print '\nODR results from the last MC trial'
mcodr.pprint()

fn = s2 + '_mcrnd_b0.txt'
np.savetxt(fn, b0, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n--- Time taken for', notr, 'MC trials', (exect/60.), 'minutes ---'
