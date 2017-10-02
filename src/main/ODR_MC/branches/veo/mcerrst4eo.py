#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 01-03-2017
Evaluate fit coefficients uncertainty via Monte Carlo (MC) trials, 
uses full error structure to generate data in each MC trial. """

import numpy as np
from datetime import datetime as dt
import readHD as rhd
import harFun as har
import errStruct as mce 

st = dt.now() # start of script run

filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]
datadir = "/home/ad6/Data" # in eoserver
ncfile = filelist[4] # netCDF file to work with 

nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
s2 = ncfile[4:7]
beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
beta[3] = inCoef[s2][3] # set a3 to the input value
#print '\nProcessing sensors', slist, 'from file', ncfile
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 
print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'
print 'Range of values of CEarth random uncertainties [',Hr[:,3].min(axis=0), ',', Hr[:,3].max(axis=0),']'


# perform odr fit and extract output
fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb)
print '\n\nODR results on Jon\'s data weighted by random uncertainty'
podr.pprint()
b0odr = podr.beta # odr fit coefficients - beta0
sd0odr = podr.sd_beta # standard error of fit coefficients - sigb0
cov0odr = podr.cov_beta # odr evaluated covariance matrix - covb0

''' Generate data for Monte Carlo run ''' 
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
sLict = Hs[0,4] # systematic error Lict
sTo = Hs[0,5] # systematic error To

# get unique scanlines, first matchup idx &number of matchup pixels per scanline
slt,midx,mcnt = np.unique(corIdx,return_index=True,return_counts=True)
rCSar = csUr[midx,:] # Cspace random uncert. per scanline: arrays of 51 slines
rCICTar = cictUr[midx,:] # Cict random uncert. per scanline: arrays of 51 slines

# calculate gaps between scanlines and create blocks that are >25 scanlines apart
cLen = int(corLen[0]) # scanlines moving average half-window 
wma = 1./(1+cLen*2) # moving average weight
sldt = 0.5 # time between consecutive scanlines; 0.5sec rounded up 
slarr,slblocks = mce.groupSln(slt,sldt,cLen)

# MC runs ODR on new data: best estimate + full correlation error draw
notr = 500 # number of MC trials
b0 = np.empty([notr, len(beta)]) # array to store beta vals from MC
print '\n\nGenerate MC data with the full error structure.'

for i in range(notr):
    ''' Run notr MC trials '''

    # compile data for the ODR run; generate errors 
    errStr = mce.genPCS(Hr,sLict,sTo,rCSar,rCICTar,wma,cLen,slarr,slblocks,mcnt) 
    # add errStr to X & Y best estimates
    Xdt = X.T + errStr[:,1:6] # X variables
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable
    
    # run ODR on new X & Y vals and Hr weights 
    modr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb)
    
    # store fit coefficients
    b0[i] = modr.beta

print '\n\nODR results from the last MC trial'
modr.pprint()

fn = s2 + '_mcerrst_b0.txt'
np.savetxt(fn, b0, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n--- Time taken for', notr, 'MC trials', (exect/60.), 'minutes ---'
