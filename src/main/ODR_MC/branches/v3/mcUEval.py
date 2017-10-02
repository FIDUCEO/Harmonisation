#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 01-02-2017
Evaluate fit coefficients uncertainty via Monte Carlo (MC) trials. 
Start with odr fit of a reference-sensor pair. """

import numpy as np
from datetime import datetime as dt
import pandas as pd
import readHD as rhd
import harFun as har
import unpFun as upf
import mcErrStruct as mce 
import mcoutVis as vmc

filelist = ["m02_n15.nc"]#,"m02_n19.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
#datadir = "/home/ad6/Data" # in eoserver
#datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS
beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
slist = rhd.sensors(filelist) # list of sensors in filelist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
ncfile = filelist[0] # netCDF file to work with 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value
print 'Processing sensors', slist, 'from file', ncfile

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 

# perform odr fit and extract output
fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb, Hs)
print 'ODR results on Jon\'s data weighted by combined random &systematic uncertainty'
podr.pprint()
b0odr = podr.beta # odr fit coefficients - beta0
sd0odr = podr.sd_beta # standard error of fit coefficients - sigb0
cov0odr = podr.cov_beta # odr evaluated covariance matrix - covb0

''' Generate data and other info for Monte Carlo runs''' 
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
sLict = np.mean(Hs[:,4]) # systematic error Lict
sTo = np.mean(Hs[:,5]) # systematic error To
cLen = int(corLen[0]) # scanlines moving average half-window 
wma = 1./(1+cLen*2) # moving average weight
szCC = 10 # sample size of calibration counts per scanline
sldt = 1. # time between consecutive scanlines; 0.5sec rounded up 

# get unique scanlines, index of 1st matchup in scanline, 
# and indices to reconstruct the original array
slt,midx,invIdx = np.unique(corIdx,return_index=True,return_inverse=True)

# uncertainty of averaged calibration counts per scanline equal to
# Earth count uncertainty of 1st matchup in the scanline 
# divided by square root of of cal.counts sample size per scanline
uCC = Hr[midx,3] / np.sqrt(szCC) # calib.counts uncertainty from Earth count uncert.

# calculate gaps between scanlines and create blocks that are >25 scanlines apart
slarr,slblocks = mce.groupSln(slt,sldt,cLen)

# MC runs ODR on new data: best estimate + full correlation error draw
noit = 3 # number of MC iterations
b0 = np.empty([noit, len(beta)]) # array to store beta vals from MC

st = dt.now() # start of MS run
for i in range(noit):
    # compile data for the ODR run; generate errors 
    errStr = mce.genPCS(Hr,sLict,sTo,uCC,wma,cLen,slarr,slblocks,invIdx) 
    # add errStr to X & Y best estimates
    Xdt = X.T + errStr[:,1:6] # X variables
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable
    
    # run ODR on new X & Y vals and Hr weights 
    podr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb)
    
    # store fit coefficients
    b0[i] = podr.beta

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print 'Time taken for', noit, 'MC trials', (exect/60.), 'minutes'

b0mc = pd.DataFrame(b0)
b0mc.to_csv("mc_errstr_b0.csv")

# print and plot odr coeffs stats from MC runs; get coeff covariance from MC stats
mcmean, mccov = vmc.mcStats(b0, b0odr, sd0odr, cov0odr) 

# evaluate AVHRR Earth radiance and uncertainty
avhrrNx = upf.avhrr(1, 2) # instance of avhrr class
inLE = avhrrNx.measEq(X.T, inCoef[s2]) # input radiance (from input coeffs)

# evaluate radiance uncertainty from coeffs uncertainty
LUodr = avhrrNx.va2ULE(X.T, b0odr, cov0odr) # at best estimate with ODR evaluated cov
LUmc = avhrrNx.va2ULE(X.T, mcmean, mccov) # at best estimate with MC sample cov

CE = (X[0,:]-X[2,:]).transpose() # Earth counts

vmc.plotLaU(CE, inLE, Y, LUodr, LUmc)
