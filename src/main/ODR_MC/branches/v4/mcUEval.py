#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 05-02-2017
Evaluate fit coefficients uncertainty via Monte Carlo (MC) trials, 
uses full error structure to generate data in each MC trial. """

import numpy as np
from datetime import datetime as dt
from os.path import join as pjoin
import readHD as rhd
import harFun as har
import unpFun as upf
import errStruct as mce 
import plot as pl
import mcStats as vmc


st = dt.now() # start of MS run

filelist = ["m02_n15.nc"]#,"m02_n19.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
#datadir = "/home/ad6/Data" # in eoserver
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results

nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors

beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
ncfile = filelist[0] # netCDF file to work with 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value
print '\nProcessing sensors', slist, 'from file', ncfile

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 

# perform odr fit and extract output
fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb, Hs)
print '\n\nODR results on Jon\'s data weighted by combined random &systematic uncertainty'
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

# get unique scanlines, index of 1st matchup, number of matchups per scanline;
# based on time ordering of matchups the 1st matchup and number of matchups  
# per scanline let rebuild matchups from scanlines.
slt,midx,mcnt = np.unique(corIdx,return_index=True,return_counts=True)

# uncertainty of averaged calibration counts per scanline equal to
# Earth count uncertainty of 1st matchup in the scanline 
# divided by square root of of cal.counts sample size per scanline
uCC = Hr[midx,3] / np.sqrt(szCC) # calib.counts uncertainty from Earth count uncert.

# calculate gaps between scanlines and create blocks that are >25 scanlines apart
slarr,slblocks = mce.groupSln(slt,sldt,cLen)

# MC runs ODR on new data: best estimate + full correlation error draw
notr = 100 # number of MC iterations
b0 = np.empty([notr, len(beta)]) # array to store beta vals from MC
print '\n\nGenerate MC data with the full error structure.'

for i in range(notr):
    ''' Run notr MC trials '''

    # compile data for the ODR run; generate errors 
    errStr = mce.genPCS(Hr,sLict,sTo,uCC,wma,cLen,slarr,slblocks,mcnt) 
    # add errStr to X & Y best estimates
    Xdt = X.T + errStr[:,1:6] # X variables
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable
    
    # run ODR on new X & Y vals and Hr weights 
    modr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb, Hs)
    
    # store fit coefficients
    b0[i] = modr.beta

print '\n\nODR results from the last MC trial'
modr.pprint()

fn = 'mcb0ueval.txt'
np.savetxt(fn, b0, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n\n\nTime taken for', notr, 'MC trials', (exect/60.), 'minutes'

## -------LOAD text file where data was stored in eoserver/CEMS-------
mcb0 = np.loadtxt(fn, delimiter=',')
# print and plot odr coeffs stats from MC runs
mcmean,mcstd,mccov,mccor = vmc.mcStats(mcb0, b0odr, sd0odr) 
#mcmean,mcstd,mccov,mccor = vmc.mcStats(b0, b0odr, sd0odr) 


# create avhrr class instance; evaluate sensor radiance and uncertainty
avhrrNx = upf.avhrr(nop, nos) # instance of class avhrr series 
Figure = pl.Plot(pltdir, Im[0,0], Im[0,1], Im[0,2], corIdx) # Plot instance

# compile data for graphs of Jon's data and last MC trial 
inL = avhrrNx.measEq(Hd[:,1:6], inCoef[s2]) # input radiance
uX = np.sqrt(Hr[:,1:6]**2 + Hs[:,1:6]**2) # uncertainty of X variables
uY = np.sqrt(Hr[:,0]**2 + Hr[:,6]**2+Hs[:,6]**2) # Lref and K uncertainties

# ODR calibration of Jon;s data: 
calL = avhrrNx.measEq(Hd[:,1:6], b0odr) # calibrated radiance from fit
cLU = avhrrNx.va2ULE(Hd[:,1:6],b0odr,podr.cov_beta) # rad.uncert. from coeffs uncert.
LU = avhrrNx.uncLE(Hd[:,1:6],b0odr,uX,podr.cov_beta) # rad uncert. from data and coeffs

# create graphs that show ODR fit results on input data (Jon's) 
Figure.plotTrueErr(Hr, podr, 'odr') 
Figure.plotTrueVar(Hd, podr, 'odr')
Figure.plotFit(inL,calL,cLU,LU,uY,podr.eps,'odr')


# graphs that last MC trial fit results (but) on input data
Figure.plotTrueErr(Hr, modr, 'mcue') 
Figure.plotTrueVar(Hd, modr, 'mcue') # input data, NOT MC trial data


# Plot graphs of evaluated radiance and uncertainty from MC mean and cov vals
mccL = avhrrNx.measEq(Hd[:,1:6], mcmean) # calibrated Earth radiance
mccLU = avhrrNx.va2ULE(Hd[:,1:6],mcmean,mccov) # uncert from calib.coeffs
mcLU = avhrrNx.uncLE(Hd[:,1:6],mcmean,uX,mccov) # uncert from data and coeffs
Figure.plotFit(inL,mccL,mccLU,mcLU,uY,modr.eps,'mcue')
