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

filelist = ["m02_n19.nc","m02_n16.nc","m02_n15.nc"]
rdir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
datadir = "D:\Projects\FIDUCEO\Data\Simulated\old_data" # data folder
#datadir = "/home/ad6/Data" # in eoserver
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results

nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
ncfile = filelist[0] # netCDF file to work with 
s2 = ncfile[4:7]

beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
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
# divided by square root of cal.counts sample size per scanline
uCC = Hr[midx,3] / np.sqrt(szCC) # calib.counts uncertainty from Earth count uncert.

# calculate gaps between scanlines and create blocks that are >25 scanlines apart
slarr,slblocks = mce.groupSln(slt,sldt,cLen)

# MC runs ODR on new data: best estimate + full correlation error draw
notr = 10 # number of MC iterations
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

fn = s2 + '_mcerrst_b0.txt'
np.savetxt(fn, b0, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n\n\nTime taken for', notr, 'MC trials', (exect/60.), 'minutes'

## -------LOAD text file where data was stored in eoserver/CEMS-------
mcb0 = np.loadtxt(fn, delimiter=',')
# print and plot odr coeffs stats from MC runs
mcmean,mcstd,mccov,mccor = vmc.mcStats(mcb0, b0odr, sd0odr) 
#mcmean,mcstd,mccov,mccor = vmc.mcStats(b0, b0odr, sd0odr) 

print '\nSample correlations of fit coefficients from MC trials'
print mccor

# create avhrr class instance; evaluate sensor radiance and uncertainty
avhrrNx = upf.avhrr(nop, nos) # instance of class avhrr series 

# compile data for graphs
inL = avhrrNx.measEq(Hd[:,1:6], inCoef[s2]) # input radiance
calL = avhrrNx.measEq(Hd[:,1:6], b0odr) # calibrated radiance from fit
# radiance uncertainty from coeffs uncert. evaluated by ODR
cLU = avhrrNx.va2ULE(Hd[:,1:6],b0odr,podr.cov_beta) 
# radiance uncertainty from coeffs uncert. evaluated via MC
mccLU = avhrrNx.va2ULE(Hd[:,1:6],b0odr,mccov) 
# radiance uncertainty from data and coeffs uncert. evaluated via MC
uX = np.sqrt(Hr[:,1:6]**2 + Hs[:,1:6]**2) # uncertainty of X variables
mcLU = avhrrNx.uncLE(Hd[:,1:6],b0odr,uX,mccov) 

# create graphs that show fit results
noMU = Im[0,2] # number of matchups
nobj = 1000 # numer of mathcups to plot
Figure = pl.Plot(pltdir, Im[0,0], Im[0,1], noMU, corIdx, nobj) # Plot instance

Figure.plotLenU(inL,calL,cLU,mccLU,mcLU,'errst')
#pl_ttl = "Input (green) and fitted radiance with 2 sigma uncertainty (gray) evaluated by ODR; full error structure"
#Figure.plotLbias(inL,calL,cLU,pl_ttl,'errst')
pl_ttl = "Input (green) and fitted radiance with 2 sigma uncertainty (gray) evaluated by MC; full error structure"
Figure.plotLbias(inL,calL,mccLU,pl_ttl,'errst')
#pl_ttl = "Input (green) and fitted radiance with 2 sigma uncertainty (gray) from data and fit coeffs; full error structure"
#Figure.plotLbias(inL,calL,mcLU,pl_ttl,'errst')
Figure.plotTrueErr(Hr, Hs, modr,'errst')


# graphs of radiance bias with 2sigma error bars
nobj = 200 # numer of mathcups to plot
plot_ttl = 'Radiance bias and ' + r'$2*\sigma$'+ ' uncertainty from MC eval. of coeffs\' covariance'
vmc.LbiasU(inL, calL, mccLU, noMU, nobj, plot_ttl) # from MC evaluation of uncertainty
plot_ttl = 'Radiance bias and ' + r'$2*\sigma$'+ ' uncertainty from data and coeffs\' uncertainty (by MC)'
vmc.LbiasU(inL, calL, mcLU, noMU, nobj, plot_ttl) # from MC evaluation of uncertainty
