#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 21-02-2017
Evaluate ODR fit coefficients uncertainty via Monte Carlo (MC) trials. 
Work with a reference-sensor pair, use the same error structure with the weights 
in ODR, the combined random and systematic. Systematic effect is constant, 
set to average of all matchups, to be comparable with Peter & Sam's results. """

import numpy as np
import random
from datetime import datetime as dt
from os.path import join as pjoin
import readHD as rhd
import harFun as har
import unpFun as upf
import mcStats as vmc


filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
oldir = "D:\Projects\FIDUCEO\Data\Simulated\old_data" # old data folder
#datadir = "/home/ad6/Data" # in eoserver
#datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results

nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors

beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
ncfile = filelist[4] # netCDF file to work with 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value
#print '\nProcessing sensors', slist, 'from file', ncfile
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(oldir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 
print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'
print 'Range of values of CEarth random uncertainties [',Hr[:,3].min(axis=0), ',', Hr[:,3].max(axis=0),']'

# perform odr fit and extract output
fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb, Hs)
print '\nODR results on Jon data, weights from combined random &systematic uncertainty'
podr.pprint()

# extract ODR output 
b0odr = podr.beta # odr fit coefficients
sd0odr = podr.sd_beta # standard error of fit coefficients 
cov0odr = podr.cov_beta # odr evaluated covariance matrix
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To


st = dt.now() # start of MC run

''' Generate errors for Monte Carlo runs''' 
Hs[:,0] = 0 # set Lref systematic uncertainty to 0
sigma = np.sqrt(Hr**2 + Hs**2) # standard uncertainty of X and Y vars
print '\nGenerate MC data with error distribution from input random and systematic uncertainty.'

# MC runs ODR on new data: best estimate + full correlation error draw
notr = 5 # number of MC trials
b0 = np.empty([notr, len(beta)]) # array to store beta vals from MC

for i in range(notr):
    ''' Run notr MC trials '''
    
    # generate errors for MC trial and add to odr best estimates
    errStr = np.random.normal(loc=0.0, scale=sigma) # errors in MC trial i
    Xdt = X.T + errStr[:,1:6] # X variables for odr in the MC trial
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable for odr 
    
    # run ODR on new X & Y vals with Hr (and Hs) weights 
    mcodr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb, Hs)
    # store fit coefficients
    b0[i] = mcodr.beta

print '\nODR results from the last MC trial'
mcodr.pprint()

#fn = pjoin(mcrdir, 'mcodrb0.txt')
fn = s2 + '_mcrnd_b0.txt'
#np.savetxt(fn, b0, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n--- Time taken for', notr, 'MC trials', (exect/60.), 'minutes ---'


## -------LOAD text file where data was stored in eoserver/CEMS-------
mcb0 = np.loadtxt(fn, delimiter=',')
mcmean,mcstd,mccov,mccor = vmc.mcStats(s2, mcb0, b0odr, sd0odr) 
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
mcLU = avhrrNx.uncLE(Hd[:,1:6],b0odr,uX,mccov) # radiance uncertainty

noMU = Im[0,2] # number of matchups
# graphs of radiance bias with 2sigma error bars
nobj = 400 # number of mathcups to plot
selMU = random.sample(range(noMU), nobj) # Select matchups for plotting
k = 2 # k value for sigma uncertainty range
plot_ttl = s2 + ' Radiance bias and ' + r'$2*\sigma$'+ ' uncertainty from ODR eval. of coeffs\' covariance'
vmc.LbiasU(inL[selMU], calL[selMU], cLU[selMU], k, plot_ttl) # from ODR evaluation of uncertainty
plot_ttl = s2 + ' Radiance bias and ' + r'$2*\sigma$'+ ' uncertainty from MC eval. of coeffs\' covariance'
vmc.LbiasU(inL[selMU], calL[selMU], mccLU[selMU], k, plot_ttl) # from MC evaluation of uncertainty
plot_ttl = s2 + ' Radiance bias and ' + r'$2*\sigma$'+ ' uncertainty from data and coeffs\' uncertainty (by MC)'
vmc.LbiasU(inL[selMU], calL[selMU], mcLU[selMU], k, plot_ttl) # from MC evaluation of uncertainty
plot_ttl = s2 + ' Radiance uncertainty bias from coeffs\' uncertainty evaluated by MC and by ODR'
vmc.LUdiff(inL[selMU], cLU[selMU], mccLU[selMU], plot_ttl)
