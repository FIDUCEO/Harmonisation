#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author:         Arta Dilo / NPL MM
    Date created:   24-01-2017
    Last update:    13-03-2017
    Version:        8.1
Evaluate ODR fit coefficients uncertainty via Monte Carlo (MC) trials. 
Work with a reference-sensor pair, use the same error structure with the weights 
in ODR, the combined random and systematic. Systematic effect is constant, 
set to average of all matchups, to be comparable with Peter & Sam's results. """

from numpy import empty, sqrt, savetxt, random, ones
from datetime import datetime as dt
from os.path import join as pjoin
import readHD as rhd
import harFun as har
import unpFun as upf
import plotMCres as mcr


st = dt.now() # start of script execution

notime = True # work with (not-) time dependant simulation data

datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
#datadir = "/home/ad6/Data" # main data folder in eoserver
#datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results
#pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]

# create instance of avhrr sensor series 
avhrrNx = upf.avhrr(datadir, filelist, notime)
p = avhrrNx.nocoefs # number of calibration parameters
m = avhrrNx.novars # # number of measured variables
nos = avhrrNx.nosensors # number of sensors in the series
slist = avhrrNx.sslab # list of sensors in the series
inCoef = avhrrNx.preHcoef # input coefficients to simulations

ncfile = filelist[4] # netCDF file to work with 
s2 = ncfile[4:7]
beta = inCoef[s2][0:p] # initial values for ODR coefficients

if notime: # work with newSim_notime data
    
    # read data from the netCDF file
    newdir = pjoin(datadir, 'newSim_notime') # data folder
    rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(newdir, ncfile)
    Hr[:,5] = 1. # change 0 uncertainty of To for ODR to work
    
    # perform odr fit, weights from random uncertainty
    podr = har.odrP(Hd, Hr, beta, avhrrNx) 
    
    # filename where MC coeffients will be stored
    fn = s2 + '_notd_mcerrst_beta.txt' 

else: # work with data in the main data folder
    
    rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
    # set systematic uncertainties equivalent to Peter&Sam GN optimisation
    Hs = rhd.resetHs(Hs, rsp) 
    
    # create ifixb array; fix a3 
    parfix = ones(p, dtype=int)
    parfix[-1] = 0
    fixb = parfix.tolist() # ifixb ODR parameter
    print '\nifixb array', fixb
    
    # create ifixx array; fix To variable
    varfix = ones(m, dtype=int)
    varfix[-1] = 0 # fix To 
    fixx = varfix.tolist() # ifixx ODR parameter
    print '\nifixx array', fixx

    # perform odr fit, weights from combined random and systematic uncertainty
    podr = har.odrP(Hd, Hr, beta, avhrrNx, fixb, fixx, Hs) # fit to adjusted ref.radiance
    
    fn = s2 + '_td_mcerrst_beta.txt' # filename for MC trials's coefficients 

print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

print '\nODR results on Jon data, weights from random uncertainty'
podr.pprint()
bodr = podr.beta # odr fit coefficients
covodr = podr.cov_beta # odr evaluated covariance matrix

# ODR output for MC trials: best estimates of X and Y variables
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To

# MC runs ODR on new data: best estimate + full correlation error draw
sigma = sqrt(Hr**2 + Hs**2) # standard uncertainty of X and Y vars
print '\n\nGenerate MC data with error distribution from input random uncertainty.'
notr = 3 # number of MC trials
mcb = empty([notr, p+1]) # array to store beta vals and ODR info of MC trials

for i in range(notr):
    ''' Run notr MC trials '''
    
    # generate errors for MC trial and add to odr best estimates
    errStr = random.normal(loc=0.0, scale=sigma) # errors in MC trial i
    if notime: # no-time dependant data: X matrix has no To column
        col = 5 # last column to read data from error matrix
    else:
        col = 6
    Xdt = X.T + errStr[:,1:col] # X variables
    ''' change To err to 0 for notime data '''
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable for odr 
    
    # run ODR on new X & Y vals and weights 
    if notime: # newSim_notime data: a3 = 0, To = 0, Hs = 0
        mcodr = har.odr4MC(Xdt, Ydt, Hr, bodr, avhrrNx)
    else: # fix a3 and To to input, weights on random & systematic uncertainty
        mcodr = har.odr4MC(Xdt, Ydt, Hr, bodr, avhrrNx, fixb, fixx, Hs)

    # store ODR fit coefficients and reason for halting
    mcb[i, 0:p] = mcodr.beta
    mcb[i,p] = mcodr.info

print '\nODR results from the last MC trial'
mcodr.pprint()

fn = pjoin(mcrdir, fn)
savetxt(fn, mcb, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n\n--- Time taken for', notr, 'MC trials', (exect/60.), 'minutes ---'


""" ===== LOAD text file with calib. coeffs in MC trials and plot ===== """

noMU = Im[0,2] # number of matchups
# graphs for ODR on simulated data, weighted residuals/errors for weight=1
nobj = 5000 # number of mathcup records to plot
mcr.plotSDfit(avhrrNx, noMU, nobj, s2, corIdx, Hr, podr, weight=1)


# heat maps for correlation of harmonisation coeffs from ODR and MC 
mcCov = mcr.mcCorr(fn, s2, bodr, covodr, 'MC with random error')

# graphs of radiance bias and MC uncertainty with 4*sigma error bars
nobj = 200 # number of mathcups to plot
mcr.plotMCres(avhrrNx,noMU,nobj,inCoef[s2],s2,Hr,Hs,podr,mcCov,4,'random error')
