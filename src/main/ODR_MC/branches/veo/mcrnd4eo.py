#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 01-03-2017
Evaluate ODR fit coefficients uncertainty via Monte Carlo (MC) trials. 
Use the same error structure with weights in ODR, i.e. combined random and 
systematic uncertainty. Systematic effect is set constant for all the matchups 
in order to have comparable results with Peter & Sam's. """

from numpy import empty, sqrt, savetxt, random, ones
from datetime import datetime as dt
from os.path import join as pjoin
import readHD as rhd
import harFun as har
import unpFun as upf

st = dt.now() # start of script execution

notime = False # work with (not-) time dependant simulation data
datadir = "/home/ad6/Data" # in eoserver
#datadir = "D:\Projects\FIDUCEO\Data\Simulated" # for test in laptop
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results

filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]
nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors

# create instance of avhrr sensor series 
avhrrNx = upf.avhrr(nop, nos) 
p = avhrrNx.nocoefs # number of calibration parameters
m = avhrrNx.novars # # number of measured variables

# simulation input coeffs, and initial ODR beta coeffs
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
beta = [-10., -4.e-3, 1.e-5, 0.0] #  initial values for ODR coefficients

ncfile = filelist[4] # netCDF file to work with 
s2 = ncfile[4:7]

if notime: # work with newSim_notime data
    
    # read data from the netCDF file
    newdir = pjoin(datadir, 'newSim_notime') # data folder
    rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(newdir, ncfile)
    Hr[:,5] = 1. # change 0 uncertainty of To for ODR to work
    
    podr = har.odrP(Hd, Hr, beta) # perform odr fit, weights from random uncert.
    
    fn = s2 + '_notd_mcrnd_beta.txt' # filename where MC coeffs will be stored

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
    beta[3] = inCoef[s2][3] # set a3 to the input value
    podr = har.odrP(Hd, Hr, beta, fixb, fixx, Hs) # fit to adjusted ref.radiance
    
    fn = s2 + '_mcrnd_beta.txt' # filename for MC trials's coefficients 

print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

print '\n\nODR results on Jon\'s simulated data'
podr.pprint()
bodr = podr.beta # odr fit coefficients
covodr = podr.cov_beta # odr evaluated covariance matrix

# ODR output for MC trials: best estimates of X and Y variables
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To


# MC runs ODR on new data: best estimate + full correlation error draw
sigma = sqrt(Hr**2 + Hs**2) # standard uncertainty of X and Y vars
print '\n\nGenerate MC data with error distribution from input random uncertainty.'
notr = 500 # number of MC trials
mcb = empty([notr, p+1]) # array to store beta vals and ODR info of MC trials

for i in range(notr):
    ''' Run notr MC trials '''
    
    # generate errors for MC trial and add to odr best estimates
    errStr = random.normal(loc=0.0, scale=sigma) # errors in MC trial i
    Xdt = X.T + errStr[:,1:6] # X variables for odr in the MC trial
    ''' change To err to 0 for notime data '''
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable for odr 
    
    # run ODR on new X & Y vals and weights 
    if notime: # newSim_notime data: a3 = 0, To = 0, Hs = 0
        mcodr = har.odr4MC(Xdt, Ydt, Hr, bodr)
    else: # fix a3 and To to input, weights on random & systematic uncertainty
        mcodr = har.odr4MC(Xdt, Ydt, Hr, bodr, fixb, fixx, Hs)

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
