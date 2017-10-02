#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 13-03-2017
    Version: 8.0
Evaluate fit coefficients uncertainty via Monte Carlo (MC) trials, 
uses full error structure to generate data in each MC trial. """

from numpy import empty, unique, savetxt
from datetime import datetime as dt
from os.path import join as pjoin
import readHD as rhd
import errStruct as mce 
import harFun as har


st = dt.now() # start of MS run

notime = False # work with (not-) time dependant simulation data
#datadir = "/home/ad6/Data" # root data folder in eoserver
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root folder in laptop
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results

filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]
nop = len(filelist) # number of sensor pairs
slist = rhd.satSens(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors

inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
beta = [-10., -4.e-3, 1.e-5, 0.0] #  initial values for ODR coefficients

ncfile = filelist[0] # netCDF file to work with 
s2 = ncfile[4:7]

if notime: # work with newSim_notime data
    
    # read data from the netCDF file
    newdir = pjoin(datadir, 'newSim_notime') # data folder
    rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(newdir, ncfile)
    Hr[:,5] = 1. # change 0 uncertainty of To for ODR to work
    
    podr = har.odrP(Hd, Hr, beta) # perform odr fit, weights from random uncert.
    
    fn = s2 + '_notd_mcerrst_beta.txt' # filename where MC coeffs will be stored

else: # work with data in the main data folder
    
    rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
    # set systematic uncertainties equivalent to Peter&Sam GN optimisation
    Hs = rhd.resetHs(Hs, rsp) 
    
    beta[3] = inCoef[s2][3] # set a3 to the input value
    
    # perform odr fit, weights from combined random and systematic uncertainty
    fixb = [1,1,1,0] # fix a3 coefficient
    fixx = [1,1,1,1,0] # fix orbit temperature, i.e. To variable
    podr = har.odrP(Hd, Hr, beta, fixb, fixx, Hs) # fit to adjusted ref.radiance
    
    fn = s2 + '_mcerrst_beta.txt' # filename for MC trials's coefficients 

print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

print '\nODR results on Jon data, weights from random uncertainty'
podr.pprint()
bodr = podr.beta # odr fit coefficients

''' Generate data for Monte Carlo run ''' 
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
sLict = Hs[0,4] # systematic error Lict
sTo = Hs[0,5] # systematic error To

# get scanlines and matchups per scanline
slt,midx,mcnt = unique(corIdx,return_index=True,return_counts=True)
rCSar = csUr[midx,:] # Cspace random uncert. per scanline: arrays of 51 slines
rCICTar = cictUr[midx,:] # Cict random uncert. per scanline: arrays of 51 slines
cLen = [25] # int(corLen[0]) # scanlines moving average half-window 

# MC runs ODR on new data: best estimate + full correlation error draw
notr = 3 # number of MC trials
p = len(beta) # number of calib. coefficients
mcb = empty([notr, p+1]) # array to store beta vals and ODR info of MC trials
print '\n\nGenerate MC data with the full error structure.'

for i in range(notr):
    ''' Run notr MC trials '''

    # compile data for the ODR run; generate errors 
    errStr = mce.genErr(Hr, sLict, sTo, rCSar, rCICTar, slt, cLen, mcnt)
    # add errStr to X & Y best estimates
    Xdt = X.T + errStr[:,1:6] # X variables
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable
    
    # run ODR on new X & Y vals and weights 
    if notime: # newSim_notime data: a3 = 0, To = 0, Hs = 0
        modr = har.odr4MC(Xdt, Ydt, Hr, bodr)
    else: # fix a3 and To to input, weights on random & systematic uncertainty
        modr = har.odr4MC(Xdt, Ydt, Hr, bodr, fixb, fixx, Hs)
    
    # store ODR fit coefficients and reason for halting
    mcb[i, 0:p] = modr.beta
    mcb[i,p] = modr.info

print '\n\nODR results from the last MC trial'
modr.pprint()

fn = pjoin(mcrdir, fn)
savetxt(fn, mcb, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n\n\n--- Time taken for', notr, 'MC trials', (exect/60.), 'minutes ---'
