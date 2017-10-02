#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 04-03-2017
    Version: 8.0
Evaluate fit coefficients uncertainty via Monte Carlo (MC) trials, 
uses full error structure to generate data in each MC trial. """

from numpy import empty, unique, savetxt
from datetime import datetime as dt
import readHD as rhd
import errStruct as mce 
import harFun as har


st = dt.now() # start of MS run

filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]
datadir = "/home/ad6/Data" # in eoserver
#datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
ncfile = filelist[4] # netCDF file to work with 

nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
s2 = ncfile[4:7]
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 
print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'
print 'Range of values of CEarth random uncertainties [',Hr[:,3].min(axis=0), ',', Hr[:,3].max(axis=0),']'


beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
beta[3] = inCoef[s2][3] # set a3 to the input value

# perform odr fit and extract output
fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb)
print '\n\nODR results on Jon\'s data weighted by random uncertainty'
podr.pprint()
bodr = podr.beta # odr fit coefficients - beta0

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
notr = 1000 # number of MC trials
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
    
    # run ODR on new X & Y vals and Hr weights 
    modr = har.odr4MC(Xdt, Ydt, Hr, bodr, fixb)
    
    # store ODR fit coefficients and reason for halting
    mcb[i, 0:p] = modr.beta
    mcb[i,p] = modr.info

print '\n\nODR results from the last MC trial'
modr.pprint()

fn = s2 + '_mcerrW_beta.txt'
savetxt(fn, mcb, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n\n\n--- Time taken for', notr, 'MC trials', (exect/60.), 'minutes ---'
