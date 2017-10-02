#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo \ NPL MM
    Date created: 12-12-2016
    Last update: 16-03-2017
Perform ODR fit on multiple pairs. """

from os.path import join as pjoin
from datetime import datetime as dt
from numpy import zeros, ones
from random import sample
import readHD as rhd
import harFun as har
import unpFun as upf
import visFun as vis

# Time the execution of the whole script
st = dt.now() # start time of script run

notime = False # work with (not-) time dependant simulation data

# set working directories
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results
#pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs

# create list of sensors from the netcdf filelist
filelist = ["m02_n19.nc","m02_n15.nc","n15_n14.nc","n14_n12.nc"]
nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors

# create instance of avhrr sensor series 
avhrrNx = upf.avhrr(nop, nos) 
p = avhrrNx.nocoefs # number of calibration parameters
m = avhrrNx.novars # # number of measured variables

# read input coefficients of sensors in the netcdf filelist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print '\n\nInput calibration coefficients for sensors', slist
print inCoef

""" Create array of initial beta values for the ODR fit """
hsCoef = zeros((nos,p))
## Arbitrary beta0 values if nothing is known; will be far from solution
#b0 = [-10., -4.e-3, 1.e-5, 0.0]
#sl = slist[0] # 1st sensor is the reference
#hsCoef[0,:] = inCoef[sl] # give input coeffs
#for sno in range(1, nos):
#    hsCoef[sno,:] = b0
#    sl = slist[sno]
#    hsCoef[sno,3] = inCoef[sl][3]

# Keep the same values as input coefficients in inCoef 
for sno in range(nos):
    sl = slist[sno]
    hsCoef[sno,:] = inCoef[sl]

beta0 = hsCoef.flatten('A') # format to ODR input for initial values
print '\n\nInitial beta values for ODR'
print beta0

if notime: # work with newSim_notime data
    
    newdir = pjoin(datadir, 'newSim_notime') # data folder
    Im,Hd,Hr,Hs,sp,mutime,corL = rhd.rHData(newdir, filelist) # read netCDF files
    
    sodr = har.seriesODR(Hd, Hr, beta0, sp) # perform odr on sensors' list
    
else: # work with data in the main data folder
    
    Im,Hd,Hr,Hs,sp,mutime,corL = rhd.rHData(datadir, filelist)    

    # create ifixb arrays; fix a3 for all series' sensors
    parfix = zeros(nos*p, dtype=int)
    for sidx in range(1,nos):
        parfix[p*sidx:p*sidx+p-1] = 1
    fixb = parfix.tolist() # ifixb ODR parameter
    print '\n\nifixb array for sensors', slist
    print fixb
    
    # create ifixx arrays; fix orbit temperature To for sensors
    varfix = ones(m*2+1, dtype=int)
    varfix[m] = 0 # fix To for 1st sensor 
    varfix[2*m] = 0 # fix To for 2nd sensor 
    fixx = varfix.tolist() # ifixx ODR parameter
    print '\nifixx array', fixx

    # perform odr fit, weights from combined random and systematic uncertainty
    sodr = har.seriesODR(Hd, Hr, beta0, sp, fixb, fixx)

print '\nIndex matrix of sensors in', filelist
print Im
print '\nODR output for sensors', slist, '\n'
sodr.pprint()

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\nTime taken for multiple-pair fit of', nos, 'sensors from', slist, (exect/60.), 'minutes'

print '\n\nRange of input K values [', min(Hd[:,11]), max(Hd[:,11]), ']'
print 'Range of estimated Y values (K*) [', min(sodr.y), max(sodr.y), ']'
print 'Range of estimated epsilon (K error) [', min(sodr.eps), max(sodr.eps), ']'
print 'Range of input Lref values [', min(Hd[:,0]), max(Hd[:,0]), ']'
print 'Range of estimated Lref values [', min(sodr.xplus[0,:]), max(sodr.xplus[0,:]), ']'
print 'Range of estimated Lref delta error [', min(sodr.delta[0,:]), max(sodr.delta[0,:]), ']'

print '\nFirst row of H data matrix'
print Hd[0,:]
print '\nLast row of H data matrix'
print Hd[-1,:]

mpbeta = sodr.beta # calibration coeffs of fitted sensors
mpcov = sodr.cov_beta # coefficients covariance
mpcor = vis.cov2cor(mpcov) # coeffs' correlation matrix
cor_ttl = 'Correlation of harmonisation coefficients for pairs\n'+', '.join(filelist)
cor_lbl = ['a0', 'a1', 'a2', 'a3'] * nos
vis.plot_corr_heatmap(mpcor, title=cor_ttl, labels=['a0'])
print 'Correlation of harmonisation coefficients for pairs '+', '.join(filelist) +'\n'
print mpcor

nobj = 200
""" Extract coefficients and covariance of each sensor, 
compute and plot radiance with 4*sigma uncertainty """
for sno in range(1,nos): # loop through fitted sensor
    sl = slist[sno] # sensor label
    slab = int(sl[1:3]) # two-digits label in Im matrix
    
    sMidx, eMidx = rhd.sliceHidx(Im, slab) # 1st and last record index
    print 'Fisrt and last record for sensor', slab, '[', sMidx, eMidx,']'
    selMU = sample(xrange(sMidx, eMidx), nobj) # Select matchups for plotting    
    
    inC = inCoef[sl] # input coeffs to simulations
    print 'Input coefficients for sensor', slab, ':', inC
    inL = avhrrNx.measEq(Hd[selMU, m+1:2*m+1], inC) # radiance from input coeffs
    
    calC = mpbeta[sno*p:(sno+1)*p] # calib. coeffs for sensor slab
    print 'Fitted coefficients for sensor', slab, ':', calC
    calL = avhrrNx.measEq(Hd[selMU, m+1:2*m+1], calC) # calibrated radiance 
    
    covCC = mpcov[sno*p:(sno+1)*p,sno*p:(sno+1)*p] # coeffs covariance from odr
    print 'Covariance of coefficients for sensor', slab
    print covCC
    # radiance uncertainty from harmonisation
    cLU = avhrrNx.va2ULE(Hd[selMU, m+1:2*m+1],calC,covCC) 

    # graphs of radiance bias with 2sigma error bars
    plot_ttl = sl + ' Radiance bias and ' + r'$4*\sigma$'+ ' uncertainty from multiple-pairs ODR covariance'
    vis.LbiasU(inL, calL, cLU, 4, plot_ttl) 
