#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo \ NPL MM
    Date created: 12-12-2016
    Last update: 16-03-2017
Perform ODR fit on multiple pairs. """

from os.path import join as pjoin
from datetime import datetime as dt
from numpy import zeros, ones, savetxt
import readHD as rhd
import harFun as har
import unpFun as upf
import visFun as vis

# Time the execution of the whole script
st = dt.now() # start time of script run

notime = False # work with (not-) time dependant simulation data

# set working directories
#datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
datadir = "/home/ad6/Data" # in eoserver
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results

# create list of sensors from the netcdf filelist
filelist = ["m02_n19.nc","n19_n17.nc","n19_n16.nc","n19_n15.nc",\
            "m02_n18.nc","n18_n17.nc","n18_n16.nc","n18_n15.nc",\
            "n17_n16.nc","n17_n15.nc","m02_n16.nc","n16_n15.nc",\
            "m02_n15.nc","n15_n14.nc","n14_n12.nc","n12_n11.nc","n11_n10.nc",\
            "n10_n09.nc","n09_n08.nc","n08_n07.nc"]
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

# store beta coefficients in a text file
fn = 'multip_beta.txt' 
fn = pjoin(mcrdir, fn)
savetxt(fn, mpbeta, delimiter=',')

# store covariance matrix of coefficients in a text file
fn = 'multip_bcov.txt' 
fn = pjoin(mcrdir, fn)
savetxt(fn, mpcov, delimiter=',')
