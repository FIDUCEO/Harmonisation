#!/usr/bin/env python

""" Perform ODR fit on multiple pairs """

from os.path import join as pjoin
from datetime import datetime as dt
from numpy import zeros
import readHD as rhd
import harFun as har

# Time the execution of the whole script
st = dt.now() # start time of script run

# set data folder and the list of netcdf files to process
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
#datadir = "/home/ad6/Data" # in eoserver
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
filelist = ["m02_n19.nc","m02_n15.nc","n15_n14.nc"]

nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist
# create list of sensors from the netcdf filelist, and their input coefficients
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print '\nInput calibration coefficients for sensors', slist
print inCoef

# create array of beta start values from input coefficients 
b0 = [-10., -4.e-3, 1.e-5, 0.0]
hsCoef = zeros((nos,4))
sl = slist[0] # 1st sensor is the reference
hsCoef[0,:] = inCoef[sl] # give input coeffs
for sno in range(1, nos):
    #hsCoef[sno,:] = b0
    sl = slist[sno]
    #hsCoef[sno,3] = inCoef[sl][3]
    hsCoef[sno,:] = inCoef[sl]
print '\nInitial values of fit coefficients for sensors', slist
print hsCoef
coef = hsCoef.flatten('A')
print '\nInitial beta values for ODR', coef

# read netdcf files from filelist and create harmonisation data matrices
Im,Hd,Hr,Hs,sp,mutime,corL = rhd.rHData(datadir, filelist)
print '\nIndex matrix of sensor pair'
print Im

# create ifixb array for beta array; fix a3 for all calibration sensors
parfix = zeros(nos*4, dtype=int)
for sidx in range(1,nos):
    parfix[4*sidx:4*sidx+3] = 1
fixb = parfix.tolist() # ifixb ODR parameter
print '\nifixb array', fixb, 'for sensors', slist

# apply odr on sensors in slist
sodr = har.seriesODR(Hd, Hr, coef, sp, fixb)
print '\nODR output for sensors', slist
sodr.pprint()

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\nTime taken for multiple-pair fit of', nos, 'sensors', slist, (exect/60.), 'minutes'

print 'Range of input K values [', min(Hd[11]), max(Hd[11]), ']'
print 'Range of estimated Y values (K*) [', min(sodr.y), max(sodr.y), ']'
print 'Range of estimated epsilon (K err) [', min(sodr.eps), max(sodr.eps), ']'
print 'Range of input Lref values [', min(Hd[0]), max(Hd[0]), ']'
print 'Range of estimated Lref values [', min(sodr.xplus[0]), max(sodr.xplus[0]), ']'
print 'Range of estimated Lref delta err [', min(sodr.delta[0]), max(sodr.delta[0]), ']'

