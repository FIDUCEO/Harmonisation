#!/usr/bin/env python

""" Perform the ODR fit on multiple pairs """

from os.path import join as pjoin
from datetime import datetime as dt
import numpy as np
import readHD as rhd
import harFun as har

# Time the execution of the whole script
st = dt.now() # start of MS run

# set data folder and the list of netcdf files to process
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
#datadir = "/home/ad6/Data" # in eoserver
#datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
filelist = ["m02_n19.nc","n19_n17.nc","n17_n15.nc","n19_n15.nc","m02_n15.nc"]

nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist
# create list of sensors from the netcdf filelist, and their input coefficients
slist = rhd.sensors(filelist) # list of sensors in filelist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
nos = len(slist) # number of sensors

# create array of beta start values from input coefficients 
hsCoef = np.zeros((nos,4))
for sno in range(nos):
    sl = slist[sno]
    hsCoef[sno,:] = inCoef[sl]
print '\nInput calibration coefficients for sensors', slist
print hsCoef
coef = hsCoef.flatten('A')
print '\nBeta values for ODR', coef

# read netdcf files from filelist and create harmonisation data matrices
Im,Hd,Hr,sp,mutime = rhd.rHData(datadir, filelist)
print '\nIndex matrix of sensor pair'
print Im

# create ifixb array for beta array; fix a3 for all calibration sensors
parfix = np.zeros(nos*4, dtype=np.int)
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

