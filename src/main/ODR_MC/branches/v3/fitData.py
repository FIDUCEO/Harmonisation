""" Perform the model fit """

import numpy as np
import readHD as rhd
import harFun as har
#import os
#import random

# set data folder and the list of netcdf files to process
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
print 'Data folder', datadir
filelist = ["m02_n19.nc","n19_n15.nc", "m02_n15.nc"]
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
print 'Input calibration coefficients for sensors', slist
print hsCoef
coef = hsCoef.flatten('F')
print 'Beta values for ODR', coef

# read netdcf files from filelist and create harmonisation data matrices
Im,Hd,Hr = rhd.rHData(datadir, filelist)    
print 'Index matrix of sensor pair'
print Im

# create ifixb array for beta array; fix a3 for all calibration sensors
parfix = np.zeros(nos*4, dtype=np.int)
for sidx in range(1,nos):
    parfix[4*sidx:4*sidx+3] = 1
fixb = parfix.tolist() # ifixb ODR parameter
print 'ifixb array', fixb, 'for sensors', slist

#selMU = random.sample(range(Hd.shape[0]), 500000)
#Hdata = Hd[selMU,:] # substruct records 
#Hrnd = Hr[selMU,:]
#print 'Dimension of subset Hdata matrix', Hdata.shape
#print 'Dimension of subset Hrnd error matrix', Hrnd.shape
#del Hd, Hr

# run odr for multiple sensor pairs
sodr = har.seriesODR2(Hd, Hr, coef, fixb)
print 'ODR output for sensors', slist
sodr.pprint()
