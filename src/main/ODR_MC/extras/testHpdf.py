""" Test functions before saving in the corresponding scripts """

import numpy as np
import readHD as rhd
from harFun import *
import random
import pandas as pd


def avhrrME(Cs,Cict,CE,Lict,To,a0,a1,a2,a3):
    # Earth radiance from calibration data and Earth counts
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) + a3*To
    
    return LE # return Earth radiance

def fcndata(coef, Hdata):
    # read data to variable names; transpose ndarrays ??
    Lr1 = Hdata[0,:] # reference radiance 1st sensor; if reference-sensor pair
    Cs1 = Hdata[1,:] # space counts 1st sensor
    Cict1 = Hdata[2,:] # ICT counts 1st sensor
    CE1 = Hdata[3,:] # Earth counts 1st sensor
    Lict1 = Hdata[4,:] # ICT radiance 1st sensor
    To1 = Hdata[5,:] # orbit temperature 1st sensor
    Cs2 = Hdata[6,:] # space counts 2nd sensor
    Cict2 = Hdata[7,:] # ICT counts 2nd sensor
    CE2 = Hdata[8,:] # Earth counts 2nd sensor
    Lict2 = Hdata[9,:] # ICT radiance 2nd sensor
    To2 = Hdata[10,:] # orbit temperature 2nd sensor
    s1 = Hdata[12,:].astype(int) # index of 1st sensor in coefficients array
    s2 = Hdata[13,:].astype(int) # index 2nd sensor
    ref = np.invert(map(bool, s1)) # boolean: reference sensor
    cs = map(bool, s1) # boolean: calibration sensor
    
    #K = avhrrME(Cs2, Cict2, CE2, Lict2, To2, a02, a12, a22, a32) - \
    #    cs * avhrrME(Cs1, Cict1, CE1, Lict1, To1, a01, a11, a21, a31) - Lr1
    coLabs = ['Lr1','Cs1','Cict1','CE1','Lict1','To1','Cs2','Cict2','CE2','Lict2','To2']
    odrData = pd.DataFrame(data=Hdata[:,0:11], columns=coLabs)
    #odrData['a01'] = a01
    odrData['s1'] = s1
    odrData['s2'] = s2
    
    return odrData

datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
print 'Data folder', datadir 
filelist = ["m02_n19.nc","n19_n15.nc", "m02_n15.nc"]
nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist

slist = rhd.sensors(filelist) # list of sensors in filelist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
nos = len(slist)
hsCoef = np.zeros((nos,4))
for sno in range(nos):
    sl = slist[sno]
    hsCoef[sno,:] = inCoef[sl]
print 'Initial values fit coefficients for sensors', slist
print hsCoef
coef = hsCoef.flatten('F')
print 'Beta values for ODR', coef

Im,Hd,Hr = rhd.rHData(datadir, filelist)    
print 'Index matrix of sensor pair'
print Im
print 'Dimensions of Hd data matrix', Hd.shape

selMU = random.sample(range(Hd.shape[0]), 500000)
# substruct records 
Hdata = Hd[selMU,:]
Hrnd = Hr[selMU,:]
print 'Dimension of subset Hdata matrix', Hdata.shape
print 'Dimension of subset Hrnd error matrix', Hrnd.shape
del Hd, Hr

parfix = np.zeros(nos*4, dtype=np.int)
for sidx in range(1,nos):
    parfix[4*sidx:4*sidx+3] = 1
fixb = parfix.tolist() # ifixb ODR parameter
print 'ifixb array', fixb, 'for sensors', slist

#sodr = seriesODR(Hdata, Hrnd, coef, nos)
#print 'ODR output for sensors', slist
#sodr.pprint()

hd2odrf = fcndata(coef, Hdata)