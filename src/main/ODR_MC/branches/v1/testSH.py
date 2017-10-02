""" Test functions before saving in the corresponding scripts"""

import numpy as np
import scipy.odr as odr
import readHD as rhd
import upFun as upf

# AVHRR measurement equation
def avhrrME(CE,Cs,Cict,Lict,To,a0,a1,a2,a3):
        
    # Earth radiance from Earth counts and  calibration data
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) + a3*To
    return LE # return Earth radiance

""" AVHRR measurement model to use for series harmonisation: two virtual sensors 
for the data matrices, a block a rows has the specific sensors. """
def fcnH(coef, Xdata):
    # read data to variable names; transpose ndarrays 
    Lr1 = Xdata[0,:] # reference radiance 1st sensor; if reference-sensor pair
    Cs1 = Xdata[1,:] # space counts 1st sensor
    Cict1 = Xdata[2,:] # ICT counts 1st sensor
    CE1 = Xdata[3,:] # Earth counts 1st sensor
    Lict1 = Xdata[4,:] # ICT radiance 1st sensor
    To1 = Xdata[5,:] # orbit temperature 1st sensor
    Cs2 = Xdata[6,:] # space counts 2nd sensor
    Cict2 = Xdata[7,:] # ICT counts 2nd sensor
    CE2 = Xdata[8,:] # Earth counts 2nd sensor
    Lict2 = Xdata[9,:] # ICT radiance 2nd sensor
    To2 = Xdata[10,:] # orbit temperature 2nd sensor
    s1 = Xdata[11,:].astype(int) # 1st sensor index in sensors list (&coeff arr)
    s2 = Xdata[12,:].astype(int) # 2nd sensor's index 
    switch = np.logical_not(s1).astype(int) # 1 if s1=0 (reference), 0 otherwise
    
    """ returns questionable results """
    p = 4 # number of measurement eq. coefficients
    a01 = coef[s1*p + 0] # fit coefficients 1st sensor
    a11 = coef[s1*p + 1]
    a21 = coef[s1*p + 2]
    a31 = coef[s1*p + 3]
    a02 = coef[s2*p + 0] # fit coefficients 2nd sensor
    a12 = coef[s2*p + 1]
    a22 = coef[s2*p + 2]
    a32 = coef[s2*p + 3]
    
    # fit model 
    K = avhrrME(CE2,Cs2,Cict2,Lict2,To2,a02,a12,a22,a32) -\
        (1-switch) * avhrrME(CE1,Cs1,Cict1,Lict1,To1,a01,a11,a21,a31) -\
        switch * Lr1
    print 'Current iteration coefficients:', \
    np.unique([a01,a11,a21,a31,a02,a12,a22,a32])
    
    return K  

""" Read harmonisation data and perform ODR fit for the whole series """
def seriesODR(Hdata, Hrnd, b0, fb):
    # extract all variables but K from Hdata matrix
    selector = [x for x in range(Hdata.shape[1]) if x != 11]
    X = Hdata[:,selector].transpose() # transpose data matrix
    Y = Hdata[:,11] # adjustment values K
    
    VX = (Hrnd[:,0:11]**2).transpose() # squared random errors in X variables
    VY = Hrnd[:,11]**2 # K squared uncertainty 
    
    # run low-level odr; no weights 
    #fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
    print 'Running ODR for multiple pairs'    
    fit = odr.odr(fcnH,b0,Y,X,ifixb=fb,maxit=20,full_output=1)
    mFit = odr.Output(fit)

    # prepare data, model and run ODR 
    #odrData = odr.Data(X, Y, wd=1./VX, we = 1./VY)
    #odrMod = odr.Model(fcnH)
    #odrF = odr.ODR(odrData, odrMod, b0, ifixb=fb, maxit=20)
    #print 'Running ODR for multiple pairs'    
    #mFit = odrF.run()   
    
    return mFit # return ODR output

datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
#filelist = ["m02_n19.nc","n19_n15.nc", "m02_n15.nc"]
filelist = ["n19_n15.nc"]
nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist

slist = rhd.sensors(filelist) # list of sensors in filelist
print 'Sensors to be calibrated', slist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print 'Input calibration coefficients', inCoef

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
print 'Dimensions of H data matrix', Hd.shape
print 'Dimensions of H random error matrix', Hr.shape

''' Calculate radiance for the first sensor '''
avhrrNx = upf.avhrr(nop, nos)
Hd[:,0] = avhrrNx.measEq(Hd[:,1:6], inCoef[slist[0]])

# create ifixb array for beta array; fix a3 for all calibration sensors
parfix = np.zeros(nos*4, dtype=np.int)
for sidx in range(1,nos):
    parfix[4*sidx:4*sidx+3] = 1
fixb = parfix.tolist() # ifixb ODR parameter
print 'ifixb array', fixb, 'for sensors', slist

# apply odr on sensors in slist
sodr = seriesODR(Hd, Hr, coef, fixb)
print 'ODR output for sensors', slist
sodr.pprint()