""" Test functions before saving in the corresponding scripts"""

import numpy as np
import scipy.odr as odr
import readHD as rhd

# AVHRR measurement equation
def avhrrME(CE,Cs,Cict,Lict,To,a0,a1,a2,a3):
        
    # Earth radiance from Earth counts and  calibration data
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) + a3*To
    return LE # return Earth radiance

def fcnH(coef, Xdata, sensors):
    # read data to variable names; transpose ndarrays 
    Cs1 = Xdata[0,:] # space counts 1st sensor
    Cict1 = Xdata[1,:] # ICT counts 1st sensor
    CE1 = Xdata[2,:] # Earth counts 1st sensor
    Lict1 = Xdata[3,:] # ICT radiance 1st sensor
    To1 = Xdata[4,:] # orbit temperature 1st sensor
    Cs2 = Xdata[5,:] # space counts 2nd sensor
    Cict2 = Xdata[6,:] # ICT counts 2nd sensor
    CE2 = Xdata[7,:] # Earth counts 2nd sensor
    Lict2 = Xdata[8,:] # ICT radiance 2nd sensor
    To2 = Xdata[9,:] # orbit temperature 2nd sensor
    K = Xdata[10,:] # adjustment values
    SeL1 = Xdata[11,:] # systematic error Lict1
    SeT1 = Xdata[12,:] # systematic error To1
    SeL2 = Xdata[13,:] # systematic error Lict2
    SeT2 = Xdata[14,:] # systematic error To2
    s1 = sensors[0,:] # 1st sensor index in sensors list (&coeff arr)
    s2 = sensors[1,:] # 2nd sensor's index 
    switch = np.logical_not(s1).astype(int) # 1 if ref.sensor, 0 otherwise
    
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
    Lref = avhrrME(CE2,Cs2,Cict2,Lict2-SeL2,To2-SeT2,a02,a12,a22,a32) -\
        (1-switch)* avhrrME(CE1,Cs1,Cict1,Lict1-SeL1,To1-SeT1,a01,a11,a21,a31)\
        - K
    #Lref = avhrrME(CE2,Cs2,Cict2,Lict2,To2,a02,a12,a22,a32) -\
    #    (1-switch)* avhrrME(CE1,Cs1,Cict1,Lict1,To1,a01,a11,a21,a31) - K
    print 'Current iteration coefficients:', \
    np.unique([a01,a11,a21,a31,a02,a12,a22,a32])
    
    return Lref  

""" Read harmonisation data and perform ODR fit for the whole series """
def seriesODR(Hdata, Hrnd, Hsys, sensors, b0, fb):
    # extract all variables but K from Hdata matrix
    selector = [4,5,9,10] # columns with systematic errors for Lict and To
    X = np.concatenate((Hdata[:,1:12], Hsys[:,selector]), axis=1) 
    X = X.transpose() # transpose matrix 
    Y = Hdata[:,0] # reference radiance
    
    nor = Y.shape[0] # number of rows
    # sigma2 errors in X variables; random error 1 for systematic corrections
    VX = np.concatenate((Hrnd[:,1:12]**2, np.ones((nor,4))), axis=1)
    VX = VX.transpose() # transpose matrix
    VY = Hrnd[:,0]**2 # reference radiance squared random error 
    
    # run low-level odr; no weights 
    print 'Running ODR for multiple pairs'    
    #fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
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
filelist = ["m02_n15.nc", "n15_n14.nc"]
nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist

slist = rhd.sensors(filelist) # list of sensors in filelist
print 'Sensors in the filelist', slist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print 'Input calibration coefficients:'
print inCoef

nos = len(slist)
hsCoef = np.zeros((nos,4))
for sno in range(nos):
    sl = slist[sno]
    hsCoef[sno,:] = inCoef[sl]
print 'Initial values fit coefficients for sensors', slist
print hsCoef
coef = hsCoef.flatten('A')
print 'Beta values for ODR', coef

Im,Hd,Hr,Hs,sensors = rhd.rHData(datadir, filelist)    
print 'Index matrix of sensor pair'
print Im
print 'Dimensions of H data matrix', Hd.shape
print 'Dimensions of H random error matrix', Hr.shape
print 'Dimensions of H systematic error matrix', Hs.shape
sensors = sensors.transpose() # transpose matrix for use in fcn
print 'Dimensions of sensors matrix', sensors.shape

# create ifixb array for beta array; fix a3 for all calibration sensors
parfix = np.zeros(nos*4, dtype=np.int)
for sidx in range(1,nos):
    parfix[4*sidx:4*sidx+3] = 1
fixb = parfix.tolist() # ifixb ODR parameter
print 'ifixb array', fixb, 'for sensors', slist

selector = [4,5,9,10] 
X = np.concatenate((Hd[:,1:12], Hs[:,selector]), axis=1) 
print 'Range of Lref values [', min(Hd[:,0]), max(Hd[:,0]), ']'
eps = fcnH(coef, X.T, sensors) - Hd[:,0] # Lr* - Lr
print 'Range of Lr* - Lr values [', min(eps), max(eps), ']'

# apply odr on sensors in slist
sodr = seriesODR(Hd, Hr, Hs, sensors, coef, fixb)
print 'ODR output for sensors', slist
sodr.pprint()