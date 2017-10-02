""" Test functions before saving in the corresponding scripts"""

import numpy as np
import scipy.odr as odr
import readHD as rhd

# AVHRR measurement equation
def avhrrME(CE,Cs,Cict,Lict,To,a0,a1,a2,a3):
        
    # Earth radiance from Earth counts and  calibration data
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) + a3*To
    return LE # return Earth radiance

""" Read harmonisation data and perform ODR fit for the whole series """
def seriesODR(Hdata, Hunc, b0, sensors, fb):
    # extract all variables but K from Hdata matrix
    selector = [x for x in range(Hdata.shape[1]) if x < 11]
    X = Hdata[:,selector].transpose() # transpose data matrix
    Y = Hdata[:,11] # adjustment values K
    bsens = sensors.transpose()
    
    VX = (Hunc[:,selector]**2).transpose() # squared random uncertainty X vars
    VY = Hunc[:,11]**2 # K squared uncertainty 
    
    def fcnH(coef, Xdata, sp=bsens):
        # read data to variable names; transpose ndarrays 
        Lr1 = Xdata[0,:] # reference radiance 1st sensor; 0 for sensor-sensor pair
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
        s1 = sp[0,:] # 1st sensor index in sensors list (&coeff arr)
        s2 = sp[1,:] # 2nd sensor's index 
        switch = (Lr1.astype(bool)).astype(int) # 1 if ref.sensor (Lr1>0), 0 otherwise
        
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

    # run low-level odr
    print 'Running ODR for multiple pairs'    
    #fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
    fit = odr.odr(fcnH,b0,Y,X,ifixb=fb,maxit=20)
    mFit = odr.Output(fit)
    
    return mFit # return ODR output

datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
filelist = ["m02_n15.nc", "n15_n14.nc"]
nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist

slist = rhd.sensors(filelist) # list of sensors in filelist
print 'Sensors to be calibrated', slist
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

Im,Hd,Hr,sp,mutime = rhd.rHData(datadir, filelist)
print 'Index matrix of sensor pair'
print Im
print 'Dimensions of H data matrix', Hd.shape
print 'Dimensions of H random error matrix', Hr.shape
print 'Dimensions of sensors', sp.shape

# create ifixb array for beta array; fix a3 for all calibration sensors
parfix = np.zeros(nos*4, dtype=np.int)
for sidx in range(1,nos):
    parfix[4*sidx:4*sidx+3] = 1
fixb = parfix.tolist() # ifixb ODR parameter
print 'ifixb array', fixb, 'for sensors', slist

# apply odr on sensors in slist
sodr = seriesODR(Hd, Hr, coef, sp, fixb)
print 'ODR output for sensors', slist
sodr.pprint()

print 'Range of input K values [', min(Hd[11]), max(Hd[11]), ']'
print 'Range of estimated Y values (K*) [', min(sodr.y), max(sodr.y), ']'
print 'Range of estimated epsilon (K err) [', min(sodr.eps), max(sodr.eps), ']'
print 'Range of input Lref values [', min(Hd[0]), max(Hd[0]), ']'
print 'Range of estimated Lref values [', min(sodr.xplus[0]), max(sodr.xplus[0]), ']'
print 'Range of estimated Lref delta err [', min(sodr.delta[0]), max(sodr.delta[0]), ']'
