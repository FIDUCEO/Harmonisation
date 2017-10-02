""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 06-12-2016
    Last update: 26-01-2017
Harmonisation functions for a pair-wise implementation and for all the sensors 
together using odr package. Functions implement different LS fits:
- ordinary leas squares (OLS)
- weighted least squres (WLS)
- unweighted orthogonal distance regression (ODR)
- weighted ODR, an error in variables (EIV) method. """

import scipy.odr as odr
import numpy as np

# AVHRR measurement equation
def avhrrME(CE,Cs,Cict,Lict,To,a0,a1,a2,a3):
        
    # Earth radiance from Earth counts and  calibration data
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) + a3*To
    return LE # return Earth radiance

# dictionary with measurement eq. function of each sensors' series 
MEfunc = {'avhrr': avhrrME}

# model function, i.e. fcn argument in ODR package
def fcnP(coef, data, slabel='avhrr'):
    a0 = coef[0] # AVHRR model coefficients
    a1 = coef[1]
    a2 = coef[2]
    a3 = coef[3]    
    # transposed ndarrays 
    CE = data[2,:] # Earth counts
    Cs = data[0,:] # space counts
    Cict = data[1,:] # ICT counts
    Lict = data[3,:] # ICT radiance
    To = data[4,:] # orbit temperature
    
    LE = MEfunc[slabel](CE,Cs,Cict,Lict,To,a0,a1,a2,a3)
    #print 'Current iteration coefficients:', [a0,a1,a2, a3]
    return LE # return Earth radiance

""" Perform LS fit for a sensor-reference pair with low-level odr function """
def odrP(Hdata, Hr, b0, fb=None):
    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    X = Hdata[:,1:6].transpose() # transpose data matrix
    # Y variable: reference radiance + adjustment values
    Y = Hdata[:,0] + Hdata[:,6]

    # cacluate weights from uncertainty matrices
    VX = (Hr[:,1:6]**2).transpose() # sigma^2 of X on random unxcertainty
    # and of Y: assume only random error (and independence of ref.rad and K)
    VY = Hr[:,0]**2 + Hr[:,6]**2 
    
    # perform odr fit; fix coefficients as set in fb if given
    if fb: # fit with some coefficients fixed by fb
        fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
    else: # odr fit for all coefficients 
        fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,full_output=1)

    odrFit = odr.Output(fit) # get odr fit output     
    return odrFit # return odr output

""" Perform ODR over data MC generated with full error structure """
def odr4MC(Xdata, Ydata, Hr, b0, fb=None):
    X = Xdata.transpose()

    # sigma^2 of X, weight on random uncertainty 
    VX = (Hr[:,1:6]**2).transpose() 
    # sigma^2 of Y (assume reference radiance and K independence)
    VY = Hr[:,0]**2 + Hr[:,6]**2
        
    #  ODR on new X,Y data from best estimates & errors with complete structure
    odrData = odr.Data(X, Ydata, wd=1./VX, we = 1./VY)
    odrMod = odr.Model(fcnP)
    # compile data and model for the odr run
    if fb: # set coefficients to be fixed
        odrF = odr.ODR(odrData, odrMod, b0, ifixb=fb)
    else: # all coefficients are free
        odrF = odr.ODR(odrData, odrMod, b0)

    mFit = odrF.run() # run ODR    
    return mFit # return ODR output

""" Perform LS fit for a sensor-reference pair with ODR function """
def pairLS(Hdata, Hr, b0, fb=None, Hs=None, rsp=1):
    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    Xdata = Hdata[:,1:6].transpose() # transpose data matrix
    # Y variable, adjusted reference radiance: Lref + K
    Ydata = Hdata[:,0] + Hdata[:,6]

    # cacluate weights from uncertainty matrices
    if Hs is not None: # weight on both random and systematic uncertainty data
        VX = (Hr[:,1:6]**2 + Hs[:,1:6]**2).transpose() # sigma^2 of X variables
        
        ''' Y = Lref+K: assume independence of ref. radiance and K
        K random: in Hr matchups uncertainty, in Hs SRF shifting uncertainty '''
        VY = Hr[:,0]**2 + (Hr[:,6]**2+Hs[:,6]**2) # sigma^2 of Y
        
    else: # weight on random uncertainty 
        VX = (Hr[:,1:6]**2).transpose() # sigma^2 of X
        VY = Hr[:,0]**2 + Hr[:,6]**2 # Y sigma^2 (no shifting uncert.)

    odrData = odr.Data(Xdata, Ydata, wd=1./VX, we = 1./VY)
    odrMod = odr.Model(fcnP)
    # compile data and model for the odr run
    if fb: # set coefficients to be fixed
        odrF = odr.ODR(odrData, odrMod, b0, ifixb=fb)
    else: # all coefficients are free
        odrF = odr.ODR(odrData, odrMod, b0)

    mFit = odrF.run() # run ODR fit
    return mFit # return ODR output


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

""" Read harmonisation data and perform ODR fit for the whole series """
def seriesODR(Hdata, Hunc2, b0, fb):
    # extract all variables but K from Hdata matrix
    selector = [x for x in range(Hdata.shape[1]) if x != 11]
    X = Hdata[:,selector].transpose() # transpose data matrix
    Y = Hdata[:,11] # adjustment values K
    
    VX = (Hunc2[:,0:11]).transpose() # squared random errors in X variables
    VY = Hunc2[:,11] # K squared uncertainty 
    
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
