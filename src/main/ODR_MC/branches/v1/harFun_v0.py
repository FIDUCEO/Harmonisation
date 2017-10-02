""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 06-12-2016
    Last update: 06-01-2017
Functions for LS fits for a pair reference-sensor and all series sensors.
-    ordinary leas squares (OLS)
-    weighted least squres (WLS)
-    unweighted orthogonal distance regression (ODR)
-    weighted ODR, an error in variables (EIV) method.
Harmonisation is implemented in a pair-wise mode and (hopefully) 
all the sensors together using odr package. """

import numpy as np
import scipy.odr as odr

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
def odrP(Hdata, Hrnd, b0, wgt=1, LStype=None):
    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    X = Hdata[:,1:6].transpose() # transpose data matrix
    # Y variable: reference radiance + adjustment values
    Y = Hdata[:,0] + Hdata[:,6]
    # sigma2 errors in X and Y variables
    VX = (Hrnd[:,1:6]**2).transpose()
    VY = Hrnd[:,0]**2 + Hrnd[:,6]**2 # assume ref.rad and K independence 
    
    fb = [1, 1, 1, 0] # fix a3 to input value
    if LStype: # run ordinary least squares
        if wgt: # weighted
            print 'Running WLS for a pair reference-sensor; a3 fixed'
            fit = odr.odr(fcnP,b0,Y,X,we=1./VY,ifixb=fb,job=00002,full_output=1)
        else: # unweighted
            print 'Running OLS for a pair reference-sensor; a3 fixed'
            fit = odr.odr(fcnP,b0,Y,X,ifixb=fb,job=00002,full_output=1)
    else: # run orthogonal least squares
        if wgt: # weighted
            print 'Running weighted ODR for a pair reference-sensor; a3 fixed'
            fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
        else: # unweighted
            print 'Running ODR for a pair reference-sensor; a3 fixed'
            fit = odr.odr(fcnP,b0,Y,X,ifixb=fb,full_output=1)
    odrFit = odr.Output(fit)
    
    return odrFit # return odr output
    
""" Perform LS fit for a sensor-reference pair with ODR function """
def pairLS(Hdata, Hrnd, beta0, wgt=1, LStype=None):
    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    Xvars = Hdata[:,1:6].transpose() # transpose data matrix
    # Y variable: reference radiance + adjustment values
    Yvar = Hdata[:,0] + Hdata[:,6]
    # sigma2 errors in X and Y variables
    VX = (Hrnd[:,1:6]**2).transpose()
    VY = Hrnd[:,0]**2 + Hrnd[:,6]**2 # assume ref.rad and K independence 

    # data and model for (weighted) ODR and OLS
    if wgt: # weighted
        odrData = odr.Data(Xvars, Yvar, wd=1./VX, we = 1./VY)
        fitT = 'weighted'
    else: # unweighted
        odrData = odr.Data(Xvars, Yvar)
        fitT = ''
    odrMod = odr.Model(fcnP)
    
    # set data, model, LS type, and print output parameters for ODR
    odrF = odr.ODR(odrData, odrMod, beta0, ifixb=[1,1,1,0])
    #odrFit = odr.ODR(odrData, odrMod, beta0, ifixb=[1,1,1,0], rptfile='pODRes')
    if LStype:
        odrF.set_job(fit_type=2) # fit for ordinary/weighted least-squares
        print 'Running', fitT, 'OLS for a pair; a3 fixed to input value'
    else:
        print 'Running', fitT, 'ODR for a pair; a3 fixed to input value'    
    #odrFit.set_iprint(init=1, iter=1, iter_step=5, final=1)
    
    mFit = odrF.run() # run ODR fit
    return mFit # return ODR output


""" AVHRR measurement model to use for series harmonisation """
def fcnH(coef, data):
    # read data to variable names; transpose ndarrays ??
    Lr1 = data[0,:] # reference radiance 1st sensor; if reference-sensor pair
    Cs1 = data[1,:] # space counts 1st sensor
    Cict1 = data[2,:] # ICT counts 1st sensor
    CE1 = data[3,:] # Earth counts 1st sensor
    Lict1 = data[4,:] # ICT radiance 1st sensor
    To1 = data[5,:] # orbit temperature 1st sensor
    Cs2 = data[6,:] # space counts 2nd sensor
    Cict2 = data[7,:] # ICT counts 2nd sensor
    CE2 = data[8,:] # Earth counts 2nd sensor
    Lict2 = data[9,:] # ICT radiance 2nd sensor
    To2 = data[10,:] # orbit temperature 2nd sensor
    s1 = data[11,:].astype(int) # index of 1st sensor in coefficients array
    s2 = data[12,:].astype(int) # index 2nd sensor
    
    """ returns questionable results; RANK DEFICIENCY """
    a01 = coef[s1+0] # fit coefficients 1st sensor
    a11 = coef[s1+1]
    a21 = coef[s1+2]
    a31 = coef[s1+3]
    a02 = coef[s2+0] # fit coefficients 2nd sensor
    a12 = coef[s2+1]
    a22 = coef[s2+2]
    a32 = coef[s2+3]
    # fit model 
    K = avhrrME(CE2,Cs2,Cict2,Lict2,To2,a02,a12,a22,a32) -\
        avhrrME(CE1,Cs1,Cict1,Lict1,To1,a01,a11,a21,a31) - Lr1
    print 'Current iteration coefficients:', \
        [a01,a11,a21,a31,a02,a12,a22,a32]    
    return K  

""" Read harmonisation data and perform ODR fit for the whole series """
def seriesODR(Hdata, Hrnd, b0, nos):
    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    selector = [x for x in range(Hdata.shape[1]) if x != 11]
    Xvars = Hdata[:,selector].transpose() # transpose data matrix
    Yvar = Hdata[:,11] # adjustment values

    # data and model for (weighted) ODR 
    odrData = odr.Data(Xvars, Yvar)
    odrMod = odr.Model(fcnH)

    # set free and fixed coefficients in fitting
    parfix = np.zeros(nos*4, dtype=np.int)
    for sidx in range(1,nos):
        parfix[4*sidx:4*sidx+3] = 1
    fixb = parfix.tolist() # values for ifixb ODR parameter

    odrF = odr.ODR(odrData, odrMod, b0, ifixb=fixb, maxit=25)
    print 'Running ODR for multiple pairs'    
    mFit = odrF.run()   
     
    return mFit # return ODR output
