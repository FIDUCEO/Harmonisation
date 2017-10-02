""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 06-12-2016
    Last update: 22-01-2017
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
    K = data[5,:] # spectral adjustemnt
    SeL = data[6,:] # systematic error Lict
    SeT = data[7,:] # systematic error To
    
    LE = MEfunc[slabel](CE,Cs,Cict,Lict-SeL,To-SeT,a0,a1,a2,a3) - K
    #print 'Current iteration coefficients:', [a0,a1,a2, a3]
    return LE # return adjusted Earth radiance

""" Perform LS fit for a sensor-reference pair with low-level odr function """
def odrP(Hdata, Hrnd, Hsys, b0, wgt=1, LStype=None):
    # X variables: Cs, Cict, CE, Lict, To & K from Hdata; Lict, To from Hsys
    X = np.concatenate((Hdata[:,1:7], Hsys[:,4:6]), axis=1) 
    Y = Hdata[:,0] # reference radiance 

    nor = Y.shape[0] # number of rows
    # sigma2 errors in X variables; random error 1 for systematic corrections
    VX = np.concatenate((Hrnd[:,1:7]**2, np.ones((nor,2))), axis=1)
    VX = VX.transpose() # transpose matrix
    VY = Hrnd[:,0]**2 # squared random error reference radiance
    
    fb = [1, 1, 1, 0] # fix a3 to input value
    if LStype: # run ordinary least squares
        if wgt: # weighted
            print 'Running WLS for a pair reference-sensor; a3 fixed'
            fit = odr.odr(fcnP,b0,Y,X.T,we=1./VY,ifixb=fb,job=00002,full_output=1)
        else: # unweighted
            print 'Running OLS for a pair reference-sensor; a3 fixed'
            fit = odr.odr(fcnP,b0,Y,X.T,ifixb=fb,job=00002,full_output=1)
    else: # run orthogonal least squares
        if wgt: # weighted ODR
            print 'Running weighted ODR for a pair reference-sensor; a3 fixed'
            fit = odr.odr(fcnP,b0,Y,X.T,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
        else: # unweighted
            print 'Running ODR for a pair reference-sensor; a3 fixed'
            fit = odr.odr(fcnP,b0,Y,X.T,ifixb=fb,full_output=1)

    odrFit = odr.Output(fit)
    return odrFit # return odr output


""" Perform LS fit for a sensor-reference pair with ODR function """
def pairLS(Hdata, Hrnd, Hsys, beta0, wgt=1, LStype=None):
    # X variables: Cs, Cict, CE, Lict, To & K from Hdata; Lict, To from Hsys
    X = np.concatenate((Hdata[:,1:7], Hsys[:,4:6]), axis=1) 
    Y = Hdata[:,0] # reference radiance 

    nor = Y.shape[0] # number of matchups/rows
    # sigma2 errors in X variables; random error 1 for systematic corrections
    VX = np.concatenate((Hrnd[:,1:7]**2, np.ones((nor,2))), axis=1)
    VX = VX.transpose() # transpose matrix
    VY = Hrnd[:,0]**2 # squared random error reference radiance

    # data and model for (weighted) ODR and OLS
    if wgt: # weighted
        odrData = odr.Data(X.T, Y, wd=1./VX, we = 1./VY)
        fitT = 'weighted'
    else: # unweighted
        odrData = odr.Data(X.T, Y)
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


""" AVHRR measurement model to use for series harmonisation: two virtual sensors 
for the data matrices, a block a rows has the specific sensors. """
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
    Lref = avhrrME(CE2,Cs2,Cict2,Lict2-SeL2,To2-SeT2,a02,a12,a22,a32) -\
        (1-switch)* avhrrME(CE1,Cs1,Cict1,Lict1-SeL1,To1-SeT1,a01,a11,a21,a31)\
        - K
    print 'Current iteration coefficients:', \
    np.unique([a01,a11,a21,a31,a02,a12,a22,a32])
    
    return Lref  

""" Read harmonisation data and perform ODR fit for the whole series """
def seriesODR(Hdata, Hrnd, Hsys, sensors, b0, fb):
    # extract all variables but K from Hdata matrix
    selector = [4,5,9,10] # columns with systematic errors for Lict and To
    X = np.concatenate((Hdata[:,1:12], Hsys[:,selector]), axis=1) 
    Y = Hdata[:,0] # reference radiance
    
    nor = Y.shape[0] # number of rows
    # sigma2 errors in X variables; random error 1 for systematic corrections
    VX = np.concatenate((Hrnd[:,1:12]**2, np.ones((nor,4))), axis=1)
    VX = VX.transpose() # transpose matrix
    VY = Hrnd[:,0]**2 # reference radiance squared random error 
    
    # run low-level odr; no weights 
    #fit = odr.odr(fcnH,b0,Y,X.T,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
    print 'Running ODR for multiple pairs'    
    fit = odr.odr(fcnH,b0,Y,X.T,ifixb=fb,maxit=20,full_output=1)
    mFit = odr.Output(fit)

    # prepare data, model and run ODR 
    #odrData = odr.Data(X, Y, wd=1./VX, we = 1./VY)
    #odrMod = odr.Model(fcnH)
    #odrF = odr.ODR(odrData, odrMod, b0, ifixb=fb, maxit=20)
    #print 'Running ODR for multiple pairs'    
    #mFit = odrF.run()   
    
    return mFit # return ODR output

""" AVHRR measurement model to use for series harmonisation: each sensor has  
its separate columns in data matrices, filled with default values if a block of 
rows has matchups from (2) other sensors. """
def fcnHL(coef, Xdata):
    nos = coef.shape[0] / 4 # number of sensors to calibrate
    

""" Read harmonisation data and perform ODR fit for the whole series """
def seriesODRL(Hdata, Hrnd, b0, fb):
    # extract variables Cs, Cict, CE, Lict, To of each sensor from Hdata matrix
    selector = [x for x in range(Hdata.shape[1]) if x > 1]
    X = Hdata[:,selector].transpose() # transpose data matrix
    Y = Hdata[:,1] + Hdata[:,0] # adjustment values + reference radiance

    # sigma2 errors in X and Y variables
    VX = (Hrnd[:,selector]**2).transpose()
    VY = Hrnd[:,1]**2 + Hrnd[:,0]**2 # combined K and Lref variance 
    
    print 'Running ODR for multiple pairs'    
    fit = odr.odr(fcnHL,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)

    mFit = odr.Output(fit)
    return mFit # return ODR output
