""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 06-12-2016
    Last update: 14-03-2017
Harmonisation functions for a pair-wise implementation and for all the sensors 
together using odr package. Functions implement weighted ODR (an EIV method)
for a pair sensor-reference and for multiple pairs of type sensor-reference and 
sensor-sensor. """

import scipy.odr as odr
from numpy import logical_not

# AVHRR measurement equation
def avhrrME(CE,Cs,Cict,Lict,To,a0,a1,a2,a3):
        
    # Earth radiance from Earth counts and  calibration data
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) + a3*To
    return LE # return Earth radiance

# dictionary with measurement eq. function of each sensors' series 
MEfunc = {'avhrr': avhrrME}

""" Model function for re-calibration, i.e. fcn argument in ODR package """
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
def odrP(Hdata, Hr, b0, fb=None, fx=None, Hs=None, rsp=1):

    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    X = Hdata[:,1:6].transpose() # transpose data matrix
    # Y is adjusted radiance: reference radiance + adjustment values
    Y = Hdata[:,0] + Hdata[:,6]

    # cacluate weights from uncertainty matrices
    if Hs is not None: # weight on both random and systematic uncertainty data
        #Hs = resetHs(Hs, rsp) # set sytematic equiv to Peter optimisation prob
        VX = (Hr[:,1:6]**2 + Hs[:,1:6]**2).transpose() # sigma^2 of X variables
        
        ''' Y = Lref+K: assume independence of ref. radiance and K
        K random: in Hr matchups uncertainty, in Hs SRF shifting uncertainty '''
        VY = Hr[:,0]**2 + (Hr[:,6]**2+Hs[:,6]**2) # sigma^2 of Y
        
    else: # weight on random uncertainty 
        VX = (Hr[:,1:6]**2).transpose() # sigma^2 of X
        VY = Hr[:,0]**2 + Hr[:,6]**2 # Y sigma^2 (no shifting uncert.)
    
    # perform odr fit (low level function)
    if fb: # keep a3 coefficient fixed (defined by fb) and To var fixed (by fx)
        fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,ifixx=fx,full_output=1)
    else: # fit all coefficients 
        fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,full_output=1)

    odrFit = odr.Output(fit) # get odr fit output     
    return odrFit # return odr output

""" Perform ODR over MC generated data with ODR best estimates from 
real or simulated data and errors """
def odr4MC(Xdata, Ydata, Hr, b0, fb=None, fx=None, Hs=None, rsp=1):
    X = Xdata.transpose()

    # cacluate weights from uncertainty matrices
    if Hs is not None: # weights from combined random & systematic uncertainty
        VX = (Hr[:,1:6]**2 + Hs[:,1:6]**2).transpose() # sigma^2 of X variables
        
        ''' Y = Lref+K: assume independence of ref. radiance and K
        K random: in Hr matchups uncertainty, in Hs SRF shifting uncertainty '''
        VY = Hr[:,0]**2 + (Hr[:,6]**2+Hs[:,6]**2) # sigma^2 of Y
        
    else: # weight on random uncertainty 
        VX = (Hr[:,1:6]**2).transpose() # sigma^2 of X
        VY = Hr[:,0]**2 + Hr[:,6]**2 # Y sigma^2 (no shifting uncert.)
        
    #  ODR on new X,Y data, perturbed best estimates  
    if fb: # keep a3 coefficient fixed (defined by fb) and To var fixed (by fx)
        fit = odr.odr(fcnP,b0,Ydata,X,we=1./VY,wd=1./VX,ifixb=fb,ifixx=fx,full_output=1)
    else: # fit all coefficients 
        fit = odr.odr(fcnP,b0,Ydata,X,we=1./VY,wd=1./VX,full_output=1)

    odrFit = odr.Output(fit) # get odr fit output     
    return odrFit # return odr output


""" Model function for series harmonisation (fcn argument in ODR package).
This setup with fcnH outside seriesODR function is not working;  
possibly the reading of sensors array sp. """
def fcnH(coef, Xdata, sensors, series):
    sp = sensors.transpose()

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
    #s1 = sensors[:,0] # 1st sensor index in sensors list (&coeff arr)
    #s2 = sensors[:,1] # 2nd sensor's index 
    switch = logical_not(s1).astype(int) 
    
    p = series.nocoefs # number of calibration coefficients
    a01 = coef[s1*p + 0] # fit coefficients 1st sensor [s*p+0 for s in s1]
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
    #print 'Current iteration coefficients:', [a01,a11,a21,a31,a02,a12,a22,a32]
    
    return K  

""" Perform ODR fit for the whole series. 
AVHRR measurement model to use for series harmonisation: two virtual sensors 
for the data matrices, a block a rows has the specific sensors. """
def seriesODR(Hdata, Hunc2, b0, sensors, series, fb=None, fx=None):
    # extract variables from Hdata matrix
    X = Hdata[:,0:11].transpose() # X vars; transpose data matrix
    Y = Hdata[:,11] # adjustment values K
    
    VX = Hunc2[:,0:11].transpose() # squared uncertainty X vars
    VY = Hunc2[:,11] # K squared uncertainty 
        
    print '\nRunning ODR for multiple pairs\n'    
    # run low-level odr
    if fb: # keep a3 coefficients fixed (fb) and To vars fixed (fx)
        fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,ifixx=fx,full_output=1)
        #fit = odr.odr(fcnH,b0,Y,X,ifixb=fb,ifixx=fx,iprint=3,rptfile='shFinal.rpt',maxit=20)
        
    else: # fit all coefficients 
        fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,full_output=1)
    
    mFit = odr.Output(fit)
    return mFit # return ODR output
