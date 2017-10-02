""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 06-12-2016
    Last update: 14-02-2017
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

""" Reset the systematic error matrix to constant values such that the ODR 
problem setting corresponds to Peter's LS optimisation problem. 
- Lref, Cspace, Cict, CEarth have 0 systematic uncertainty
- Lict and To have a constant systematic uncertainty 
- K random uncertainty from SRF shifting is stored in the Hs matrix """
def resetHs(Hs, rspair):
    Hsys = np.zeros(Hs.shape)
    nor = Hs.shape[0]

    # set Lict and To systematic to mean of corresponding Hs column
    # columns [4] and [5] are respectively Lict and To of series sensor in 
    # ref-sensor pair, and Lict and To of the 1st sensor in sensor-sensor pair
    sLict = np.mean(Hs[:,4]) # mean Lict through all matchups
    sTo = np.mean(Hs[:,5])  # mean To 
    Hsys[:,4] = sLict * np.ones(nor)
    Hsys[:,5] = sTo  * np.ones(nor)
        
    if rspair: # reference sensor pair, i.e. Hs has 7 columns     
        # keep K random uncertainty of SRF shifting   
        Hsys[:,6] = Hs[:,6] 

    else: # sensor-sensor pair, Hs has 12 columns
        Hsys[:,11] = Hs[:,11] # keep K random uncertainty from SRF shifting
        
        # set 2nd sensor' Lict and To systematic to mean of Hs values
        sLict = np.mean(Hs[:,9]) 
        sTo = np.mean(Hs[:,10]) 
        Hsys[:,9] = sLict * np.ones(nor)
        Hsys[:,10] = sTo * np.ones(nor)
       
    return Hsys # return the new set of sytematic uncertainties

""" Model function, i.e. fcn argument in ODR package """
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
def odrP(Hdata, Hr, b0, fb=None, Hs=None, rsp=1):

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
    if fb: # keep some coefficients fixed, defined by fb
        fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
    else: # fit all coefficients 
        fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,full_output=1)

    odrFit = odr.Output(fit) # get odr fit output     
    return odrFit # return odr output

""" Perform ODR over data MC generated with full error structure """
def odr4MC(Xdata, Ydata, Hr, b0, fb=None, Hs=None, rsp=1):
    X = Xdata.transpose()

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
        
    #  ODR on new X,Y data, perturbed best estimates  
    if fb: # keep some coefficients fixed, defined by fb
        fit = odr.odr(fcnP,b0,Ydata,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
    else: # fit all coefficients 
        fit = odr.odr(fcnP,b0,Ydata,X,we=1./VY,wd=1./VX,full_output=1)

    odrFit = odr.Output(fit) # get odr fit output     
    return odrFit # return odr output


""" Perform ODR fit for the whole series (gives runtime error). 
AVHRR measurement model to use for series harmonisation: two virtual sensors 
for the data matrices, a block a rows has the specific sensors. """
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
        #print 'Current iteration coefficients:', \
        #np.unique([a01,a11,a21,a31,a02,a12,a22,a32])
        
        return K  

    # run low-level odr
    print '\nRunning ODR for multiple pairs'    
    fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
    #fit = odr.odr(fcnH,b0,Y,X,ifixb=fb,maxit=20)
    mFit = odr.Output(fit)

    # prepare data, model and run ODR 
    #odrData = odr.Data(X, Y, wd=1./VX, we = 1./VY)
    #odrMod = odr.Model(fcnH)
    #odrF = odr.ODR(odrData, odrMod, b0, ifixb=fb, maxit=20)
    #print 'Running ODR for multiple pairs'    
    #mFit = odrF.run()   
    
    return mFit # return ODR output
