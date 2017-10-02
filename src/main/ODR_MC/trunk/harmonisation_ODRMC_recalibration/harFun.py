""" FCDR sensor recalibration modules
    Project:        H2020 FIDUCEO
    Author:         Arta Dilo \NPL MM
    Reviewer:       Peter Harris \NPL MM, Jon Mittaz, Sam Hunt \NPL ECO
    Date created:   06-12-2016
    Last update:    14-03-2017
    Version:        10.0
Harmonisation functions for a pair-wise implementation and for all the sensors 
together using odr package. Functions implement weighted ODR (an EIV method)
for a pair sensor-reference and for multiple pairs of type sensor-reference and 
sensor-sensor. """

import scipy.odr as odr

# AVHRR measurement equation
def avhrrME(X, a, notime):
        
    a0 = a[0] # AVHRR model coefficients
    a1 = a[1]
    a2 = a[2]
    Cs = X[0,:] # space counts 
    Cict = X[1,:] # ICT counts 
    CE = X[2,:] # Earth counts
    Lict = X[3,:] # ICT radiance
    
    # Earth radiance from Earth counts and calibration data
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) 
    
    if not notime:
        a3 = a[3]    
        To = X[4,:] # orbit temperature
        LE += a3*To # Earth radiance
        
    return LE # return Earth radiance

# dictionary with measurement equation function for each sensors' series 
MEfunc = {'avhrr': avhrrME}


""" Perform LS fit for a sensor-reference pair with low-level odr function """
def odrP(Hdata, Hr, b0, series, fb=None, fx=None, Hs=None, rsp=1):
    
    if series.notime: # work with no-time dependant data
        selector = [x for x in range(Hdata.shape[1]) if (x>0 and x<5)]
    else:
        selector = [x for x in range(Hdata.shape[1]) if (x>0 and x<6)]
        
    # extract variables Cs, Cict, CE, Lict (and To) from Hdata matrix
    X = Hdata[:,selector].transpose() # transpose data matrix
    # Y is adjusted radiance: reference radiance + adjustment values
    Y = Hdata[:,0] + Hdata[:,6]

    # cacluate weights from uncertainty matrices
    if Hs is not None: # weight on both random and systematic uncertainty data
        #Hs = resetHs(Hs, rsp) # set sytematic equiv to Peter optimisation prob
        VX = (Hr[:,selector]**2 + Hs[:,selector]**2).transpose() # sigma^2 of X variables
        
        ''' Y = Lref+K: assume independence of ref. radiance and K
        K random: in Hr matchups uncertainty, in Hs SRF shifting uncertainty '''
        VY = Hr[:,0]**2 + (Hr[:,6]**2+Hs[:,6]**2) # sigma^2 of Y
        
    else: # weight on random uncertainty 
        VX = (Hr[:,selector]**2).transpose() # sigma^2 of X
        VY = Hr[:,0]**2 + Hr[:,6]**2 # Y sigma^2 (no shifting uncert.)
    
    # ODR model
    def fcnP(coef, Xdata, ss=series):
        slabel = ss.slabel # series label to pick measurement model
        LE = MEfunc[slabel](Xdata, coef, ss.notime)
        return LE # return Earth radiance

    # perform odr fit (low level function)
    if fb: # keep a3 coefficient fixed (defined by fb) and To var fixed (by fx)
        fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,ifixx=fx,full_output=1)
    else: # fit all coefficients 
        fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,full_output=1)

    odrFit = odr.Output(fit) # get odr fit output     
    return odrFit # return odr output


""" Perform ODR over MC generated data with ODR best estimates from 
real or simulated data and errors """
def odr4MC(Xdata, Ydata, Hr, b0, series, fb=None, fx=None, Hs=None, rsp=1):
    X = Xdata.transpose()
    if series.notime: # work with no-time dependant data
        selector = [x for x in range(Hr.shape[1]) if (x>0 and x<5)]
    else:
        selector = [x for x in range(Hr.shape[1]) if (x>0 and x<6)]

    # cacluate weights from uncertainty matrices
    if Hs is not None: # weights from combined random & systematic uncertainty
        VX = (Hr[:,selector]**2 + Hs[:,selector]**2).transpose() # sigma^2 of X variables
        
        ''' Y = Lref+K: assume independence of ref. radiance and K
        K uncert: in Hr matchups uncertainty, in Hs SRF shifting uncertainty '''
        VY = Hr[:,0]**2 + (Hr[:,6]**2+Hs[:,6]**2) # sigma^2 of Y
        
    else: # weight on random uncertainty 
        VX = (Hr[:,selector]**2).transpose() # sigma^2 of X
        VY = Hr[:,0]**2 + Hr[:,6]**2 # Y sigma^2 (no shifting uncert.)
        
    # ODR model
    def fcnP(coef, Xdata, ss=series):
        slabel = ss.slabel # series label to pick measurement model
        LE = MEfunc[slabel](Xdata, coef, ss.notime)
        return LE # return Earth radiance

    #  ODR on new X,Y data, perturbed best estimates  
    if fb: # keep a3 coefficient fixed (defined by fb) and To var fixed (by fx)
        fit = odr.odr(fcnP,b0,Ydata,X,we=1./VY,wd=1./VX,ifixb=fb,ifixx=fx,full_output=1)
    else: # fit all coefficients 
        fit = odr.odr(fcnP,b0,Ydata,X,we=1./VY,wd=1./VX,full_output=1)

    odrFit = odr.Output(fit) # get odr fit output     
    return odrFit # return odr output
