""" FCDR harmonisation modules
    Project:        H2020 FIDUCEO
    Author:         Arta Dilo \NPL MM
    Reviewer:       Jon Mitazz, Peter Harris \NPL MM, Sam Hunt \NPL ECO
    Date created:   06-12-2016
    Last update:    20-03-2017
    Version:        12.0
Harmonisation functions for a pair-wise implementation and for all the sensors 
together using odr package. Functions implement weighted ODR (an EIV method)
for a pair sensor-reference and for multiple pairs of type sensor-reference and 
sensor-sensor. """

import scipy.odr as odr
from numpy import logical_not

# AVHRR measurement equation
def avhrrME(X, a0,a1,a2, notime, a3=0):
        
    Cs = X[0,:] # space counts 
    Cict = X[1,:] # ICT counts 
    CE = X[2,:] # Earth counts
    Lict = X[3,:] # ICT radiance
    
    # Earth radiance from Earth counts and calibration data
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) 
    
    if not notime: # time dependent measurement model
        To = X[4,:] # orbit temperature
        LE += a3*To # add time-dependant component to Earth radiance
        
    return LE # return Earth radiance


# dictionary with measurement equation function for each sensors' series 
MEfunc = {'avhrr': avhrrME}


""" Perform ODR fit for the whole series. 
AVHRR measurement model to use for series harmonisation: two virtual sensors 
for the data matrices, a block a rows has the specific sensors. """
def seriesODR(Xdata,Y,Xrnd,Yrnd,b0,sensors,series,fb=None,fx=None):
    notd = series.notime
    # use slab to choose meas. model: avhrrME -> MEfunc[slab](Xdata, coef, notd)
    slab = series.slabel # series label; not yet used

    
    bsens = sensors.transpose()
    X = Xdata.transpose() # X vars; transpose data matrix
    VX = (Xrnd**2).transpose() # squared uncertainty X vars
    VY = Yrnd**2 # squared Y uncertainty 

    def fcnH(coef, Xdata, sp=bsens):
        # read data to variable names; transpose ndarrays 
        Lr1 = Xdata[0,:] # reference radiance 1st sensor; 0 for sensor-sensor pair
        s1 = sp[0,:] # 1st sensor index in sensors list (&coeff arr)
        s2 = sp[1,:] # 2nd sensor's index 
        switch = logical_not(s1).astype(int) 
        
        m = series.novars # number of measured variables
        Xs1 = Xdata[1:1+m,:] # model variables for the 1st sensor
        Xs2 = Xdata[1+m:1+2*m,:] # model variables for the 2nd sensor

        p = series.nocoefs # number of calibration coefficients        
        a01 = coef[s1*p + 0] # fit coefficients 1st sensor [s*p+0 for s in s1]
        a11 = coef[s1*p + 1]
        a21 = coef[s1*p + 2]
        a02 = coef[s2*p + 0] # fit coefficients 2nd sensor
        a12 = coef[s2*p + 1]
        a22 = coef[s2*p + 2]

        if not notd: # time dependent measurement model
            a31 = coef[s1*p + 3]
            a32 = coef[s2*p + 3]
            # fit model 
            K = avhrrME(Xs2,a02,a12,a22,notd,a32) - \
            (1-switch) * avhrrME(Xs1,a01,a11,a21,notd,a31) - switch * Lr1
            
        else: # not time-dependant model 
            K = avhrrME(Xs2,a02,a12,a22,notd) - \
                (1-switch) * avhrrME(Xs1,a01,a11,a21,notd) - switch * Lr1
        
        return K  
    
    # run low-level odr
    if fb is not None: # keep a3 coefficients fixed (fb) and To vars fixed (fx)
        if fx is not None: 
            fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,ifixx=fx,full_output=1)
        else: # fix reference sensor cal.coeffs for non time-dependent model
            fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,ifixb=fb,full_output=1)
    else: # fit all coefficients 
        fit = odr.odr(fcnH,b0,Y,X,we=1./VY,wd=1./VX,full_output=1)
    
    mFit = odr.Output(fit)
    return mFit # return ODR output
