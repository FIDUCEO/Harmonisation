""" Generate data for Monte Carlo simulation of beta coefficients 
uncertainty evaluation. """

import numpy as np

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

''' Generate the matrix of errors respecting the correlation structure '''
def genpCS(Hr,Lsys,Tsys,Crnd,mawlen):
    err = np.empty(Hr.shape) # matrix of errors
    nor = err.shape[0] # number of matchups
    
    # draw a value of systematic error for Lict and To
    errL = np.random.normal(scale=Lsys) # Lict systematic error
    errT = np.random.normal(scale=Tsys) # To systematic error
    
    # constant random uncertainty for moving average on calib. counts; temporary
    CCrma = Crnd/np.sqrt(Crnd)
    
    # create error matrix for the whole dataset
    for i in range(nor):
        # Lref, K, CE: random error from Gaussian with sigma from Hr data &mu=0
        err[i, 0] = np.random.normal(scale=Hr[i,0]) # Lref random error
        err[i, 3] = np.random.normal(scale=Hr[i,3]) # CE random error
        err[i, 6] = np.random.normal(scale=Hr[i,6]) # K random error
        
        # Lict, To: Gaussian random with sigma from Hr & mu=0, + systematic err
        err[i, 4] = np.random.normal(scale=Hr[i,3]) + errL # Lict error
        err[i, 5] = np.random.normal(scale=Hr[i,3]) + errT # To error
        
        # Cs, Cict: Gaussian random with correlations NOT YET finished
        err[i, 1] = np.random.normal(scale=CCrma) # Cs error
        err[i, 2] = np.random.normal(scale=CCrma) # Cict error
        
    return err
