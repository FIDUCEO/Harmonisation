""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 01-02-2017
    Last update: 21-02-2017
Create graphs to display results of harmonisation uncertainty via MC. """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def mcStats(slab, mcb0, beta0, sb0):

    # sample mean of MC trials coeffs
    mcmb0 = np.mean(mcb0, axis=0)     
    # sample stdev of MC trials calibration coefficients
    mcsdb0 = np.std(mcb0, axis=0)    
    # sample covariance of MC trials calibration coefficients
    mccovb0 = np.cov(mcb0.T)
    # Sample correlations of calib. coefficients from MC trials
    mccorb0 = np.corrcoef(mcb0.T)
    
    for i in range(3):
        # Plot histogram of beta values from MC runs         
        plt.figure()

        n,bins,patches = plt.hist(mcb0[:,i], 40, normed=1, fc='blue', alpha=0.6)
        # add a 'best fit' Gaussian with mu and sigma evaluated by ODR
        y = mlab.normpdf(bins, beta0[i], sb0[i])
        l = plt.plot(bins, y, 'r--', linewidth=1.5)        
        title = 'Histogram of '+slab+' a[' +str(i)+ '] coeff from MC runs, Gaussian w ' +r'$\mu, \sigma$' +' from ODR eval'
        plt.title(title)
        plt.xlabel('a coeff values')
        plt.ylabel('Probability')
        plt.grid(True)

        plt.show()
        

        # plot histogram with ODR beta and sigma as vertical lines
        plt.figure()
        
        plt.hist(mcb0[:,i], bins, fc='blue', alpha=0.6)
        plt.axvline(x=beta0[i], color='g', linestyle='solid', linewidth=4)
        plt.axvline(x=beta0[i]-sb0[i], color='g', ls='dashed', lw=2)
        plt.axvline(x=beta0[i]+sb0[i], color='g', ls='dashed', lw=2)
        title = 'Histogram of '+slab+' a[' +str(i)+ '] coefficient from MC runs with ' +r'$\mu, \sigma$' +' from ODR eval'
        plt.title(title)
        plt.xlabel('a values')
        plt.ylabel('Frequency')

        plt.show()
    
    return mcmb0, mcsdb0, mccovb0, mccorb0 # return coeffs stats from MC eval


# graph to show bias of fitted radiance to input with error bars
def LbiasU(inL, calL, uL, k, ttl):
    """ Create errorbar graph for radiance bias. Arguments to the function 
    - inL: radiance evaluated from input coefficients
    - calL: radiance evaluated from ODR fitted coefficients 
    - uL: radiance uncertainty from cal. coeffs/coeffs and data uncertainty """
    
    Lbias = calL - inL # bias: fitted minus input radiance
    sigma = k * uL # k sigma uncertainy
    
    plt.figure()
    
    plt.errorbar(inL, Lbias, yerr=sigma, fmt='o', color='green')
    plt.title(ttl)
    plt.xlabel('Radiance (mW/m2/sr/cm-1)')
    plt.ylabel('Radiance bias (mW/m2/sr/cm-1)')
    
    plt.show()    
    return 0


# graph of differences between radiance uncertainties evaluated by ODR and MC
def LUdiff(inL, odrLU, mcLU, ttl):
    """ Show the difference between uncertainty evaluation from ODR and MC. 
    Arguments to the function 
    - inL: radiance evaluated from input coefficients
    - odrLU: radiance uncertainty evaluated by ODR covariance matrix
    - mcLU: radiance uncertainty evaluated by MC covariance matrix """
    
    plt.figure()
    
    LUbias = mcLU - odrLU # bias: MC eval minus ODR eval 
    
    plt.scatter(inL, LUbias, s=35, color='maroon')
    plt.title(ttl)
    plt.ylim(min(LUbias), max(LUbias))
    plt.xlabel('Radiance (mW/m2/sr/cm-1)')
    plt.ylabel('Radiance uncertainty (mW/m2/sr/cm-1)')
    
    plt.show()    
    return 0

