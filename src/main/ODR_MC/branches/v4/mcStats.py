import numpy as np
import matplotlib.pyplot as plt

def mcStats(mcb0, beta0, sb0):

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
        no_bin = 20 # number of bins for the histograms
        plt.figure()
        plt.hist(mcb0[:,i], no_bin, fc='blue', alpha=0.6)
        plt.axvline(x=beta0[i], color='g', linestyle='solid', linewidth=4)
        plt.axvline(x=beta0[i]-sb0[i], color='g', ls='dashed', lw=2)
        plt.axvline(x=beta0[i]+sb0[i], color='g', ls='dashed', lw=2)
        title = 'Histogram of a[' + str(i) + '] coefficient from MC runs'
        plt.title(title)
        plt.xlabel('a values')
        plt.ylabel('Frequency')
        plt.show()
    
    return mcmb0, mcsdb0, mccovb0, mccorb0 # return coeffs stats from MC eval

# graph to show bias of fitted radiance to input with error bars
def LbiasU(inL, calL, ucL):
    """ Create errorbar graph for radiance bias. Arguments to the function 
    - inL: radiance evaluated from input coefficients
    - calL: radiance evaluated from ODR fitted coefficients 
    - ucL:  radiance uncertainty from cal. coeffs uncertainty """
    
    Lbias = calL - inL # bias: fitted minus input radiance
    sigma = 2 * ucL # 2 sigma uncertainy
    
    plt.figure()
    
    plt.errorbar(inL, Lbias, yerr=sigma)
    ttl = 'Radiance bias with' + r'$2*\sigma$' + 'uncertainty'
    plt.title(ttl)
    plt.xlabel('Radiance (mW/m2/sr/cm-1)')
    plt.ylabel('Radiance bias (mW/m2/sr/cm-1)')
    
    plt.show()
    
    return 0
