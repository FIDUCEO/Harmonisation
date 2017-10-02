""" Function to print results of different fit methods & functions. """

import matplotlib.pyplot as plt
import numpy as np

def mcStats(mcb0, b0, sb0, cb0):

    # sample mean of MC trials coeffs
    mcmb0 = np.mean(mcb0.T, axis=0)     
    # sample stdev of MC trials calibration coefficients
    mcsdb0 = np.std(mcb0.T, axis=0)    
    # sample covariance of MC trials calibration coefficients
    mccovb0 = np.cov(mcb0.T)
    # Sample correlations of calib. coefficients from MC trials
    mccorb0 = np.corrcoef(mcb0.T)
    
    for i in range(3):
        # Plot histogram of beta values from MC runs 
        no_bin = 20 # number of bins for the histograms
        plt.figure()
        plt.hist(mcb0[:,i], no_bin, fc='blue', alpha=0.6)
        plt.axvline(b0[i], color='g', linestyle='solid', linewidth=4)
        plt.axvline(b0[i]-sb0[i], color='g', ls='dashed', lw=2)
        plt.axvline(b0[i]+sb0[i], color='g', ls='dashed', lw=2)
        title = 'Histogram of a[' + str(i) + '] coefficient from MC runs'
        plt.title(title)
        plt.xlabel('a values')
        plt.ylabel('Frequency')
        plt.show()

        ## plot distribution of beta values
        #plt.figure() 
        #plt.plot(mcb0[:,i], 'go')#, label='Measured')
        #plt.title('Spread of a2 coefficients in MC runs')
        #plt.xlabel('MC iterations')
        #plt.ylabel('a2 coefficient')
        ## plot beta0 +/- standard deviation
        #sbnd = plt.axhspan(b0[i] - sb0[i], b0[i] + sb0[i], fc='c', alpha=0.8)
        #sm = plt.axhline(b0[i], linewidth=1, color='b') # plot beta0 
        #plt.show()
    
    return mcmb0, mcsdb0, mccovb0, mccorb0 # return coeffs covariance from MC eval

def plotLaU(CE, inL, odrL, LUodr, LUmc):
    
    plt.figure()    
    lodr = plt.errorbar(CE, odrL, yerr=LUodr, fmt='o', color='green')
    lin = plt.scatter(CE, inL, color='red')
    plt.title('Offset of odr fitted to input radiance: sigma range from odr uncertainty eval.')
    # add legend
    plt.show()

    plt.figure()    
    lodr = plt.errorbar(CE, odrL, yerr=LUmc, fmt='o', color='green')
    lin = plt.scatter(CE, inL, color='red')
    plt.title('Offset of odr fitted to input radiance: sigma range from MC uncertainty eval.')    
    plt.show()

def podrRes(fitM):
    #print s2, 'Betta coefficients:', \
    #    np.array_str(fitM.beta, precision=5, suppress_small=True)
    #print 'Betta Covariance:'
    #print np.array_str(fitM.cov_beta, precision=2)
    #print 'Reason for returning:', fitM.stopreason
    print 'Sum of squares error:', fitM.sum_square
    print 'Sum of squares of delta error:', fitM.sum_square_delta
    print 'Sum of squares of eps error:', fitM.sum_square_eps
    print 'Relative error in function values computed within fcn:', \
        fitM.rel_error
    #print 'Indices into work for drawing out values', fitM.work_ind
    print 'Array of estimated errors in input variables:' 
    print fitM.delta
    print 'Array of estimated errors in response variables:'
    print fitM.eps
    print 'Array of estimated true X variables (x + delta):'
    print fitM.xplus
    print 'Array of estimated true Y values fcn(x+delta):'
    print fitM.y
    
    return 0
