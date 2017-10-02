""" Function to print results of different fit methods & functions. """

import matplotlib.pyplot as plt
import numpy as np

def MCwODRstats(mcb0, b0, sb0, cb0):
    # print statistics of fit coefficients from odr and estimation from MC runs 
    print 'Fit coefficients from odr run on data weighted by random uncertainties'
    print b0
    print 'Standard error of fit coefficients evaluated by odr'
    print sb0
    print 'Covariance matrix of fit coefficients evaluated by odr'
    print cb0
    print 'Beta coefficients from MC ODRs on simulated data with full error structure'
    print mcb0
    mcsb0 = np.zeros((4))
    for i in range(4):
        mcsb0[i] = np.std(mcb0[i])
    print 'Standard error of fit coefficients evaluated from statistics on MC runs'
    print mcsb0
    mccb0 = np.cov(mcb0.T)
    print 'Covariance matrix evaluated statistically from fit coefficients in MC runs'
    print mccb0
    
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

        # plot distribution of beta values
        #plt.figure() 
        #plt.plot(mcb0[:,i], 'go')#, label='Measured')
        #plt.title('Spread of a2 coefficients in MC runs')
        #plt.xlabel('MC iterations')
        #plt.ylabel('a2 coefficient')
        ## plot beta0 +/- standard deviation
        #sbnd = plt.axhspan(b0[i] - sb0[i], b0[i] + sb0[i], fc='c', alpha=0.8)
        #sm = plt.axhline(b0[i], linewidth=1, color='b') # plot beta0 
        #plt.show()
    
    return mccb0 # return coeffs covariance from MC eval

def plotLaU(CE, inL, odrL, LUodr, LUmc):
    Lbias = odrL - inL # offset odr fitted to input radiance 
    
    plt.figure()
    
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
    ax = axs[0,0]
    ax.errorbar(CE, Lbias, yerr=LUodr, fmt='o')
    ax.set_title('Offset of odr fitted to input radiance: sigma range from odr uncertainty eval.')
    
    # With 4 subplots, reduce the number of axis ticks to avoid crowding.
    ax.locator_params(nbins=4)
    
    ax = axs[1,0]
    ax.errorbar(CE, Lbias, yerr=LUmc, fmt='o')
    ax.set_title('Offset of odr fitted to input radiance: sigma range from MC uncertainty eval.')
    
    fig.suptitle('Earth radiance uncertainty from coefficients uncertainty')
    
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
