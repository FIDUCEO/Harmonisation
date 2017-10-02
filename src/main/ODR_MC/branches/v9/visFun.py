""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo NPL\MM and Sam Hunt NPL\ENV
    Date created: 07-02-2016
    Last update: 14-03-2017
Functions for visualisation. """

from numpy import mean, cov, corrcoef, sqrt, zeros
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from numpy import ndenumerate


""" Compute correlation from covariance matrix """ 
def cov2cor(cov):
    dim = cov.shape[0] # dimention of squared cov matrix
    cor = zeros((dim,dim))
    
    # fill the upper triangle of the corr matrix
    for i in range(dim):
        for j in range(i, dim):
            if cov[i,i]>0 and cov[j,j]>0:
                cor[i, j] = cov[i, j] / sqrt(cov[i,i] * cov[j,j])
            elif cov[i, j] == 0:
                cor[i, j] = 0.

    # fill the lower triangle with the value of opposite to diagonal cell          
    for i in range(1, dim):
        for j in range(i):
            cor[i, j] = cor[j, i]
            
    return cor # return correlation matrix

def plot_corr_heatmap(A, title=None, labels=None, save_path=None):
    """
    Plot heatmap of input array

    :param A: numpy.ndarray
        input array to plot heatmap
    :param save_path: str
        path to save chart image to
    """

    fig1, (ax1) = plt.subplots(1)
    ax1.imshow(A, cmap="bwr", interpolation='nearest', vmin=-1, vmax=1)
    ax1.set_title(title)
    for (j, i), label in ndenumerate(A):
        label = round(label, 2)
        ax1.text(i, j, label, ha='center', va='center')

    ax1.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='on', right='off', left='off',
                    labelleft='on')

    tick_labels = []
    for l in labels:
        tick_labels.append("")
        tick_labels.append(l)
    ax1.set_xticklabels(tick_labels)
    ax1.set_yticklabels(tick_labels)

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)

    return 0

""" Plot heat map of correlations between fit coefficients """
def corrMap(slab, corr, tlab):
    
    print '\nCorrelations of fit coefficients from', tlab
    print corr
    # plot map of correlation of ODR fit and MC coefficients 
    pttl = slab + ' correlation of harmonisation coefficients in ' + tlab
    plot_corr_heatmap(corr, title=pttl, labels=['a0', 'a1', 'a2', 'a3'])


""" Calculates statitics of fit coefficients in MC trials: mean, covariance, 
correlation matrix; plots histograms of coefficients. """
def mcStats(slab, mcb, beta, covb):
    # sample mean of MC trials coeffs
    mcmb0 = mean(mcb, axis=0)     
    # sample covariance of MC trials calibration coefficients
    mccovb0 = cov(mcb.T)
    # Sample correlations of calib. coefficients from MC trials
    mccorb0 = corrcoef(mcb.T)
    stdb = zeros(beta.shape[0])
    
    for i in range(3):
        stdb[i] = sqrt(covb[i ,i])
        # Plot histogram of beta values from MC runs         
        plt.figure()
        n,bins,patches = plt.hist(mcb[:,i], 40, normed=1, fc='blue', alpha=0.6)
        # add a 'best fit' Gaussian with mu and sigma evaluated by ODR
        y = mlab.normpdf(bins, beta[i], stdb[i])
        plt.plot(bins, y, 'r--', linewidth=1.5)        
        title = 'Histogram of '+slab+' a[' +str(i)+ '] coeff from MC runs, Gaussian w '\
                +r'$\mu, \sigma$' +' from ODR eval'
        plt.title(title)
        plt.xlabel('a coeff values')
        plt.ylabel('Probability')
        plt.grid(True)
        # show graph in plot window
        plt.show()
        
        # plot histogram with ODR beta and sigma as vertical lines
        plt.figure()
        plt.hist(mcb[:,i], bins, fc='blue', alpha=0.6)
        plt.axvline(x=beta[i], color='g', linestyle='solid', linewidth=4)
        plt.axvline(x=beta[i]-stdb[i], color='g', ls='dashed', lw=2)
        plt.axvline(x=beta[i]+stdb[i], color='g', ls='dashed', lw=2)
        title = 'Histogram of '+slab+' a[' +str(i)+ '] coefficient from MC runs with ' \
                +r'$\mu, \sigma$' +' from ODR eval'
        plt.title(title)
        plt.xlabel('a values')
        plt.ylabel('Frequency')
        plt.show()
    
    return mcmb0, mccovb0, mccorb0 # return coeffs stats from MC eval


""" Create errorbar graph for radiance bias. Arguments to the function 
    - inL: radiance evaluated from input coefficients
    - calL: radiance evaluated from ODR fitted coefficients 
    - uL: radiance uncertainty from calib. coefficients (& data uncertainty) """
def LbiasU(inL, calL, uL, k, ttl):
    
    Lbias = calL - inL # bias: fitted minus input radiance
    sigma = k * uL # k sigma uncertainy
    
    plt.figure() # open graphical window for plotting
    
    plt.errorbar(inL, Lbias, yerr=sigma, fmt='o', color='green')
    plt.title(ttl)
    plt.xlabel('Radiance')
    plt.ylabel('Radiance bias')
    
    plt.show()    
    return 0


""" Show the difference between uncertainty evaluation from ODR and MC. 
    Arguments to the function 
    - inL: radiance evaluated from input coefficients
    - odrLU: radiance uncertainty evaluated by ODR covariance matrix
    - mcLU: radiance uncertainty evaluated by MC covariance matrix """
def LUdiff(inL, odrLU, mcLU, ttl):
    
    plt.figure()
    
    LUbias = mcLU - odrLU # bias: MC eval minus ODR eval 
    
    plt.scatter(inL, LUbias, s=35, color='maroon')
    plt.title(ttl)
    plt.ylim(min(LUbias), max(LUbias))
    plt.xlabel('Earth radiance')
    plt.ylabel('Radiance uncertainty')
    
    plt.show()    
    return 0


""" Plot radiance fit residuals against different varaibles in the x-axis. 
    Arguments to the function 
    - resid: ODR evaluated true error of Y and X variables
    - xvar: covariate to plot residuals against: radiance, Lict, To, OR time
    - xlab: label for the x-axis with units (add)
    - ylab: label for the y-axis with units (add)
    - ttl: title of the graph """
def resLfit(resid, xvar, xlab, ylab, ttl):
    
    plt.figure()
        
    plt.scatter(xvar, resid, s=15, color='brown')
    plt.title(ttl)
    #plt.ylim(min(resid), max(resid))
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    
    plt.show()    
    return 0
