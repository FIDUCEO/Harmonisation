#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author:       Arta Dilo \ NPL MM
    Date created: 12-12-2016
    Last update:  20-03-2017
    Version:      10.0

Perform harmonisation of a satellite series using matchup records of pairs of 
series sensors and reference sensor. Harmonisation runs an ODR regression on 
matchups and is based on two models: 
    - measurement model of a series sensor; series are AVHRR, HIRS and MW
    - adjustment model for a pair of measurements, i.e. matchup pixels 
    from two (sensor) instruments of a series; the adjustment accounts for the
    spectral differences between two sensors, and the matchup process.

Harmonisation returns calibration coefficients for each sensor in the series and 
their covariance matrix, it propagates coefficients uncertainty to the measured 
quantity (radiance for the considered series), i.e. evaluates the harmonisation 
uncertainty to the series FCDR. """


from numpy import zeros, ones
from random import sample
from os.path import join as pjoin
from optparse import OptionParser
from datetime import datetime as dt
import readHD as rhd
import harFun as har
import unpFun as upf
import visFun as vis


# Set GLOBAL variables 
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # main data folder
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results
#pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
newdir = pjoin(datadir, 'newSim_notime') # no-time dependant simulation data

""" Perform multiple-pair regression with ODR """
def multipH(filelist, series):
    
    p = series.nocoefs # number of calibration parameters
    m = series.novars # # number of measured variables
    nos = series.nosensors # number of sensors in the series
    slist = series.sslab # list of sensors in the series
    inCoef = series.preHcoef # input coefficients to simulations
    
    # Create array of initial beta values for the ODR fit
    hsCoef = zeros((nos,p))    

    # Keep the same values as input coefficients in inCoef 
    for sno in range(nos):
        sl = slist[sno]
        hsCoef[sno,:] = inCoef[sl][0:p]
    
    beta0 = hsCoef.flatten('A') # format to ODR input for initial values
    print '\n\nInitial beta values for ODR'
    print beta0
    
    if series.notime: # work with no-time dependant dataset
    
        Im,Hd,Hr,Hs,sp,mutime,corL = rhd.rHData(newdir, filelist) # read netCDF files 
        series.setIm(Im) # set the series index matrix
          
        # perform odr fit with weights from H uncertainties
        HU2 = Hr**2 + Hs**2 # combined random and systematic
        sodr = har.seriesODR(Hd,HU2,beta0,sp,series) # perform odr on sensors' list
        
    else: # work with data in the main data folder
        
        Im,Hd,Hr,Hs,sp,mutime,corL = rhd.rHData(datadir, filelist)    
        series.setIm(Im) # set the series index matrix      
    
        # create ifixb arrays; fix a3 for all series' sensors
        parfix = zeros(nos*p, dtype=int)
        for sidx in range(1,nos):
            parfix[p*sidx:p*sidx+p-1] = 1
        fixb = parfix.tolist() # ifixb ODR parameter
        print '\n\nifixb array for sensors', slist
        print fixb
        
        # create ifixx arrays; fix orbit temperature To for sensors
        varfix = ones(m*2+1, dtype=int)
        varfix[m] = 0 # fix To for 1st sensor 
        varfix[2*m] = 0 # fix To for 2nd sensor 
        fixx = varfix.tolist() # ifixx ODR parameter
        print '\nifixx array', fixx
    
        # perform odr fit with weights from H uncertainties
        HU2 = Hr**2 + Hs**2 # combined random and systematic
        sodr = har.seriesODR(Hd, HU2, beta0, sp, series, fixb, fixx)
        
    # print summary of odr results and diagnostics from data
    print '\nIndex matrix of sensors in', filelist
    print Im
    print '\nODR output for sensors', slist, '\n'
    sodr.pprint()
    
    print '\n\nRange of input K values [', min(Hd[:,11]), max(Hd[:,11]), ']'
    print 'Range of estimated K values (ODR y) [', min(sodr.y), max(sodr.y), ']'
    print 'Range of estimated K error (ODR epsilon) [', min(sodr.eps), max(sodr.eps), ']'
    print 'Range of input Lref values [', min(Hd[:,0]), max(Hd[:,0]), ']'
    print 'Range of estimated Lref values (from ODR xplus) [', min(sodr.xplus[0,:]), max(sodr.xplus[0,:]), ']'
    print 'Range of estimated Lref error (from ODR delta) [', min(sodr.delta[0,:]), max(sodr.delta[0,:]), ']'
    
    print '\nFirst row of H data matrix'
    print Hd[0,:]
    print '\nLast row of H data matrix'
    print Hd[-1,:]

    return sodr, Hd, Hr, Hs, mutime

# Plot harmonisation results for series sensors    
def plotSSH(sodr, Hd, series, nobj):
    
    nos = series.nosensors # number of sensors in the series
    p = series.nocoefs # number of calibration parameters
    m = series.novars # # number of measured variables
    slist = series.sslab # list of sensors in the series
    inCoef = series.preHcoef # input coefficients to simulations
    Im = series.im # index matrix for series matchups
    
    mpbeta = sodr.beta # calibration coeffs of fitted sensors
    mpcov = sodr.cov_beta # coefficients covariance
    mpcor = vis.cov2cor(mpcov) # coeffs' correlation matrix
    
    cor_ttl = 'Correlation of harmonisation coefficients for pairs\n'+', '.join(filelist)
    #cor_lbl = ['a0', 'a1', 'a2', 'a3'] * nos
    vis.plot_corr_heatmap(mpcor, title=cor_ttl, labels=['a0'])
    print '\nCorrelation of harmonisation coefficients for pairs '+', '.join(filelist) +'\n'
    print mpcor
    
    """ Extract coefficients and covariance of each sensor, 
    compute and plot radiance with 4*sigma uncertainty """
    for sno in range(1,nos): # loop through fitted sensor
        sl = slist[sno] # sensor label
        slab = int(sl[1:3]) # two-digits label in Im matrix
        
        sMidx, eMidx = rhd.sliceHidx(Im, slab) # 1st and last record index
        print '\nFirst and last record for sensor', slab, '[', sMidx, eMidx,']'
        selMU = sample(xrange(sMidx, eMidx), nobj) # Select matchups for plotting    
        
        inC = inCoef[sl] # input coeffs to simulations
        print 'Input coefficients for sensor', slab, ':', inC
        inL = avhrrNx.measEq(Hd[selMU, m+1:2*m+1], inC) # radiance from input coeffs
        
        calC = mpbeta[sno*p:(sno+1)*p] # calib. coeffs for sensor slab
        print 'Fitted coefficients for sensor', slab, ':', calC
        calL = avhrrNx.measEq(Hd[selMU, m+1:2*m+1], calC) # calibrated radiance 
        
        covCC = mpcov[sno*p:(sno+1)*p,sno*p:(sno+1)*p] # coeffs covariance from odr
        print 'Covariance of coefficients for sensor', slab
        print covCC
        # radiance uncertainty from harmonisation
        cLU = avhrrNx.va2ULE(Hd[selMU, m+1:2*m+1],calC,covCC) 
    
        # graphs of radiance bias with 2sigma error bars
        plot_ttl = sl + ' Radiance bias and ' + r'$4*\sigma$'+ ' uncertainty from multiple-pairs ODR covariance'
        vis.LbiasU(inL, calL, cLU, 4, plot_ttl) 
        
    return mpcor

    
if __name__ == "__main__":

    usage = "usage: %prog filelist time-flag series-label"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if 3 != len(args):
        parser.error("incorrect number of arguments")

    # 1st argument: list of netCDF files with matchup harmonisation data
    filelist = args[0] # probably input from a text file
    filelist = ["m02_n19.nc","m02_n15.nc","n15_n14.nc","n14_n12.nc"] # temp solution
    
    # 2nd argument boolean: work with no-/ time dependant simulation data
    notime = args[1] # if False work with time dependant data
    if not isinstance(notime, bool): # input is not boolean type
        notime = str(args[1]).lower()
        if notime in ("yes", "true", "t", "1"):
            notime = True
        else:
            notime = False # preferred
    
    # 3rd argument: series label, not yet used
    series = args[2] 
    
    # Time the execution of the script
    st = dt.now() # start time of script run
    
    # create instance of series class, currently assumed 'avhrr' only
    # TODO: change for different series (label)
    avhrrNx = upf.avhrr(datadir, filelist, notime) 
    
    # perform regression on multiple pairs
    sodr, Hd, Hr, Hs, mutime = multipH(filelist, avhrrNx)
    
    # plot results of harmonisation 
    cCorr = plotSSH(sodr, Hd, avhrrNx, 200)
    
    et = dt.now() # end of script run
    exect = (et-st).total_seconds()
    print '\nTime taken for fitting pairs from', filelist
    print (exect/60.), 'minutes\n'
