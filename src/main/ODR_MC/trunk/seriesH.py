#!/usr/bin/env python

""" FCDR harmonisation modules
    Project:        H2020 FIDUCEO
    Author:         Arta Dilo \NPL MM
    Reviewer:       Peter Harris \NPL MM, Sam Hunt \NPL ECO
    Date created:   12-12-2016
    Last update:    02-05-2017
    Version:        12.0

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


from numpy import zeros, ones, savetxt
from random import sample
from os.path import join as pjoin
from optparse import OptionParser
from datetime import datetime as dt
import readHD as rhd
import harFun as har
import unpFun as upf
import visFun as vis
from plotHaR import plotSDfit as pltErrs

# Set GLOBAL variables 
datadir = "D:\Projects\FIDUCEO\Data" # main data folder in laptop
#datadir = "/home/ad6/Data" # main data folder in eoserver
#datadir = "/group_workspaces/cems2/fiduceo/Data/Matchup_Simulated/Data" # in CEMS
mcrdir = pjoin(datadir, 'Results') # folder for MC trials results
#pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
hvars = 12 # number of columns in the H data matrices of a series


""" Perform multiple-pair regression with ODR """
def multipH(filelist, series, dstype):
    
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

    b0 = hsCoef.flatten('A') # format to ODR input for initial values
    print '\n\nInitial beta values for ODR'
    print b0
    
    if series.notime: # work with not time-dependant dataset
        folder = pjoin(datadir, 'newSim_notime') # data folder
        # columns for X variables in the H matrices
        cols = [x for x in range(hvars) if (x!=5 and x<10)]
        
        ## create ifixb arrays; fix a coeffs for reference sensor
        #parfix = zeros(nos*p, dtype=int)
        #parfix[0:p] = 1
        #fixb = parfix.tolist() # ifixb ODR parameter
        #print '\n\nifixb array for sensors', slist
        #print fixb
        fixb = None
        fixx = None
        
    else: # work with data in the main/time-dependent data folder        
        folder = pjoin(datadir, 'newSim') # data folder
        # columns for X variables in the H matrices
        cols = [x for x in range(hvars) if x < 11]
    
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
        
    #if Hd.shape[1] != hvars:
    #    sys.exit('Incorrect shape of harmonisation matrices')
    
    # work with real datasets; currently a different folder for the data
    if dstype == 'r': 
        folder = pjoin(datadir, 'Harm_RealData') # real data folder

    # read data from the list of netCDF files   
    Im,Hd,Hr,Hs,sp,mutime,corL,Is,sl1,mxsl1,sl2,mxsl2,CsU1,CictU1,CsU2,CictU2 = rhd.rHData(folder, filelist) 
    series.setIm(Im) # set the series index matrix    

    # perform odr on all sensors from the list
    print '\nRunning ODR for multiple pairs\n'    
    sodr = har.seriesODR(Hd[:,cols],Hd[:,11],Hr[:,cols],Hr[:,11],b0,sp,series,fixb,fixx)

    print '\nODR output for sensors', slist, '\n'
    sodr.pprint() # print summary of odr results 
    print '\nCost function in final iteration (Sum of squares):', sodr.sum_square
    print '\nSum of squares of epsilon error (K):', sodr.sum_square_eps
    print '\nSum of squares of delat error (H variables):', sodr.sum_square_delta

    #print '\n odr iwork array'
    #print sodr.iwork
    
    print '\nIndex matrix of sensors in', filelist
    print Im
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
    for sno in range(1,nos): # loop through fitted sensors
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
        cLU = avhrrNx.harUnc(Hd[selMU, m+1:2*m+1],calC,covCC) 
    
        # graphs of radiance bias with 2sigma error bars
        plot_ttl = sl + ' Radiance bias and ' + r'$4*\sigma$'+ ' uncertainty from multiple-pairs ODR covariance'
        vis.LbiasU(inL, calL, cLU, 4, plot_ttl) 
        
    return mpcor

    
if __name__ == "__main__":

    usage = "usage: %prog time-flag data-type series-label"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if len(args) < 3:
        parser.error("Insufficient number of arguments")

    # 1st argument: boolean, defines what dataset to work with 
    notime = args[0] 
    if not isinstance(notime, bool): # if input is not True/False
        notime = str(args[0]).lower()
        if notime in ("yes", "true", "t", "1"):
            notime = True # work with not-time dependent dataset
        else:
            notime = False # work with time dependent dataset
    
    dtype = args[1] # 2nd argument: type of dataset, r for real, s for simulated
    sslab = args[2] # 3rd argument: series label e.g. avhrr; currently not used
    
    # TODO: add 4th argument the list of netCDFs with harmonisation data
    #filelist = args[3] # input from a text file ?
    #filelist = ["m02_n15.nc","n15_n14.nc","n14_n12.nc","n12_n11.nc"] #,"n11_n10.nc","n10_n09.nc"] 
    #filelist = ["m02_n19.nc","m02_n17.nc","m02_n15.nc","n19_n15.nc","n17_n15.nc"] 
    filelist = ["m02_n19.nc","m02_n15.nc","n19_n15.nc","n15_n14.nc"] 
   
    # Time the execution of harmonisation
    st = dt.now() # start time of script run
    
    # create instance of series class, currently assumed 'avhrr' only
    # TODO: change for different series label
    avhrrNx = upf.avhrr(datadir, filelist, notime) 
    
    # perform regression on multiple pairs
    sodr, Hd, Hr, Hs, mutime = multipH(filelist, avhrrNx, dtype)

    # Store beta and covariance to text files
    fn = 'sh_' # compile filename
    if dtype =='r':
        fn += 'rd_'
    if notime:
        fn += 'notd_'
    else:
        fn += 'td_'    
    # store calibration coefficients for sensors in the avhrr/.. class
    fnb = fn+'beta.txt' # filename for beta coefficients
    fnb = pjoin(mcrdir, fnb) # path & filename
    savetxt(fnb, sodr.beta, delimiter=',')
    # store coefficients' covariance of sensors in the avhrr/.. class
    fnc = fn+'bcov.txt' # filename for covariance matrix
    fnc = pjoin(mcrdir, fnc) # path & filename
    savetxt(fnc, sodr.cov_beta, delimiter=',')


    et = dt.now() # end of harmonisation run
    exect = (et-st).total_seconds()
    print '\nTime taken for fitting pairs', filelist
    print (exect/60.), 'minutes\n'
    
    # plot results of harmonisation 
    cCorr = plotSSH(sodr, Hd, avhrrNx, 250)
    
    # plot weighted residuals
    noMU = Hd.shape[0] # number of matchups in multiple pairs
    nobj = 10000 # number of mathcup records to plot
    pltErrs(avhrrNx, noMU, nobj, '', mutime, Hr, sodr, weight=1)
