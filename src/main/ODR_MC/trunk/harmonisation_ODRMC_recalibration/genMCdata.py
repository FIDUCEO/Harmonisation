#!/usr/bin/env python

import netCDF4 as nc
from numpy import ones, unique
from os.path import join as pjoin
from optparse import OptionParser
import readHD as rhd
import harFun as har
import unpFun as upf


# Set GLOBAL variables 
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # main data folder
#datadir = "/home/ad6/Data" # main data folder in eoserver
#datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS
filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]  


""" Create netCDF file with data that is used in MC trials """
def genMCdata(fn,M,Xbe,Ybe,Hr,CsUr,CictUr,slTime,mCnt,cLen,beta0,series,Hs=None):
    
    # create netCDF file with data for MC trials
    ncid = nc.Dataset(fn,mode='w',format='NETCDF4',clobber=True)

    M = ncid.createDimension('M',M) # number of matchups
    slno = len(slTime)
    S = ncid.createDimension('S',slno) # number of scanlines
    m = ncid.createDimension('m',7) # number of columns in Hr and Hs
    n = ncid.createDimension('n',series.novars) # number of columns in Xbe data
    p = ncid.createDimension('p',series.nocoefs) # number of calib.coefficients   
    scol = ncid.createDimension('scol',1) 
    wlen = 2*cLen[0] + 1 # length of moving average window for cal.counts
    w = ncid.createDimension('w',wlen) # no cols calib.count uncert CsUr & CictUr

    X = ncid.createVariable('X','f4',('M','n',),zlib=True,complevel=9)
    X.Description='Best estimates of X vars (M,n): Cs, Cict, CE, Lict, To'
    Y = ncid.createVariable('Y','f4',('M','scol',),zlib=True,complevel=9)
    Y.Description='Best estimates of Y (M,1) per matchup'
    Ur = ncid.createVariable('Ur','f4',('M','m',),zlib=True,complevel=9)
    Ur.Description='Random uncertainties for Hdata: Lref, Cs, Cict, CE, Lict, To, K'
    UCs = ncid.createVariable('UCs','f4',('S','w',),zlib=True,complevel=9)
    UCs.Description='Uncertainty per scanline for count calibration data (Space)'
    UCict = ncid.createVariable('UCict','f4',('S','w',),zlib=True,complevel=9)
    UCict.Description='Uncertainty per scanline for count calibration data (ICT)'
    tidx = ncid.createVariable('tidx','f8',('S'),zlib=True,complevel=9)
    tidx.Description='Input for pixel-to-pixel correlations (scanline times)'
    mxSl = ncid.createVariable('mxSl','f8',('S'),zlib=True,complevel=9)
    mxSl.Description='Number of matchups per scanline'
    corrLen = ncid.createVariable('corrLen','f4',('scol'))        
    corrLen.Description='Half-length of the moving average window'
    beta = ncid.createVariable('beta','f4',('p'))
    beta.Description='Initial coefficient estimates for the fit in MC trials'
    notdep = ncid.createVariable('notdep','f4',('scol'))
    notdep.Description='Flag for (no-) time depenedent data'
    
    X[:,:] = Xbe[:,:] # best estimates for explanatory variables
    Y[:] = Ybe[:] # best estimate of dependant variable
    Ur[:,:] = Hr[:,:] # random uncertainties of all variables
    UCs[:,:] = CsUr[:,:] # array of space calibration counts uncertainty
    UCict[:,:] = CictUr[:,:] # array of ICT calibration counts uncertainty
    tidx[:] = slTime[:] # scanline times
    mxSl[:] = mCnt[:] # matchups per scanline
    corrLen[:] = cLen[:] # half-length of moving average window
    beta[:] = beta0[:] # beta coefficients from ODR run on Jon's data
    notdep[:] = [int(series.notime)]
    
    if not series.notime: # add systematic uncertainty matrix for time dependant data
        Us = ncid.createVariable('Us','f4',('M','m',),zlib=True,complevel=9)
        Us.Description='Systematic uncertainties for Hdata: Lict and To columns'
        Us[:,:] = Hs[:,:] # systematic uncertainties

    ncid.close()   
   

""" Run ODR on Jon's data and create a netCDF file with ODR results and 
unecrtainties from simulated data to use in MC trials for error structure. """
if __name__ == "__main__":

    usage = "usage: %prog time-flag"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if 1 != len(args):
        parser.error("Expected boolean argument")

    # 2nd argument boolean: work with no-/ time dependant simulation data
    notime = args[0] # if False work with time dependant data
    
    # create instance of avhrr sensor series 
    avhrrNx = upf.avhrr(datadir, filelist, notime)
    p = avhrrNx.nocoefs # number of calibration parameters
    m = avhrrNx.novars # # number of measured variables
    nos = avhrrNx.nosensors # number of sensors in the series
    slist = avhrrNx.sslab # list of sensors in the series
    inCoef = avhrrNx.preHcoef # input coefficients to simulations
    
    ncfile = filelist[4] # netCDF file to work with 
    s2 = ncfile[4:7]
    beta = inCoef[s2][0:p] # initial values for ODR coefficients

    if notime: # work with newSim_notime data
        
        # read data from the netCDF file
        newdir = pjoin(datadir, 'newSim_notime') # no-time dependent simulation data
        rsp,Im,Hd,Hr,Hs,corIdx,cLen,csUr,cictUr = rhd.rHDpair(newdir, ncfile)
        Hr[:,5] = 1. # change 0 uncertainty of To for ODR to work
        
        # perform odr fit, weights from random uncertainty
        podr = har.odrP(Hd, Hr, beta, avhrrNx) 
        
        tfn = s2 + '_ntd_mcdata.nc' # filename where ODR output will be stored
        
    else:      
        newdir = pjoin(datadir, 'newSim') # time dependent simulation data
        rsp,Im,Hd,Hr,Hs,corIdx,cLen,csUr,cictUr = rhd.rHDpair(newdir, ncfile)
        # data in the main data folder from 12-Mar-2017
        #rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
        # set systematic uncertainties equivalent to Peter&Sam GN optimisation
        Hs = rhd.resetHs(Hs, rsp) 
        
        # create ifixb array; fix a3 
        parfix = ones(p, dtype=int)
        parfix[-1] = 0
        fixb = parfix.tolist() # ifixb ODR parameter
        print '\nifixb array', fixb
        
        # create ifixx array; fix To variable
        varfix = ones(m, dtype=int)
        varfix[-1] = 0 # fix To 
        fixx = varfix.tolist() # ifixx ODR parameter
        print '\nifixx array', fixx
    
        # perform odr fit, weights from combined random and systematic uncertainty
        podr = har.odrP(Hd, Hr, beta, avhrrNx, fixb, fixx, Hs) # fit to adjusted ref.radiance
        
        tfn = s2 + '_td_mcdata.nc'
   
    ''' Generate data and other info for Monte Carlo runs''' 
    bodr = podr.beta # odr fit coefficients
    covodr = podr.cov_beta # odr evaluated covariance matrix
    Y = podr.y # best est.of adjusted reference radiance: Lref + K
    X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
    nor = Im[0,2] #  number of matchups
    
    # get unique scanlines, first matchup idx &number of matchup pixels per scanline
    slt,midx,mcnt = unique(corIdx,return_index=True,return_counts=True)
    rCSar = csUr[midx,:] # Cspace random uncert. per scanline: arrays of 51 slines
    rCICTar = cictUr[midx,:] # Cict random uncert. per scanline: arrays of 51 slines
    
    if notime:
        genMCdata(tfn,nor,X.T,Y,Hr,rCSar,rCICTar,slt,mcnt,cLen,bodr,avhrrNx)
    else:
        genMCdata(tfn,nor,X.T,Y,Hr,rCSar,rCICTar,slt,mcnt,cLen,bodr,avhrrNx,Hs)
    