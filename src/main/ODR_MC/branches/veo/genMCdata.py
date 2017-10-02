#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
import readHD as rhd
import harFun as har

""" Create netCDF file with data that is used in MC trials """
def genMCdata(fn,M,Xbe,Ybe,Hr,Hs,CsUr,CictUr,corIdx,cLen,beta0):
    
    # create netCDF file with data for MC trials
    ncid = nc.Dataset(fn,mode='w',format='NETCDF4',clobber=True)

    M = ncid.createDimension('M',M) # number of matchups
    m = ncid.createDimension('m',7) # number of columns in Hr and Hs
    n = ncid.createDimension('n',5) # number of columns Xbe data
    p = ncid.createDimension('p',4) # number of calib.coefficients   
    scol = ncid.createDimension('scol',1) 
    wlen = 2*cLen[0] + 1 # length of moving average window for cal.counts
    w = ncid.createDimension('w',wlen) # no cols calib.count uncert CsUr & CictUr

    X = ncid.createVariable('X','f4',('M','n',),zlib=True,complevel=9)
    X.Description='Best estimates of X vars (M,n): Cs, Cict, CE, Lict, To'
    Y = ncid.createVariable('Y','f4',('M','scol',),zlib=True,complevel=9)
    Y.Description='Best estimates of Y (M,1) per matchup'
    Ur = ncid.createVariable('Ur','f4',('M','m',),zlib=True,complevel=9)
    Ur.Description='Random uncertainties for Hdata: Lref, Cs, Cict, CE, Lict, To, K'
    Us = ncid.createVariable('Us','f4',('M','m',),zlib=True,complevel=9)
    Us.Description='Systematic uncertainties for Hdata: Lref, Cs, Cict, CE, Lict, To, K'
    UCs = ncid.createVariable('UCs','f4',('M','w',),zlib=True,complevel=9)
    UCs.Description='Uncertainty per scanline for count calibration data (Space)'
    UCict = ncid.createVariable('UCict','f4',('M','w',),zlib=True,complevel=9)
    UCict.Description='Uncertainty per scanline for count calibration data (ICT)'
    tidx = ncid.createVariable('tidx','f8',('M'),zlib=True,complevel=9)
    tidx.Description='Input for pixel-to-pixel correlations (scanline times)'
    corrLen = ncid.createVariable('corrLen','f4',('scol'))        
    corrLen.Description='Half-length of the moving average window'
    beta = ncid.createVariable('beta','f4',('p'))
    beta.Description='Initial coefficient estimates for the fit in MC trials'
    
    X[:,:] = Xbe[:,:] # best estimates for explanatory variables
    Y[:] = Ybe[:] # best estimate of dependant variable
    Ur[:,:] = Hr[:,:] # random uncertainties of all variables
    Us[:,:] = Hs[:,:] # systematic uncertainties
    UCs[:,:] = CsUr[:,:] # array of space calibration counts uncertainty
    UCict[:,:] = CictUr[:,:] # array of ICT calibration counts uncertainty
    tidx[:] = corIdx[:] # scanline times
    corrLen[:] = cLen[:] # half-length of moving average window
    beta[:] = beta0[:] # beta coefficients from ODR run on Jon's data

    ncid.close()   
   

""" Run ODR on Jon's data and create a netCDF file with ODR results and 
unecrtainties from simulated data to use in MC trials for error structure. """
if __name__ == "__main__":

    # data folder and netCDF files
    datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS
    datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder in laptop
    filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]    
    ncfile = filelist[4] # netCDF file to work with 

    s2 = ncfile[4:7]
    print '\nProcessing sensor', s2, 'from file', ncfile
    inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
    beta = [-10., -4.e-3, 1.e-5, 0.0] # initialise fit coefficients
    beta[3] = inCoef[s2][3] # set a3 to the input value

    # read simulation data from Jon's netCDF file
    rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
    # set systematic uncertainties equivalent to Peter&Sam LS optimisation
    Hs = har.resetHs(Hs, rsp) 
    print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'
    print 'Range of values of CEarth random uncertainties [',Hr[:,3].min(axis=0), ',', Hr[:,3].max(axis=0),']'
    Hr[:,3] = np.absolute(Hr[:,3])
    print 'Range of values of CEarth random uncertainties [',Hr[:,3].min(axis=0), ',', Hr[:,3].max(axis=0),']'
    
    # perform odr fit and extract output
    fixb = [1,1,1,0] # fix a3 coefficient
    podr = har.odrP(Hd, Hr, beta, fixb, Hs)
    print '\n\nODR results on Jon\'s data weighted by combined random &systematic uncertainty'
    podr.pprint()
    
    ''' Generate data and other info for Monte Carlo runs''' 
    b0odr = podr.beta # odr fit coefficients - beta0
    Y = podr.y # best est.of adjusted reference radiance: Lref + K
    X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
    cLen = [25] # int(corLen[0]) # scanlines moving average half-window 
    nor = Im[0,2] #  number of matchups
    
    tfn = s2 + '_mcdata.nc'
    genMCdata(tfn,nor,X.T,Y,Hr,Hs,csUr,cictUr,corIdx,cLen,b0odr)
    