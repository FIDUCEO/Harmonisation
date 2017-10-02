""" Generate data for Monte Carlo uncertainty evaluation of odr fit coefficients """

import numpy as np
import netCDF4 as nc
from optparse import OptionParser
import harFun as har
import calc_CC_err as cer

def readMCdt(filename):     
    #print 'Opening netCDF file', filename
    ncid = nc.Dataset(filename,'r')
    
    X = ncid.variables['X'][:,:] # best estimates for explanatory variables
    Y = ncid.variables['Y'][:] # best estimate of dependant variable
    Ur = ncid.variables['Ur'][:,:] # random uncertainty for H vars
    Us = ncid.variables['Us'][:,:] # systematic uncertainty for H vars
    CsUr = ncid.variables['UCs'][:,:] # array for space counts
    CictUr = ncid.variables['UCict'][:,:] # array for ICT counts
    corIdx = ncid.variables['tidx'][:] # matchup time; internal format
    corLen = ncid.variables['corrLen'][:] # length of averaging window
    beta = ncid.variables['beta'][:] # initial vals for beta coefficients

    ncid.close()   
    
    return X,Y,Ur,Us,corIdx,corLen,CsUr,CictUr,beta
    
''' Generate the matrix of errors using W matrix from Sam & Peter '''
def genErr(Hr, Lsys, Tsys, uCs, uCict, corIdx, clen):
    err = np.empty(Hr.shape) # matrix of errors
    nor = err.shape[0] # number of matchups
    v1 = np.ones(nor) # array of ones with size no. of matchups
    
    # Lref, K, CE: random error from Gaussian with sigma from Hr data &mu=0
    err[:,0] = np.random.normal(scale=Hr[:,0]) # Lref random error
    err[:,6] = np.random.normal(scale=Hr[:,6]) # K random error
    err[:,3] = np.random.normal(scale=Hr[:,3]) # CE random error
    
    # Run Sam's function calc_CC_err to generate Cspace averaged errors
    Cs_err, Cs_raw = cer.calc_CC_err(uCs, corIdx, clen)
    err[:,1] = Cs_err # averaged Space count error
    
    # Run Sam's function calc_CC_err to generate Cict averaged errors
    Cict_err, Cict_raw = cer.calc_CC_err(uCict, corIdx, clen)
    err[:,2] = Cict_err # averaged ICT count error
    
    # draw a value of systematic error for Lict and To
    errL = np.random.normal(scale=Lsys) # Lict systematic error
    errT = np.random.normal(scale=Tsys) # To systematic error
    # Lict, To: Gaussian random with sigma from Hr & mu=0, + systematic err
    err[:,4] = np.random.normal(scale=Hr[:,4]) + errL*v1 # Lict error
    err[:,5] = np.random.normal(scale=Hr[:,5]) + errT*v1 # To error
            
    return err

if __name__ == "__main__":

    usage = "usage: %prog sensor_label"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    if 1 != len(args):
        parser.error("No data file given")
   
    ncfile = args[0] # netCDF filename / path +filename 
    X,Y,Hr,Hs,corIdx,corLen,CsUr,CictUr,b0odr = readMCdt(ncfile)     
    #wma = 1./(1+corLen[0]*2) # moving average weight
    #win = np.ones(CsUr.shape[0]) * wma
    
    ''' compile data for ODR run in the MC trial'''
    sLict = Hs[0,4] # systematic error Lict
    sTo = Hs[0,5] # systematic error To
    # generate error structure 
    errStr = genErr(Hr, sLict, sTo, CsUr, CictUr, corIdx, corLen)
    # add error structure to X & Y best estimates
    Xdt = X + errStr[:,1:6] # X variables
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable
    
    # run ODR on new X & Y vals and Hr weights 
    fixb = [1,1,1,0] # fix a3 coefficient
    modr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb, Hs)
    
    b0 = modr.beta # store fit coefficients
    print b0[0], b0[1], b0[2], b0[3]
