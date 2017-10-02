""" Generate data for Monte Carlo uncertainty evaluation of odr fit coefficients """

from numpy import zeros, ones
import numpy.random as random
import netCDF4 as nc
from optparse import OptionParser
import harFun as har
import errStruct as mce 
import unpFun as upf

# Set GLOBAL variables 
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # main data folder
#datadir = "/home/ad6/Data" # main data folder in eoserver
#datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS
filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]  

def readMCdata(filename):     
    #print 'Opening netCDF file', filename
    ncid = nc.Dataset(filename,'r')
    
    X = ncid.variables['X'][:,:] # best estimates for explanatory variables
    Y = ncid.variables['Y'][:] # best estimate of dependant variable
    Ur = ncid.variables['Ur'][:,:] # random uncertainty for H vars
    CsUr = ncid.variables['UCs'][:,:] # array for space counts
    CictUr = ncid.variables['UCict'][:,:] # array for ICT counts
    corIdx = ncid.variables['tidx'][:] # matchup time; internal format
    muCnt = ncid.variables['mxSl'][:] # matchup time; internal format
    corLen = ncid.variables['corrLen'][:] # length of averaging window
    beta = ncid.variables['beta'][:] # initial vals for beta coefficients
    notdep = ncid.variables['notdep'][:] # initial vals for beta coefficients
    notd = notdep[0] # time dependancy flag
    if not notd:
        Us = ncid.variables['Us'][:,:] # systematic uncertainty for H vars
    else:
        Us = zeros(Ur.shape) # zero systematic uncertainty

    ncid.close()   
    
    return X,Y,Ur,Us,corIdx,muCnt,corLen,CsUr,CictUr,beta, notd
    
''' Generate the matrix of errors using W matrix from Sam & Peter '''
def genErr(Hr, Lsys, Tsys, uCs, uCict, corIdx, clen, notime):
    err = zeros(Hr.shape) # matrix of errors
    nor = err.shape[0] # number of matchups
    v1 = ones(nor) # array of ones with size no. of matchups
    
    # Lref, K, CE: random error from Gaussian with sigma from Hr data &mu=0
    err[:,0] = random.normal(scale=Hr[:,0]) # Lref random error
    err[:,6] = random.normal(scale=Hr[:,6]) # K random error
    err[:,3] = random.normal(scale=Hr[:,3]) # CE random error
    
    # Run Sam's function calc_CC_err to generate Cspace averaged errors
    Cs_err = mce.calc_CC_err(uCs, corIdx, clen)
    err[:,1] = Cs_err # averaged Space count error
    
    # Run Sam's function calc_CC_err to generate Cict averaged errors
    Cict_err = mce.calc_CC_err(uCict, corIdx, clen)
    err[:,2] = Cict_err # averaged ICT count error
    
    # combined systematic and random errors for Lict and To
    err[:,4] = random.normal(scale=Hr[:,4]) # random error Lict
    if notime:
        err[:,5] = 0 # To error 0; this column will not be used
    else:   # time dependant data
        
        # add systematic error to Lict random error
        errL = random.normal(scale=Lsys) # Lict systematic error
        err[:,4] += errL*v1 # Lict error
    
        # combined systematic and random To error
        errT = random.normal(scale=Tsys) # To systematic error
        err[:,5] = random.normal(scale=Hr[:,5]) + errT*v1 # To error
            
    return err

if __name__ == "__main__":

    usage = "usage: %prog netCDF-filename"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    if 1 != len(args):
        parser.error("No data file given")
   
    ncfile = args[0] # netCDF filename / path +filename 
    X,Y,Hr,Hs,corIdx,muCnt,cLen,CsUr,CictUr,bodr,notd = readMCdata(ncfile)     
    avhrrNx = upf.avhrr(datadir, filelist, notd)   
    
    ''' compile data for ODR run in the MC trial'''
    if notd: # no-time dependant data: X matrix has no To column
        col = 5 # last column to read data from error matrix
        sLict = 0
        sTo = 0
    else:
        col = 6    
        sLict = Hs[0,4] # systematic error Lict
        sTo = Hs[0,5] # systematic error To
        fixb = [1,1,1,0] # fix a3 coefficient
        fixx = [1,1,1,1,0] # fix orbit temperature variable

    # Generate errors with the weight matrix W from Peter & Sam
    #errStr = genErr(Hr, sLict, sTo, CsUr, CictUr, corIdx, cLen, notd)
    errStr = mce.genErr(Hr, sLict, sTo, CsUr, CictUr, corIdx, muCnt, cLen, notd)

    # add errStr to X & Y best estimates
    Xdt = X + errStr[:,1:col] # X variables
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable
    
    # run ODR on new X & Y vals and weights 
    if notd: # newSim_notime data: a3 = 0, To = 0, Hs = 0
        mcodr = har.odr4MC(Xdt, Ydt, Hr, bodr, avhrrNx)
        b0 = mcodr.beta # store fit coefficients
        print b0[0], b0[1], b0[2]
    
    else: # fix a3 and To to input, weights on random & systematic uncertainty
        mcodr = har.odr4MC(Xdt, Ydt, Hr, bodr, avhrrNx, fixb, fixx, Hs)
        print b0[0], b0[1], b0[2], b0[3]
