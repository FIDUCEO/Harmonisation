""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 01-02-2017
Evaluate ODR fit coefficients uncertainty via Monte Carlo (MC) trials. 
Work with a reference-sensor pair, use the same error structure with the weights 
in ODR, the combined random and systematic. Systematic effect is constant, 
set to average of all matchups, to be comparable with Peter & Sam's results. """

import numpy as np
from datetime import datetime as dt
import pandas as pd
import readHD as rhd
import harFun as har
import unpFun as upf
import mcoutVis as vmc

filelist = ["m02_n15.nc"]#,"m02_n19.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
#datadir = "/home/ad6/Data" # in eoserver
#datadir = "/group_workspaces/cems2/fiduceo/Users/adilo/Data" # in CEMS
beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
slist = rhd.sensors(filelist) # list of sensors in filelist
ncfile = filelist[0] # netCDF file to work with 
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value
print 'Processing sensors', slist, 'from file', ncfile

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 

# perform odr fit and extract output
fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb, Hs)
print 'ODR results on Jon data, weights from combined random &systematic uncertainty'
podr.pprint()
b0odr = podr.beta # odr fit coefficients
sd0odr = podr.sd_beta # standard error of fit coefficients 
cov0odr = podr.cov_beta # odr evaluated covariance matrix

''' Generate data and other info for Monte Carlo runs''' 
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
sigma = np.sqrt(Hr**2 + Hs**2) # standard uncertainty of X and Y vars

# MC runs ODR on new data: best estimate + full correlation error draw
noit = 10 # number of MC iterations
b0 = np.empty([noit, len(beta)]) # array to store beta vals from MC

st = dt.now() # start of MS run
for i in range(noit):
    
    errStr = np.random.normal(sigma) # generate errors for MC trial
    # add error to X & Y best estimates
    Xdt = X.T + errStr[:,1:6] # X variables for odr in the MC trial
    Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable for odr 
    
    # run ODR on new X & Y vals and Hr weights 
    fixb = [1,1,1,0] # fix a3 coefficient
    podr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb)
    
    # store fit coefficients
    b0[i] = podr.beta

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print 'Time taken for', noit, 'MC trials', (exect/60.), 'minutes'

b0mc = pd.DataFrame(b0)
b0mc.to_csv("mcb0.csv")
#
mcb0 = np.loadtxt(open("mcb0.csv", "rb"), delimiter=",", skiprows=1)

''' CHECK the stats and graphs below '''
# print and plot odr coeffs stats from MC runs; get coeff covariance from MC stats
mcmean, mccov = vmc.mcStats(mcb0, b0odr, sd0odr, cov0odr) 

# evaluate AVHRR Earth radiance and uncertainty
avhrrNx = upf.avhrr(1, 2) # instance of avhrr class
inLE = avhrrNx.measEq(X.T, inCoef[s2]) # input radiance (from input coeffs)

# evaluate radiance uncertainty from coeffs uncertainty
LUodr = avhrrNx.va2ULE(X.T, b0odr, cov0odr) # at best estimate with ODR evaluated cov
LUmc = avhrrNx.va2ULE(X.T, mcmean, mccov) # at best estimate with MC sample cov

CE = (X[0,:]-X[2,:]).transpose() # Earth counts

#vmc.plotLaU(CE, inLE, Y, LUodr, LUmc)
