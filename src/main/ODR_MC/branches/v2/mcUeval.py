""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 26-01-2017
Evaluation of ODR fit coefficients uncertainty via Monte-Carlo. 
Start with odr fit of a reference-sensor pair. """

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as dt
import readHD as rhd
import harFun as har
import upFun as upf
import genMCdata as gdt
import visMCr as vmc

filelist = ["m02_n15.nc","m02_n19.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
slist = rhd.sensors(filelist) # list of sensors in filelist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 

ncfile = filelist[0] # netCDF file to work with 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value

# read data from the netCDF file
st = dt.now() # start of reading
rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
et = dt.now() # end of reading
exect = (et-st).total_seconds()
print 'Time taken to read data (in seconds)', exect

# perform odr fit and extract output
st = dt.now() # start of odr fit
podr = har.odrP(Hd, Hr, beta)
beta0 = podr.beta # fit coefficients
sigb0 = podr.sd_beta # standard error of fit coefficients 
covb0 = podr.cov_beta
print 'ODR results from the run on simulated data and random uncertainty'
podr.pprint()
et = dt.now() # end of odr fit
exect = (et-st).total_seconds()
print 'Time taken from odr fit (in seconds)', exect

# Extract X,Y best estimates from ODR output, and fixed value uncertainties
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
sLict = np.mean(Hs[:,4]) # systematic uncertainty Lict
sTo = np.mean(Hs[:,5]) # systematic uncertainty To
szCC = 10 # sample size of calib. counts per scanline
rrCC = np.mean(Hr[:,3])/np.sqrt(szCC) # constant calib. count random uncertainty
maLen = 51 # width of average moving window for calibration counts

# Monte Carlo for ODR on new draw of data: best estimate + err respecting correlations
noit = 10 # number of MC iterations
b0 = np.empty([noit, len(beta)]) # array to store beta vals from MC

st = dt.now() # start of MS run
for i in range(noit):
    # compile data for the ODR run
    errStr = gdt.genpCS(Hr,sLict,sTo,rrCC,maLen) # generated errors for the data
    # add errStr to X & Y vbest estimates
    Xdt = X.T + errStr[:,1:6]
    Ydt = Y + errStr[:,0] + errStr[:,6]
    
    # run ODR on new X & Y vals and Hr weights 
    fixb = [1,1,1,0] # fix a3 coefficient
    podr = har.odr4MC(Xdt, Ydt, Hr, beta0, fixb)
    
    # store fit coefficients
    b0[i] = podr.beta
print 'Beta coefficients from ODR run on MC simulated data with error structure'
print b0
et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print 'Time taken from MC run (in minutes)', (exect/60.)

# print and plot odr coeffs stats from MC runs; get coeff covariance from MC stats
mccov = vmc.MCwODRstats(b0, beta0, sigb0, covb0) 

# evaluate AVHRR Earth radiance and uncertainty
avhrrNx = upf.avhrr(1, 2) # instance of avhrr class
beta = inCoef[s2] # input calibration coeffs 
inLE = avhrrNx.measEq(X.T, beta) # input radiance (from input coeffs)

# evaluate radiance uncertainty from coeffs uncertainty
LUodr = avhrrNx.va2ULE(X.T, beta0, covb0) # at best estimate with ODR evaluated cov
LUmc = avhrrNx.va2ULE(X.T, beta0, mccov) # at best estimate with MC sample cov

CE = (X[0,:]-X[2,:]).transpose() # Earth counts

vmc.plotLaU(CE, inLE, Y, LUodr, LUmc)
