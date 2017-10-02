#!/usr/bin/env python

from os.path import join as pjoin
import numpy as np
import readHD as rhd
import harFun as har
import unpFun as upf
import plot as pl

filelist = ["m02_n15.nc"]#,"m02_n19.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs

beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
slist = rhd.sensors(filelist) # list of sensors in filelist
ncfile = filelist[0] # netCDF file to work with 
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value
nop = len(filelist) # number of sensor pairs
nos = len(slist) # number of sensors

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 

fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb, Hs) # run odr with random & syst weights
b0odr = podr.beta # odr fit coefficients
Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
print 'ODR output for', s2
podr.pprint()

# plot results
avhrrNx = upf.avhrr(nop, nos) # instance of class avhrr series 
inL = avhrrNx.measEq(Hd[:,1:6], inCoef[s2])
# Y uncertainty: combined random uncert. of ref.radiance and adjustment vals
Yu = np.sqrt(Hr[:,0]**2+Hr[:,6]**2 + Hs[:,0]**2+Hs[:,6]**2)
calL = avhrrNx.measEq(Hd[:,1:6], b0odr)
#cLU = avhrrNx.uncLE(Hd[:,1:6],b0odr,Hr[:,1:6],podr.cov_beta)
cLU = avhrrNx.va2ULE(Hd[:,1:6],b0odr,Hr[:,1:6],podr.cov_beta)

# graphs showing fit results
Figure = pl.Plot(pltdir, Im[0,0], Im[0,1], Im[0,2], corIdx) # Plot instance
Figure.plotTrueErr(Hr, podr, 'odr') 
Figure.plotTrueVar(Hd, podr, 'odr')
# space clamped Earth counts -> Hd[:,1]-Hd[:,3]
Figure.plotFit(inL,calL,Yu,podr.eps,cLU,Hd[:,1]-Hd[:,3],'odr')

# generate errors for MC trial and add to X & Y best estimates
sigma = np.sqrt(Hr**2 + Hs**2) # standard uncertainty of X and Y vars
errStr = np.random.normal(sigma) 
Xdt = X.T + errStr[:,1:6] # X variables for odr in the MC trial
Ydt = Y + errStr[:,0] + errStr[:,6] # Y var: Lref + K

# run ODR on new X & Y vals and Hr weights
mcodr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb, Hs) # ODR result
# store fit coefficients
b0 = mcodr.beta # beta coeffs for the MC trial
print 'ODR output for perturbed', s2
mcodr.pprint()
mcL = avhrrNx.measEq(Hd[:,1:6], b0)
mcLU = avhrrNx.va2ULE(Hd[:,1:6],b0,Hr[:,1:6],mcodr.cov_beta)

# Plot results
Figure.plotTrueErr(Hr, mcodr, 'mct') 
Figure.plotTrueVar(Hd, mcodr, 'mct')
# space clamped Earth counts -> Hd[:,1]-Hd[:,3]
Figure.plotFit(inL,mcL,Yu,mcodr.eps,mcLU,Hd[:,1]-Hd[:,3],'mct')

# OUTPUT i need from each CEMS job: Calibration Coefficients 
print b0
