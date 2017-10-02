#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author:         Arta Dilo / NPL MM
    Date created:   03-03-2017
    Last update:    03-03-2017
    Version:        1.0
Perform harmonisation of a reference-sensor pair using ODR. Re-parameterise for 
orbit temperature: test by moving the term a3*To on the left side. """

import scipy.odr as odr
from numpy import zeros, append, insert
from random import sample
import readHD as rhd
import harFun as har
import unpFun as upf
import visFun as vis


""" AVHRR measurement equation not including orbit temperature """
def avhrrME(CE,Cs,Cict,Lict,a0,a1,a2):        
    # Earth radiance from Earth counts and  calibration data
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE)
    return LE # return Earth radiance

""" Model function, i.e. fcn argument in ODR package """
def fcnP(coef, data):
    a0 = coef[0] # AVHRR model coefficients
    a1 = coef[1]
    a2 = coef[2]
    # transposed ndarrays 
    CE = data[2,:] # Earth counts
    Cs = data[0,:] # space counts
    Cict = data[1,:] # ICT counts
    Lict = data[3,:] # ICT radiance
    
    LE = avhrrME(CE,Cs,Cict,Lict,a0,a1,a2)
    #print 'Current iteration coefficients:', [a0,a1,a2]
    return LE # return Earth radiance

""" Perform LS fit for a sensor-reference pair with low-level odr function """
def harP(Hdata, Hr, b0, a3):
    # extract variables Cs, Cict, CE, Lict from Hdata matrix
    X = Hdata[:,1:5].transpose() # transpose data matrix
    ''' Y is adjusted radiance, i.e. reference radiance with spectral adjustment
    and the time effect modelled by orbit temperature: Y = Lref + K '''
    Y = Hdata[:,0] + Hdata[:,6] #- a3 * Hdata[:,5]

    VX = (Hr[:,1:5]**2).transpose() # sigma^2 of X
    VY = Hr[:,0]**2 + Hr[:,6]**2 #+ a3**2 * Hr[:,5]**2 # sigma^2 of Y 
    
    # perform odr fit (low level function)
    fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,full_output=1)
    odrFit = odr.Output(fit) # get regression output     
    
    return odrFit # return odr output


filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
ncfile = filelist[4] # netCDF file to work with 

s2 = ncfile[4:7] # label of 2nd sensor
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices\n'

# perform odr fit without To effect and print fit result
print '\nODR fit on Jon\'s data, weights from random uncertainty, no To effect\n'
beta = [-10., -4.e-3, 1.e-5] #  [a0, a1, a2] coefficients
a3 = inCoef[s2][3] # set a3 to the input value
rpc = harP(Hd, Hr, beta, a3) # weights from random uncertainties
print '\nODR results on Jon\'s data no To effect\n'
rpc.pprint()

# perform odr fit for avhrr eq. keeping a3 fixed
fixb = [1,1,1,0] # fix a3 coefficient
beta0 = [-10., -4.e-3, 1.e-5, a3]
print '\nODR fit on Jon\'s data, weights from random uncertainty, keep a3 fixed\n'
podr = har.odrP(Hd, Hr, beta0, fixb) # weights from random uncertainties
print '\nODR results on Jon data with a3 kept fixed'
podr.pprint()


""" Plot fit results """
nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
noMU = Im[0,2] # number of matchups
nobj = 5000 # number of mathcups to plot

avhrrNx = upf.avhrr(nop, nos) # instance of class avhrr series 
inL = avhrrNx.measEq(Hd[:,1:6], inCoef[s2]) # input radiance

selMU = sample(range(noMU), nobj) # Select matchups for plotting

# fit coeffs, covariance and corr matrix from odr w/o To
b0 = append(rpc.beta, 0)#a3) # fit coefficients + a3
cov0 = insert(rpc.cov_beta, 3, zeros(3), axis=0) # odr evaluation matrix
cov0 = insert(cov0, 3, zeros(4), axis=1) # extended with 0 row & column for a3
cor0 = vis.cov2cor(cov0) # correlation matrix of odr fit coeffs & a3
vis.corrMap(s2, cor0, 'ODR w/o To') # heat map of correlations of fit coeffs

# plot graphs of residuals and radiance bias for odr w/o To 
LEres = rpc.eps # odr estimated error of Y var: Lref + K 
plot_ttl = s2 + ' radiance residuals from ODR fit w/o To'
vis.resLfit(LEres[selMU], inL[selMU], 'Radiance (mW/m2/sr/cm-1)', plot_ttl)
sLict = Hd[selMU, 4] # Lict for the selected matchups
vis.resLfit(LEres[selMU], sLict, 'ICT radiance (mW/m2/sr/cm-1)', plot_ttl)

cLwoTo = avhrrNx.measEq(Hd[:,1:6], b0)
cLU = avhrrNx.va2ULE(Hd[:,1:6],b0,cov0) # radiance uncertainty from coeffs uncert.

# coeffs, covariance and corr matrix from odr with fixed a3
b0odr = podr.beta # odr fit coefficients
cov0odr = podr.cov_beta # odr evaluated covariance matrix
cor0odr = vis.cov2cor(cov0odr) # correlation matrix of odr fit coeffs & a3
vis.corrMap(s2, cor0odr, 'ODR with fixed a3') # heat map of correlations of fit coeffs
cLwTo = avhrrNx.measEq(Hd[:,1:6], b0odr) # calibrated radiance from fit
cLUwTo = avhrrNx.va2ULE(Hd[:,1:6],b0odr,cov0odr) # radiance uncert. from coeffs 

plot_ttl = s2 + ' radiance residuals from ODR fit with fixed a3'
vis.resLfit(LEres[selMU], cLwTo[selMU], 'Radiance (mW/m2/sr/cm-1)', plot_ttl)
vis.resLfit(LEres[selMU], sLict, 'ICT radiance (mW/m2/sr/cm-1)', plot_ttl)


nobj = 200 # number of mathcups to plot
selMU = sample(range(noMU), nobj) # Select matchups for plotting
k = 4 # k value for sigma uncertainty range

plot_ttl = s2 + ' Radiance bias and ' + r'$4*\sigma$'+ ' uncertainty from ODR w/o To'
vis.LbiasU(inL[selMU], cLwoTo[selMU], cLU[selMU], k, plot_ttl) # ODR eval. of uncertainty
plot_ttl = s2 + ' Radiance bias and ' + r'$4*\sigma$'+ ' uncertainty from ODR with fixed a3'
vis.LbiasU(inL[selMU], cLwTo[selMU], cLUwTo[selMU], k, plot_ttl) # ODR eval. of uncertainty
