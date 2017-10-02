#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 24-01-2017
    Last update: 04-03-2017
    Version: 8.0
Evaluate fit coefficients uncertainty via Monte Carlo (MC) trials, 
uses full error structure to generate data in each MC trial. """

import scipy.odr as odr
from numpy import empty, unique, savetxt
from datetime import datetime as dt
import readHD as rhd
import errStruct as mce 
import harFun as har

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
    print 'Current iteration coefficients:', [a0,a1,a2]
    return LE # return Earth radiance

""" Perform LS fit for a sensor-reference pair with low-level odr function """
def harP(Hdata, Hr, b0, a3):
    # extract variables Cs, Cict, CE, Lict from Hdata matrix
    X = Hdata[:,1:5].transpose() # transpose data matrix
    ''' Y is adjusted radiance, i.e. reference radiance with spectral adjustment
    and the time effect modelled by orbit temperature: Y = Lref + K - a3 * To '''
    Y = Hdata[:,0] + Hdata[:,6] - a3 * Hdata[:,5]

    VX = (Hr[:,1:5]**2).transpose() # sigma^2 of X
    VY = Hr[:,0]**2 + Hr[:,6]**2 + Hr[:,5]**2 # sigma^2 of Y 
    
    # perform odr fit (low level function)
    fit = odr.odr(fcnP,b0,Y,X,we=1./VY,wd=1./VX,full_output=1)
    odrFit = odr.Output(fit) # get regression output     
    
    return odrFit # return odr output

""" Perform ODR over data MC generated with full error structure """
def odr4MC(Xdata, Ydata, Hr, b0):

    X = Xdata.transpose()
    VX = (Hr[:,1:5]**2).transpose() # sigma^2 of X
    ''' Y = Lref + K - a3 * To '''
    # cacluate weights from uncertainty matrices
    VY = Hr[:,0]**2 + Hr[:,6]**2 + Hr[:,5]**2 # sigma^2 of Y 
        
    #  ODR on new X,Y data, perturbed best estimates  
    fit = odr.odr(fcnP,b0,Ydata,X,we=1./VY,wd=1./VX,full_output=1)
    odrFit = odr.Output(fit) # get odr fit output   
      
    return odrFit # return odr output

st = dt.now() # start of MS run

filelist = ["m02_n19.nc","m02_n18.nc","m02_n17.nc","m02_n16.nc","m02_n15.nc"]
datadir = "/home/ad6/Data" # in eoserver
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
ncfile = filelist[4] # netCDF file to work with 

nop = len(filelist) # number of sensor pairs
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
s2 = ncfile[4:7]
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print '\nCalibrating sensor', s2, 'against the reference'
print '\nInput coefficients for', s2, ':', inCoef[s2]

# read data from the netCDF file
rsp,Im,Hd,Hr,Hs,corIdx,corLen,csUr,cictUr = rhd.rHDpair(datadir, ncfile)
# set systematic uncertainties equivalent to Peter's LS optimisation
Hs = har.resetHs(Hs, rsp) 
print Im[0,2], 'matchup data from', ncfile, 'passed to harmonisation matrices'


# perform odr fit without To effect and print fit result
print '\nODR fit on Jon\'s data, weights from random uncertainty, no To effect\n'
beta = [-10., -4.e-3, 1.e-5] #  [a0, a1, a2] coefficients
a3 = inCoef[s2][3] # set a3 to the input value
rpc = harP(Hd, Hr, beta, a3) # weights from random uncertainties
print '\nODR results on Jon\'s data no To effect\n'
rpc.pprint()

''' Generate data for Monte Carlo run ''' 
bodr = rpc.beta # odr fit coefficients
Y = rpc.y # best est.of adjusted reference radiance: Lref + K - a3 * To
X = rpc.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict
sLict = Hs[0,4] # systematic error Lict
sTo = Hs[0,5] # systematic error To

# get scanlines and matchups per scanline
slt,midx,mcnt = unique(corIdx,return_index=True,return_counts=True)
rCSar = csUr[midx,:] # Cspace random uncert. per scanline: arrays of 51 slines
rCICTar = cictUr[midx,:] # Cict random uncert. per scanline: arrays of 51 slines

# MC runs ODR on new data: best estimate + full correlation error draw
notr = 500 # number of MC trials
p = len(beta) # number of calib. coefficients
mcb = empty([notr, p+1]) # array to store beta vals and ODR info of MC trials
print '\n\nRun MC with the full error structure.'

for i in range(notr):
    ''' Run notr MC trials '''

    # compile data for the ODR run; generate errors 
    errStr = mce.genErr(Hr, sLict, sTo, rCSar, rCICTar, slt, corLen, mcnt)
    # add errStr to X & Y best estimates
    Xdt = X.T + errStr[:,1:5] # X variables
    Ydt = Y + errStr[:,0] + errStr[:,5] + errStr[:,6] # Y variable
    
    # run ODR on new X & Y vals and Hr weights 
    modr = odr4MC(Xdt, Ydt, Hr, bodr)
    
    # store ODR fit coefficients and reason for halting
    mcb[i, 0:p] = modr.beta
    mcb[i,p] = modr.info

print '\n\nODR results from the last MC trial'
modr.pprint()

fn = s2 + '_mcerrPS_beta.txt'
savetxt(fn, mcb, delimiter=',')

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print '\n\n\n--- Time taken for', notr, 'MC trials', (exect/60.), 'minutes ---'
