#!/usr/bin/env python

""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 11-01-2017
    Last update: 05-02-2017
Perform pairwise harmonisation of sensors in a series, i.e. apply regression in
a sequential manner. Use pairLS function to perform a weighted ODR fit of 
2nd sensor to the 1st of a pair. When the 1st sensor is not a reference sensor, 
get calibration data from a previous fit in the sequence. """

from datetime import datetime as dt
import readHD as rhd
import harFun as har
import unpFun as upf

# Time the execution of the whole script
st = dt.now() # start of pair-wise harmonisation exec

datadir = "/home/ad6/Data" # in eoserver
filelist = ["m02_n19.nc","n19_n17.nc","n17_n15.nc"]
#filelist = ["n19_n15.nc","m02_n15.nc"]
slabel = 'avhrr' # series label

nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
print '\n', nos, 'sensors in the file list', slist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print '\nInput calibration coefficients:\n'
print inCoef

b0 = [-10., -4.e-3, 1.e-5, 0.0] # initialise array of cal.coeff; a3 to fix later
fixb = [1, 1, 1, 0] # keep a3 fixed for all pairs
avhrrNx = upf.avhrr(nop, nos) # instance of avhrr class
calC = {} # dictionary of sensors' calib. coefficients from ODR fits
covC = {} # dictionary of calib. coefficients covariance of sensors

for ncfile in filelist:
    # read harmonisation data: variables and uncertainties
    rsp,Im,Hd,Hr,Hs,corIdx,cLen = rhd.rHDpair(datadir, ncfile)
    print '\nNetCDF data from', ncfile, 'passed to harmonisation variables.'
        
    s2l  = 'n' + str(Im[0,1]) # label of 2nd sensor in the pair
    b0[3] = inCoef[s2l][3] # set a3 to the input value
    
    if s2l not in calC: # 2nd sensor of the pair not yet calibrated 
        if rsp: # if a reference-sensor pair
            s1l = 'm02' # label of 1st sensor in the pair

            # apply weighted ODR to data 
            podr = har.odrP(Hd, Hr, b0, fixb)
            print '\nODR output for', s2l
            podr.pprint()
            
            # add results to calC and covC dictionaries of calibration info
            calC[s2l] = podr.beta
            covC[s2l] = podr.cov_beta
        else: # if a sensor-sensor pair
            
            s1l = 'n' + str(Im[0,0]) # label of 1st sensor in the pair
            if s1l in calC: # 1st sensor calibration data in dictionary
                
                # compute Earth radiance for 1st sensor data and cal.coeff.
                LE = avhrrNx.measEq(Hd[:,0:5], calC[s1l])
                # compute Earth radiance uncertainty
                uLE = avhrrNx.uncLE(Hd[:,0:5],calC[s1l],Hr[:,0:5],covC[s1l])
                
                # fill 1st column of Hd and Hr matrix with Earth radiance data
                Hd[:,0] = LE # Earth radiance data
                Hr[:,0] = uLE # radiance uncertainty
                didx = [0, 5, 6, 7, 8, 9, 10] # data columns to use in the LS fit
                
                # perform weighted ODR fit for 2nd sensor in the pair
                ''' DEBUG: ODR takes unusual long time for the number of records
                and results are far from the input '''
                podr = har.odrP(Hd[:,didx], Hr[:,didx], b0, fixb)
                print '\nODR output for', s2l
                podr.pprint()
                print '\nInput coefficients for', s2l, '\n', inCoef[s2l]
                
                # add results to calC and covC dictionaries of calibration info
                calC[s2l] = podr.beta
                covC[s2l] = podr.cov_beta

        print '\n', s2l, "calibrated from netCDF file", ncfile

print '\nCalibration coefficients for fitted pairs'
for keys,values in calC.items():
    print(keys)
    print(values)

print '\nCovariance of calibration coefficients of fitted pairs'
for keys,values in covC.items():
    print(keys)
    print(values)


et = dt.now() # end of pair-wise harmonisation
exect = (et-st).total_seconds()
print '\nTime taken for pairwise fitting of', nos, 'sensors', slist, (exect/60.), 'minutes'
