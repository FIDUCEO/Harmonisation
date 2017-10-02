import numpy as np
import scipy.odr as odr
import netCDF4 as nc
import csv
import os

def rHDpair(folder, filename):     
    pfn = os.path.join(folder, filename) # filename with path  
    print 'Opening netCDF file', pfn
    ncid = nc.Dataset(pfn,'r')
    
    Im = ncid.variables['lm'][:] # matchup index array
    H = ncid.variables['H'][:,:] # harmonisation variables; empty vars included
    Ur = ncid.variables['Ur'][:,:] # random uncertainty for H vars
    Us = ncid.variables['Us'][:,:] # systematic uncertainty for H vars
    K = ncid.variables['K'][:] # evaluated K adjustment values
    Kr = ncid.variables['Kr'][:] # matchup random uncertainty
    Ks = ncid.variables['Ks'][:] # systematic uncertainty for K values
    corIdx = ncid.variables['CorrIndexArray'][:] # matchup time; internal format
    corLen = ncid.variables['corrData'][:] # length of averaging window

    ncid.close()   
    # Extract non-empty columns
    if Im[0,0] == -1: 
        rspair = 1 # reference-sensor pair
        didx = [0, 5, 6, 7, 8, 9] # non-empty columns in H, Ur, Us
        H = H[:,didx] 
        Ur = Ur[:,didx] 
        Us = Us[:,didx] 
    else:
        rspair = 0 # sensor-sensor pair
    # compile ndarrays of harmonisation data
    nor = Im[0,2]
    noc = H.shape[1] + 1 # plus one column for K's
    # harmonisation variables
    Hdata = np.zeros((nor, noc)) 
    Hdata[:,:-1] = H
    Hdata[:,noc-1] = K # adjustment values in last column
    # random uncertainties 
    Hrnd = np.zeros((nor, noc)) 
    Hrnd[:,:-1] = Ur
    Hrnd[:,noc-1] = Kr
    # systematic uncertainty
    Hsys = np.zeros((nor, noc)) 
    Hsys[:,:-1] = Us
    Hsys[:,noc-1] = Ks
    
    return rspair,Im,Hdata,Hrnd,Hsys,corIdx,corLen

def avhrr(coef, data):
    a0 = coef[0] # AVHRR model coefficients
    a1 = coef[1]
    a2 = coef[2]
    a3 = coef[3]    
    # transposed ndarrays 
    Cs = data[0,:] # space counts
    Cict = data[1,:] # ICT counts
    CE = data[2,:] # Earth counts
    Lict = data[3,:] # ICT radiance
    To = data[4,:] # orbit temperature * input a3
    
    LE = a0 + (0.98514+a1)*Lict*(Cs-CE)/(Cs-Cict) + a2*(Cict-CE)*(Cs-CE) + a3*To
    #print 'Current iteration coefficients:', [a0,a1,a2, a3]
    return LE  
    
""" read harmonisation data and feed to ODR """
def pairODR(Hdata, Hrnd, civ):
   
    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    Xvars = Hdata[:,1:6].transpose() # transpose data matrix
    Yvar = Hdata[:,0] + Hdata[:,6] # reference radiance + adjustment values
    # uncertainties in X and Y variables; to use as initail error values
    UX = Hrnd[:,1:6].transpose()
    UY = np.sqrt(Hrnd[:,0]**2 + Hrnd[:,6]**2)

    # weights to be defined differently
    odrData = odr.Data(Xvars, Yvar, wd=1./UX, we = 1./UY)
    odrMod = odr.Model(avhrr)
    odrFit = odr.ODR(odrData, odrMod, beta0=civ, ifixb=[1,1,1,0])
    mFit = odrFit.run()
    
    return mFit # return ODR output

# sensors in the netCDF file from the list 
def sensors(nclist):
   fsl = [] # list of sensor labels in nclist
   for netcdf in nclist:
       s1 = netcdf[0:3]
       s2 = netcdf[4:7]
       if s1 in fsl:
           pass
       else:
           fsl.append(s1)
       if s2 in fsl:
           pass
       else:
           fsl.append(s2)
   return fsl

# input calibration coefficients of sensors from nclist
def sInCoeff(csvfolder, nclist):
   fsl = sensors(nclist) # list of sensors in nc files 
   
   inCc = {} # dictionary of fsl sensors input coefficients
   cfn = os.path.join(csvfolder, 'CalCoeff.csv')
   with open(cfn, 'rb') as f:
       reader = csv.reader(f)
       reader.next() # skip header
       for row in reader:
           sl = row[0] # sensor label
           coefs = np.array(row[1:5]).astype('float')
           if sl in fsl:
               inCc[sl] = coefs
   return inCc

datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
filelist = ['m02_n15.nc']#, 'm02_n19.nc', 'm02_n16.nc']
slist = sensors(filelist) # list of sensors in filelist
inCoef = sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print 'Sensors in the file list', slist
print 'Input calibration coefficients', inCoef

initC = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
for ncfile in filelist:
    s2 = ncfile[4:7]
    initC[3] = inCoef[s2][3] # set a3 to the input value
    
    # read harmonisation data: variables and uncertainties
    rsp,Im,Hd,Hr,Hs,corIdx,corLen = rHDpair(datadir, ncfile)
    print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
    # apply odr to data using random uncertainty only, fix a3
    podr = pairODR(Hd, Hr, initC)
    print 'ODR output for', s2
    podr.pprint()
