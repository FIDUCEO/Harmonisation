""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 06-10-2016
    Last update: 12-12-2016
Read matchup data from a netCDF file, create a class Hdata for har. vars, 
store harmonisation data in an object of class Hdata. 
Storing function versions for potentail use later. """

import netCDF4 as nc
import numpy as np
# import datetime
import os
import scipy.odr as odr


""" This class contains the harmonisation data as it is to be fed to the 
least squares optimisation: one row for the stochastic solution, a subset or 
the full dataset for a TLS solution."""
class Hdata(object):
    
    def __init__(self, nRow):
        """ Initialises the class with the specified number of rows """
        
        self.nor = nRow
        
""" ported from readHD with K input for possibly later use with HIRS data """
def rHDpair(folder, filename):
    print 'Setting the working directory to', folder
    os.chdir(folder) # change working directory to folder

    print 'Opening netCDF file', filename
    ncid = nc.Dataset(filename,'r')
    #print 'NetCDF variables for', filename
    #print ncid.variables
    #print 'NetCDF dimensions for', filename
    #print ncid.dimensions
    
    Im = ncid.variables['lm'][:] # matchup index array
    H = ncid.variables['H'][:,:] # harmonisation variables; empty vars included
    Ur = ncid.variables['Ur'][:,:] # random uncertainty for H vars
    Us = ncid.variables['Us'][:,:] # systematic uncertainty for H vars
    Kin = ncid.variables['K_InputData'][:] # Input to compute K values
    Kaux = ncid.variables['K_Auxilliary'][:,:] # auxiliary data to compute K 
    K = ncid.variables['K'][:] # evaluated K adjustment values
    Kr = ncid.variables['Kr'][:] # matchup random uncertainty
    Ks = ncid.variables['Ks'][:] # systematic uncertainty for K values
    corIdx = ncid.variables['CorrIndexArray'][:] # matchup time; internal format
    corLen = ncid.variables['corrData'][:] # length of averaging window

    ncid.close()
   
    """ Extract and/or return ndarrays with harmonisation data """
    #K = np.asarray(Adj) # K values and uncertainties as nparrays
    if Im[0,0] == -1: 
        rspair = 1 # reference-sensor pair
        didx = [0, 5, 6, 7, 8, 9] # non-empty columns in H, Ur, Us
        Hdata = H[:,didx] # harmonisation variables from H matrix
        Hrnd = Ur[:,didx] # random uncertainties from Ur
        Hsys = Us[:,didx] # systematic uncertainty from Us
        return rspair,Im,Hdata,Hrnd,Hsys,Kin,Kaux,K,Kr,Ks,corIdx,corLen
    else:
        rspair = 0 
        return rspair, Im, H, Ur, Us, Kin, Kaux, K, Kr, Ks, corIdx, corLen

""" ODR with total uncertatinty as vars variance, i.e. using systematic """
def pairODR(Hdata, Hrnd, Hsys, K, Krnd, civ):   
    # extract variables Cs, Cict, CE, Lict, To from Hdata
    Xvars = (Hdata[:,1:6]).transpose() # transpose data marix
    Yvar = Hdata[:,0] + K # reference radiance + adjustment values
    
    # total variance, random & systematic, for both sides of meas. eq.
    VX = (Hrnd[:,1:6] **2 + Hsys[:,1:6] **2).transpose()
    VY = Hrnd[:,0]**2 + Hsys[:,0]**2 + Krnd**2
    
    odrData = odr.Data(Xvars, Yvar, wd=1./np.sqrt(VX), we = 1./np.sqrt(VY))
    odrMod = odr.Model(avhrr4odr)
    odrFit = odr.ODR(odrData, odrMod, beta0=civ)
    mFit = odrFit.run()
    mFit.pprint()
    
    return mFit # return ODR output

""" read harmonisation data as are & feed to ODR using avhrr model function """
def pairODRnot(Hdata, Hrnd, K, Krnd, civ):
   
    # extract variables Cs, Cict, CE, Lict, To from Hdata matrix
    Xvars = Hdata[:,1:6] # transpose data marix
    Yvar = Hdata[:,0] + K # reference radiance + adjustment values
    
    # total variance, random & systematic, for both sides of meas. eq.
    VX = Hrnd[:,1:6] **2
    VY = Hrnd[:,0]**2  + Krnd**2
    
    odrData = odr.Data(Xvars, Yvar, wd=1./np.sqrt(VX), we = 1./np.sqrt(VY))
    odrMod = odr.Model(avhrr)
    odrFit = odr.ODR(odrData, odrMod, beta0=civ)
    mFit = odrFit.run()
    
    return mFit # return ODR output