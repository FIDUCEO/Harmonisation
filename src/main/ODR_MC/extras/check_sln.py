import netCDF4 as nc
import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime as dt

""" Read netcdf of a sensors' pair to harmonisation data arrays """
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
    Ks = ncid.variables['Ks'][:] # SRF uncertainty for K values
    corIdx = ncid.variables['CorrIndexArray'][:] # matchup time; internal format
    corLen = ncid.variables['corrData'][:] # length of averaging window
    # arrays of Cspace and Cict random uncertainty for 51 mav scanlines per each matchup
    CsUr = ncid.variables['cal_Sp_Ur'][:,:] # array for space counts
    CictUr = ncid.variables['cal_BB_Ur'][:,:] # array for ICT counts
    
    print '\ncorrData value for calculating pixel-to-pixel correlation', corLen[0]

    ncid.close()   

    ''' Compile ndarrays of harmonisation data '''
    nor = Im[0,2] # number of matchups in the pair
    
    if Im[0,0] == -1: # reference-sensor pair
        rspair = 1 
        
        # Extract non-empty columns: 0 [Lref], 
        # 5 [Cspace], 6 [Cict], 7 [CEarth], 8 [Lict], 9 [To]
        didx = [0, 5, 6, 7, 8, 9] # non-empty columns in H, Ur, Us
        H = H[:,didx] 
        Ur = Ur[:,didx] 
        Us = Us[:,didx] 
        
        # create data and uncertainty arrays
        noc = H.shape[1] + 1 # plus one column for K data

        # data variables
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
        
    else: # sensor-sensor pair
        rspair = 0 
        noc = H.shape[1] + 2 # plus two columns for K data and Lref
        
        Hdata = np.zeros((nor, noc)) 
        Hdata[:,1:11] = H
        Hdata[:,noc-1] = K # adjustment values in last column       

        # random uncertainties 
        Hrnd = np.zeros((nor, noc)) 
        Hrnd[:,1:11] = Ur
        Hrnd[:,noc-1] = Kr

         # systematic uncertainty
        Hsys = np.zeros((nor, noc)) 
        Hsys[:,1:11] = Us
        Hsys[:,noc-1] = Ks
       
    return rspair,Im,Hdata,Hrnd,Hsys,corIdx,corLen, CsUr, CictUr

def conv2date(inTime):
    """Return unix time in seconds (from 1970) from AVHRR time in seconds (from 1975)"""

    # Calculate difference from AVHRR start time and unix start time in seconds
    start_time_AVHRR = dt(1975, 1, 1)
    start_time_unix = dt(1970, 1, 1)
    time_diff = (start_time_AVHRR - start_time_unix).total_seconds()

    # Convert to time from 1975 to date
    outTime = [dt.fromtimestamp(time+time_diff) for time in inTime]

    return outTime

# Function to give error correlation between scan lines
# Operated on a case by case basis (does not assume inputs are arrays)
# USAGE:
#    corr_val = return_correlation(CorrIndexArray,corrData,i,j)
# where
#    CorrIndexArray  - Data from file (this name)
#    corrData        - Auxil data from file (this name)
#    i               - central scanline of averaging
#    j               - outer scanline of interest
def return_correlation(index,corr_array,cent_pos,req_pos):

    diff = abs(index[cent_pos]-index[req_pos])
    if diff > corr_array[0]:
        return 0.
    else:
        return 1.-(diff/corr_array[0])


if __name__ == '__main__':
    def main():

        filelist = ["m02_n19.nc","m02_n16.nc","m02_n15.nc"]
        datadir = "D:\Projects\FIDUCEO\Data\Simulated" # root data folder
        ncfile = filelist[2] # netCDF file to work with 

        # read data from the netCDF file
        rsp,Im,Hd,Hr,Hs,corIdx,corLen,CsUr, CictUr = rHDpair(datadir, ncfile)

        # get unique scanlines, 1st matchup (idx) and no.of matchups per scanline
        uslt,midx = np.unique(corIdx,return_index=True)
        nol = len(uslt) # number of scanlines
        sltime = conv2date(uslt) # convert to time format
        
        print nol, 'unique times in matchups of', ncfile
        
        UCs = CsUr[:,25] # Cs random uncertainty of matchup's scanline 
        UCict = CictUr[:,25] # Cict random uncertainty of matchup's scanline 
        
        # get unique values for UCs and UCict; should be one value per scanline
        # check that number of unique values is the same as unique times
        ucsu, midxcs = np.unique(UCs,return_index=True)
        ucictu, midxcict = np.unique(UCict,return_index=True)
        
        #check that both matchup index arrays for unique values are the same
        print 'Same unique indices for Cs and Cict uncert:', np.array_equal(midxcs, midxcict)
        # check that the unique indices for Cs values are within scanline times indices
        for idx in midxcs:
            if idx not in midx:
                print 'Matchup index', idx, 'not in matchup indices of unique scanlines'

        print 'Unique matchup indices of Cs unc are within the array of unique scanline times', np.in1d(midxcs, midx)
        
        # 
        print len(ucsu), 'unique values in space count uncertainty per matchup'
        # get unique arrays of 51 cal uncertainties per matchup
        uacs = np.unique(CsUr)
        print 'unique arrays of 51 scanlines Cspace uncert', len(uacs)
        
        print len(ucictu), 'unique values in ICT count uncertainty per matchup'
        uacict = np.unique(CictUr)
        print 'unique arrays of 51 scanlines Cict uncert', len(uacict)
        
        # plot scanline times 
        plt.figure() 
        plt.plot(sltime, 'bo')
        plt.title('Unique scanline times')
        plt.ylabel('Time')
        plt.show()

        # zoom in the first 30 scanlines
        plt.figure() 
        plt.plot(sltime[0:30], 'bo')
        plt.title('Times of the first 30 scanlines')
        plt.ylabel('Time (sec)')
        plt.show()
        
        # check the time difference between consecutive scanlines
        slnTdiff = list()
        for i in range(1,nol):
            slnTdiff.append((sltime[i] - sltime[i-1]).total_seconds())
        
        # plot time difference of consecutive scanlines
        plt.figure() 
        plt.plot(slnTdiff, 'go')
        plt.title('Time difference between consecutive scanlines')
        plt.ylabel('Time (sec)')
        plt.show()
        
        # get unique time differences
        slnTDs = np.unique(slnTdiff)
        
        # plot unique values of time differences
        plt.figure() 
        plt.plot(slnTDs, 'go')
        plt.title('Unique values of time difference of consecutive scanlines')
        plt.ylabel('Time (sec)')
        plt.show()
         
        # zoom in the first 20 entries
        plt.figure() 
        plt.plot(slnTDs[0:20], 'go')
        plt.title('First 20 values of time difference of consecutive scanlines')
        plt.ylabel('Time (sec)')
        plt.show()
        print 'First 20 values of unique time differences of consecutive scanlines'
        print slnTDs[0:20]
        
        return 0

    main()