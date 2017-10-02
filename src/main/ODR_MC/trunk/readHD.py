""" FCDR harmonisation modules
    Project:        H2020 FIDUCEO
    Author:         Arta Dilo \NPL MM
    Reviewer:       Peter Harris \NPL MM, Sam Hunt \NPL ECO
    Date created:   06-10-2016
    Last update:    20-04-2017
    Version:        12.0
Read matchup data from netCDF files of sensors pairs, 
create matrices of variables to give to the harmonisation. 
The read functions use netCDF4 package to work with netCDF files, 
and numpy ndarrays to store harmonisation data. """

from netCDF4 import Dataset
from numpy import append, array, empty, ones, zeros, where, unique
from numpy import mean as npmean
from os.path import join as pjoin
from csv import reader as csvread


""" Extract sensors from a list of netCDF files """
def satSens(nclist):
   fsl = [] # list of satellite sensor labels in nclist
   for netcdf in nclist:
       s1 = netcdf[0:3]  # label of 1st sensor 
       s2 = netcdf[4:7] # label of 2nd sensor 
       
       if s1 not in fsl: # if not yet in the list of satellites
           fsl.append(s1)  # add to the list    
       if s2 not in fsl:
           fsl.append(s2)
   return fsl

""" Get input calibration coefficients of sensors in netCDF filelist """
def sInCoeff(csvfolder, nclist, notd):
   fsl = satSens(nclist) # list of sensors in nc files 
   
   if notd: # known calib. coeffs for non-/ time dependency
       csvname = 'CalCoeff_notd.csv'
   else:
       csvname = 'CalCoeff.csv'
   cfn = pjoin(csvfolder, csvname)
   
   inCc = {} # dictionary of fsl sensors input coefficients
   with open(cfn, 'rb') as f:
       reader = csvread(f)
       reader.next() # skip header
       for row in reader:
           sl = row[0] # sensor label
           coefs = array(row[1:5]).astype('float')
           if sl in fsl:
               inCc[sl] = coefs
   return inCc


""" Reads netCDF file of a sensor pair to harmonisation data arrays. 
Function arguments: folder where netCDF file is; netCDF filename. 
TODO: simplify code by keeping the same matrix size for both types of pairs.  
The function returns:
    - rspair: boolean; 1 for reference-sensor pair
    - ...    """
def rHDpair(folder, filename):     
    pfn = pjoin(folder, filename) # filename with path  
    print '\nOpening netCDF file', pfn
    ncid = Dataset(pfn,'r')
    
    Im = ncid.variables['lm'][:] # matchup index array
    H = ncid.variables['H'][:,:] # harmonisation variables; empty vars included
    Ur = ncid.variables['Ur'][:,:] # random uncertainty for H vars
    Us = ncid.variables['Us'][:,:] # systematic uncertainty for H vars
    K = ncid.variables['K'][:] # evaluated K adjustment values
    Kr = ncid.variables['Kr'][:] # matchup random uncertainty
    Ks = ncid.variables['Ks'][:] # SRF uncertainty for K values
    # scanline time 2nd sensor in internal format
    corIdx = ncid.variables['CorrIndexArray'][:]
    rcIdx = ncid.variables['ref_CorrIndexArray'][:] # scanline time 1st sensor
    corLen = ncid.variables['corrData'][:] # length of averaging window
    # arrays of Cspace and Cict random uncertainty of surrounding 51 scanlines 
    RCs1 = ncid.variables['ref_cal_Sp_Ur'][:,:] # 1st sensor Cs random uncert's
    RCict1 = ncid.variables['ref_cal_BB_Ur'][:,:] # 1st sensor Cict random uncert's
    RCs2 = ncid.variables['cal_Sp_Ur'][:,:] # 2nd sensor Cs random uncert's
    RCict2 = ncid.variables['cal_BB_Ur'][:,:] # 2nd sensor Cict random uncert's
    
    #print '\ncorrData value for calculating pixel-to-pixel correlation', corLen[0]
    ncid.close()   

    ''' Compile ndarrays of harmonisation data '''
    nor = Im[0,2] # number of matchups in the pair
    
    if Im[0,0] == -1: # reference-sensor pair
        rspair = 1 
        
        # Extract non-empty columns in data matrices, H, Ur, Us
        # 0 [Lref], 5 [Cspace], 6 [Cict], 7 [CEarth], 8 [Lict], 9 [To]
        didx = [0, 5, 6, 7, 8, 9] # non-empty columns in H, Ur, Us
        H = H[:,didx] 
        Ur = Ur[:,didx] 
        Us = Us[:,didx] 
        
        # create data and uncertainty arrays
        noc = H.shape[1] + 1 # plus one column for K data
        # data variables
        Hdata = zeros((nor, noc)) 
        Hdata[:,:-1] = H
        Hdata[:,noc-1] = K # adjustment values in last column        
        # random uncertainties 
        Hrnd = zeros((nor, noc)) 
        Hrnd[:,:-1] = Ur
        Hrnd[:,noc-1] = Kr
        # systematic uncertainty
        Hsys = zeros((nor, noc)) 
        Hsys[:,:-1] = Us
        Hsys[:,noc-1] = Ks
        
    else: # sensor-sensor pair
        rspair = 0 
        noc = H.shape[1] + 2 # plus two columns for K and Lref
        
        Hdata = zeros((nor, noc)) 
        Hdata[:,1:11] = H
        Hdata[:,noc-1] = K # adjustment values in last column       
        # random uncertainties 
        Hrnd = zeros((nor, noc)) 
        Hrnd[:,1:11] = Ur
        Hrnd[:,noc-1] = Kr
         # systematic uncertainty
        Hsys = zeros((nor, noc)) 
        Hsys[:,1:11] = Us
        Hsys[:,noc-1] = Ks
        
    # set systematic uncert. conform Peter and Sam's method
    Hsys = resetHs(Hsys, rspair) 
        
    # get unique scanlines, idx of 1st matchup pixel, no of matchups per scanline
    # for 2nd sensor
    slt,midx,mcnt = unique(corIdx,return_index=True,return_counts=True) 
    # for 1st sensor
    slt1,midx1,mcnt1 = unique(rcIdx,return_index=True,return_counts=True)     

    ''' Use midx to extract unique scanline rows from the large arrays of 
    random calibration count uncertainties of 51 sorrunding scanlines, e.g.
    RCs[midx,:] - 51 Cspace random uncertainties per scanline
    RCict[midx,:] - 51 Cict random uncertainties per scanline.
    slt array needed for computing errors, corIdx returned potentially for 
    plotting on time and further checks. '''

    Is = zeros((1,4), dtype=int) # index matrix of scanlines per sensor in a pair
    Is[0,0:2] = Im[0,0:2] # copy sensor indices from Im matrix
    Is[0,3] = len(slt)  # number of scanlines for 2nd sensor

    if rspair: 
        Is[0,2] = 0 # number of scanlines for 1st sensor
        slt1= zeros(0) # number of scanlines 1st sensor
        mcnt1 = zeros(0, dtype=int) # number of matchups 1st sensor
        UrCs1 = zeros([0,51], order = 'C') # Cs uncert. arrays
        UrCict1 = zeros([0,51], order = 'C') # Cict uncert. arrays
        
    else: # TODO: change for corresponding scanlines of 1st sensor
        Is[0,2] = len(slt1) # 1st sensor scanlines
        UrCs1 = RCs1[midx1,:] # change unique indices for corresponding scanlines
        UrCict1 = RCict1[midx1,:]

    #print 'Im matrix of matchups per pair', Im
    #print 'Is matrix of scanlines per pair', Is
    
    return rspair,Im,Hdata,Hrnd,Hsys,corIdx,corLen,Is,slt1,mcnt1,slt,mcnt,UrCs1,UrCict1,RCs2[midx,:],RCict2[midx,:]


""" Create harmonisation data from netcdf files in the list.
Data matrix with fixed number of columns (12) and number of rows equal to the 
number of matchup observations for multiple pairs of a series. """          
def rHData(folder, nclist):
    fsl = satSens(nclist) # list of sat.sensors; same order as coeffs array

    # initialise harmonisation data arrays 
    Im = empty([0,3], dtype=int) # matrix index of matchups per pair
    Is = empty([0,4], dtype=int) # index of scanlines per pair
    Hdata = empty([0,12], order = 'C')
    Hrnd = empty([0,12], order = 'C')
    Hsys = empty([0,4], order = 'C')
    psens = empty([0,2], dtype=int) # sensor labels in a pair
    UCs1 = empty([0,51], order = 'C') # Cs random uncert. arrays 1st sensor
    UCict1 = empty([0,51], order = 'C') # Cict random uncert. arrays 1st sensor
    UCs2 = empty([0,51], order = 'C') # Cs random uncert. arrays 1st sensor
    UCict2 = empty([0,51], order = 'C') # Cict random uncert. arrays 1st sensor
    slTs1 = empty([0]) # unique scanlines 1st sensor
    mxSl1 = empty([0], dtype=int) # number of matchups per scanline 1st sensor
    slines = empty([0]) # unique scanlines 2nd sensor
    muXsl = empty([0], dtype=int) # number of matchups per scanline 2nd sensor
    tidx = empty([0]) # array of scanline time for each matchup

    # loop through the netCDFs filelist 
    for ncfile in nclist:
        # read harmonisation data, variables and uncertainties
        rsp,lm,Hd,Hr,Hs,corIdx,cLen,ls,slt1,mxsl1,slt2,mxsl2,CsU1,CbbU1,CsU2,CbbU2 = rHDpair(folder, ncfile)
        print 'NetCDF data from', ncfile, 'passed to harmonisation arrays.'
        
        nor = lm[0,2] # no of matchups for the current pair
        tH = zeros((nor, 12))
        tHr = ones((nor, 12)) 
        tps = zeros((nor, 2), dtype=int) 
        if rsp: # reference-sensor pair
            tH[:,0] = Hd[:,0] # reference radiance
            tH[:,1].fill(1.) # space counts 1st sensor
            tH[:,3].fill(1.) # Earth counts 1st sensor
            tH[:,6:12] = Hd[:,1:7] # 2nd sensor calibration data
            tHr[:,0] = Hr[:,0] # reference radiance random uncertainty
            tHr[:,6:12] = Hr[:,1:7] # 2nd sensor random uncertainty data
            sl1 = 'm02'# 1st sensor label
        else: # sensor-sensor pair
            tH = Hd # data variables 
            tHr[:,1:12] = Hr[:,1:12] # 1st and 2nd sensor random uncertainty 
            if lm[0,0] < 10:
                sl1 = 'n0' + str(lm[0,0]) # 1st sensor label
            else:
                sl1 = 'n' + str(lm[0,0]) # 1st sensor label            

        sidx1 = fsl.index(sl1) # 1st sensor index in the list of sensors
        if lm[0,1] < 10:
            sl2 = 'n0' + str(lm[0,1]) # 2nd sensor label
        else:
            sl2 = 'n' + str(lm[0,1])    
        sidx2 = fsl.index(sl2) # 2nd sensor idx in list of sensors
        
        tps[:,0] = sidx1 # 1st sensor's index in coeffs array 
        tps[:,1] = sidx2 # 2nd sensor's index
        
        # add data of the current pair to harmonisation arrays
        Im = append(Im, lm, axis=0)
        Is = append(Is, ls, axis=0)
        Hdata = append(Hdata, tH, axis=0)
        Hrnd = append(Hrnd, tHr, axis=0)
        Hsys = append(Hsys, Hs, axis=0)
        psens = append(psens, tps, axis=0)
        UCs1 = append(UCs1, CsU1, axis=0)
        UCict1 = append(UCict1, CbbU1, axis=0)
        UCs2 = append(UCs2, CsU2, axis=0)
        UCict2 = append(UCict2, CbbU2, axis=0)
        slTs1 = append(slTs1, slt1)
        mxSl1 = append(mxSl1, mxsl1)
        slines = append(slines, slt2)
        muXsl = append(muXsl, mxsl2)
        tidx = append(tidx, corIdx)
        
    return Im,Hdata,Hrnd,Hsys,psens,tidx,cLen,Is,slTs1,mxSl1,slines,muXsl,UCs1,UCict1,UCs2,UCict2


""" Reset the systematic error matrix to constant values so that  
ODR problem setting corresponds to Peter's LS optimisation problem. 
- Lref, Cspace, Cict, CEarth have 0 systematic uncertainty
- Lict and To have a constant systematic uncertainty 
K random uncertainty from SRF shifting is stored in the Hs matrix; ignored! """
def resetHs(Hs, rspair):
    Hsys = zeros([1,4]) # 2 sensors * 2 vars having systematic uncertainty
        
    if rspair: # reference sensor pair, i.e. Hs has 7 columns     
        # set Lict and To systematic to mean of corresponding Hs column
        # columns [4] and [5] are respectively Lict and To of series sensor in 
        sLict = npmean(Hs[:,4]) # mean Lict through all matchups
        sTo = npmean(Hs[:,5])  # mean To 
        Hsys[0,2] = sLict
        Hsys[0,3] = sTo

    else: # sensor-sensor pair, Hs has 12 columns
        # set Lict and To systematic for 1st sensor
        sLict = npmean(Hs[:,4]) # mean Lict through all matchups
        sTo = npmean(Hs[:,5])  # mean To 
        Hsys[0,0] = sLict
        Hsys[0,1] = sTo
        
        # set 2nd sensor' Lict and To systematic to mean of Hs values
        sLict = npmean(Hs[:,9]) 
        sTo = npmean(Hs[:,10]) 
        Hsys[0,2] = sLict
        Hsys[0,3] = sTo
    # TODO: set the same value for a sensor in each pair it appears
       
    return Hsys # return the new set of sytematic uncertainties


""" Extract start and end row indices of the 1st block of matchups for a sensor. 
Arguments 
- Im: matchups index array
- slab: 2 digits label of a sensor
- cidx: column in Im, i.e. 1st (cidx=0) or 2nd (cidx=1) sensor in a pair """
def sliceHidx(Im, slab, cidx=1):
    
    bidx = where(Im[:,cidx]==slab) # entries in Im for sensor slab
    sidx = bidx[0][0] # first occurance of slab in Im 
    
    startMU = 0 # start index of sensor matchups
    for i in range(sidx): # loop through Im rows until 1st slab occurance
        startMU += Im[i,2]
    endMU = startMU + Im[sidx,2]
    
    return startMU, endMU # return indices of first and last matchup record 

""" Extract start and end row indices for the block of a sensor pair. 
Arguments:
- IM: index matrix of matchups or scanlines for sensor pairs
- s1: 2 digits label of sensor 1
- s2: 2 digits label of sensor 2 """
def HblockIdx(IM, s1, s2):
    
    bidx = where((IM[:,0]==s1) & (IM[:,1]==s2)) # entry in IM for the pair
    sidx = bidx[0][0] # read index value in the 2D array
    
    startR = 0 # start index of matchups/scanlines block
    for i in range(sidx): # loop through IM rows until the entry of the pair
        startR += IM[i,2]
    endR = startR + IM[sidx,2]
    
    return startR, endR # return indices of first and last record 
