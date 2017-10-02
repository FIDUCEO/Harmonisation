""" FCDR harmonisation modules
    Project:      H2020 FIDUCEO
    Authors:      Arta Dilo, Peter Harris /NPL MM, Sam Hunt, Jon Mitazz /NPL ECO
    Date created: 02-02-2017
    Last update:  20-04-2017
    Version:      12.0
Generate errors respecting the full correlation structure for a sensor series. 
An error matrix E is generated from genErr function which has the same dimension 
of the H matrix i.e. harmonisation variables through all matchups. 

An error matrix E (computed for each MC trial) is added to best estimates of 
explanatory and dependent (harmonisation) variables; an MC trial from harODRMC
performs ODR regression over these generated data and returns the calibration 
coefficients. """

from scipy.sparse import csr_matrix
from numpy import zeros, ones, trim_zeros, arange, repeat
import numpy.random as random
#import numpy.ma as ma
from readHD import HblockIdx as block


""" Peter & Sam's computation of moving average error with the weight matrix """
def calc_CC_err(u, times, corrData):
    """ :param u: float
        standard uncertainties
    :param times: numpy.ndarray
        match-up times for match-ups in match-up series
    :param corrData: numpy.ndarray
        match-up time data

    :return:
        :CC_err: numpy.ndarray
            error for averaged calibration counts  """

    n_var = len(times)  # number of match_ups
    N_W = len(u[0])     # length of maximum averaging kernel

    # initialise sparse matrix index and values arrays (of maximum size, i.e. if all windows are n_w)
    ir = zeros(n_var * N_W)
    jc = zeros(n_var * N_W)
    ws = zeros(n_var * N_W)

    col = 0  # column number
    iend = 0
    for i in xrange(n_var):
        ui = u[i][u[i] != 0]  # scanline non-zero uncertainties
        n_w = len(ui)  # width of match-up specific averaging window

        # find col_step of match-up compared to last match-up (if not first match-up)
        col_step = 0
        if i > 0:
            corr_val = return_correlation(times, corrData, i, i - 1)
            col_step = int(round(n_w * (1 - corr_val)))
        col += col_step

        # fill sparse matrix index and value arrays
        istart = iend
        iend = istart + n_w

        ir[istart:iend] = ones(n_w) * i
        jc[istart:iend] = arange(n_w) + col
        ws[istart:iend] = ui * ones(n_w) / n_w

    # trim off trailing zeros if maximum size not required, i.e. if all windows are not n_w in length
    ir = trim_zeros(ir, trim='b')
    jc = trim_zeros(jc, trim='b')
    ws = trim_zeros(ws, trim='b')

    # build sparse matrix
    W = csr_matrix((ws, (ir, jc)))

    # generate raw scanline errors (uncertainy normalised to 1)
    CC_raw_err = random.normal(loc=zeros(W.indices[-1]+1))

    # average raw errors to generate CC_err (scaling by raw CC uncertainty)
    CC_err = W.dot(CC_raw_err)

    return CC_err

""" Jon's function for correlation in moving averages. """
def return_correlation(index, corr_array, cent_pos, req_pos):
    """ Function to give error correlation between scan lines
    Operated on a case by case basis (does not assume inputs are arrays)

    :param index: numpy.ndarray
        CorrIndexArray data from file (this name)
    :param corr_array: numpy.ndarray
        corrData auxiliary data from file (this name)
    :param cent_pos: int
        central scanline of averaging
    :param req_pos: int
        outer scanline of interest

    :return
        :corr_val:
            correlation between scanlines  """

    diff = abs(index[cent_pos] - index[req_pos])
    if diff > corr_array[0]:
        return 0.
    else:
        return 1. - (diff / corr_array[0])


""" Generate matrix of errors for multiple pairs of the series to be harmonised.
The error matrix has the same shape as the harmonisation data matrix i.e. 
12 columns, and number of rows equal to the number of matchups for all pairs. 
The columns are 
- errors for reference radiance (Lref); 
- five measured variables for each sensor in the pair: 
    - space count (Cs)
    - ICT count (Cict)
    - Earth count (CE)
    - ICT radiance (Lict) 
    - orbit temperature (To), when working with time-dependent meas. model
- spectral adjustment K. 

The error matrix consists of blocks (much in the same way as all the other
harmonisation matrices) according to the type of pair, reference-sensor or 
sensor-sensor. A reference-sensor pair has 0 values for all variables of the 1st 
sensor, a sensor-sensor pair has o values in the Lref column for all the matchups
of that pair. 
The current code reads each block of matchups for a pair and processes according 
to the pair type. Probably a better way would be to use masked arrays. """
def genErr(Hr,Hsys, CsUr,CbbUr,CsU,CbbU, Is,slT1,muCnt1,slT2,muCnt2,clen, series):
    
    notd = series.notime # True/False indicating time-dependency of meas.model
    nop = series.nopairs # number of sensor pairs 
    noM = series.im[:,2] # number of matchups per pair
    err = zeros(Hr.shape) # matrix of errors
    
    # 2nd sensor calibration counts errors from moving average, 
    # i.e. Cs and Cict; use weight matrix W to generate errors 
    Cs_err = calc_CC_err(CsU, slT2, clen) # Space count error per scanline 
    err[:,6] = repeat(Cs_err, muCnt2)  # Cs errors per matchup        
    Cict_err = calc_CC_err(CbbU, slT2, clen) # ICT count error per scanline
    err[:,7] = repeat(Cict_err, muCnt2) # Cict error per matchup
    
    # CE 2nd sensor: random error from Gaussian with sigma from Hr &mu=0
    err[:,8] = random.normal(scale=Hr[:,8]) # CE random error of 2nd sensor

    # Lict 2nd sensor: random error; to combine with systematic (if not 0)  
    err[:,9] = random.normal(scale=Hr[:,9]) # random error Lict 2nd sensor        
    
    # To 2nd sensor: systematic and random error for time-dependent model/data
    if not notd: 
        errT = random.normal(scale=Hsys[:,3]) # To systematic error per pair
        err[:,10] = random.normal(scale=Hr[:,10]) + repeat(errT, noM)       

    # K error: random from Gaussian with sigma from Hr & mu=0
    err[:,11] = random.normal(scale=Hr[:,11]) # K random error    

    # compute other errors /error components per block 
    for i in range(nop): # loop through sensor pairs
        
        # get start and end matchup index of a block
        sMidx, eMidx = block(series.im, series.im[i,0], series.im[i,1]) 
        
        if series.im[i,0] == -1: # reference-sensor pair
            # generate Lref random error for the block
            err[sMidx:eMidx,0] = random.normal(scale=Hr[sMidx:eMidx,0])
            
        else: # sensor-sensor pair; compute errors for 1st sensor variables
            
            # get start and end scanlines for the block
            sSidx, eSidx = block(Is, Is[i,0], Is[i,1]) 
            
            # 1st sensor Cs errors per matchup
            Cs_err = calc_CC_err(CsUr[sSidx:eSidx,:], slT1[sSidx:eSidx], clen)
            err[sMidx:eMidx,1] = repeat(Cs_err, muCnt1[sSidx:eSidx])  
            # 1st sensor Cict error per matchup
            Cict_err = calc_CC_err(CbbUr[sSidx:eSidx,:], slT1[sSidx:eSidx], clen)
            err[sMidx:eMidx,2] = repeat(Cict_err, muCnt1[sSidx:eSidx]) 

             # CE random error for 1st sensor
            err[sMidx:eMidx,3] = random.normal(scale=Hr[sMidx:eMidx,3])
            
            # Lict random error for 1st sensor
            err[sMidx:eMidx,4] = random.normal(scale=Hr[sMidx:eMidx,4]) 
            if Hsys[i,0]: # add Lict systematic component if not 0
                errL = random.normal(scale=Hsys[i,0]) 
                err[sMidx:eMidx,4] += repeat(errL, noM[i]) 
                
            if not notd: # time-dependent data & model; add To error
                err[sMidx:eMidx,5] = random.normal(scale=Hr[sMidx:eMidx,5]) 
                errT = random.normal(scale=Hsys[i,1]) # systematic error
                err[sMidx:eMidx,5] += repeat(errT, noM[i]) # combined To error
        
        # 2nd sensor Lict systematic error component; add to random errors
        if Hsys[i,2]: # if not zero
            errL = random.normal(scale=Hsys[i,2]) # Lict systematic error
            err[sMidx:eMidx,9] += repeat(errL, noM[i]) # combined Lict error
                        
    return err
