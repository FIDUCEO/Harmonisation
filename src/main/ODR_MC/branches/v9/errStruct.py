""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo, Peter Harris / NPL MM, Sam Hunt, Jon Mittaz / NPL ECO
    Date created: 02-02-2017
    Last update: 16-03-2017
Generate data for Monte Carlo uncertainty evaluation of odr fit coefficients. """

from numpy import empty, repeat, array, concatenate
from datetime import datetime as dt
from scipy.sparse import csr_matrix
from numpy import zeros, ones, trim_zeros, arange
import numpy.random as random


def calc_CC_err(u, times, corrData):
    """ Peter & Sam's weight matrix:
    Return weighting matrix: sparse representation.

    :param u: float
        standard uncertainties
    :param times: numpy.ndarray
        match-up times for match-ups in match-up series
    :param corrData: numpy.ndarray
        match-up time data

    :return:
        :CC_err: numpy.ndarray
            error for averaged calibration counts
    """

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


def return_correlation(index, corr_array, cent_pos, req_pos):
    """ Jon's function for correlation in moving averages. 
    Function to give error correlation between scan lines
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
            correlation between scanlines

    """

    diff = abs(index[cent_pos] - index[req_pos])
    if diff > corr_array[0]:
        return 0.
    else:
        return 1. - (diff / corr_array[0])


''' Generate the matrix of errors using W matrix from Sam & Peter '''
def genErr(Hr, Lsys, Tsys, uCs, uCict, slTidx, clen, muCnt, notime):
    err = zeros(Hr.shape) # matrix of errors
    nor = err.shape[0] # number of matchups
    v1 = ones(nor) # array of ones with size no. of matchups
    
    # Lref, K, CE: random error from Gaussian with sigma from Hr data &mu=0
    err[:,0] = random.normal(scale=Hr[:,0]) # Lref random error
    err[:,6] = random.normal(scale=Hr[:,6]) # K random error
    err[:,3] = random.normal(scale=Hr[:,3]) # CE random error
    
    # Run Sam's function calc_CC_err to generate Cspace averaged errors per sline
    Cs_err = calc_CC_err(uCs, slTidx, clen)
    err[:,1] = repeat(Cs_err, muCnt)  # averaged Space count error per matchup
    
    # Sam's function calc_CC_err to generate Cict averaged errors per scanline
    Cict_err = calc_CC_err(uCict, slTidx, clen)
    err[:,2] = repeat(Cict_err, muCnt) # averaged ICT count error per matchup
    
    # combined systematic and random errors for Lict and To
    err[:,4] = random.normal(scale=Hr[:,4]) # random error Lict
    if notime:
        err[:,5] = 0 # To error 0; this column will not be used
    else:   # time dependant data
        
        # add systematic error to Lict random error
        errL = random.normal(scale=Lsys) # Lict systematic error
        err[:,4] += errL*v1 # Lict error
    
        # combined systematic and random To error
        errT = random.normal(scale=Tsys) # To systematic error
        err[:,5] = random.normal(scale=Hr[:,5]) + errT*v1 # To error
            
    return err

    
""" Group contiguous scanlines (possibly with gaps) in blocks such that two
blocks are away from each other more than half the averaging window, i.e. 
25 scanlines. Moving average is applied within a block. """
def groupSln(sltidx, sltd, cLen):
    nol = len(sltidx) # number of scanlines
    sltime = conv2date(sltidx) # convert scanline time idx to time format
    
    # create array of scanline idx and distance to next scanline in the array
    # gap = the number of lines that current scanline is away from the previous 
    # scanline, e.g. gap = 1 for consecutive scanlines in the array sltidx,
    # gap = 1 at the start of a new block.
    gaps = ones((nol,2)) # array to store scanline time idx and gaps 
    gaps[:,0] = sltidx # scanline time idx in first column
    
    # create blocks of scanlines that are >25 lines apart
    # first column stores the block number, start at 0 increase by 1;
    # second column stores the index of the 1st scanline in the block
    blocks = zeros((1,2), dtype=long) # 1st element: group=0, scanline idx=0
    bno = blocks[0, 0] # block number of 1st block = 0
    
    for i in range(1, nol):
        gaps[i,1] = (sltime[i] - sltime[i-1]).total_seconds()/sltd
        #print 'gap[',i,'] =',gaps[i,1]
        if gaps[i,1] > cLen: # 25: half window size in scanlines number
            gaps[i,1] = 1 # start of new block 
            #print 'gap[',i,'] =',gaps[i,1]
            bno += 1 # increase block number 
            brow = array([[bno, i]]) # new row block
            #print 'New block row [',bno,i,']'
            blocks = concatenate((blocks, brow), axis = 0) # add in blocks array 

    # add an empty block at the end
    bno += 1
    brow = array([[bno, i]])
    #print 'Last row [',bno,i,']'
    blocks = concatenate((blocks, brow), axis = 0)
    
    return gaps, blocks

""" Generate moving average of the array errCC of calibration count errors per 
scanline, with window length clen and (constant) weight. Moving average is 
performed within a block from blocks array filling the gaps between scanlines 
in a block stored in scanlines array. """    
def genMAerr(errCC,weight,clen,scanlines,blocks):
    nob = blocks.shape[0] - 1 # number of scanline blocks
    nol = scanlines.shape[0] # number of scanlines
    maErr = zeros(nol) # moving average error per scanline
    
    # scanlines[:,1] contains the gaps
    
    for j in range(nob): # loop throught the blocks
        # moving average through scanlines blocks[j,1] to blocks[j+1,1]
        start = blocks[j,1] # first scanline (index) in the block
        end = blocks[j+1,1] # last scanline in the block
        
        # loop through scanlines in the block and calculate moving average
        if (end - start) > 2*clen: # the block has at list 51 scanlines
            for i in range(start, end):
                if i < (start+clen):
                    # fill the start of window with the error of first scanline  
                    maErr[i] += errCC[start]*(start+clen-i)*weight 
                    for k in range(start, i+clen+1): # build up weighted sum
                        maErr[i] += errCC[k]*scanlines[k,1]*weight
                elif i > (end-clen+1):
                    # fill the end of the window with the error of last scanline
                    maErr[i] += errCC[end-1]*(clen+i-end-1)*weight 
                    for k in range(i-clen, end): # build up weighted sum
                        maErr[i] += errCC[k]*scanlines[k,1]*weight
                else:
                    for k in range(i-clen,i+clen): # build up weighted sum
                        maErr[i] += errCC[k]*scanlines[k,1]*weight
        else: # block is shorter than averaging window
            for i in range(start, end):
                # fill missing scanlines with the start scanline error weighted
                maErr[i] = errCC[start]*(2*clen+1-i)*weight
                for k in range(start, i): # add weighted err of block scanlines
                    maErr[i] += errCC[k]*weight
            
    return maErr


''' Generate the matrix of errors respecting the correlation structure '''
def genPCS(Hr,Lsys,Tsys,Ccr,Cictr,maWgt,clen,scanlines,blocks,mcounts):
    err = empty(Hr.shape) # matrix of errors
    nor = err.shape[0] # number of matchups
    v1 = ones(nor) # array of ones with size no. of matchups
    
    # Lref, K, CE: random error from Gaussian with sigma from Hr data &mu=0
    err[:,0] = random.normal(scale=Hr[:,0]) # Lref random error
    err[:,6] = random.normal(scale=Hr[:,6]) # K random error
    err[:,3] = random.normal(scale=Hr[:,3]) # CE random error
    
    # Run moving average on scanlines to generate Cs error 
    errCs = random.normal(scale=Ccr) # Cs count error per scanline
    maeCs = genMAerr(errCs,maWgt,clen,scanlines,blocks) # moving average on scanlines
    # reconstruct err on matchups from errors in scanlines
    err[:,1] = repeat(maeCs, mcounts)
    
     # Run moving average on scanlines to generate Cict error 
    errCict = random.normal(scale=Cictr) # generate Cict error per scanline
    maeCict = genMAerr(errCict,maWgt,clen,scanlines,blocks) # moving average on scanlines
    err[:,2] = repeat(maeCict, mcounts) # recunstruct full err array
    
    # combined systematic and random Lict error
    errL = random.normal(scale=Lsys) # Lict systematic error
    err[:,4] = random.normal(scale=Hr[:,4]) + errL*v1 # Lict error
    
    # combined systematic and random To error
    errT = random.normal(scale=Tsys) # To systematic error
    err[:,5] = random.normal(scale=Hr[:,5]) + errT*v1 # To error
            
    return err

    
""" Return unix time in seconds (from 1970) from AVHRR time in sec (from 1975)
from SEH script readHD_SH.py  """
def conv2date(inTime):
    # Calculate difference from AVHRR start time and unix start time in seconds
    start_time_AVHRR = dt(1975, 1, 1)
    start_time_unix = dt(1970, 1, 1)
    time_diff = (start_time_AVHRR - start_time_unix).total_seconds()

    # Convert to time from 1975 to date
    outTime = [dt.fromtimestamp(time+time_diff) for time in inTime]

    return outTime
