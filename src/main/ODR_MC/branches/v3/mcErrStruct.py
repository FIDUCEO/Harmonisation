""" Generate data for Monte Carlo uncertainty evaluation of odr fit coefficients """

import numpy as np
from datetime import datetime as dt
    
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

""" Group contiguous scanlines (possibly with gaps) in blocks such that two
blocks are away from each other more than half the averaging window, i.e. 
25 scanlines. Moving average is applied within a block. """
def groupSln(sltidx, sltd, cLen):
    nol = len(sltidx) # number of scanlines
    sltime = conv2date(sltidx) # convert scanline time idx to time format
    
    # create array of scanline idx and distance to next scanline in the array
    # gap = the number of lines that current scanline is away from the previous 
    # in the usl array, e.g. gap = 1 for consecutive scanlines in the array usl,
    # gap = 1 at the start of a new block.
    gaps = np.ones((nol,2)) # array to store scanline time idx and gaps 
    gaps[:,0] = sltidx # scanline time idx in first column
    
    # create blocks of scanlines that are >25 lines apart
    # first column stores the block number, start at 0 increase by 1;
    # second column stores the index of the 1st scanline in the block
    blocks = np.zeros((1,2), dtype=long) # 1st element: group=0, scanline idx=0
    bno = blocks[0, 0] # block number of 1st block = 0
    
    for i in range(1, nol):
        gaps[i,1] = (sltime[i] - sltime[i-1]).total_seconds()/sltd
        #print 'gap[',i,'] =',gaps[i,1]
        if gaps[i,1] > cLen: # 25: half window size in scanlines number
            gaps[i,1] = 1 # start of new block 
            #print 'gap[',i,'] =',gaps[i,1]
            bno += 1 # increase block number 
            brow = np.array([[bno, i]]) # new row block
            #print 'New block row [',bno,i,']'
            blocks = np.concatenate((blocks, brow), axis = 0) # add in blocks array 

    # add an empty block at the end
    bno += 1
    brow = np.array([[bno, i]])
    #print 'Last row [',bno,i,']'
    blocks = np.concatenate((blocks, brow), axis = 0)
    
    return gaps, blocks

""" Generate moving average of the array errCC of calibration count errors per 
scanline, with window length clen and (constant) weight. Moving average is 
performed within a block from blocks array filling the gaps between scanlines 
in a block stored in scanlines array. """    
def genMAerr(errCC,weight,clen,scanlines,blocks):
    nob = blocks.shape[0] - 1 # number of scanline blocks
    nol = scanlines.shape[0] # number of scanlines
    maErr = np.zeros(nol) # moving average error per scanline
    
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
def genPCS(Hr,Lsys,Tsys,CCrnd,maWgt,clen,scanlines,blocks,invIdx):
    err = np.empty(Hr.shape) # matrix of errors
    nor = err.shape[0] # number of matchups
    
    # draw a value of systematic error for Lict and To
    errL = np.random.normal(scale=Lsys) # Lict systematic error
    errT = np.random.normal(scale=Tsys) # To systematic error
        
    # Run moving average on scanlines; fill err cols of Cs and Cict, all matchups
    errCs = np.random.normal(scale=CCrnd) # Cs count error per scanline
    maeCs = genMAerr(errCs,maWgt,clen,scanlines,blocks) # apply moving average on scanlines
    # reconstruct err on matchups from indices of unique scanlines; to CHECK !!
    err[:1] = maeCs[invIdx] 
    
    errCict = np.random.normal(scale=CCrnd) # generate Cict error per scanline
    maeCict = genMAerr(errCict,maWgt,clen,scanlines,blocks) # moving average of Cict errs
    err[:,2] = maeCict[invIdx] # recunstruct full err array
    
    # create error matrix for the whole dataset
    for i in range(nor):
        # Lref, K, CE: random error from Gaussian with sigma from Hr data &mu=0
        err[i, 0] = np.random.normal(scale=Hr[i,0]) # Lref random error
        err[i, 3] = np.random.normal(scale=Hr[i,3]) # CE random error
        err[i, 6] = np.random.normal(scale=Hr[i,6]) # K random error
        
        # Lict, To: Gaussian random with sigma from Hr & mu=0, + systematic err
        err[i, 4] = np.random.normal(scale=Hr[i,3]) + errL # Lict error
        err[i, 5] = np.random.normal(scale=Hr[i,3]) + errT # To error
        
    return err

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma