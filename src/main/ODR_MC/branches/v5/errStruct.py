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
    # scanline, e.g. gap = 1 for consecutive scanlines in the array sltidx,
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
def genPCS(Hr,Lsys,Tsys,CCrnd,maWgt,clen,scanlines,blocks,mcounts):
    err = np.empty(Hr.shape) # matrix of errors
    nor = err.shape[0] # number of matchups
    v1 = np.ones(nor) # array of ones with size no. of matchups
    
    # Lref, K, CE: random error from Gaussian with sigma from Hr data &mu=0
    err[:,0] = np.random.normal(scale=Hr[:,0]) # Lref random error
    err[:,6] = np.random.normal(scale=Hr[:,6]) # K random error
    err[:,3] = np.random.normal(scale=Hr[:,3]) # CE random error
    
    # Run moving average on scanlines to generate Cs error 
    errCs = np.random.normal(scale=CCrnd) # Cs count error per scanline
    maeCs = genMAerr(errCs,maWgt,clen,scanlines,blocks) # moving average on scanlines
    # reconstruct err on matchups from errors in scanlines
    err[:,1] = np.repeat(maeCs, mcounts)
    
     # Run moving average on scanlines to generate Cict error 
    errCict = np.random.normal(scale=CCrnd) # generate Cict error per scanline
    maeCict = genMAerr(errCict,maWgt,clen,scanlines,blocks) # moving average on scanlines
    err[:,2] = np.repeat(maeCict, mcounts) # recunstruct full err array
    
    # draw a value of systematic error for Lict and To
    errL = np.random.normal(scale=Lsys) # Lict systematic error
    errT = np.random.normal(scale=Tsys) # To systematic error
    # Lict, To: Gaussian random with sigma from Hr & mu=0, + systematic err
    err[:,4] = np.random.normal(scale=Hr[:,4]) + errL*v1 # Lict error
    err[:,5] = np.random.normal(scale=Hr[:,5]) + errT*v1 # To error
            
    return err

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma