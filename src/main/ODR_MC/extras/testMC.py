import numpy as np
from datetime import datetime as dt
import readHD as rhd
import harFun as har

def conv2date(inTime):
    # Calculate difference from AVHRR start time and unix start time in seconds
    start_time_AVHRR = dt(1975, 1, 1)
    start_time_unix = dt(1970, 1, 1)
    time_diff = (start_time_AVHRR - start_time_unix).total_seconds()

    # Convert to time from 1975 to date
    outTime = [dt.fromtimestamp(time+time_diff) for time in inTime]

    return outTime

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

def genPCS(Hr,Lsys,Tsys,CCrnd,maWgt,clen,scanlines,blocks,mcounts):
    err = np.empty(Hr.shape) # matrix of errors
    nor = err.shape[0] # number of matchups
    v1 = np.ones(nor) # array of ones with size no. of matchups
    
    # Lref, K, CE: random error from Gaussian with sigma from Hr data &mu=0
    err[:,0] = np.random.normal(scale=Hr[:,0]) # Lref random error
    err[:,6] = np.random.normal(scale=Hr[:,6]) # K random error
    err[:,3] = np.random.normal(scale=Hr[:,3]) # CE random error
    
    # Run moving average on scanlines; fill err cols of Cs and Cict, all matchups
    errCs = np.random.normal(scale=CCrnd) # Cs count error per scanline
    maeCs = genMAerr(errCs,maWgt,clen,scanlines,blocks) # apply moving average on scanlines
    # reconstruct err on matchups from indices of unique scanlines; to CHECK !!
    err[:,1] = np.repeat(maeCs, mcounts)
    
    errCict = np.random.normal(scale=CCrnd) # generate Cict error per scanline
    maeCict = genMAerr(errCict,maWgt,clen,scanlines,blocks) # moving average of Cict errs
    err[:,2] = np.repeat(maeCict, mcounts) # recunstruct full err array
    
    # draw a value of systematic error for Lict and To
    errL = np.random.normal(scale=Lsys) # Lict systematic error
    errT = np.random.normal(scale=Tsys) # To systematic error
    # Lict, To: Gaussian random with sigma from Hr & mu=0, + systematic err
    err[:,4] = np.random.normal(scale=Hr[:,4]) + errL*v1 # Lict error
    err[:,5] = np.random.normal(scale=Hr[:,5]) + errT*v1 # To error
            
    return err

# Start execution !
filelist = ["m02_n15.nc"]#,"m02_n19.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
beta = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
slist = rhd.sensors(filelist) # list of sensors in filelist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
ncfile = filelist[0] # netCDF file to work with 
s2 = ncfile[4:7]
beta[3] = inCoef[s2][3] # set a3 to the input value

rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
Hs = har.resetHs(Hs, rsp) 

fixb = [1,1,1,0] # fix a3 coefficient
podr = har.odrP(Hd, Hr, beta, fixb, Hs)
podr.pprint()
b0odr = podr.beta # odr fit coefficients - beta0
sd0odr = podr.sd_beta # standard error of fit coefficients - sigb0
cov0odr = podr.cov_beta # odr evaluated covariance matrix - covb0

Y = podr.y # best est.of adjusted reference radiance: Lref + K
X = podr.xplus # best est. of explanatory variables: Cs,Cict,CE,Lict,To
sLict = np.mean(Hs[:,4]) # systematic error Lict
sTo = np.mean(Hs[:,5]) # systematic error To
cLen = int(corLen[0]) # scanlines moving average half-window 
wma = 1./(1+cLen*2) # moving average weight
szCC = 10 # sample size of calibration counts per scanline
sldt = 1. # time between consecutive scanlines; 0.5sec rounded up 

slt,midx,mcnt = np.unique(corIdx,return_index=True,return_counts=True)
uCC = Hr[midx,3] / np.sqrt(szCC) # calib.counts uncertainty from Earth count uncert.
slarr,slblocks = groupSln(slt,sldt,cLen)

st = dt.now() # start of MS run
# compile data for the ODR run; generate errors 
errStr = genPCS(Hr,sLict,sTo,uCC,wma,cLen,slarr,slblocks,mcnt) 
# add errStr to X & Y best estimates
Xdt = X.T + errStr[:,1:6] # X variables
Ydt = Y + errStr[:,0] + errStr[:,6] # Y variable
# run ODR on new X & Y vals and Hr weights 
podr = har.odr4MC(Xdt, Ydt, Hr, b0odr, fixb)
podr.pprint()

et = dt.now() # end of MC run
exect = (et-st).total_seconds()
print 'Time taken for one MC trials', (exect/60.), 'minutes'
