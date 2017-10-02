from os.path import join as pjoin
from numpy import sqrt
import readHD as rhd
import harFun as har
#import py_compile
#py_compile.compile('plot_AD.py')
import plot_AD as pl
import upFun as upf

filelist = ["m02_n15.nc"]#,"m02_n19.nc"]
datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
#slabel = 'avhrr' # series label

slist = rhd.sensors(filelist) # list of sensors in filelist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print 'Sensors in the file list', slist
print 'Input calibration coefficients', inCoef
inC = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input
nop = len(filelist) # number of pairs
nos = len(slist) # number of sensors

ncfile = "m02_n15.nc"
for ncfile in filelist:
    s2 = ncfile[4:7]
    inC[3] = inCoef[s2][3] # set a3 to the input value
    avhrrNx = upf.avhrr(nop, nos) # instance of avhrr class
    
    # read harmonisation data: variables and uncertainties
    #rsp,Im,Hd,Hr,Hs,corIdx,corLen = rhd.rHDpair(datadir, ncfile)
    rsp,Im,Hd,Hr,corIdx = rhd.rHDpair(datadir, ncfile)
    print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
        
    # perform odr fit and display results
    podr = har.odrP(Hd, Hr, inC)
    podr.pprint()    
    
    # compute LE from fitted coefficients and uncertainty by GUM law of propag.
    fL = avhrrNx.measEq(Hd[:,1:6], podr.beta)
    fLu = avhrrNx.uncLE(Hd[:,1:6],podr.beta,Hr[:,1:6],podr.cov_beta)
    print 'Min and max ratio of Y var from fit with L from meas.eq.:'
    print min(podr.y/fL), max(podr.y/fL)
    print 'Min and max ratio of eps from fit with uncertainty by GUM law:'
    print min(abs(podr.eps)/fLu), max(abs(podr.eps)/fLu)
    
    """ Plot input and estimated errors and values of harmonisation vars """
    Figure = pl.Plot(pltdir, Im[0,0], Im[0,1], Im[0,2], corIdx) # Plot instance
    
    # plot graphs of estimated true values for X, Y vars and errors
    Figure.plotTrueErr(Hr, podr) 
    Figure.plotTrueVar(Hd, podr)
    
    # input radiance in simulations of data generation
    inL = avhrrNx.measEq(Hd[:,1:6], inCoef[s2]) # calculated from input cal.coeff

    # Y uncertainty: combined random uncert. of ref.radiance and adjustment vals
    Yu = sqrt(Hr[:,0]**2 + Hr[:,6]**2)
    
    # Y values: reference radiance + adjustment values -> Hd[:,0] + Hd[:,6]
    # space clamped Earth counts -> Hd[:,1]-Hd[:,3]
    Figure.plotFit(inL,Hd[:,0]+Hd[:,6],fL,Yu,podr.eps,fLu,Hd[:,1]-Hd[:,3])
