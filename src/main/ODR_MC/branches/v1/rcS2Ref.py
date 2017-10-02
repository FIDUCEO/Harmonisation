""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 19-12-2016
    Last update: 12-01-2017
Perform different model fitting for reference-sensor pairs: 
-    ordinary leas squares (OLS)
-    weighted least squres (WLS)
-    unweighted orthogonal distance regression (ODR)
-    weighted ODR, an error in variables (EIV) method.
Compare results of LS fits with the curve of input calibration coefficients. """

from os.path import join as pjoin
import readHD as rhd
import harFun as har
import plot_AD as pl
import upFun as upf

datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
filelist = ["m02_n15.nc"]#,"m02_n19.nc"]
slabel = 'avhrr' # series label

nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
print nos, 'sensors in the file list', slist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print 'Input calibration coefficients', inCoef
b0 = [-10., -4.e-3, 1.e-5, 0.0] #  a3 value to fix to input

#ncfile = "m02_n19.nc"
for ncfile in filelist:
    s2 = ncfile[4:7]
    b0[3] = inCoef[s2][3] # set a3 to the input value
    avhrrNx = upf.avhrr(nop, nos) # instance of class avhrr series 
    
    # read harmonisation data: variables and uncertainties
    rsp,Im,Hd,Hr,corIdx = rhd.rHDpair(datadir, ncfile)
    print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
    inL = avhrrNx.measEq(Hd[:,1:6], inCoef[s2])
    
    if rsp:
        # apply OLS fit
        #pols = har.odrP(Hd, Hr, b0, wgt=0, LStype=2)
        pols = har.pairLS(Hd, Hr, b0, wgt=0, LStype=2)
        print 'OLS output for', s2
        pols.pprint()
        # radiance and uncertainty from OLS fit
        oL = avhrrNx.measEq(Hd[:,1:6], pols.beta)
        oLU = avhrrNx.uncLE(Hd[:,1:6],pols.beta,Hr[:,1:6],pols.cov_beta)
        print 'Min and max ratio of Y var from fit to L from meas.eq.:'
        print min(pols.y/oL), max(pols.y/oL)
        print 'Min and max ratio of eps from fit to GUM uncertainty:'
        print min(abs(pols.eps)/oLU), max(abs(pols.eps)/oLU)
        
        # apply WLS fit
        #pwls = har.odrP(Hd, Hr, b0, LStype=2)
        pwls = har.pairLS(Hd, Hr, b0, wgt=1, LStype=2)
        print 'WLS output for', s2
        pwls.pprint()
        # radiance and uncertainty from WLS fit
        wL = avhrrNx.measEq(Hd[:,1:6], pwls.beta)
        wLU = avhrrNx.uncLE(Hd[:,1:6],pwls.beta,Hr[:,1:6],pwls.cov_beta)
        print 'Min and max ratio of Y var from fit to L from meas.eq.:'
        print min(pwls.y/wL), max(pwls.y/wL)
        print 'Min and max ratio of eps from fit to GUM uncertainty:'
        print min(abs(pwls.eps)/wLU), max(abs(pwls.eps)/wLU)
        
        # apply ODR to data without error' weighting 
        #podr = har.odrP(Hd, Hr, b0, wgt=0)
        podr = har.pairLS(Hd, Hr, b0, wgt=0)
        print 'Unweighted ODR output for', s2
        podr.pprint()
        # radiance and uncertainty from ODR fit
        rL = avhrrNx.measEq(Hd[:,1:6], podr.beta)
        rLU = avhrrNx.uncLE(Hd[:,1:6],podr.beta,Hr[:,1:6],podr.cov_beta)
        print 'Min and max ratio of Y var from fit to L from meas.eq.:'
        print min(podr.y/rL), max(podr.y/rL)
        print 'Min and max ratio of eps from fit to GUM uncertainty:'
        print min(abs(podr.eps)/rLU), max(abs(podr.eps)/rLU)
        
        # apply weighted ODR using X & Y errors
        #peiv = har.odrP(Hd, Hr, b0)
        peiv = har.pairLS(Hd, Hr, b0)
        print 'Weighted ODR output for', s2
        peiv.pprint()
        # radiance and uncertainty from weighted ODR fit
        eL = avhrrNx.measEq(Hd[:,1:6], peiv.beta)
        eLU = avhrrNx.uncLE(Hd[:,1:6],peiv.beta,Hr[:,1:6],peiv.cov_beta)
        print 'Min and max ratio of Y var from fit to L from meas.eq.:'
        print min(peiv.y/eL), max(peiv.y/eL)
        print 'Min and max ratio of eps from fit to GUM uncertainty:'
        print min(abs(peiv.eps)/eLU), max(abs(peiv.eps)/eLU)
        
        # graphs comparing different LS fits results
        Figure = pl.Plot(pltdir, Im[0,0], Im[0,1], Im[0,2], corIdx) # Plot instance
        Figure.plotLSfits(inL,oL,wL,rL,eL,oLU,wLU,rLU,eLU,Hd[:,1]-Hd[:,3])
        