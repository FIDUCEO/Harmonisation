""" FIDUCEO FCDR harmonisation 
    Author: Arta Dilo / NPL MM
    Date created: 11-01-2017
    Last update: 13-01-2017
Perform pairwise harmonisation of sensors in a series, i.e. apply regression in
a sequential manner. Use pairLS function to perform a weighted ODR fit of 
2nd sensor to the 1st of a pair. When the 1st sensor is not a reference sensor, 
get calibration data from a previous fit in the sequence. """

from os.path import join as pjoin
import readHD as rhd
import harFun as har
import unpFun as upf

datadir = "D:\Projects\FIDUCEO\Data\Simulated" # data folder
pltdir = pjoin(datadir, 'Graphs') # folder for png images of graphs
filelist = ["m02_n19.nc","n19_n15.nc", "m02_n15.nc"]
slabel = 'avhrr' # series label

nop = len(filelist) # number of sensor pairs
print nop, 'pairs in filelist', filelist
slist = rhd.sensors(filelist) # list of sensors in filelist
nos = len(slist) # number of sensors
print nos, 'sensors in the file list', slist
inCoef = rhd.sInCoeff(datadir, filelist) # dictionary with input cal. coeffs 
print 'Input calibration coefficients:'
print inCoef

b0 = [-10., -4.e-3, 1.e-5, 0.0] # initialise array of cal.coeff; a3 to set later
avhrrNx = upf.avhrr(nop, nos) # instance of avhrr class
calC = {} # dictionary for cal. info from the fit

# for testing only; work within the loop
ncfile = filelist[1]

for ncfile in filelist:
    # read harmonisation data: variables and uncertainties
    rsp,Im,Hd,Hr,corIdx = rhd.rHDpair(datadir, ncfile)
    print 'NetCDF data from', ncfile, 'passed to harmonisation variables.'
        
    s2l  = 'n' + str(Im[0,1]) # label of 2nd sensor in the pair
    b0[3] = inCoef[s2l][3] # set a3 to the input value
    
    if s2l not in calC: # 2nd sensor of the pair not yet calibrated 
        if rsp: # if a reference-sensor pair
            s1l = 'm02' # label of 1st sensor in the pair

            # apply weighted ODR to data 
            podr = har.pairLS(Hd, Hr, b0)
            print 'ODR output for', s2l
            podr.pprint()
            
            # add fit results to calC dictionary of calibration info
            calC[s2l] = podr
        else:
            s1l = 'n' + str(Im[0,0]) # label of 1st sensor in the pair
            if s1l in calC: # 1st sensor calibration data in dictionary
                
                # compute Earth radiance for 1st sensor data and cal.coeff.
                LE = avhrrNx.measEq(Hd[:,0:5], calC[s1l].beta)
                # compute Earth radiance uncertainty
                uLE = avhrrNx.uncLE(Hd[:,0:5],calC[s1l].beta,Hr[:,0:5],calC[s1l].cov_beta)
                
                # fill 1st column of Hd and Hr matrix with Earth radiance data
                Hd[:,0] = LE # Earth radiance data
                Hr[:,0] = uLE # radiance uncertainty
                didx = [0, 5, 6, 7, 8, 9, 10] # data columns to use in the LS fit
                
                # perform weighted ODR fit for 2nd sensor in the pair
                podr = har.pairLS(Hd[:,didx], Hr[:,didx], b0)
                print 'ODR output for', s2l
                podr.pprint()
                
                # add fit results to calC dictionary of calibration info
                calC[s2l] = podr
        print s2l, "calibrated from netCDF file", ncfile
    kot = raw_input("Press enter to continue ...")
