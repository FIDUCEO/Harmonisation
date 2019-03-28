import numpy as np

#
# Convert temperature to radiance using band coefficients
#
class convertBT(object):

    def return_band_coefs(self,instr):

        if 'm02' == instr:
            in_nuc = [2685.11,926.605,837.805]
            in_aval = [1.99126,0.448099,0.382986]
            in_bval = [0.996287,0.998631,0.998879]
        elif 'n19' == instr:
            in_nuc = [2678.51,928.179,832.010]
            in_aval = [1.90192,0.406106,0.391425]
            in_bval = [0.999063,0.998698,0.998849]
        elif 'n18' == instr:
            in_nuc = [2675.31,928.680,834.073]
            in_aval = [2.10086,0.348003,0.408531]
            in_bval = [1.00006,0.998861,0.998805]
        elif 'n17' == instr:
            in_nuc = [2678.19,927.486,840.187]
            in_aval = [1.93286,0.431676,0.370802]
            in_bval = [0.999144,0.998637,0.998917]
        elif 'n16' == instr:
            in_nuc = [2688.43,919.600,836.710]
            in_aval = [1.85548,0.693914,0.440396]
            in_bval = [0.995179,0.998170,0.998707]
        elif 'n15' == instr:
            in_nuc = [2687.15,926.258,840.231]
            in_aval = [1.40674,0.451169,0.365465]
            in_bval = [0.996299,0.998636,0.998931]
        elif'n14' == instr:
            in_nuc = [2674.39,930.521,835.544]
            in_aval = [2.39486,0.621590,0.419563]
            in_bval = [1.00004,0.998294,0.998769]
        elif 'n12' == instr:
            in_nuc = [2671.81,921.925,838.003]
            in_aval = [2.42853,0.560212,0.398581]
            in_bval = [1.00102,0.998415,0.998832]
        elif 'n11' == instr:
            in_nuc = [2682.85,928.047,842.869]
            in_aval = [1.68218,0.393074,0.406868]
            in_bval = [0.997630,0.998759,0.998803]
        elif 'n10' == instr:
            in_nuc = [2679.11,910.515,-1.]
            in_aval = [1.93677,0.457244,-1.]
            in_bval = [0.998787,0.998768,-1.]
        elif 'n09' == instr:
            in_nuc = [2685.84,929.574,846.292]
            in_aval = [1.74231,0.347473,0.493176]
            in_bval = [0.996351,0.998839,0.998698]
        elif 'n08' == instr:
            in_nuc = [2671.12,915.281,-1.]
            in_aval = [2.31640,0.497544,-1.]
            in_bval = [1.00146,0.998654,-1.]
        elif 'n07' == instr:
            in_nuc = [2683.18,927.508,841.403]
            in_aval = [1.73290,0.405709,0.392737]
            in_bval = [0.997424,0.998741,0.998845]
        elif 'n06' == instr:
            in_nuc = [2678.59,913.461,-1.]
            in_aval = [1.94528,0.505261,-1.]
            in_bval = [0.998976,0.998636,-1.]
        elif 'tn' == instr:
            in_nuc = [2672.50,913.041,-1.]
            in_aval = [2.11244,0.531096,-1.]
            in_bval = [1.00119,0.998563,-1.]
        elif 'aatsr' == instr:
            in_nuc = [2681.59,923.695,832.715]
            in_aval = [1.77887,0.454752,0.498667]
            in_bval = [0.998005,0.998814,0.998525]
        elif 'atsr1' == instr:
            in_nuc = [2686.35,918.255,841.035]
            in_aval = [2.13333,0.532039,0.432575]
            in_bval = [0.995620,0.998695,0.998729]
        elif 'atsr2' == instr:
            in_nuc = [2699.00,917.236,831.569]
            in_aval = [2.04002,0.597643,0.467731]
            in_bval = [0.990740,0.998578,0.998620]
        else:
            raise Exception,"Cannot find instrument IR for SRF coeficients"
        
        self.nuc = in_nuc
        self.aval = in_aval
        self.bval = in_bval

    #
    # Convert temperature to radiance
    #
    def temp2rad(self,temperature,filterval):
    
        if self.nuc[filterval-1] < 0:
            raise Exception,"This AVHRR does not have a 12 micron channel"

        C1 = 1.191062E-5
        C2 = 1.4387863
        
        if isinstance(temperature,np.ndarray):
            outR = np.zeros(temperature.shape)
            if 3 == len(temperature.shape):
                outR[:,:,:] = float('nan')
            elif 2 == len(temperature.shape):
                outR[:,:] = float('nan')
            elif 1 == len(temperature.shape):
                outR[:] = float('nan')
                    
            gd = (temperature > 0.) & np.isfinite(temperature)
            tstar = self.aval[filterval-1]+\
                self.bval[filterval-1]*temperature[gd]
        
            outR[gd] = C1*(self.nuc[filterval-1]**3)/\
                (np.exp(C2*self.nuc[filterval-1]/tstar)-1.)
        else:
            tstar = self.aval[filterval-1]+\
                self.bval[filterval-1]*temperature
        
            outR = C1*(self.nuc[filterval-1]**3)/\
                (np.exp(C2*self.nuc[filterval-1]/tstar)-1.)

        return outR

    #
    # Convert radiance to temperature using band coefficients
    #
    def rad2temp(self,rad,filterval):
        
        if self.nuc[filterval-1] < 0:
            raise Exception,"This AVHRR does not have a 12 micron channel"

        C1 = 1.191062E-5
        C2 = 1.4387863
        
        if isinstance(temperature,np.ndarray):
            outBT = np.zeros(rad.shape)
            if 3 == len(rad.shape):
                outBT[:,:,:] = float('nan')
            elif 2 == len(rad.shape):
                outBT[:,:] = float('nan')
            elif 2 == len(rad.shape):
                outBT[:] = float('nan')

            gd = (rad > 0.) & np.isfinite(rad)     
            tstar = C2*self.nuc[filterval-1]/\
                np.log((C1*(self.nuc[filterval-1]**3)/rad[gd])+1.)        

            outBT[gd] = (tstar-self.aval[filterval-1])/self.bval[filterval-1]

        else:
            tstar = C2*self.nuc[filterval-1]/\
                np.log((C1*(self.nuc[filterval-1]**3)/rad)+1.)        

            outBT = (tstar-self.aval[filterval-1])/self.bval[filterval-1]

        return outBT

    def __init__(self,instr):

        self.return_band_coefs(instr)
        
