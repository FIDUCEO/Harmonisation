"""
Created on Wed Nov 30 2016
@author: seh2

Functions for visualising FIDUCEO harmonisation data
(Ported to Python from Arta Dilo's R scripts)
AD: added functions for plotting fit results and uncertainty, 
removed make_sifn
"""

'''---Python Modules---'''
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from os.path import join as pjoin
import sys

'''---Harmonisation Modules---'''
sys.path.insert(0, '.\extras')
from readHD_SH import HData


class Plot:
    def __init__(self, saveDirectory, sensor1, sensor2, n_mu, timedata, sub):
        """Set variables used for the creation of graphs and image files
        from all the scripts for plotting"""

        # select a subset of data for plotting, e.g. max of 10,000 records
        if n_mu > 10000:
            n_mu_sample = sub
        else:
            n_mu_sample = n_mu

        self.n_mu = n_mu
        self.selMU = random.sample(range(n_mu), n_mu_sample) # Select matchups
        self.timedata = timedata
        self.seltimedata = [self.timedata[idx] for idx in self.selMU]

        # Set names
        self.refsensor_name = "MetopA"
        self.sensor1 = sensor1
        self.sensor2 = sensor2
        
        if self.sensor1 == -1:
            self.sensor1_name = self.refsensor_name
            sensor1_accronym = 'm02'
        else:
            self.sensor1_name = "NOAA" + str(sensor1)
            sensor1_accronym = 'n' + str(sensor1)
        self.sensor2_name = "NOAA" + str(sensor2)
        sensor2_accronym = 'n' + str(sensor2)
        
        self.sifn = sensor1_accronym + '_' + sensor2_accronym
        self.saveDirectory = saveDirectory

        # Set master title for all the graphs
        self.ttl = " ".join(("Matchups between AVHRR in", self.sensor1_name, "and", self.sensor2_name,
                             ": subset of", str(n_mu_sample), "matchups"))

    def set_xvar(self, idx):
        """Return x axis, label and limits for plotting and image filename suffix"""

        # If timedata available, use this as xvar
        if idx is not None:
            xvar = self.seltimedata                       # chart x axis data
            xlbl = 'Matchup Time'                         # chart x axis label
            xlims = [self.timedata[0], self.timedata[-1]] # chart x range
            fname_sfx = 'time.png'                        # image filename suffix

        # If timedata not available default to selMU as xvar
        else:
            xvar = self.selMU      # chart x axis data
            xlbl = 'Matchup Index' # chart x axis label
            xlims = [0, self.n_mu] # chart x range
            fname_sfx = 'mno.png'  # image filename suffix

        return xvar, xlbl, xlims, fname_sfx

    def plot_scatter(self, xvar, yvars, loc=((1,1),(0,0)), colspan=1, colors=['blue'],
                s=[25], alpha=1, ttl=None, xlbl=None, ylbl=None, 
                xlims=None, ylims=None,y_errorbar=None):
        """Plot scatter graph for input parameters"""
        
        ax = plt.subplot2grid(loc[0], loc[1], colspan=colspan)
	for yvar, color, size in zip(yvars, colors, s):
		if y_errorbar is None:
			ax.scatter(xvar, yvar, s=size, facecolors='none', edgecolors=color)
		else:
			ax.errorbars(xvar, yvar, yerr=y_errorbar, s=size, facecolors='none', edgecolors=color)
		ax.set_title(ttl)
		ax.set_xlabel(xlbl)
		ax.set_ylabel(ylbl)
		ax.set_xlim(xlims)
		if str(type(xvar[0])) == "<class 'datetime.datetime'>":
			ax.xaxis.set_major_locator(mdates.YearLocator())
			ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
			ax.xaxis.set_minor_locator(mdates.MonthLocator())

		if str(type(yvars[0][0])) == "<class 'datetime.datetime'>":
			ax.yaxis.set_major_locator(mdates.YearLocator())
			ax.yaxis.set_major_formatter(mdates.DateFormatter('%Y'))
			ax.yaxis.set_minor_locator(mdates.MonthLocator())

        return 0

    def picksensor(self, sensor_num):
        """Return sensor name and data frame colID suffix for given sensor"""

        # Get sensor name and suffix to data frame keys
        if sensor_num == 1:
            sensor_name = self.sensor1_name
            key_sfx = ""
        else:
            sensor_name = self.sensor2_name
            key_sfx = "2"

        return sensor_name, key_sfx

    def plotLbias(self, inL, calL, uL, plot_ttl, flab, idx=None):
        """ Create errorbar graph for radiance. Arguments to the function 
        - inL: radiance evaluated from input coefficients
        - calL: radiance evaluated from ODR fitted coefficients 
        - uL: radiance uncertainty from cal. coeffs/data and coeffs uncertainty 
        - plot_ttl: title of the plot 
        - flab: file naming label """
        
        # Set plotting x variable and label; Set output image path
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        imgFName = "_".join((flab, "fit", self.sifn))
        savePath = pjoin(self.saveDirectory, "".join((imgFName, ".png")))
        
        # range of +/- 2sigma uncertainties around fitted radiance
        up2sLE = calL + 2 * uL 
        low2sLE = calL - 2 * uL  

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Bias between calibrated radiances from fitted coeffs and input coefs 
        self.plot_scatter([inL[idx] for idx in self.selMU], 
                                [[inL[idx] for idx in self.selMU],
                                [up2sLE[idx] for idx in self.selMU],
                                [low2sLE[idx] for idx in self.selMU]],
                          colors=['limegreen','dimgray','dimgray'], s=[30,20,20],
                          loc=((1, 1), (0, 0)), colspan=1, alpha=0.5,
                          ttl=plot_ttl, xlbl="Radiance (mW/m2/sr/cm-1)",
                          ylbl="Radiance w.2sigma bounds (mW/m2/sr/cm-1)")
        
        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        if flab=='rnd':
            fig_ttl = 'Radiance calibration with random errors'
        elif flab=='errst':
            fig_ttl = 'Radiance calibration with full error structure'
        else:
            fig_ttl = ''
        fig.suptitle(fig_ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')
        
        return 0
    
    
    def plotLenU(self, inL, calL, ucLodr, ucLmc, ucdL, flab, idx=None):
        """ Creates PNG image for fit results. Arguments to the function 
        - inL: radiance evaluated from input coefficients
        - calL: radiance evaluated from ODR fitted coefficients
        - ucLodr: radiance uncertainty from cal.coeffs uncertainty evaluated by ODR
        - ucLmc: radiance uncertainty from cal.coeffs uncertainty evaluated by MC
        - ucdL: radiance uncertainty from calib. coeffs and data
        - flab: label 'rnd' | 'errst' for random and full error structure
        """

        # Set plotting x variable and label; Set output image path
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        imgFName = "_".join((flab, "fit", self.sifn))
        savePath = pjoin(self.saveDirectory, "".join((imgFName, ".png")))
        
        LElim = [min(inL), max(inL)]
        ULEodr = [min(ucLodr), max(ucLodr)]
        ULEmc = [min(ucLmc), max(ucLmc)]
        ULE = [min(ucdL), max(ucdL)]

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Bias between calibrated radiances from fitted coeffs and input coefs 
        self.plot_scatter(xvar, [[calL[idx] for idx in self.selMU],
                                [inL[idx] for idx in self.selMU]],
                          colors=['black','limegreen'], s=[30,20],
                          loc=((4, 1), (0, 0)), colspan=1, alpha=0.5,
                          ttl="Radiance from input coeffs (green) and ODR fitted (black) ",
                          xlbl="Matchups", ylbl="Radiance (mW/m2/sr/cm-1)", 
                          xlims=xlims)

        # Plot #2 - Radiance uncertainty from cal.coeffs uncertainty eval. by ODR
        self.plot_scatter([inL[idx] for idx in self.selMU], 
                            [[ucLodr[idx] for idx in self.selMU]],
                          colors=['maroon'], loc=((4, 1), (1, 0)), colspan=1,
                          ttl='Radiance uncertainty from ODR fitcoefficisnts, with uncertainty evaluated by ODR', 
                          xlbl="Radiance (mW/m2/sr/cm-1)",
                          ylbl="Uncertainty (mW/m2/sr/cm-1)",
                          xlims=LElim,ylims=ULEodr)

         # Plot #3 - Radiance uncertainty from coeffs uncertainty eval. by MC
        self.plot_scatter([inL[idx] for idx in self.selMU], 
                            [[ucLmc[idx] for idx in self.selMU]],
                          colors=['maroon'], loc=((4, 1), (2, 0)), colspan=1,
                          ttl='Radiance uncertainty from ODR fit coefficients, with uncertainty evaluated via MC', 
                          xlbl="Radiance (mW/m2/sr/cm-1)", 
                          ylbl="Uncertainty (mW/m2/sr/cm-1)",xlims=LElim,ylims=ULEmc)
                          
         # Plot #4 - Radiance uncertainty from cal.coeffs and data uncertainty
        self.plot_scatter([inL[idx] for idx in self.selMU], 
                            [[ucdL[idx] for idx in self.selMU]],
                          colors=['maroon'], loc=((4, 1), (3, 0)), colspan=1,
                          ttl='Radiance uncertainty from data and fit coefficients', 
                          xlbl="Radiance (mW/m2/sr/cm-1)", 
                          ylbl="Uncertainty (mW/m2/sr/cm-1)",xlims=LElim,ylims=ULE)
        

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        if flab=='rnd':
            fig_ttl = 'Radiance calibration with random errors'
        elif flab=='errst':
            fig_ttl = 'Radiance calibration with full error structure'
        else:
            fig_ttl = ''
        fig.suptitle(fig_ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotTrueErr(self, Hr, Hs, pfR, flab, idx=None):
        """ Creates PNG image with input and estimated errors from pair fitting: 
        random errors in harmonisation variables and 
        estimated true errors from odr fit. """

        # Set plotting x variable and label; Set output image path
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join((flab, "err", self.sifn, fname_sfx)))

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))
        Yr = np.sqrt(Hr[:,0]**2 + Hr[:,6]**2 + Hs[:,6]**2)
        sigLict = np.sqrt(Hr[:,4]**2 + Hs[:,4]**2)
        sigTo = np.sqrt(Hr[:,5]**2 + Hs[:,5]**2)

        # Plot #1 - Reference radiance and spectral adjustment, Y var in odr
        self.plot_scatter(xvar, [[pfR.eps[idx] for idx in self.selMU],
                                [Yr[idx] for idx in self.selMU]],
                          colors=['dimgray','magenta'], s=[25,15],
                          loc=((2, 3), (0, 0)), colspan=1, alpha=0.5, 
                          ttl="Radiance: input uncert. (magenta), ODR err (gray)",
                          xlbl=xlbl, ylbl="Radiance err (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #2 - Space count error
        self.plot_scatter(xvar, [[pfR.delta[0][idx] for idx in self.selMU],
                                [Hr[idx][1] for idx in self.selMU]],
                          colors=['dimgray','magenta'], s=[25,15],
                          loc=((2, 3), (0, 1)), colspan=1,
                          ttl="Space count: input uncert.(m), ODR err (g)",
                          xlbl=xlbl, ylbl="Count err", xlims=xlims)

        # Plot #3 - ICT count error
        self.plot_scatter(xvar, [[pfR.delta[1][idx] for idx in self.selMU],
                                [Hr[idx][2] for idx in self.selMU]],
                          colors=['dimgray','magenta'], s=[25,15],
                          loc=((2, 3), (0, 2)), colspan=1,
                          ttl="ICT count: input uncert.(m), ODR err (g)",
                          xlbl=xlbl, ylbl="Count err", xlims=xlims)

        # Plot #4 - Earth count error
        self.plot_scatter(xvar, [[pfR.delta[2][idx] for idx in self.selMU],
                                [Hr[idx][3] for idx in self.selMU]],
                          colors=['dimgray','magenta'], s=[25,15],
                          loc=((2, 3), (1, 0)), colspan=1,
                          ttl="Earth count: input uncert.(m), ODR err (g)",
                          xlbl=xlbl, ylbl="Count err", xlims=xlims)

        # Plot #5 - ICT radiance error
        self.plot_scatter(xvar, [[pfR.delta[3][idx] for idx in self.selMU],
                                [sigLict[idx][4] for idx in self.selMU]],
                          colors=['dimgray','magenta'], s=[25,15],
                          loc=((2, 3), (1, 1)), colspan=1,
                          ttl="ICT radiance: input uncert.(m), ODR err (g)",
                          xlbl=xlbl,ylbl="Radiance err (mW/m2/sr/cm-1)",xlims=xlims)

        # Plot #6 - orbit temperature error
        self.plot_scatter(xvar, [[pfR.delta[4][idx] for idx in self.selMU],
                                [sigTo[idx][5] for idx in self.selMU]],
                          colors=['dimgray','magenta'], s=[25,15],
                          loc=((2, 3), (1, 2)), colspan=1, 
                          ttl="Orbit temperature: input uncert.(m), ODR err (g)",
                          xlbl=xlbl, ylbl="Temperature err (K)", xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        if flab=='rnd':
            fig_ttl = 'Radiance calibration with random error structure'
        elif flab=='errst':
            fig_ttl = 'Radiance calibration with full error structure'
        else:
            fig_ttl = ''
        fig.suptitle(fig_ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')
        
        return 0


if __name__ == "__main__":
    def sample_code():
        """Example code illustrate how to use Plot class """

        inputPath = r"D:\Projects\FIDUCEO\Data\Simulated\m02_n15.nc"
        saveDirectory = r"D:\Projects\FIDUCEO\Data\Simulated\Graphs"

        hdata = HData(inputPath)
        if hdata.sensor1 == -1:
            sensor1 = 02
        else:
            sensor1 = hdata.sensor1

        # Create Plot class object
        Figure = Plot(saveDirectory, sensor1, hdata.sensor2, hdata.n_mu, hdata.timedata)

        # Plot data
        Figure.plotAll(hdata.H, hdata.Ur, hdata.Us, idx='time')
        Figure.plotAll(hdata.H, hdata.Ur, hdata.Us)

        return 0

    sample_code()