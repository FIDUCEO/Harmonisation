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
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from os.path import join as pjoin
import sys

'''---Harmonisation Modules---'''
sys.path.insert(0, '.\extras')
from readHD_SH import HData


class Plot:
    def __init__(self, saveDirectory, sensor1, sensor2, n_mu, timedata):
        """Set variables used for the creation of graphs and image files
        from all the scripts for plotting"""

        # select a subset of data for plotting, e.g. max of 10,000 records
        if n_mu > 10000:
            n_mu_sample = 10000
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

    def plot_scatter(self, xvar, yvars, loc=((1,1),(0,0)), colspan=1, colors=['blue'], \
                s=[25], alpha=1, ttl=None, xlbl=None, ylbl=None, xlims=None, y_errorbar=None):
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

    def plotHDrefL(self, H, Ur, Us, idx=None):
        """Returns chart of reference radiance and uncertainties plots"""

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(("RefL ", self.sensor1_name, self.sifn, fname_sfx)))

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Reference Radiance
        self.plot_scatter(xvar, [[H['Lref'][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 3), (0, 0)), colspan=2,
                          ttl=self.sensor1_name+" Reference Radiance", 
                          xlbl=xlbl, ylbl="Radiance", xlims=xlims)

        # Plot #2 - Correlation index (as time)
        self.plot_scatter(self.selMU, [self.seltimedata],
                          colors=['grey'], loc=((2, 3), (0, 2)), colspan=1,
                          ttl="Correlation Index (as time)", 
                          xlbl="Matchup Index", ylbl="Time", xlims=[0, self.n_mu])

        # Plot #3 - Random error
        self.plot_scatter(xvar, [[Ur['Lref'][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 3), (1, 0)), colspan=1,
                          ttl=self.sensor1_name+" Radiance error - Random", 
                          xlbl=xlbl, ylbl="Error (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #4 - Random error histogram
        ax = plt.subplot2grid((2, 3), (1, 1))
        ax.hist([Ur['Lref'][idx] for idx in self.selMU],
                 bins=10, color="grey")
        ax.set_title("Hist. of "+self.sensor1_name+" Radiance error - Random")
        ax.set_xlabel("Error (mW/m2/sr/cm-1)")
        ax.set_ylabel("Frequency")

        # Plot #5 - Systematic error
        self.plot_scatter(xvar, [[Us['Lref'][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 3), (1, 2)), colspan=1,
                          ttl=self.sensor1_name+" Radiance error - Systematic", 
                          xlbl=xlbl, ylbl="Error (mW/m2/sr/cm-1)", xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotCnts(self, sensor_num, H, Ur, Us, idx=None):
        """Returns charts of for all counts: space, ICT and Earth"""

        sensor_name, key_sfx = self.picksensor(sensor_num)
        if sensor_name == self.refsensor_name:
            return 0 # If sensor is reference sensor no count data so exit

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(('Cnts', sensor_name, self.sifn, fname_sfx)))

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Plot Counts
        self.plot_scatter(xvar, [[H['Cspace'+key_sfx][idx] for idx in self.selMU],
                                 [H['CICT' + key_sfx][idx] for idx in self.selMU],
                                 [H['CEarth' + key_sfx][idx] for idx in self.selMU]],
                          colors=['grey', 'cyan', 'limegreen'], loc=((3, 3), (0, 0)), colspan=2,
                          ttl=sensor_name+" Counts: Earth (green), ICT (cyan), space (grey)", 
                          xlbl=xlbl, ylbl="Counts", xlims=xlims)

        # Plot #2 - Earth Count error
        self.plot_scatter(xvar, [[Ur['CEarth'+key_sfx][idx] for idx in self.selMU],
                                 [Us['CEarth'+key_sfx][idx] for idx in self.selMU]],
                          colors=['magenta', 'blue'], loc=((3, 3), (0, 2)), colspan=1,
                          ttl="$C_{Earth}$ Error components", 
                          xlbl=xlbl, ylbl="Count error", xlims=xlims)

        # Plot #3 - Space Counts
        self.plot_scatter(xvar, [[H['Cspace'+key_sfx][idx] for idx in self.selMU]],
                          colors=['grey'], loc=((3, 3), (1, 0)), colspan=1,
                          ttl=sensor_name+" Space Counts", xlbl=xlbl, ylbl="Counts", xlims=xlims)

        # Plot #4 - Histogram Space Counts
        ax = plt.subplot2grid((3, 3), (1, 1))
        ax.hist([H['Cspace'+key_sfx][idx] for idx in self.selMU],
                 bins=10, color="grey")
        ax.set_title("Hist. "+sensor_name+" Space Counts")
        ax.set_xlabel("Space Counts")
        ax.set_ylabel("Frequency")

        # Plot #5 - Space Counts error
        self.plot_scatter(xvar, [[Ur['Cspace'+key_sfx][idx] for idx in self.selMU],
                                 [Us['Cspace'+key_sfx][idx] for idx in self.selMU]],
                          colors=['magenta', 'blue'], loc=((3, 3), (1, 2)), colspan=1,
                          ttl="$C_{Space}$ Error Components", 
                          xlbl=xlbl, ylbl="Count error", xlims=xlims)

        # Plot #6 - ICT Counts
        self.plot_scatter(xvar, [[H['CICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['cyan'], loc=((3, 3), (2, 0)), colspan=1,
                          ttl=sensor_name+" ICT Counts", xlbl=xlbl, ylbl="Counts", xlims=xlims)

        # Plot #7 - ICT Counts Histogram
        ax = plt.subplot2grid((3, 3), (2, 1))
        ax.hist([H['CICT'+key_sfx][idx] for idx in self.selMU],
                 bins=10, color="grey")
        ax.set_title("Hist. "+sensor_name+" ICT Counts")
        ax.set_xlabel("Space Counts")
        ax.set_ylabel("Frequency")

        # Plot #8 - ICT Counts Error
        self.plot_scatter(xvar, [[Ur['CICT'+key_sfx][idx] for idx in self.selMU],
                                 [Us['CICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['magenta', 'blue'], loc=((3, 3), (2, 2)), colspan=1,
                          ttl="ICT Count Error Components", 
                          xlbl=xlbl, ylbl="Count error", xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotLICT(self, sensor_num, H, Ur, Us, idx=None):
        """Returns charts of ICT Radiance"""

        sensor_name, key_sfx = self.picksensor(sensor_num)
        if sensor_name == self.refsensor_name:
            return 0 # If sensor is reference sensor no count data so exit

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(("Lict", sensor_name, self.sifn, fname_sfx)))

        # Calculate variables to plot
        Cict = H["CICT"+key_sfx] - H["CEarth"+key_sfx]
        gain = H["LICT"+key_sfx] / Cict

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - ICT Radiance
        self.plot_scatter(xvar, [[H['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 3), (0, 0)), colspan=2,
                          ttl=sensor_name+" ICT radiance", 
                          xlbl=xlbl, ylbl="Radiance", xlims=xlims)

        # Plot #2 - Plot Instrument Gain
        self.plot_scatter(xvar, [[gain[idx] for idx in self.selMU]],
                          colors=['yellow'], loc=((2, 3), (0, 2)), colspan=1,
                          ttl=sensor_name+" Instrument Gain", 
                          xlbl=xlbl, ylbl="Gain (rad / count)", xlims=xlims)

        # Plot #3 - LICT Random Error
        self.plot_scatter(xvar, [[Ur['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 3), (1, 0)), colspan=1,
                          ttl=sensor_name+" $L_{ICT}$ Error - Random", 
                          xlbl=xlbl, ylbl="Error (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #4 - LICT Systematic Error
        self.plot_scatter(xvar, [[Us['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['blue'], loc=((2, 3), (1, 1)), colspan=1,
                          ttl=sensor_name+" $L_{ICT}$ Error - Systematic", 
                          xlbl=xlbl, ylbl="Error (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #5 - Space clamped ICT Counts
        self.plot_scatter(xvar, [[Cict[idx] for idx in self.selMU]],
                          colors=['cyan'], loc=((2, 3), (1, 2)), colspan=1,
                          ttl=sensor_name+" Space Clamped ICT Counts", 
                          xlbl=xlbl, ylbl="Counts", xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotTorb(self, sensor_num, H, Ur, Us, idx=None):
        """Return chart of orbital temperature"""

        sensor_name, key_sfx = self.picksensor(sensor_num)
        if sensor_name == self.refsensor_name:
            return 0 # If sensor is reference sensor no count data so exit

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(("Torb", sensor_name, self.sifn, fname_sfx)))

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Orbital Temperature
        self.plot_scatter(xvar, [[H['Torb'+key_sfx][idx] for idx in self.selMU]],
                          colors=['red'], loc=((2, 3), (0, 0)), colspan=2,
                          ttl=sensor_name+" Orbit Temperature", 
                          xlbl=xlbl, ylbl="Temperature (K)", xlims=xlims)

        # Plot #2 - ICT Radiance
        self.plot_scatter(xvar, [[H['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['green'], loc=((2, 3), (0, 2)), colspan=1,
                          ttl=sensor_name+" ICT Radiance", 
                          xlbl=xlbl, ylbl="Radiance", xlims=xlims)

        # Plot #3 - Orbital Temperature Error
        self.plot_scatter(xvar, [[Us['Torb'+key_sfx][idx] for idx in self.selMU],
                                 [Ur['Torb'+key_sfx][idx] for idx in self.selMU]],
                          colors=['blue', 'magenta'], loc=((2, 3), (1, 0)), colspan=1,
                          ttl=sensor_name+" Orbit Temperature Error Components", 
                          xlbl=xlbl, ylbl="Temperature error (K)", xlims=xlims)

        # Plot #4 - ICT radiance vs  Orbit Temperature
        self.plot_scatter([H['Torb'+key_sfx][idx] for idx in self.selMU],
                          [[H['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['grey'], loc=((2, 3), (1, 1)), colspan=1,
                          ttl=sensor_name+" ICT radiance vs  Orbit Temperature", 
                          xlbl="Orbit Temperature", ylbl="Radiance")

        # Plot #5 - ICT radiance vs  Orbit Temperature
        self.plot_scatter(xvar, [[Us['LICT'+key_sfx][idx] for idx in self.selMU],
                                 [Ur['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['blue', 'magenta'], loc=((2, 3), (1, 2)), colspan=1,
                          ttl=sensor_name+" L$_{ICT}$ Error Components", 
                          xlbl=xlbl, ylbl="Error (mW/m2/sr/cm-1)", xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotHDsAdj(self, H, Ur, Us, idx=None):
        """Returns chart of spectral adjustment values K and uncertainties"""

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(("SpctAdj", self.sifn, fname_sfx)))

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Spectral Adjustment
        self.plot_scatter(xvar, [[H['sAdj'][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 3), (0, 0)), colspan=1,
                          ttl="Spectral Adjustment between "+self.sensor1_name+
                          " and "+self.sensor2_name, 
                          xlbl=xlbl, ylbl="K Values", xlims=xlims)

        # Plot #2 - Spectral Adjustment Systematic Error
        self.plot_scatter(xvar, [[Us['sAdj'][idx] for idx in self.selMU]],
                          colors=['blue'], loc=((2, 3), (0, 1)), colspan=1,
                          ttl="Adjustment Error - Systematic", xlbl=xlbl, 
                          ylbl="Error (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #3 - Match Random Error
        self.plot_scatter(xvar, [[Ur['sAdj'][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 3), (0, 2)), colspan=1,
                          ttl="Matchup Error - Random", xlbl=xlbl, 
                          ylbl="Error (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #4 - K input
        self.plot_scatter(xvar, [[H['Kin'][idx] for idx in self.selMU]],
                          colors=['yellow'], loc=((2, 3), (1, 0)), colspan=1,
                          ttl="K Input Values", xlbl=xlbl, ylbl="K Input", xlims=xlims)

        # If a reference sensor pair
        if self.sensor1 == -1:
            # Plot #5 - Reference Radiance
            self.plot_scatter(xvar, [[H['Lref'][idx] for idx in self.selMU]],
                              colors=['limegreen'], loc=((2, 3), (1, 1)), colspan=1,
                              ttl=self.sensor1_name+" Reference Radiance", 
                              xlbl=xlbl, ylbl="Radiance", xlims=xlims)

            # Plot #6 - Time
            self.plot_scatter(self.selMU, [[self.seltimedata[idx] for idx in self.selMU]],
                              colors=["grey"], loc=((2, 3), (1, 2)), colspan=1,
                              ttl="Correlation index (as time)",
                              xlbl=xlbl, ylbl="Time", xlims=xlims)

        else:
            # Plot #5 - Hist. of Matchup Error - Random
            ax = plt.subplot2grid((2, 3), (1, 1))
            ax.hist([Ur['sAdj'][idx] for idx in self.selMU],
                    bins=10, color="grey")
            ax.set_title("Hist. of Matchup Error - Random")
            ax.set_xlabel("Error (mW/m2/sr/cm-1)")
            ax.set_ylabel("Frequency")

            # Plot #6 - Hist. of Adjustment Error - Systematic
            ax = plt.subplot2grid((2, 3), (1, 2))
            ax.hist([Us['sAdj'][idx] for idx in self.selMU],
                    bins=10, color="grey")
            ax.set_title("Hist. of Spectral Adjustment Error - Systematic")
            ax.set_xlabel("Error (mW/m2/sr/cm-1)")
            ax.set_ylabel("Frequency")

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotSGvrel(self, SGData, idx=None):
        """Returns chart of """

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(('SGdt', self.sifn, fname_sfx)))

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        if self.sensor1 == -1:
            # Plot #1 - Spectral Adjustment
            self.plot_scatter(xvar, [[SGData['K'][idx] for idx in self.selMU]],
                              colors=["limegreen"], loc=((2, 3), (0, 0)), colspan=2,
                              ttl="Spectral Adjustment, K", xlbl=xlbl, 
                              ylbl="K: radiance", xlims=xlims)

            # Plot #2 - Xlr
            self.plot_scatter(xvar, [[SGData['Xlr'][idx] for idx in self.selMU]],
                              colors=["black"], loc=((2, 3), (0, 2)), colspan=1,
                              ttl="Reference Radiance", xlbl=xlbl, 
                              ylbl="XLr radiance", xlims=xlims)

            # Plot #3 - Xlin
            self.plot_scatter(xvar, [[SGData['Xlin'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((2, 3), (1, 0)), colspan=1,
                              ttl="Linear component: gain * Earth counts", 
                              xlbl=xlbl, ylbl="XLin radiance", xlims=xlims)

            # Plot #4 - Xquad
            self.plot_scatter(xvar, [[SGData['Xquad'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((2, 3), (1, 1)), colspan=1,
                              ttl="Quadratic component: Earth counts^2 +...", 
                              xlbl=xlbl, ylbl="Xquad", xlims=xlims)

            # Plot #5 - Xtd
            self.plot_scatter(xvar, [[SGData['Xtd'][idx] for idx in self.selMU]],
                              colors=["grey"],  loc=((2, 3), (1, 1)), colspan=1,
                              ttl="Orbit Temperature", xlbl=xlbl, 
                              ylbl="Xtd Temperature", xlims=xlims)

        else:
            # Plot #1 - Spectral Adjustment
            self.plot_scatter(xvar, [[SGData['K'][idx] for idx in self.selMU]],
                              colors=["black"], loc=((3, 3), (0, 0)), colspan=3,
                              ttl="Spectral Adjustment, K", xlbl=xlbl, 
                              ylbl="K: radiance", xlims=xlims)

            # Plot #2 - Linear Component 1
            self.plot_scatter(xvar, [[SGData['Xlin1'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((3, 3), (1, 0)), colspan=1,
                              ttl=self.sensor1_name+" Linear Component", 
                              xlbl=xlbl, ylbl="Xlin1", xlims=xlims)

            # Plot #3 - Quadratic Component 1
            self.plot_scatter(xvar, [[SGData['Xquad1'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((3, 3), (1, 1)), colspan=1,
                              ttl=self.sensor1_name + " Quadratic Component", 
                              xlbl=xlbl, ylbl="Xquad1", xlims=xlims)

            # Plot #4 - Orbital Temperature 1
            self.plot_scatter(xvar, [[SGData['Xtd1'][idx] for idx in self.selMU]],
                              colors=["limegreen"], loc=((3, 3), (1, 2)), colspan=1,
                              ttl=self.sensor1_name + " Orbital Temperature", 
                              xlbl=xlbl, ylbl="Xtd1 temperature", xlims=xlims)

            # Plot #5 - Linear Component 2
            self.plot_scatter(xvar, [[SGData['Xlin2'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((3, 3), (1, 0)), colspan=1,
                              ttl=self.sensor1_name + " Linear Component", 
                              xlbl=xlbl, ylbl="Xlin2", xlims=xlims)

            # Plot #6 - Quadratic Component 2
            self.plot_scatter(xvar, [[SGData['Xquad2'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((3, 3), (1, 1)), colspan=1,
                              ttl=self.sensor1_name + " Quadratic Component", 
                              xlbl=xlbl, ylbl="Xquad2", xlims=xlims)

            # Plot #7 - Orbital Temperature 2
            self.plot_scatter(xvar, [[SGData['Xtd2'][idx] for idx in self.selMU]],
                              colors=["limegreen"], loc=((3, 3), (1, 2)), colspan=1,
                              ttl=self.sensor1_name + " Orbital Temperature", 
                              xlbl=xlbl, ylbl="Xtd2 temperature",
                              xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle("Variables for the stochastic gradient from matchups of "+ 
        self.sensor1_name+" and "+self.sensor2_name, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotFit(self, simLE, yL, fitLE, yU, fitRes, gumLU, CE, idx=None):
        """Returns PNG display of fit results: graph of fitted radiances,
        graphs of Y input error, ODR estimate of Y error, 
        radiance uncertainty propagated by GUM law. """

        # Set plotting x variable and label; Set output image path
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        imgFName = "_".join((self.sensor2_name, "FitL", self.sifn))
        savePath = pjoin(self.saveDirectory, "".join((imgFName, ".png")))
        CElim = [min(CE),max(CE)]
        LElim = [min(simLE), max(simLE)]

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Reference (sensor 1), fitted and simulations input radiance (sensor 2)
        self.plot_scatter([CE[idx] for idx in self.selMU], 
                            [[yL[idx] for idx in self.selMU],
                            [simLE[idx] for idx in self.selMU],
                            [fitLE[idx] for idx in self.selMU]],
                          colors=['limegreen', 'red', 'black'], s=[25,20,15],
                          loc=((2, 2), (0, 0)), colspan=1, alpha=0.5,
                          ttl="Adjusted "+self.sensor1_name+" reference (green), "+
                          self.sensor2_name+" input (red) and fitted (black)",
                          xlbl="Earth counts (space clamped)", 
                          ylbl="Radiance (mW/m2/sr/cm-1)",xlims=CElim)

        # Plot #2 - Y input error to ODR: combined random components of Lref & K
        self.plot_scatter([simLE[idx] for idx in self.selMU], 
                            [[yU[idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 2), (0, 1)), colspan=1,
                          ttl="Y input uncertainty to LS fit: combined random uncert. Lref and K",
                          xlbl="Radiance (mW/m2/sr/cm-1)", 
                          ylbl="Y error (mW/m2/sr/cm-1)",xlims=LElim)

        # Plot #3 - GUM uncertainty of Earth radiance in radiance scale
        if self.sensor1 == -1:  # if a reference-sensor pair
            ittl = 'GUM uncertainty of '+self.sensor2_name+' radiance'
            ulab = 'Uncertainty (mW/m2/sr/cm-1)'
        else:
            ittl = self.sensor1_name+' radiance fit residuals'
            ulab = 'Residuals (mW/m2/sr/cm-1)'
        self.plot_scatter([CE[idx] for idx in self.selMU], 
                            [[gumLU[idx] for idx in self.selMU]],
                          colors=['black'], loc=((2, 2), (1, 0)), colspan=1,
                          ttl=ittl, xlbl="Earth counts", ylbl=ulab,xlims=CElim)

        # Plot #4 - Radiance Fit Residuals in Radiance Scale
        self.plot_scatter([simLE[idx] for idx in self.selMU], 
                        [[fitRes[idx] for idx in self.selMU]],
                          colors=['dimgray'], loc=((2, 2), (1, 1)), colspan=1,
                          ttl='Fit residuals for '+self.sensor2_name+' radiance', 
                          xlbl="Radiance (mW/m2/sr/cm-1)", 
                          ylbl="Residuals (mW/m2/sr/cm-1)",xlims=LElim)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0
        
    def plotLSfits(self,incL,olsL,wlsL,odrL,eivL,olsLU,wlsLU,odrLU,eivLU,CE,idx=None):
        """ Compares results of LS fits: OLS, WLS, unweighted & weighted ODR.
        Plots the difference of LS calibrated curve with the input cal. curve to
        the simulation; plots uncertainty of calibrated radiance from each LS fit.
        """

        # Set plotting x variable and label; Set output image path
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        imgFName = "_".join((self.sensor2_name, "cmpLSfits", self.sifn))
        savePath = pjoin(self.saveDirectory, "".join((imgFName, ".png")))
        
        # calculate arrays for plots 
        OLSbias = olsL-incL # offset of ols fitted radiance to sim. input rad
        WLSbias = wlsL-incL # bias of WLS fitted radiance to input radiance
        ODRbias = odrL-incL # bias of unweighted ODR fitted radiance to input
        EIVbias = eivL-incL # bias of weighted ODR fitted radiance to input
        CElim = [min(CE),max(CE)]

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - radiance bias: fitted to input 
        self.plot_scatter([CE[i] for i in self.selMU], 
                            [[OLSbias[i] for i in self.selMU],
                            [WLSbias[i] for i in self.selMU],
                            [ODRbias[i] for i in self.selMU],
                            [EIVbias[i] for i in self.selMU]],
                          colors=['gold', 'orange', 'sienna','maroon'],alpha=0.5,
                          loc=((3, 2), (0, 0)), colspan=2, s=[35,30,25,15],
                          ttl=self.sensor2_name+" bias of fitted radiance to input:"
                          +" OLS (gold), WLS (orange), ODR (sienna), weighted ODR (maroon)",
                          xlbl="Earth counts (space clamped)",
                          ylbl="Radiance bias", xlims=CElim)

        # Plot #2 - radiance uncertainty from OLS fit
        self.plot_scatter([CE[idx] for idx in self.selMU], 
                            [[olsLU[idx] for idx in self.selMU]],
                          colors=['gold'], loc=((3, 2), (1, 0)), colspan=1,
                          ttl=self.sensor2_name+" OLS fitted radiance uncertainty",
                          xlbl="Earth counts", ylbl="Radiance err",xlims=CElim)

        # Plot #3 - radiance uncertainty from WLS fit
        self.plot_scatter([CE[idx] for idx in self.selMU], 
                            [[wlsLU[idx] for idx in self.selMU]],
                          colors=['orange'], loc=((3, 2), (1, 1)), colspan=1,
                          ttl=self.sensor2_name+" WLS fitted radiance uncertainty",
                          xlbl="Earth counts", ylbl="Radiance err",xlims=CElim)

        # Plot #4 - radiance uncertainty from unweighted ODR fit
        self.plot_scatter([CE[idx] for idx in self.selMU], 
                            [[odrLU[idx] for idx in self.selMU]],
                          colors=['sienna'], loc=((3, 2), (2, 0)), colspan=1,
                          ttl=self.sensor2_name+" radiance uncertainty unweighted ODR fit",
                          xlbl="Earth counts", ylbl="Radiance err",xlims=CElim)

        # Plot #5 - radiance uncertainty from weighted ODR fit
        self.plot_scatter([CE[idx] for idx in self.selMU], 
                            [[eivLU[idx] for idx in self.selMU]],
                          colors=['maroon'], loc=((3, 2), (2, 1)), colspan=1,
                          ttl=self.sensor2_name+" radiance uncertainty weighted ODR fit",
                          xlbl="Earth counts", ylbl="Radiance err",xlims=CElim)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotTrueErr(self, Hr, pfR, idx=None):
        """ Creates PNG image with input and estimated errors from pair fitting: 
        random errors in harmonisation variables and 
        estimated true errors from odr fit. """

        # Set plotting x variable and label; Set output image path
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(("TrueErr", self.sifn, fname_sfx)))

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Reference radiance and spectral adjustment, Y var in odr
        self.plot_scatter(xvar, [[pfR.eps[idx] for idx in self.selMU],
                                [Hr[idx][0] for idx in self.selMU],
                                [Hr[idx][6] for idx in self.selMU]],
                          colors=['dimgray','magenta','violet'], s=[15,25,15],
                          loc=((2, 3), (0, 0)), colspan=1, alpha=0.5, 
                          ttl="Y error: ref.radiance (magenta), K (violet), ODR estimate (gray)",
                          xlbl=xlbl, ylbl="Radiance err (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #2 - Space count error
        self.plot_scatter(xvar, [[Hr[idx][1] for idx in self.selMU],
                                [pfR.delta[0][idx] for idx in self.selMU]],
                          colors=['magenta','dimgray'], s=[25,15],
                          loc=((2, 3), (0, 1)), colspan=1,
                          ttl=self.sensor2_name+" space count error: input (m), est.truth (g)",
                          xlbl=xlbl, ylbl="Count error", xlims=xlims)

        # Plot #3 - ICT count error
        self.plot_scatter(xvar, [[Hr[idx][2] for idx in self.selMU],
                                [pfR.delta[1][idx] for idx in self.selMU]],
                          colors=['magenta','dimgray'], s=[25,15],
                          loc=((2, 3), (0, 2)), colspan=1,
                          ttl="ICT count error: input (m), est.truth (g)",
                          xlbl=xlbl, ylbl="Count error", xlims=xlims)

        # Plot #4 - Earth count error
        self.plot_scatter(xvar, [[Hr[idx][3] for idx in self.selMU],
                                [pfR.delta[2][idx] for idx in self.selMU]],
                          colors=['magenta','dimgray'], s=[25,15],
                          loc=((2, 3), (1, 0)), colspan=1,
                          ttl="Earth count error: input (m), est.truth (g)",
                          xlbl=xlbl, ylbl="Count error", xlims=xlims)

        # Plot #5 - ICT radiance error
        self.plot_scatter(xvar, [[Hr[idx][4] for idx in self.selMU],
                                [pfR.delta[3][idx] for idx in self.selMU]],
                          colors=['magenta','dimgray'], s=[25,15],
                          loc=((2, 3), (1, 1)), colspan=1,
                          ttl="ICT radiance error: input (m), est.truth (g)",
                          xlbl=xlbl,ylbl="Radiance err (mW/m2/sr/cm-1)",xlims=xlims)

        # Plot #6 - orbit temperature error
        self.plot_scatter(xvar, [[Hr[idx][5] for idx in self.selMU],
                                [pfR.delta[4][idx] for idx in self.selMU]],
                          colors=['magenta','dimgray'], s=[25,15],
                          loc=((2, 3), (1, 2)), colspan=1, 
                          ttl="Orbit temperature error: input (m), est.truth (g)",
                          xlbl=xlbl, ylbl="Temperature err (K)", xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')
        
        return 0

    def plotTrueVar(self, Hd, pfR, idx=None):
        """ Creates PNG image with input and estimated variable values from 
        pair fitting: harmonisation variables & estimated true values from odr fit. """

        # Set plotting x variable and label; Set output image path
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(("TrueVal", self.sifn, fname_sfx)))

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Reference radiance and spectral adjustment, Y var in odr
        self.plot_scatter(xvar, [[pfR.y[idx] for idx in self.selMU],
                                [Hd[idx,0]+Hd[idx,6] for idx in self.selMU]],
                          colors=['black','limegreen'], s=[30,20],
                          loc=((2, 3), (0, 0)), colspan=1, alpha=0.5,
                          ttl=self.sensor1_name+" ref.radiance + K value, and estimated Y (black)",
                          xlbl=xlbl, ylbl="Radiance (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #2 - Space counts
        self.plot_scatter(xvar, [[pfR.xplus[0,idx] for idx in self.selMU],
                                [Hd[idx,1] for idx in self.selMU]],
                          colors=['black','grey'], s=[30,20],
                          loc=((2, 3), (0, 1)), colspan=1, alpha=0.5,
                          ttl=self.sensor2_name+" space counts: input and est.truth (b)",
                          xlbl=xlbl, ylbl="Space count", xlims=xlims)

        # Plot #3 - ICT counts
        self.plot_scatter(xvar, [[pfR.xplus[1,idx] for idx in self.selMU],
                                [Hd[idx,2] for idx in self.selMU]],
                          colors=['black','cyan'], s=[15,30],
                          loc=((2, 3), (0, 2)), colspan=1, alpha=0.5,
                          ttl=self.sensor2_name+" ICT counts: input and est.truth (b)",
                          xlbl=xlbl, ylbl="ICT count", xlims=xlims)

        # Plot #4 - Earth count
        self.plot_scatter(xvar, [[pfR.xplus[2,idx] for idx in self.selMU],
                                [Hd[idx,3] for idx in self.selMU]],
                          colors=['black','limegreen'], s=[30,20],
                          loc=((2, 3), (1, 0)), colspan=1, alpha=0.5,
                          ttl=self.sensor2_name+" Earth counts: input and est.truth (b)",
                          xlbl=xlbl, ylbl="Count", xlims=xlims)

        # Plot #5 - ICT radiance
        self.plot_scatter(xvar, [[Hd[idx,4] for idx in self.selMU],
                                [pfR.xplus[3,idx] for idx in self.selMU]],
                          colors=['limegreen','black'], s=[30,20],
                          loc=((2, 3), (1, 1)), colspan=1, alpha=0.5,
                          ttl=self.sensor2_name+" ICT radiance: input and est.truth (b)",
                          xlbl=xlbl, ylbl="Radiance (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #6 - orbit temperature
        self.plot_scatter(xvar, [[pfR.xplus[4,idx] for idx in self.selMU],
                                [Hd[idx,5] for idx in self.selMU]],
                          colors=['black','red'], s=[30,20],
                          loc=((2, 3), (1, 2)), colspan=1, alpha=0.5,
                          ttl=self.sensor2_name+" orbit temperature: input and est.truth (b)",
                          xlbl=xlbl, ylbl="Temperature (K)", xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')
        
        return 0

    def plotAdjFit(self, newK, H, Ur, Us, idx=None):
        """Returns chart of """

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = pjoin(self.saveDirectory, "_".join(("SpctAdj", self.sifn, fname_sfx)))

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Spectral Adjustment
        self.plot_scatter(xvar, [[H['sAdj'][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((3, 3), (0, 0)), colspan=1,
                          ttl="Spectral Adjustment between " + self.sensor1_name + " and " + self.sensor2_name,
                          xlbl=xlbl, ylbl="K Values", xlims=xlims)

        # Plot #2 - Spectral Adjustment Systematic Error
        self.plot_scatter(xvar, [[Us['sAdj'][idx] for idx in self.selMU]],
                          colors=['blue'], loc=((3, 3), (0, 1)), colspan=1,
                          ttl="Adjustment Error - Systematic", 
                          xlbl=xlbl, ylbl="Error (mW/m2/sr/cm-1)", xlims=xlims)

        # Plot #3 - Match Random Error
        self.plot_scatter(xvar, [[Ur['sAdj'][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((3, 3), (0, 2)), colspan=1,
                          ttl="Matchup Error - Random", 
                          xlbl=xlbl, ylbl="Matchup error", xlims=xlims)

        # Plot #4 - K input
        self.plot_scatter(xvar, [[H['Kin'][idx] for idx in self.selMU]],
                          colors=['yellow'], loc=((3, 3), (1, 0)), colspan=1,
                          ttl="K Input Values", xlbl=xlbl, ylbl="K Input", xlims=xlims)

        # Plot #5 - Hist. of Matchup Error - Random
        ax = plt.subplot2grid((3, 3), (1, 1))
        ax.hist([Ur['sAdj'][idx] for idx in self.selMU],
                bins=10, color="grey")
        ax.set_title("Hist. of Matchup Error - Random")
        ax.set_xlabel("Matchup error")
        ax.set_ylabel("Frequency")

        # Plot #6 - Hist. of Adjustment Error - Systematic
        ax = plt.subplot2grid((3, 3), (1, 2))
        ax.hist([Us['sAdj'][idx] for idx in self.selMU],
                bins=10, color="grey")
        ax.set_title("Hist. of Spectral Adjustment Error - Systematic")
        ax.set_xlabel("Error (mW/m2/sr/cm-1)")
        ax.set_ylabel("Frequency")

        # Plot #7 - Reference Radiance
        self.plot_scatter(xvar, [[H['Lref'][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((3, 3), (2, 0)), colspan=1,
                          ttl=self.sensor1_name + " Reference Radiance", xlbl=xlbl, ylbl="Radiance", xlims=xlims)

        # Plot #8 - New K
        self.plot_scatter(xvar, [[newK[idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((3, 3), (2, 1)), colspan=1,
                          ttl="New K values from model fitting", xlbl=xlbl, ylbl="K Values", xlims=xlims)

        # Plot #9 - Hist. of new values of K
        ax = plt.subplot2grid((3, 3), (2, 2))
        ax.hist([newK[idx] for idx in self.selMU],
                bins=10, color="grey")
        ax.set_title("Hist. of new K Values")
        ax.set_xlabel("K Values")
        ax.set_ylabel("Frequency")

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotHDsensV(self, H, Ur, Us, idx=None):
        """Returns charts of counts, ICT radiance and orbital temperatures"""

        sensor_nums = [1, 2]
        for sensor_num in sensor_nums:
            self.plotCnts(sensor_num, H, Ur, Us, idx=idx)
            self.plotLICT(sensor_num, H, Ur, Us, idx=idx)
            self.plotTorb(sensor_num, H, Ur, Us, idx=idx)

        return 0

    def plotAll(self, H, Ur, Us, idx=None):

        self.plotHDrefL(H, Ur, Us, idx=idx)
        self.plotHDsensV(H, Ur, Us, idx=idx)

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