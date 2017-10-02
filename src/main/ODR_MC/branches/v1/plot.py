"""
Created on Wed Nov 30 2016
@author: seh2

Functions for visualising FIDUCEO harmonisation data
(Ported to Python from Arta Dillo's R scripts)
"""

'''---Python Modules---'''
import random
import os.path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

'''---Harmonisation Modules---'''
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
        self.sifn = self.make_sifn(sensor1, sensor2)
        self.sensor1 = sensor1

        self.refsensor_name = "AVHRR02"
        self.sensor1_name = "AVHRR" + str(sensor1)
        if self.sensor1 == -1:
            self.sensor1_name = self.refsensor_name
        self.sensor2 = sensor2
        self.sensor2_name = "AVHRR" + str(sensor2)
        self.saveDirectory = saveDirectory

        # Set master title for all the graphs
        self.ttl = " ".join(("Matchups between", self.sensor1_name, "and", self.sensor2_name,
                             ": subset of", str(n_mu_sample), "matchups"))

    def make_sifn(self, sensor1, sensor2):
        def get_mission_letter(sensor):
            if int(sensor) == 2 or int(sensor) == -1:
                missionLetter = "m"
            else:
                missionLetter = "n"
            return missionLetter

        letter1 = get_mission_letter(sensor1)
        letter2 = get_mission_letter(sensor2)
        return "".join((letter1, str(sensor1), "_", letter2, str(sensor2)))

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
                     ttl=None, xlbl=None, ylbl=None, xlims=None):
        """Plot scatter graph for input parameters"""

        ax = plt.subplot2grid(loc[0], loc[1], colspan=colspan)
        for yvar, color in zip(yvars, colors):
            ax.scatter(xvar, yvar, s=35, facecolors='none', edgecolors=color)
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
        savePath = os.path.join(self.saveDirectory, "_".join(("RefL ", self.sensor1_name, self.sifn, fname_sfx)))

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Reference Radiance
        self.plot_scatter(xvar, [[H['Lref'][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 3), (0, 0)), colspan=2,
                          ttl=self.sensor1_name+" Reference Radiance", xlbl=xlbl, ylbl="Radiance", xlims=xlims)

        # Plot #2 - Correlation index (as time)
        self.plot_scatter(self.selMU, [self.seltimedata],
                          colors=['grey'], loc=((2, 3), (0, 2)), colspan=1,
                          ttl="Correlation Index (as time)", xlbl="Matchup Index", ylbl="Time", xlims=[0, self.n_mu])

        # Plot #3 - Random uncertainty
        self.plot_scatter(xvar, [[Ur['Lref'][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 3), (1, 0)), colspan=1,
                          ttl=self.sensor1_name+" Radiance Uncertainty - Random", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Plot #4 - Random uncertainty histogram
        ax = plt.subplot2grid((2, 3), (1, 1))
        ax.hist([Ur['Lref'][idx] for idx in self.selMU],
                 bins=10, color="grey")
        ax.set_title("Hist. of "+self.sensor1_name+" Radiance Uncertainty - Random")
        ax.set_xlabel("Standard Uncertainty")
        ax.set_ylabel("Frequency")

        # Plot #5 - Systematic uncertainty
        self.plot_scatter(xvar, [[Us['Lref'][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 3), (1, 2)), colspan=1,
                          ttl=self.sensor1_name+" Radiance Uncertainty - Systematic", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

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
        savePath = os.path.join(self.saveDirectory, "_".join(('Cnts', sensor_name, self.sifn, fname_sfx)))

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Plot Counts
        self.plot_scatter(xvar, [[H['Cspace'+key_sfx][idx] for idx in self.selMU],
                                 [H['CICT' + key_sfx][idx] for idx in self.selMU],
                                 [H['CEarth' + key_sfx][idx] for idx in self.selMU]],
                          colors=['grey', 'cyan', 'limegreen'], loc=((3, 3), (0, 0)), colspan=2,
                          ttl=sensor_name+" Counts: Earth (green), ICT (cyan), space (grey)", xlbl=xlbl, ylbl="Counts", xlims=xlims)

        # Plot #2 - Earth Count Uncertainty
        self.plot_scatter(xvar, [[Ur['CEarth'+key_sfx][idx] for idx in self.selMU],
                                 [Us['CEarth'+key_sfx][idx] for idx in self.selMU]],
                          colors=['magenta', 'blue'], loc=((3, 3), (0, 2)), colspan=1,
                          ttl="$C_{Earth}$ Uncertainty Components", xlbl=xlbl, ylbl="Standard uncertainty", xlims=xlims)

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

        # Plot #5 - Space Counts Uncertainty
        self.plot_scatter(xvar, [[Ur['Cspace'+key_sfx][idx] for idx in self.selMU],
                                 [Us['Cspace'+key_sfx][idx] for idx in self.selMU]],
                          colors=['magenta', 'blue'], loc=((3, 3), (1, 2)), colspan=1,
                          ttl="$C_{Space}$ Uncertainty Components", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

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

        # Plot #8 - ICT Counts Uncertainty
        self.plot_scatter(xvar, [[Ur['CICT'+key_sfx][idx] for idx in self.selMU],
                                 [Us['CICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['magenta', 'blue'], loc=((3, 3), (2, 2)), colspan=1,
                          ttl="ICT Count Uncertainty Components", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

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
        savePath = os.path.join(self.saveDirectory, "_".join(("Lict", sensor_name, self.sifn, fname_sfx)))

        # Calculate variables to plot
        Cict = H["CICT"+key_sfx] - H["CEarth"+key_sfx]
        gain = H["LICT"+key_sfx] / Cict

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - ICT Radiance
        self.plot_scatter(xvar, [[H['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 3), (0, 0)), colspan=2,
                          ttl=sensor_name+" ICT radiance", xlbl=xlbl, ylbl="Radiance", xlims=xlims)

        # Plot #2 - Plot Instrument Gain
        self.plot_scatter(xvar, [[gain[idx] for idx in self.selMU]],
                          colors=['yellow'], loc=((2, 3), (0, 2)), colspan=1,
                          ttl=sensor_name+" Instrument Gain", xlbl=xlbl, ylbl="Gain (rad / count)", xlims=xlims)

        # Plot #3 - LICT Random Uncertainty
        self.plot_scatter(xvar, [[Ur['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 3), (1, 0)), colspan=1,
                          ttl=sensor_name+" $L_{ICT}$ Uncertainty - Random", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Plot #4 - LICT Systematic Uncertainty
        self.plot_scatter(xvar, [[Us['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['blue'], loc=((2, 3), (1, 1)), colspan=1,
                          ttl=sensor_name+" $L_{ICT}$ Uncertainty - Systematic", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Plot #5 - Space clamped ICT Counts
        self.plot_scatter(xvar, [[Cict[idx] for idx in self.selMU]],
                          colors=['cyan'], loc=((2, 3), (1, 2)), colspan=1,
                          ttl=sensor_name+" Space Clamped ICT Counts", xlbl=xlbl, ylbl="Counts", xlims=xlims)

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
        savePath = os.path.join(self.saveDirectory, "_".join(("Torb", sensor_name, self.sifn, fname_sfx)))

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Orbital Temperature
        self.plot_scatter(xvar, [[H['Torb'+key_sfx][idx] for idx in self.selMU]],
                          colors=['red'], loc=((2, 3), (0, 0)), colspan=2,
                          ttl=sensor_name+" Orbit Temperature", xlbl=xlbl, ylbl="Temperature (K)", xlims=xlims)

        # Plot #2 - ICT Radiance
        self.plot_scatter(xvar, [[H['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['green'], loc=((2, 3), (0, 2)), colspan=1,
                          ttl=sensor_name+" ICT Radiance", xlbl=xlbl, ylbl="Radiance", xlims=xlims)

        # Plot #3 - Orbital Temperature Uncertainty
        self.plot_scatter(xvar, [[Us['Torb'+key_sfx][idx] for idx in self.selMU],
                                 [Ur['Torb'+key_sfx][idx] for idx in self.selMU]],
                          colors=['blue', 'magenta'], loc=((2, 3), (1, 0)), colspan=1,
                          ttl=sensor_name+" Orbit Temperature Uncertainty Components", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Plot #4 - ICT radiance vs  Orbit Temperature
        self.plot_scatter([H['Torb'+key_sfx][idx] for idx in self.selMU],
                          [[H['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['grey'], loc=((2, 3), (1, 1)), colspan=1,
                          ttl=sensor_name+" ICT radiance vs  Orbit Temperature", xlbl="Orbit Temperature", ylbl="Radiance")

        # Plot #5 - ICT radiance vs  Orbit Temperature
        self.plot_scatter(xvar, [[Us['LICT'+key_sfx][idx] for idx in self.selMU],
                                 [Ur['LICT'+key_sfx][idx] for idx in self.selMU]],
                          colors=['blue', 'magenta'], loc=((2, 3), (1, 2)), colspan=1,
                          ttl=sensor_name+" L$_{ICT}$ Uncertainty Components", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotHDsAdj(self, H, Ur, Us, idx=None):
        """Returns chart of spectral adjustment values K and uncertainties"""

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = os.path.join(self.saveDirectory, "_".join(("SpctAdj", self.sifn, fname_sfx)))

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Spectral Adjustment
        self.plot_scatter(xvar, [[H['sAdj'][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 3), (0, 0)), colspan=1,
                          ttl="Spectral Adjustment between "+self.sensor1_name+" and "+self.sensor2_name, xlbl=xlbl, ylbl="K Values", xlims=xlims)

        # Plot #2 - Spectral Adjustment Systematic Uncertainty
        self.plot_scatter(xvar, [[Us['sAdj'][idx] for idx in self.selMU]],
                          colors=['blue'], loc=((2, 3), (0, 1)), colspan=1,
                          ttl="Adjustment Uncertainty - Systematic", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Plot #3 - Match Random Uncertainty
        self.plot_scatter(xvar, [[Ur['sAdj'][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((2, 3), (0, 2)), colspan=1,
                          ttl="Matchup Uncertainty - Random", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Plot #4 - K input
        self.plot_scatter(xvar, [[H['Kin'][idx] for idx in self.selMU]],
                          colors=['yellow'], loc=((2, 3), (1, 0)), colspan=1,
                          ttl="K Input Values", xlbl=xlbl, ylbl="K Input", xlims=xlims)

        # If a reference sensor pair
        if self.sensor1 == -1:
            # Plot #5 - Reference Radiance
            self.plot_scatter(xvar, [[H['Lref'][idx] for idx in self.selMU]],
                              colors=['limegreen'], loc=((2, 3), (1, 1)), colspan=1,
                              ttl=self.sensor1_name+" Reference Radiance", xlbl=xlbl, ylbl="Radiance", xlims=xlims)

            # Plot #6 - Time
            self.plot_scatter(self.selMU, [[self.seltimedata[idx] for idx in self.selMU]],
                              colors=["grey"], loc=((2, 3), (1, 2)), colspan=1,
                              ttl="Correlation index (as time)",xlbl=xlbl, ylbl="Time", xlims=xlims)

        else:
            # Plot #5 - Hist. of Matchup Uncertainty - Random
            ax = plt.subplot2grid((2, 3), (1, 1))
            ax.hist([Ur['sAdj'][idx] for idx in self.selMU],
                    bins=10, color="grey")
            ax.set_title("Hist. of Matchup Uncertainty - Random")
            ax.set_xlabel("Standard uncertainty")
            ax.set_ylabel("Frequency")

            # Plot #6 - Hist. of Adjustment Uncertainty - Systematic
            ax = plt.subplot2grid((2, 3), (1, 2))
            ax.hist([Us['sAdj'][idx] for idx in self.selMU],
                    bins=10, color="grey")
            ax.set_title("Hist. of Spectral Adjustment Uncertainty - Systematic")
            ax.set_xlabel("Standard uncertainty")
            ax.set_ylabel("Frequency")

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotSGvrel(self, SGData):
        """Returns chart of """

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = os.path.join(self.saveDirectory, "_".join(('SGdt', self.sifn, fname_sfx)))

        # Initialise plot
        fig = plt.figure(figsize=(20, 12))

        if self.sensor1 == -1:
            # Plot #1 - Spectral Adjustment
            self.plot_scatter(xvar, [[SGData['K'][idx] for idx in self.selMU]],
                              colors=["limegreen"], loc=((2, 3), (0, 0)), colspan=2,
                              ttl="Spectral Adjustment, K", xlbl=xlbl, ylbl="K: radiance", xlims=xlims)

            # Plot #2 - Xlr
            self.plot_scatter(xvar, [[SGData['Xlr'][idx] for idx in self.selMU]],
                              colors=["black"], loc=((2, 3), (0, 2)), colspan=1,
                              ttl="Reference Radiance", xlbl=xlbl, ylbl="XLr radiance", xlims=xlims)

            # Plot #3 - Xlin
            self.plot_scatter(xvar, [[SGData['Xlin'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((2, 3), (1, 0)), colspan=1,
                              ttl="Linear component: gain * Earth counts", xlbl=xlbl, ylbl="XLin radiance", xlims=xlims)

            # Plot #4 - Xquad
            self.plot_scatter(xvar, [[SGData['Xquad'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((2, 3), (1, 1)), colspan=1,
                              ttl="Quadratic component: Earth counts^2 +...", xlbl=xlbl, ylbl="Xquad", xlims=xlims)

            # Plot #5 - Xtd
            self.plot_scatter(xvar, [[SGData['Xtd'][idx] for idx in self.selMU]],
                              colors=["grey"],  loc=((2, 3), (1, 1)), colspan=1,
                              ttl="Orbit Temperature", xlbl=xlbl, ylbl="Xtd Temperature", xlims=xlims)

        else:
            # Plot #1 - Spectral Adjustment
            self.plot_scatter(xvar, [[SGData['K'][idx] for idx in self.selMU]],
                              colors=["black"], loc=((3, 3), (0, 0)), colspan=3,
                              ttl="Spectral Adjustment, K", xlbl=xlbl, ylbl="K: radiance", xlims=xlims)

            # Plot #2 - Linear Component 1
            self.plot_scatter(xvar, [[SGData['Xlin1'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((3, 3), (1, 0)), colspan=1,
                              ttl=self.sensor1_name+" Linear Component", xlbl=xlbl, ylbl="Xlin1", xlims=xlims)

            # Plot #3 - Quadratic Component 1
            self.plot_scatter(xvar, [[SGData['Xquad1'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((3, 3), (1, 1)), colspan=1,
                              ttl=self.sensor1_name + " Quadratic Component", xlbl=xlbl, ylbl="Xquad1", xlims=xlims)

            # Plot #4 - Orbital Temperature 1
            self.plot_scatter(xvar, [[SGData['Xtd1'][idx] for idx in self.selMU]],
                              colors=["limegreen"], loc=((3, 3), (1, 2)), colspan=1,
                              ttl=self.sensor1_name + " Orbital Temperature", xlbl=xlbl, ylbl="Xtd1 temperature", xlims=xlims)

            # Plot #5 - Linear Component 2
            self.plot_scatter(xvar, [[SGData['Xlin2'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((3, 3), (1, 0)), colspan=1,
                              ttl=self.sensor1_name + " Linear Component", xlbl=xlbl, ylbl="Xlin2", xlims=xlims)

            # Plot #6 - Quadratic Component 2
            self.plot_scatter(xvar, [[SGData['Xquad2'][idx] for idx in self.selMU]],
                              colors=["grey"], loc=((3, 3), (1, 1)), colspan=1,
                              ttl=self.sensor1_name + " Quadratic Component", xlbl=xlbl, ylbl="Xquad2", xlims=xlims)

            # Plot #7 - Orbital Temperature 2
            self.plot_scatter(xvar, [[SGDat['Xtd2'][idx] for idx in self.selMU]],
                              colors=["limegreen"], loc=((3, 3), (1, 2)), colspan=1,
                              ttl=self.sensor1_name + " Orbital Temperature", xlbl=xlbl, ylbl="Xtd2 temperature",
                              xlims=xlims)

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle("Variables for the stochastic gradient from matchups of "+self.sensor1_name+" and "+self.sensor2_name,
                     fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotFit(self, H, s2fL, s2cL, UrL, VL):
        """Returns PNG display results of fitting: plot graph of fitted radiances.
        Print calibration coefficients & covariance matrix"""

        # Set output image path
        savePath = os.path.join(self.saveDirectory, "_".join(("FitL", self.sensor2_name, self.sifn, ".png")))

        # Calculate variables to plot
        CE = H["Cspace2"] - H["CEarth2"]  # space clamped Earth counts

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Reference (sensor 1), fitted and simulations input radiance (sensor 2)
        self.plot_scatter([CE[idx] for idx in self.selMU], [[H['Lref'][idx] for idx in self.selMU],
                                                            [s2fL[idx] for idx in self.selMU],
                                                            [s2cL[idx] for idx in self.selMU]],
                          colors=['limegreen', 'black', 'red'], loc=((2, 2), (0, 0)), colspan=1,
                          ttl=self.sensor1_name+" reference (green) and "+self.sensor2_name+" fitted (black) and input (red)",
                          xlbl="Earth counts (space clamped)", ylbl="Radiance")

        # Plot #2 - Spectral Adjustment
        self.plot_scatter([CE[idx] for idx in self.selMU], [[H['sAdj'][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 2), (0, 1)), colspan=1,
                          ttl="Spectral Adjustment "+self.sensor1_name+" to "+self.sensor2_name,
                          xlbl="Earth counts (space clamped)", ylbl="K Values")

        # Plot #3 - plot total uncertainty for reference radiance in radiance scale
        if self.sensor1 == -1:  # if a reference-sensor pair
            ittl = 'Total uncertainty of '+self.sensor1_name+'reference radiance'
            UrL = sqrt(VL)  # total standard uncertainty
            ulab = 'Standard uncertainty'
        else:
            ittl = self.sensor1_name+' radiance fit residuals'
            ulab = 'Residuals'
        self.plot_scatter([H['Lref'][idx] for idx in self.selMU], [[UrL[idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 2), (1, 0)), colspan=1,
                          ttl=ittl, xlbl="Radiance", ylbl=ulab)

        # Plot #4 - Radiance Fit Residuals in Radiance Scale
        self.plot_scatter([H['Lref'][idx] for idx in self.selMU], [[fErrL[idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((2, 2), (1, 1)), colspan=1,
                          ttl='Fit residuals for '+self.sensor2_name+' radiance', xlbl="Radiance", ylbl="Residuals")

        # Resize plot, add title and save
        plt.tight_layout(pad=1.08, w_pad=0.5, h_pad=5.0, rect=(0, 0, 1, 1))
        fig.suptitle(self.ttl, fontsize=20, y=1.08)
        plt.savefig(savePath, bbox_inches='tight')

        return 0

    def plotAdjFit(self, newK, H, Ur, Us, idx=None):
        """Returns chart of """

        # Set plotting x variable and label
        xvar, xlbl, xlims, fname_sfx = self.set_xvar(idx)
        savePath = os.path.join(self.saveDirectory, "_".join(("SpctAdj", self.sifn, fname_sfx)))

        # Initialise Plot
        fig = plt.figure(figsize=(20, 12))

        # Plot #1 - Spectral Adjustment
        self.plot_scatter(xvar, [[H['sAdj'][idx] for idx in self.selMU]],
                          colors=['limegreen'], loc=((3, 3), (0, 0)), colspan=1,
                          ttl="Spectral Adjustment between " + self.sensor1_name + " and " + self.sensor2_name,
                          xlbl=xlbl, ylbl="K Values", xlims=xlims)

        # Plot #2 - Spectral Adjustment Systematic Uncertainty
        self.plot_scatter(xvar, [[Us['sAdj'][idx] for idx in self.selMU]],
                          colors=['blue'], loc=((3, 3), (0, 1)), colspan=1,
                          ttl="Adjustment Uncertainty - Systematic", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Plot #3 - Match Random Uncertainty
        self.plot_scatter(xvar, [[Ur['sAdj'][idx] for idx in self.selMU]],
                          colors=['magenta'], loc=((3, 3), (0, 2)), colspan=1,
                          ttl="Matchup Uncertainty - Random", xlbl=xlbl, ylbl="Standard Uncertainty", xlims=xlims)

        # Plot #4 - K input
        self.plot_scatter(xvar, [[H['Kin'][idx] for idx in self.selMU]],
                          colors=['yellow'], loc=((3, 3), (1, 0)), colspan=1,
                          ttl="K Input Values", xlbl=xlbl, ylbl="K Input", xlims=xlims)

        # Plot #5 - Hist. of Matchup Uncertainty - Random
        ax = plt.subplot2grid((3, 3), (1, 1))
        ax.hist([Ur['sAdj'][idx] for idx in self.selMU],
                bins=10, color="grey")
        ax.set_title("Hist. of Matchup Uncertainty - Random")
        ax.set_xlabel("Standard Uncertainty")
        ax.set_ylabel("Frequency")

        # Plot #6 - Hist. of Adjustment Uncertainty - Systematic
        ax = plt.subplot2grid((3, 3), (1, 2))
        ax.hist([Us['sAdj'][idx] for idx in self.selMU],
                bins=10, color="grey")
        ax.set_title("Hist. of Spectral Adjustment Uncertainty - Systematic")
        ax.set_xlabel("Standard Uncertainty")
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

        inputPath = r"C:\Users\seh2\data\Harmonisation\m02_n15.nc"
        saveDirectory = r"C:\Users\seh2\data\Harmonisation\Graphs"

        hdata = HData(inputPath)

        # Create Plot class object
        Figure = Plot(saveDirectory, hdata.sensor1, hdata.sensor2, hdata.n_mu, hdata.timedata)

        # Plot data
        Figure.plotAll(hdata.H, hdata.Ur, hdata.Us, idx='time')
        Figure.plotAll(hdata.H, hdata.Ur, hdata.Us)

        return 0

    sample_code()