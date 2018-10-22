"""
Module containing harmonisation result plotting functions
"""

'''___Python Modules____'''
import sys
from os.path import dirname
from os.path import join as pjoin
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from numpy import isnan, sort, arange, sort
from numpy.random import rand

above_nplcore_directory = dirname(dirname(dirname(dirname(__file__))))
from core.utils.plots import plot_scatter, plot_density
from core.utils.stats import time_average, bin

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "12/12/2017"
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

MAX_SCATTER_SAMPLES = 10000


class MatchUpVis:

    def __init__(self, MatchUpData = None):

        self.MatchUpData = MatchUpData

        self.Xmaxs = None
        self.Xmins = None

        self.L1_nom = None
        self.L2_nom = None
        self.kres_nom_nom = None
        self.times = None

        self.i_musamples = None
        self.musamples_total = None
        self.mu_total = None

        self.mu_labels = None
        self.mu_fname_suffix = None

        self.kres_nom_nom_musamples = None
        self.times_musamples = None
        self.L1_nom_musamples = None
        self.L2_nom_musamples = None

        self.kres_nom_nom_monthlymean = None
        self.kres_nom_nom_monthlymean_months = None
        self.months = None

    def get_L1_nom(self):
        if self.L1_nom is None:
            if self.MatchUpData is not None:
                try:
                    self.L1_nom = self.MatchUpData.getAdditionalData("nominal_BT1")
                except:
                    raise ValueError("Nominal Measurand 1 data not available in Matchup data")
            else:
                raise ValueError("Matchup data not available")
        return self.L1_nom

    def get_L2_nom(self):
        if self.L2_nom is None:
            if self.MatchUpData is not None:
                try:
                    self.L2_nom = self.MatchUpData.getAdditionalData("nominal_BT2")
                except:
                    raise ValueError("Nominal Measurand 2 data not available in Matchup data")
            else:
                raise ValueError("Matchup data not available")
        return self.L2_nom

    def get_kres_nom_nom(self):

        if self.kres_nom_nom is None:
            self.kres_nom_nom = self.get_L2_nom() - self.get_L1_nom() - self.MatchUpData.ks
            # try:
            #     self.kres_nom_nom = self.get_L2_nom() - self.get_L1_nom() - self.MatchUpData.ks
            # except ValueError:
            #     return None
        return self.kres_nom_nom

    def get_times(self):
        if self.times is None:
            try:
                self.times = self.MatchUpData.time2
            except AttributeError:
                raise ValueError("matchup time not available")
        return self.times

    def get_Xmax(self, cov):
        if self.Xmaxs is None:
            self.set_Xmax_Xmin()
        return self.Xmaxs[cov-1]

    def get_Xmin(self, cov):
        if self.Xmins is None:
            self.set_Xmax_Xmin()
        return self.Xmins[cov-1]

    def get_X_label(self, cov, sensor_name):

        if "sensor_model_variable_labels" in self.MatchUpData.sensor_data[sensor_name]:
            return self.MatchUpData.sensor_data[sensor_name]["sensor_model_variable_labels"][cov - 1]
        return "X" + str(cov)

    def get_X_name(self, cov, sensor_name):

        if "sensor_model_variable_names" in self.MatchUpData.sensor_data[sensor_name]:
            return self.MatchUpData.sensor_data[sensor_name]["sensor_model_variable_names"][cov - 1]
        return "X" + str(cov)

    def set_Xmax_Xmin(self):

        Xmaxs = []
        Xmins = []
        for cov in range(1, max(self.MatchUpData.idx['sensor_ms']) + 1):
            Xmax_i = []
            Xmin_i = []
            for i, mu in enumerate(range(1, max(self.MatchUpData.idx['n_mu']) + 1)):
                sensor_2_name = self.MatchUpData.idx['sensors'][self.MatchUpData.idx['Im'][mu - 1][1]]
                sensor_2_m = self.MatchUpData.idx['sensor_ms'][self.MatchUpData.idx['Im'][mu - 1][1]]

                if cov <= sensor_2_m:
                    Xmax_i.append(self.MatchUpData.getVariableDataMax(cov, sensor_2_name, mu))
                    Xmin_i.append(self.MatchUpData.getVariableDataMin(cov, sensor_2_name, mu))

            Xmaxs.append(max(Xmax_i))
            Xmins.append(min(Xmin_i))

        self.Xmaxs = Xmaxs
        self.Xmins = Xmins
        return 0

    def get_L1_nom_musamples(self):
        if self.L1_nom_musamples is None:
            self.L1_nom_musamples = self.musample_values(self.get_L1_nom())
        return self.kres_nom_nom_musamples

    def get_L2_nom_musamples(self):
        if self.L2_nom_musamples is None:
            self.L2_nom_musamples = self.musample_values(self.get_L2_nom())
        return self.L2_nom_musamples

    def get_kres_nom_nom_musamples(self):
        if self.kres_nom_nom_musamples is None:
            self.kres_nom_nom_musamples = self.musample_values(self.get_kres_nom_nom())
        return self.kres_nom_nom_musamples

    def get_times_musamples(self):
        if self.times_musamples is None:
            self.times_musamples = self.musample_values(self.get_times())
        return self.times_musamples

    def get_kres_nom_nom_monthlymean(self):
        if self.kres_nom_nom_monthlymean is None:
            self.kres_nom_nom_monthlymean,\
            self.kres_nom_nom_monthlymean_months = self.monthlymean_values(self.get_kres_nom_nom())
        return self.kres_nom_nom_monthlymean, self.kres_nom_nom_monthlymean_months

    def get_i_musamples(self):

        if self.i_musamples is None:
            if self.MatchUpData is not None:
                total = len(self.get_kres_nom_nom())
                sample = False
                if total > MAX_SCATTER_SAMPLES:
                    sample = True
                    sample_ratio = float(MAX_SCATTER_SAMPLES) / total

                self.i_musamples = []
                for i, (i_start, i_end) in enumerate(zip(self.MatchUpData.idx["cNm"][:-1], self.MatchUpData.idx["cNm"][1:])):

                    total_i = i_end - i_start
                    if sample:
                        n_samples = int(total_i * sample_ratio)
                        self.i_musamples.append(sort((total_i * rand(n_samples)).astype(int)))
                    else:
                        self.i_musamples.append(arange(total_i))

            else:
                return None
        return self.i_musamples

    def get_musamples_total(self):
        if self.musamples_total is None:
            if self.i_musamples is None:
                return None
            else:
                self.musamples_total = int(sum([len(i_musample) for i_musample in self.get_i_musamples()]))
        return self.musamples_total

    def get_mu_total(self):
        if self.mu_total is None:
            if self.MatchUpData is None:
                return None
            else:
                self.mu_total = self.MatchUpData.idx['cNm'][-1]
        return self.mu_total

    def get_mu_labels(self):
        if self.mu_labels is None:
            if self.MatchUpData is not None:
                self.mu_labels = []
                for i in range(len(self.MatchUpData.idx['Nm'])):
                    sensor_1_name = str(self.MatchUpData.idx['sensors'][self.MatchUpData.idx['Im'][i][0]])
                    sensor_2_name = str(self.MatchUpData.idx['sensors'][self.MatchUpData.idx['Im'][i][1]])
                    matchup_label = sensor_1_name + " - " + sensor_2_name
                    self.mu_labels.append(matchup_label)
            else:
                raise ValueError("Match Up data required to generate plotting labels")
        return self.mu_labels

    def get_mu_fname_suffix(self):
        if self.mu_fname_suffix is None:
            if self.MatchUpData is not None:
                self.mu_fname_suffix = []
                for i in range(len(self.MatchUpData.idx['Nm'])):
                    sensor_1_name = str(self.MatchUpData.idx['sensors'][self.MatchUpData.idx['Im'][i][0]])
                    sensor_2_name = str(self.MatchUpData.idx['sensors'][self.MatchUpData.idx['Im'][i][1]])
                    matchup_label = sensor_1_name + "_" + sensor_2_name
                    self.mu_fname_suffix.append(matchup_label)
            else:
                raise ValueError("Match Up data required to generate plotting labels")
        return self.mu_fname_suffix

    def musample_values(self, values):

        values_musamples = []
        i_musamples = self.get_i_musamples()
        for i_start, i_end, i_musample in zip(self.MatchUpData.idx["cNm"][:-1], self.MatchUpData.idx["cNm"][1:], i_musamples):
            values_musamples.append(values[i_start:i_end][i_musample])
        return values_musamples

    def monthlymean_values(self, values):

        averages = []
        average_times = []
        # 2. Plot by match-up series
        for i_start, i_end in zip(self.MatchUpData.idx["cNm"][:-1], self.MatchUpData.idx["cNm"][1:]):

            # Bin monthly means for match-up series (match-up series data range)
            averages_i, average_times_i = time_average(values[i_start:i_end],
                                                       self.get_times()[i_start:i_end],
                                                       time_period="month")

            averages.append(averages_i[~isnan(averages_i)])
            average_times.append(average_times_i[~isnan(averages_i)])

        return averages, average_times

    def plot_kres_nom_nom_scatter(self, directory, ylim=None):
        plot_scatter(pjoin(directory, "k_res_nominal_scatter.pdf"),
                     self.get_kres_nom_nom_musamples(),
                     self.get_times_musamples(),
                     xlbl="Time", ylbl= "Nominal K Residual (" + str(int(self.get_musamples_total())) + " Samples)",
                     title="", legend=self.get_mu_labels(), legend_loc="outside",
                     solid_ylines=[0], fig_size=(12, 8), ylim=ylim)
        return 0

    def plot_kres_nom_nom_scatter_monthlymean(self, directory, ylim=None):
        kres_nom_nom_monthlymean, months = self.get_kres_nom_nom_monthlymean()
        plot_scatter(pjoin(directory, "k_res_nominal_monthlymean.pdf"),
                     kres_nom_nom_monthlymean, months,
                     xlbl="Time", ylbl="Nominal K Residual (Monthly Mean)",
                     title="", legend=self.get_mu_labels(), legend_loc="outside",
                     solid_ylines=[0], fig_size=(12, 8), ylim=ylim, line="-")
        return 0

    def plot_kres_nom_nom_density(self, directory, ylim=None):
        plot_density(pjoin(directory, "k_res_nominal_density.pdf"),
                     self.get_kres_nom_nom_musamples(),
                     self.get_times_musamples(),
                     bins=200, title="", ylbl="Nominal K Residual", xlbl="Time",
                     txt=["M = " + str(self.get_mu_total()) + " "], solid_ylines=[0], ylim=ylim)
        return 0

    def plot_L1_nom_v_L2_nom_scatter(self, directory):
        plot_scatter(pjoin(self, directory, "L1_vs_L2_nominal_scatter.pdf"),
                     self.get_L1_nom_musamples(),
                     self.get_L2_nom_musamples(),
                     xlbl="Nominal Sensor 2 Measurand",
                     ylbl="Nominal Sensor 1 Measurand (" + str(self.get_musamples_total()) + " Samples)", title="",
                     yerr=None, xerr=None, legend=self.get_mu_labels(), legend_loc="outside",
                     solid_ylines=[0], fig_size=(12, 8))
        return 0

    def plot_L1_nom_v_L2_nom_density(self, directory):
        plot_density(pjoin(directory, "L1_vs_L2_nominal_density.pdf"),
                     self.get_L1_nom(),
                     self.get_L2_nom(),
                     bins=400,
                     title="", ylbl="Nominal Sensor 1 Measurand", xlbl="Nominal Sensor 2 Measurand",
                     txt=["M = " + str(self.get_mu_total()) + " "], solid_ylines=[0])
        return 0

    def plot_Y_X_binned_line(self, directory, Y, title, fname, ylim=None, mus=None):

        MatchUpData = self.MatchUpData

        # 1. Initialise
        matchup_labels = self.get_mu_labels()
        matchup_fname_suffix = self.get_mu_fname_suffix()

        if mus is None:
            mus = list(range(1, max(MatchUpData.idx['n_mu']) + 1))

        # 2. Plot by X by match-up series
        # Loop through X
        for cov in range(1, max(MatchUpData.idx['sensor_ms']) + 1):

            Xmin = self.get_Xmin(cov)
            Xmax = self.get_Xmax(cov)

            # Initialise data lists
            averages = []
            Xaverages = []

            # Plot per match-up series
            for i, mu in enumerate(mus):

                # Indices and labels
                i1 = MatchUpData.idx["cNm"][mu - 1]
                i2 = MatchUpData.idx["cNm"][mu]
                sensor_2_name = MatchUpData.idx['sensors'][MatchUpData.idx['Im'][mu - 1][1]]
                cov_name = self.get_X_name(cov, sensor_2_name)
                cov_label = "$"+self.get_X_label(cov, sensor_2_name) + "$ (Binned)"
                mu_fn = matchup_fname_suffix[mu-1]

                # Plot if sensor 2 of given match-up series has an cov-th X
                sensor_2_m = MatchUpData.idx['sensor_ms'][MatchUpData.idx['Im'][mu - 1][1]]
                if cov <= sensor_2_m:
                    X = MatchUpData.getVariableData(cov, sensor_2_name, mu)

                    # Bin data by match-up series limits
                    averages_i, Xaverages_i = bin(Y[i1:i2], X, bins=50, X_range=[min(X), max(X)])

                    # Plot match-up series
                    fname_i = pjoin(directory, fname+"." + cov_name + "_" + sensor_2_name + "_binned." + mu_fn + ".pdf")
                    plot_scatter(fname_i, averages_i[~isnan(averages_i)], Xaverages_i[~isnan(averages_i)],
                                 xlbl=cov_label+" - "+sensor_2_name, ylbl=title+" ("+matchup_labels[i]+")",
                                 title="", yerr=None, xerr=None, txt=None,
                                 line='-', legend=None, solid_ylines=[0], ylim=ylim)

                    # Bin data by full limits
                    averages_i, Xaverages_i = bin(Y[i1:i2], X, bins=50, X_range=[Xmin, Xmax])
                    averages.append(averages_i[~isnan(averages_i)])
                    Xaverages.append(Xaverages_i[~isnan(averages_i)])

            # Plot X for all match-up series
            plot_scatter(pjoin(directory, fname + "." + cov_name + "_sensor_2_binned.pdf"), averages, Xaverages,
                         xlbl=cov_label + " - Sensor 2", ylbl=title, title="",
                         yerr=None, xerr=None, txt=None, line='-', legend=matchup_labels, legend_loc="outside",
                         solid_ylines=[0], fig_size=(12, 8), ylim=ylim)

    def plot_Y_x_binned_line(self, directory, Y, x, x_name, x_label, title, fname, ylim=None, mus=None):
        """
        Output plot to given directory for the difference between data and estimated k against variable x

        :type directory: str
        :param directory: directory to write plots to

        :type MatchUpData: eopy.matchup.matchupIO.MatchUp
        :param MatchUpData: MatchUp data object to plot

        :type k_est: numpy.ndarray
        :param k_est: estimate of k to plot

        :type x: numpy.ndarray
        :param x: variable to plot against

        :type x_name: str
        :param x_name: name for file for x variable

        :type x_label: str
        :param x_label: label for plotting x variable

        :type mus: list
        :param mus: match-up series to plot for
        """

        # 1. Preamble
        # a. Requried data
        # - Full range of x-variable
        MatchUpData = self.MatchUpData

        # 1. Initialise
        matchup_labels = self.get_mu_labels()
        matchup_fname_suffix = self.get_mu_fname_suffix()

        xmin = min(x)
        xmax = max(x)

        # b. Initialise data lists
        averages = []
        averages_x = []

        if mus is None:
            mus = list(range(1, max(MatchUpData.idx['n_mu'])+1))

        # 2. Plot by match-up series
        for i, mu in enumerate(mus):

            # Indices and labels
            i1 = MatchUpData.idx["cNm"][mu - 1]
            i2 = MatchUpData.idx["cNm"][mu]
            sensor_2_name = MatchUpData.idx['sensors'][MatchUpData.idx['Im'][mu - 1][1]]
            mu_fn = matchup_fname_suffix[mu - 1]

            # Bin by x for match-up series (match-up series data range)
            x_i = x[i1:i2]
            x_i_range = [min(x_i), max(x_i)]
            averages_i, averages_x_i = bin(Y[i1:i2], x_i, bins=50, X_range=x_i_range)

            # Plot match-up series data
            fname_i = pjoin(directory, fname + "." + x_name + "_" + sensor_2_name + "_binned." + mu_fn + ".pdf")
            plot_scatter(fname_i, averages_i[~isnan(averages_i)], averages_x_i[~isnan(averages_i)],
                         xlbl=x_label + " - " + sensor_2_name, ylbl=title + " (" + matchup_labels[i] + ")",
                         title="", yerr=None, xerr=None, txt=None,
                         line='-', legend=None, solid_ylines=[0], ylim=ylim)

            # Bin by x for match-up series (all data range)
            averages_i, averages_x_i = bin(Y[i1:i2], x_i, bins=50, X_range=[xmin, xmax])
            averages.append(averages_i[~isnan(averages_i)])
            averages_x.append(averages_x_i[~isnan(averages_i)])

        # 3. Plot all match-up series together
        plot_scatter(pjoin(directory, fname + "." + x_name + "_sensor_2_binned.pdf"), averages, averages_x,
                     xlbl=x_label + " - Sensor 2", ylbl=title, title="",
                     yerr=None, xerr=None, txt=None, line='-', legend=matchup_labels, legend_loc="outside",
                     solid_ylines=[0], fig_size=(12,8), ylim=ylim)

        return 0

if __name__ == "__main__":
    pass
