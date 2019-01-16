"""
Module containing harmonisation result plotting functions
"""

'''___Python Modules____'''
from os.path import join as pjoin
from copy import deepcopy

'''__Third Party Modules___'''
import matplotlib
matplotlib.use('Agg')  # Force matplotlib to not use any Xwindows backend.

'''___harmonisation Modules___'''
from harmonisation.core.matchup.matchupToolbox.utils.matchup_maths import evaluate_K, evaluate_adjusted_measurand, evaluate_measurand
from harmonisation.core.matchup.matchupProcessing.simulate_matchup.SimulateMatchUp import SimulateMatchUp
from harmonisation.core.matchup.matchupVis.MatchUpVis import MatchUpVis
from harmonisation.core.utils.plots import plot_scatter, plot_2dhist, plot_density


'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "12/12/2017"
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class HarmonisationVis(MatchUpVis):

    def __init__(self, MatchUpData=None, HarmResult=None):

        self.MatchUpData = MatchUpData
        self.HarmResult = HarmResult

        self.Xmaxs = None
        self.Xmins = None

        self.MatchUpSim = None

        self.L1_nom = None
        self.L2_nom = None
        self.kres_nom_nom = None
        self.times = None

        self.L1_harm = None
        self.L2_harm = None
        self.kres_harm_harm = None

        self.kres_nom_harm = None

        self.i_musamples = None
        self.musamples_total = None
        self.mu_total = None
        self.mu_labels = None
        self.mu_fname_suffix = None

        self.kres_nom_nom_musamples = None
        self.times_musamples = None
        self.L1_nom_musamples = None
        self.L2_nom_musamples = None

        self.kres_harm_harm_musamples = None
        self.L1_harm_musamples = None
        self.L2_harm_musamples = None

        self.kres_nom_nom_monthlymean = None
        self.kres_nom_nom_monthlymean_months = None

        self.kres_harm_harm_monthlymean = None
        self.kres_harm_harm_monthlymean_months = None

        self.kres_nom_harm_musamples = None
        self.kres_nom_harm_monthlymean = None
        self.kres_nom_harm_monthlymean_months = None

    def get_L1_harm(self):
        if self.L1_harm is None:
            if self.MatchUpData is not None:
                try:
                    old_a = deepcopy(self.MatchUpData.a)
                    old_parameter_sensors = deepcopy(self.MatchUpData.idx['parameter_sensor'])
                    self.MatchUpData.a = deepcopy(self.HarmResult.parameter)
                    self.MatchUpData.idx['parameter_sensor'] = deepcopy(self.HarmResult.parameter_sensors)
                    L_est = evaluate_adjusted_measurand(self.MatchUpData)
                    self.L1_harm = L_est[:, 0]
                    self.L2_harm = L_est[:, 1]
                    self.MatchUpData.a = old_a
                    self.MatchUpData.idx['parameter_sensor'] = old_parameter_sensors
                except:
                    raise ValueError("Unable to evaluate harmonised radiances")
            else:
                raise ValueError("Matchup data not available")
        return self.L1_harm

    def get_L2_harm(self):
        if self.L2_harm is None:
            if self.MatchUpData is not None:
                try:
                    L_est = evaluate_adjusted_measurand(self.MatchUpData)
                    self.L1_harm = L_est[:, 0]
                    self.L2_harm = L_est[:, 1]
                except:
                    raise ValueError("Unable to evaluate harmonised radiances")
            else:
                raise ValueError("Matchup data not available")
        return self.L2_harm

    def get_kres_harm_harm(self):

        if self.kres_harm_harm is None:
            try:
                self.kres_harm_harm = self.get_L2_harm() - self.get_L1_harm() - self.MatchUpData.ks
            except:
                return ValueError("Couldn't evaluate harmonised kres")
        return self.kres_harm_harm

    def get_kres_nom_harm(self):

        if self.kres_nom_harm is None:
            try:
                self.kres_nom_harm = self.get_L2_nom() - self.get_L1_harm() - self.MatchUpData.ks
            except:
                return ValueError("Couldn't evaluate harmonised kres")
        return self.kres_nom_harm

    def get_L1_harm_musamples(self):
        if self.L1_harm_musamples is None:
            self.L1_harm_musamples = self.musample_values(self.get_L1_harm())
        return self.kres_harm_harm_musamples

    def get_L2_harm_musamples(self):
        if self.L2_harm_musamples is None:
            self.L2_harm_musamples = self.musample_values(self.get_L2_harm())
        return self.L2_harm_musamples

    def get_kres_harm_harm_musamples(self):
        if self.kres_harm_harm_musamples is None:
            self.kres_harm_harm_musamples = self.musample_values(self.get_kres_harm_harm())
        return self.kres_harm_harm_musamples

    def get_kres_harm_harm_monthlymean(self):
        if self.kres_harm_harm_monthlymean is None:
            self.kres_harm_harm_monthlymean, \
            self.kres_harm_harm_monthlymean_months = self.monthlymean_values(self.get_kres_harm_harm())
        return self.kres_harm_harm_monthlymean, self.kres_harm_harm_monthlymean_months

    def get_kres_nom_harm_musamples(self):
        if self.kres_nom_harm_musamples is None:
            self.kres_nom_harm_musamples = self.musample_values(self.get_kres_nom_harm())
        return self.kres_nom_harm_musamples

    def get_kres_nom_harm_monthlymean(self):
        if self.kres_nom_harm_monthlymean is None:
            self.kres_nom_harm_monthlymean, \
            self.kres_nom_harm_monthlymean_months = self.monthlymean_values(self.get_kres_nom_harm())
        return self.kres_nom_harm_monthlymean, self.kres_nom_harm_monthlymean_months

    def plot_kres_harm_harm_scatter(self, directory, ylim=None):
        plot_scatter(pjoin(directory, "k_res_harmonised_scatter.pdf"),
                     self.get_kres_harm_harm_musamples(),
                     self.get_times_musamples(),
                     xlbl="Time", ylbl="Harmonised K Residual (" + str(int(self.get_musamples_total())) + " Samples)",
                     title="", legend=self.get_mu_labels(), legend_loc="outside",
                     solid_ylines=[0], fig_size=(12, 8), ylim=ylim)
        return 0

    def plot_kres_harm_harm_scatter_monthlymean(self, directory, ylim=None):
        kres_harm_harm_monthlymean, months = self.get_kres_harm_harm_monthlymean()
        plot_scatter(pjoin(directory, "k_res_harmonised_monthlymean.pdf"),
                     kres_harm_harm_monthlymean, months,
                     xlbl="Time", ylbl="Harmonised K Residual (Monthly Mean)",
                     title="", legend=self.get_mu_labels(), legend_loc="outside",
                     solid_ylines=[0], fig_size=(12, 8), ylim=ylim, line="-")
        return 0

    def plot_kres_harm_harm_density(self, directory, ylim=None):
        plot_density(pjoin(directory, "k_res_harmonised_density.pdf"),
                     self.get_kres_harm_harm_musamples(),
                     self.get_times_musamples(),
                     bins=200, title="", ylbl="Harmonised K Residual", xlbl="Time",
                     txt=["M = " + str(self.get_mu_total()) + " "], solid_ylines=[0], ylim=ylim)
        return 0

    def plot_kres_nom_harm_scatter(self, directory, ylim=None):
        plot_scatter(pjoin(directory, "k_res_nom_harm_scatter.pdf"),
                     self.get_kres_nom_harm_musamples(),
                     self.get_times_musamples(),
                     xlbl="Time", ylbl="K = K - (nominal-harmonised) (" + str(int(self.get_musamples_total())) + " Samples)",
                     title="", legend=self.get_mu_labels(), legend_loc="outside",
                     solid_ylines=[0], fig_size=(12, 8), ylim=ylim)
        return 0

    def plot_kres_nom_harm_scatter_monthlymean(self, directory, ylim=None):
        kres_nom_harm_monthlymean, months = self.get_kres_nom_harm_monthlymean()
        plot_scatter(pjoin(directory, "k_res_nom_harm_monthlymean.pdf"),
                     kres_nom_harm_monthlymean, months,
                     xlbl="Time", ylbl="K = K - (nominal-harmonised) (Monthly Mean)",
                     title="", legend=self.get_mu_labels(), legend_loc="outside",
                     solid_ylines=[0], fig_size=(12, 8), ylim=ylim, line="-")
        return 0

    def plot_kres_nom_harm_density(self, directory, ylim=None):
        plot_density(pjoin(directory, "k_res_nom_harm_density.pdf"),
                     self.get_kres_nom_harm_musamples(),
                     self.get_times_musamples(),
                     bins=200, title="", ylbl="K = K - (nominal-harmonised)", xlbl="Time",
                     txt=["M = " + str(self.get_mu_total()) + " "], solid_ylines=[0], ylim=ylim)
        return 0

    def plot_L1_harm_v_L2_harm_scatter(self, directory):
        plot_scatter(pjoin(directory, "L1_vs_L2_harmonised_scatter.pdf"),
                     self.get_L1_harm_musamples(),
                     self.get_L2_harm_musamples(),
                     xlbl="Harmonised Sensor 2 Measurand",
                     ylbl="Harmonised/Reference Sensor 1 Measurand (" + str(self.get_musamples_total()) + " Samples)",
                     title="", yerr=None, xerr=None, legend=self.get_mu_labels(), legend_loc="outside",
                     solid_ylines=[0], fig_size=(12, 8))
        return 0

    def plot_L1_harm_v_L2_harm_density(self, directory):
        plot_density(pjoin(directory, "L1_vs_L2_harmonised_density.pdf"),
                     self.get_L1_harm(),
                     self.get_L2_harm(),
                     bins=400,
                     title="", ylbl="Harmonised/Reference Sensor 1 Measurand", xlbl="Harmonised Sensor 2 Measurand",
                     txt=["M = " + str(self.get_mu_total()) + " "], solid_ylines=[0])
        return 0

    def plot_kres_harm_harm_X_binned_line(self, directory, ylim=None, mus=None):
        self.plot_Y_X_binned_line(directory, self.get_kres_harm_harm(),
                                  "Harmonised K Residual",  "k_res_harmonised", ylim, mus)

    def plot_kres_nom_harm_X_binned_line(self, directory, ylim=None, mus=None):
        self.plot_Y_X_binned_line(directory, self.get_kres_nom_harm(),
                                  "K Residual = K - (Nominal - Harmonised)", "k_res_nom_harm", ylim, mus)

    def plot_kres_harm_harm_additional_variables_binned_line(self, directory, ylim=None, mus=None):
        for x_name in self.MatchUpData.getAdditionalDataNames():
            self.plot_Y_x_binned_line(directory,
                                      self.get_kres_harm_harm(), self.MatchUpData.getAdditionalData(x_name),
                                      x_name, x_name, "Harmonised K Residual",  "k_res_harmonised",
                                      ylim, mus)

    def plot_kres_nom_harm_additional_variables_binned_line(self, directory, ylim=None, mus=None):
        for x_name in self.MatchUpData.getAdditionalDataNames():
            self.plot_Y_x_binned_line(directory,
                                      self.get_kres_nom_harm(), self.MatchUpData.getAdditionalData(x_name),
                                      x_name, x_name, "K Residual = K - (Nominal - Harmonised)", "k_res_nom_harm",
                                      ylim, mus)

    def get_MatchUpSim(self, n_sim=20, verbose=False):
        MatchUpSimOp = SimulateMatchUp()
        MatchUpSim = MatchUpSimOp.run(self.MatchUpData, n_samples_mu=n_sim, verbose=verbose)
        return MatchUpSim

    def plot_compare_calibration(self, directory, parameter_comparison, parameter_covariance_matrix_comparison, verbose=False):

        MatchUpSim = self.get_MatchUpSim(n_sim=20, verbose=verbose)
        measurand_1 = evaluate_measurand(MatchUpSim)

        MatchUpSim.a = parameter_comparison
        measurand_2, u_measurand_cov_2 = evaluate_measurand(MatchUpSim,
                                                  parameter_covariance_matrix=parameter_covariance_matrix_comparison)

        measurand_diff = measurand_1 - measurand_2

        # Just make one plot per sensor
        sensors = MatchUpSim.idx['sensors']
        plotted_sensors = []

        i_mu = 0
        while set(plotted_sensors) != set(sensors):
            istart = MatchUpSim.idx['cNm'][i_mu]
            iend = MatchUpSim.idx['cNm'][i_mu+1]

            pair = MatchUpSim.idx['Im'][i_mu]
            for i_pair, n_sensor in enumerate(pair):

                sensor_name = sensors[n_sensor]

                measurand_1_i = measurand_1[istart:iend, i_pair]
                measurand_diff_i = measurand_diff[istart:iend, i_pair]
                u_measurand_cov_2_i = u_measurand_cov_2[istart:iend, i_pair]

                print measurand_2[istart, i_pair]

                path = pjoin(directory, "compare_calibration." + sensor_name + ".pdf")
                plot_scatter(path, measurand_diff_i, measurand_1_i, yerr=2*u_measurand_cov_2_i,
                             title="Sampled Radiance Residual - " + str(sensor_name),
                             ylbl="Radiance Residual, $\Delta L = L_{Data} - L_{Est}$",
                             xlbl="Data Radiance, $L_{DATA}$",
                             dash_ylines=[0])
                plotted_sensors.append(sensor_name)
            i_mu += 1

if __name__ == "__main__":
    pass
