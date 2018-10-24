"""
Harmonsation Output Plotting Implementation
"""

'''___Built-In Modules___'''
from os import makedirs

'''___Third Party Modules___'''
from numpy import concatenate

'''___Harmonisation Modules___'''
from common import *
from harmonisation import MatchUp, HarmonisationResult, HarmonisationVis


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "12/12/2017"
__credits__ = ["Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class HarmonisationPlottingOp:
    """
    This class runs the harmonisation process of Peter Harris for satellite sensor match-up data contained in the
    input directory

    Sample Code:

    .. code-block:: python

        H = HarmOp(directory)
        H.run()

    :Attributes:
        .. py:attribute:: dataset_dir

            *list:str*

            Directory of matchup series files

        .. py:attribute:: parameter_path

            *str*

            Path of parameters to be used as starting point for solver

        .. py:attribute:: output_dir

            *str*

            Path of directory to store output data files in

        .. py:attribute:: sensor_model

            *func*

            Python function to calculate radiance and derivatives given input sensor state data

        .. py:attribute:: adjustment_model

            *func*

            Python function to calculate spectral adjustment factor between two sensors, k, and derivatives given the
            two sensor radiances

        .. py:attribute:: software

            *str*

            software implementation name

        .. py:attribute:: software_version

            *str*

            software implementation version

        .. py:attribute:: software_tag

            *str*

            software implementation vcs

        .. py:attribute:: job_id

            *str*

            job configuratoin ID

        .. py:attribute:: matchup_dataset

            *str*

            harmonisation dataset set name

        .. py:attribute:: data_reader

            *cls*

            Harmonisation data reader

        .. py:attribute:: hout_path

            *str*

            path to store harmonisation output file

        .. py:attribute:: hres_paths

            *str*

            path to store harmonisation residual files

    :Methods:
        .. py:method:: run(...):

            This function runs the harmonisation of satellite instrument calibration parameters for group of sensors
            with a reference sensor from the match-up data located in the input directory

        . py:method:: calculate_parameter_covariance_ij(...):

            Return element of harmonisation output parameter covariance matrix

        .. py:method:: combine_parameter_covariance_ij(...):

            Combine parameter covariance elements to form full matrix and save to file

    """

    def __init__(self, dataset_dir, sensor_data_path, output_dir, software_cfg, hout_path, hres_paths):
        """
        Initialise harmonisation algorithm class

        :type dataset_paths: list:str
        :param dataset_paths: Paths of matchup series files in matchup dataset directory

        :type sensor_data_path: str
        :param sensor_data_path: Path of sensor data for of matchup sensors

        :type output_dir: str
        :param output_dir: Path of directory to store output data files in

        :type software_cfg: dict:str
        :param software_cfg: dictionary of software configuration information

        :type hout_path: str
        :param hout_path: path of harmonisation output file

        :type hres_paths: str
        :param hres_paths: path of harmonisation residual files
        """

        self.dataset_dir = None
        self.sensor_data_path = None
        self.output_dir = None
        self.plot_dir = None
        self.software = None
        self.software_version = None
        self.software_tag = None
        self.job_id = None
        self.matchup_dataset = None
        self.HarmData = None
        self.hout_path = None
        self.hres_paths = None

        if (dataset_dir is not None) and (sensor_data_path is not None):
            self.dataset_dir = dataset_dir
            self.sensor_data_path = sensor_data_path
            self.output_dir = output_dir
            self.plots_dir = pjoin(output_dir, "plots")
            try:
                makedirs(self.plots_dir)
            except OSError:
                pass

        if software_cfg is not None:
            #self.software = software_cfg['software']
            #self.software_version = software_cfg['version']
            #self.software_tag = software_cfg['tag']
            self.job_id = software_cfg["job_id"]
            self.matchup_dataset = software_cfg['matchup_dataset']

        else:
            self.software = "MM"
            self.software_version = "V.V"
            self.software_tag = "TTTTTTT"
            self.job_id = "CC"
            self.matchup_dataset = "TEST"

        if hout_path is not None:
            self.hout_path = hout_path
            self.hres_paths = hres_paths

    def run(self):
        """
        Generates plots for dataset specified by paths defined object attributes
        """

        # Initialise
        dataset_dir = self.dataset_dir
        sensor_data_path = self.sensor_data_path
        plots_dir = self.plots_dir
        hout_path = self.hout_path
        hres_paths = self.hres_paths

        ################################################################################################################
        # 1.	Read Match-up Data and Harmonisation Result Data
        ################################################################################################################

        print"Match-up Dataset:"
        print ">", dataset_dir

        print"\nSensor Data:"
        print ">", sensor_data_path

        print "\nHarmonisation Output Data:"
        print ">", hout_path
        for path in hres_paths:
            print ">", path

        print("\nOpening Files...")
        MatchUpData = MatchUp(dataset_dir, sensor_data_path, open_uncertainty=False)
        HarmResult = HarmonisationResult(hout_path) # , hres_paths)

        ################################################################################################################
        # 2.	Make plots
        ################################################################################################################

        print("Generating Plots...")

        HarmVisOp = HarmonisationVis(MatchUpData, HarmResult)

        # Define plot directories
        overview_plots_dir = pjoin(plots_dir, "overview")
        diagnostic_plots_dir = pjoin(plots_dir, "diagnostics")

        diagnostic_kres_harmonised_plots_dir = pjoin(diagnostic_plots_dir, "kres_harmonised")
        diagnostic_kres_harmonised_X_plots_dir = pjoin(diagnostic_kres_harmonised_plots_dir, "sensor_state_variables")
        diagnostic_kres_harmonised_x_plots_dir = pjoin(diagnostic_kres_harmonised_plots_dir, "additional_variables")

        diagnostic_kres_nom_harm_plots_dir = pjoin(diagnostic_plots_dir, "kres_nom_harm")
        diagnostic_kres_nom_harm_X_plots_dir = pjoin(diagnostic_kres_nom_harm_plots_dir, "sensor_state_variables")
        diagnostic_kres_nom_harm_x_plots_dir = pjoin(diagnostic_kres_nom_harm_plots_dir, "additional_variables")

        # Make plot directories
        try_makedirs(overview_plots_dir)
        try_makedirs(diagnostic_plots_dir)
        try_makedirs(diagnostic_kres_harmonised_plots_dir)
        try_makedirs(diagnostic_kres_harmonised_X_plots_dir)
        try_makedirs(diagnostic_kres_harmonised_x_plots_dir)
        try_makedirs(diagnostic_kres_nom_harm_plots_dir)
        try_makedirs(diagnostic_kres_nom_harm_X_plots_dir)
        try_makedirs(diagnostic_kres_nom_harm_x_plots_dir)

        # 1. Overview plots
        # a. Kres vs. time

        # Determine plot limits
        kres_nom = concatenate(HarmVisOp.get_kres_nom_nom_musamples())
        kres = concatenate(HarmVisOp.get_kres_harm_harm_musamples())

        kres_max = max([max(kres), max(kres_nom), abs(min(kres_nom)), abs(min(kres))])
        kres_ylim = [-kres_max * 1.2, kres_max * 1.2]

        # Plot
        HarmVisOp.plot_kres_nom_nom_scatter(overview_plots_dir) #, ylim=kres_ylim)
        HarmVisOp.plot_kres_nom_nom_scatter_monthlymean(overview_plots_dir) #, ylim=kres_ylim)
        HarmVisOp.plot_kres_harm_harm_scatter(overview_plots_dir) #, ylim=kres_ylim)
        HarmVisOp.plot_kres_harm_harm_scatter_monthlymean(overview_plots_dir) #, ylim=kres_ylim)
        HarmVisOp.plot_kres_nom_harm_scatter(overview_plots_dir)
        HarmVisOp.plot_kres_nom_harm_scatter_monthlymean(overview_plots_dir)

        # b. Nominal L vs. Nominal L1
        HarmVisOp.plot_L1_nom_v_L2_nom_scatter(overview_plots_dir)
        HarmVisOp.plot_L1_harm_v_L2_harm_scatter(overview_plots_dir)

        # c. Kres vs. sensor state variables, X
        HarmVisOp.plot_kres_harm_harm_X_binned_line(diagnostic_kres_harmonised_X_plots_dir)
        HarmVisOp.plot_kres_nom_harm_X_binned_line(diagnostic_kres_nom_harm_X_plots_dir)
        HarmVisOp.plot_kres_harm_harm_additional_variables_binned_line(diagnostic_kres_harmonised_x_plots_dir)
        HarmVisOp.plot_kres_nom_harm_additional_variables_binned_line(diagnostic_kres_nom_harm_x_plots_dir)

        print "\nPlots written to:", plots_dir

if __name__ == "__main__":
    pass
