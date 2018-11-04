"""
Harmonsation Output Plotting
"""

'''___Built-In Modules___'''

'''___Third Party Modules___'''
from numpy import concatenate

'''___Harmonisation Modules___'''
from harmonisation.version import __version__, __tag__
from common import *
from harmonisation import MatchUp, HarmonisationResult, HarmonisationVis


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "12/12/2017"
__version__ = __version__
__tag__ = __tag__
__credits__ = ["Jon Mittaz"]
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class HarmonisationPlottingOp:
    """
    Class to run harmonisation result plotting routine

    Sample Code:

    .. code-block:: python

        H = HarmOp(dataset_dir, sensor_data, output_dir, software_cfg)
        H.run()

    :Attributes:
        .. py:attribute:: dataset_dir

            *list:str*

            Match-up dataset directory

        .. py:attribute:: sensor_data

            *dict*

            Sensor data dictionary

        .. py:attribute:: output_dir

            *str*

            Harmonisation output directory path

        .. py:attribute:: software_cfg

            *dict*

            Software configuration directory


    :Methods:
        .. py:method:: run(...):

            Run dataset plotting

    """

    def __init__(self, dataset_dir, sensor_data, output_dir, software_cfg):
        """
        Initialise harmonisation algorithm class

        :type dataset_dir: list
        :param dataset_dir: Match-up dataset directory

        :type sensor_data: dict
        :param sensor_data: Sensor data dictionary

        :type output_dir: str
        :param output_dir: Harmonisation output directory path

        :type software_cfg: dict
        :param software_cfg: Software configuration directory
        """

        self.dataset_dir = None
        self.sensor_data = None
        self.output_dir = None
        self.plot_dir = None
        self.software = None
        self.software_version = None
        self.software_tag = None
        self.job_id = None
        self.matchup_dataset = None
        self.HarmData = None

        if (dataset_dir is not None) and (sensor_data is not None):
            self.dataset_dir = dataset_dir
            self.sensor_data = sensor_data
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

    def run(self):
        """
        Generates plots for dataset specified by paths defined object attributes
        """

        # Initialise
        dataset_dir = self.dataset_dir
        sensor_data = self.sensor_data
        output_dir = self.output_dir
        plots_dir = self.plots_dir

        ################################################################################################################
        # 1.	Read Match-up Data and Harmonisation Result Data
        ################################################################################################################

        print"Match-up Dataset:"
        print ">", dataset_dir

        print"\nSensor Data:"
        #print ">", sensor_data_path

        print("\nOpening Files...")
        MatchUpData = MatchUp(dataset_dir, sensor_data, open_uncertainty=False)
        HarmResult = HarmonisationResult(output_dir, open_residual=False)

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
