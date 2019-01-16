"""
Harmonsation Output Plotting
"""

'''___Built-In Modules___'''

'''___Third Party Modules___'''
from numpy import concatenate

'''___Harmonisation Modules___'''
from harmonisation.version import __version__, __tag__
from common import *
from harmonisation import open_matchup
from harmonisation.core.matchup.matchupToolbox.harmonisation.harmonisationIO.HarmonisationResult import HarmonisationResult
from harmonisation.core.matchup.matchupToolbox.harmonisation.harmonisationVis.HarmonisationVis import HarmonisationVis


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

        .. py:attribute:: plot_dir

            *str*

            Plot directory path

    :Methods:
        .. py:method:: run(...):

            Run dataset plotting

    """

    def __init__(self, dataset_dir, sensor_data, output_dir, logger=None):
        """
        Initialise harmonisation algorithm class

        :type dataset_dir: str
        :param dataset_dir: Match-up dataset directory

        :type sensor_data: dict
        :param sensor_data: Sensor data dictionary

        :type output_dir: str
        :param output_dir: Harmonisation output directory path

        :type logger: logging.Logger
        :param logger: Logger
        """

        self.dataset_dir = dataset_dir
        self.sensor_data = sensor_data
        self.output_dir = output_dir
        self.plots_dir = pjoin(output_dir, "plots")
        self.logger = logger

        try:
            makedirs(self.plots_dir)
        except OSError:
            pass

    def run(self):
        """
        Generates plots for dataset specified by paths defined object attributes
        """

        # 1. Read Match-up Data and Harmonisation Result Data
        self.logger.info("Match-up Dataset: "+self.dataset_dir)

        self.logger.debug("Opening Match-up Dataset...")
        MatchUpData = open_matchup(self.dataset_dir, open_uncertainty=False)
        MatchUpData.setSensorData(self.sensor_data)
        self.logger.debug("Complete")

        self.logger.info("Harmonisation Output Dataset: " + self.output_dir)

        self.logger.debug("Opening Harmonisation Result Dataset...")
        HarmResult = HarmonisationResult(self.output_dir, open_residuals=False)
        self.logger.debug("Complete")

        # 2.	Make plots
        self.logger.debug("Generating Plots...")
        HarmVisOp = HarmonisationVis(MatchUpData, HarmResult)

        # Define plot directories
        overview_plots_dir = pjoin(self.plots_dir, "overview")
        diagnostic_plots_dir = pjoin(self.plots_dir, "diagnostics")

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
        self.logger.debug("Complete")

        self.logger.info("Plots written to: "+self.plots_dir)


if __name__ == "__main__":
    pass
