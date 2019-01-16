"""
Harmonsation Comparison Output Plotting Implementation
"""

'''___Built-In Modules___'''

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from harmonisation.version import __version__, __tag__
from common import *
from harmonisation import open_matchup, HarmonisationResult, HarmonisationVis


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "12/12/2017"
__version__ = __version__
__tag__ = __tag__
__credits__ = ["Jon Mittaz"]
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def try_makedirs(directory):
    try:
        makedirs(directory)
    except OSError:
        pass
    return 0


class HarmonisationComparisonPlottingOp:
    """
    Plotting routines to compare the result of different harmonisation processes
    (or with true parameters in the case of simulated data)

    :Attributes:
        .. py:attribute:: dataset_paths

            *list:str*

            Paths of matchup series files in matchup dataset directory

        .. py:attribute:: parameter_path

            *str*

            Path of parameters to be used as starting point for solver

        .. py:attribute:: output_dir

            *str*

            Path of directory to store output data files in

        .. py:attribute:: sensor_data_path

            *func*

            Sensor data file path

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

            Run comparison plotting routine
    """

    def __init__(self, dataset_dir=None, sensor_data=None, output_dir=None):
        """
        Initialise harmonisation comparison plotting class

        :type dataset_dir: list:str
        :param dataset_dir: Match-up dataset directory

        :type sensor_data: str
        :param sensor_data: Path of sensor data for of match-up sensors, contains reference set of parameters

        :type output_dir: str
        :param output_dir: Path of directory to store output data files in
        """

        self.dataset_dir = None
        self.sensor_data = None
        self.output_dir = None
        self.plot_dir = None
        self.matchup_dataset = None
        self.HarmData = None
        self.hout_path = None
        self.hres_paths = None

        if (dataset_dir is not None) and (sensor_data is not None):
            self.dataset_dir = dataset_dir
            self.sensor_data = sensor_data
            self.output_dir = output_dir
            self.plots_dir = pjoin(output_dir, "plots")
            try:
                makedirs(self.plots_dir)
            except OSError:
                pass

    def run(self, verbose=True):
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
        print ">", sensor_data

        print "\nHarmonisation Output:"
        print ">", output_dir

        print("\nOpening Files...")
        MatchUpData = open_matchup(dataset_dir, open_uncertainty=False)
        MatchUpData.setSensorData(sensor_data)
        HarmResult = HarmonisationResult(output_dir, open_residuals=False)

        if verbose:
            print HarmResult.parameter_sensors
            print HarmResult.parameter
            print HarmResult.parameter_covariance_matrix
            print "\n"

        ################################################################################################################
        # 2.	Make plots
        ################################################################################################################

        print("Generating Plots...")

        HarmVisOp = HarmonisationVis(MatchUpData, HarmResult)
        HarmVisOp.plot_compare_calibration(plots_dir, HarmResult.parameter,
                                           HarmResult.parameter_covariance_matrix, verbose=verbose)

        print "\nPlots written to:", plots_dir


if __name__ == "__main__":
    pass

