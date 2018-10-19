"""
Harmonsation Comparison Output Plotting Implementation
"""

'''___Python Modules___'''
from os import makedirs
from os.path import basename, dirname
from os.path import join as pjoin
import sys
from sys import argv

'''___Third Party Modules___'''
from numpy import concatenate

'''___Harmonisation Modules___'''
from common import *

main_directory = dirname(dirname(__file__))
sys.path.append(main_directory)
from nplcore import HarmonisationResult
from nplcore import HarmonisationVis

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "12/12/2017"
__credits__ = ["Jon Mittaz"]
__version__ = "0.0"
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

    def __init__(self, dataset_paths=None, sensor_data_path=None, output_dir=None,
                 data_reader=None, hout_path=None, hres_paths=None):
        """
        Initialise harmonisation comparison plotting class

        :type dataset_paths: list:str
        :param dataset_paths: Paths of matchup series files in matchup dataset directory

        :type sensor_data_path: str
        :param sensor_data_path: Path of sensor data for of matchup sensors, contains reference set of parameters

        :type output_dir: str
        :param output_dir: Path of directory to store output data files in

        :type data_reader: cls
        :param data_reader: Python class to open harmonisation data. If none given default data reader used.

        :type hout_path: str
        :param hout_path: path of harmonisation output file

        :type hres_paths: str
        :param hres_paths: path of harmonisation residual files
        """

        self.dataset_paths = None
        self.sensor_data_path = None
        self.output_dir = None
        self.plot_dir = None
        self.matchup_dataset = None
        self.HarmData = None
        self.hout_path = None
        self.hres_paths = None

        if (dataset_paths is not None) and (sensor_data_path is not None):
            self.dataset_paths = dataset_paths
            self.sensor_data_path = sensor_data_path
            self.output_dir = output_dir
            self.plots_dir = pjoin(output_dir, "plots")
            try:
                makedirs(self.plots_dir)
            except OSError:
                pass

        if data_reader is not None:
            self.data_reader = data_reader

        else:
            from nplcore import MatchUp
            self.data_reader = MatchUp

        if hout_path is not None:
            self.hout_path = hout_path
            self.hres_paths = hres_paths

    def run(self, verbose=True):
        """
        Generates plots for dataset specified by paths defined object attributes
        """

        # Initialise
        dataset_paths = self.dataset_paths
        sensor_data_path = self.sensor_data_path
        plots_dir = self.plots_dir
        hout_path = self.hout_path
        hres_paths = self.hres_paths

        ################################################################################################################
        # 1.	Read Match-up Data and Harmonisation Result Data
        ################################################################################################################

        print"Match-up Dataset:"
        for path in dataset_paths:
            print ">", path

        print"\nSensor Data:"
        print ">", sensor_data_path

        print "\nHarmonisation Output File:"
        print ">", hout_path

        print("\nOpening Files...")
        MatchUpData = self.data_reader(dataset_paths, sensor_data_path, open_uncertainty=False)
        HarmResult = HarmonisationResult(hout_path)

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
        HarmVisOp.plot_compare_calibration(plots_dir, HarmResult.parameter, HarmResult.parameter_covariance_matrix, verbose=verbose)

        print "\nPlots written to:", plots_dir


if __name__ == "__main__":
    pass

