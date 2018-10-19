"""
Functions to read harmonisation job configuration files

Created on Tues May 9  2017 15:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
import ConfigParser
import logging
import argparse
import importlib
import os.path
from os import listdir
from os.path import isfile, abspath, split, dirname
from os.path import join as pjoin
import sys

'''___Harmonisation Modules___'''
sys.path.append(dirname(dirname(__file__)))

'''___Constants___'''
software_short_name = "EV"


def parse_cmdline():
    parser = argparse.ArgumentParser(
        description="Run harmonisation of match-up dataset",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("config_file", action="store",
                        help="Path of harmonisation configuration file")

    parser.add_argument("--mode", action="store", choices=["normal", "setup_pc", "setup_trans2ind"], default="normal",
                        help="Mode to run harmonisation solver in")

    parser.add_argument("--no_covariance", action="store_true",
                        help="Option to not compute covariance matrix of harmonised calibration parameters")

    parser.add_argument("--pc_input",
                        help="Path pre-computed harmonisation output file for dataset, to be used as "
                             "preconditioner for this run")

    parser.add_argument("--save_pc", action="store",
                        help="Path to save preconditioner generated in (or used by if opened from elsewhere) "
                             "the harmonisation processing to")

    parser.add_argument("--gn_input",
                        help="Path pre-computed Gauss Newton solver state to start harmonisation processing at")

    parser.add_argument("--save_gn", action="store",
                        help="Path to write Gauss Newton solver state to at end of the harmonisation processing")

    log_options = parser.add_mutually_exclusive_group()
    log_options.add_argument("--verbose", action="store_true",
                             help="Option for verbose output")

    log_options.add_argument("--quiet", action="store_true",
                             help="Option for quiet output")

    parser.add_argument("--log", action="store", type=str,
                        help="Log file to write to. Leave out for stdout.")

    parser.add_argument("--version", action="version", version='v%s' % __version__)

    return parser.parse_args()

def configure_logging(fname=None, verbose=False, quiet=False):
    """
    Configure logger

    :param fname: str
    :param fname: path to directory to write log file to (None for not to write file)
    """

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)

    # File logging
    if fname is not None:
        file_formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(message)s')
        file_handler = logging.FileHandler(fname)

        if quiet:
            file_handler.setLevel(logging.INFO)
        else:
            file_handler.setLevel(logging.DEBUG)

        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    # Stream logging
    stream_formatter = logging.Formatter('%(message)s')
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(stream_formatter)

    if verbose:
        stream_handler.setLevel(logging.DEBUG)
    elif quiet:
        stream_handler.setLevel(logging.DEBUG)
    else:
        stream_handler.setLevel(logging.INFO)

    logger.addHandler(stream_handler)

    return logger


def read_job_cfg(filename):
    """
    Return data from harmonisation job configuration file

    :type filename: str
    :param filename: path of job configuration file

    :return:
        :job_id: *str*

        ID of harmonisation job

        :matchup_dataset: *str*

        Abbreviated harmonisation dataset name, of the form:
        INSTR_TYPE_M_CORREL_MCN_YYYYMMDD_YYYYMMDD

        where,
        INSTR - instrument name
        TYPE - dataset type
        M - number of harmonised parameters in sensor model
        CORREL - correlation structures present in data
        MCN - Monte Carlo trial number (filled with ___ if not a MC dataset)

        :dataset_dir:  *str*

        Path of directory containing matchup data files only (as named in matchup_dataset)

        :sensor_data_path: *str*

        Path of sensor data for of matchup sensors

        :output_dir: *str*

        Path of directory to store output data files in

        :data_reader_path: str

        Path to sensor functions
    """

    # First take string of whole configuration file
    with open(filename) as f:
        job_text = f.read()

    # Open file
    config = ConfigParser.RawConfigParser()
    config.read(filename)

    # Get naming info
    job_id = config.get('DEFAULT', 'job_id')
    matchup_dataset = config.get('DEFAULT', 'matchup_dataset')

    # Get data directories and paths
    dataset_dir = abspath(config.get('DATA', 'dataset_dir'))
    sensor_data_path = abspath(config.get('DATA', 'sensor_data_path'))
    output_dir = abspath(config.get('DATA', 'output_dir'))

    return job_id, matchup_dataset, dataset_dir, sensor_data_path, output_dir, job_text


def get_harm_paths(output_dir):
    """
    Return path of harmonisation output file and list of paths of harmonisation residuals files in directory

    :param output_dir: str
        path of directory containing only harmonisation output and residual data files

    :return:
        :dataset_paths: list: str
            Paths of matchup series files in matchup dataset directory
    """

    harm_output_path = None
    harm_res_paths = []
    for f in listdir(output_dir):
        if (isfile(pjoin(output_dir, f))) and (f[-3:] == ".nc") and ("res" in f):
            harm_res_paths.append(pjoin(output_dir, f))
        elif (isfile(pjoin(output_dir, f))) and (f[-3:] == ".nc"):
            harm_output_path = pjoin(output_dir, f)

    return harm_output_path, harm_res_paths


def import_file(path):
    """
    Return imported python module

    :param path: str
        path of python module to import

    :return:
        :mod: module
            imported module
    """

    directory = split(path)[0]
    fname = os.path.splitext(os.path.basename(path))[0]

    sys.path.insert(0, directory)

    mod = importlib.import_module(fname)

    return mod

if __name__ == "__main__":
    pass

