"""
Functions to read harmonisation job configuration files

Created on Tues May 9  2017 15:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
import ConfigParser
import importlib
import os.path
from os import listdir
from os.path import isfile, abspath, split, dirname
from os.path import join as pjoin
import sys

'''___Harmonisation Modules___'''
sys.path.append(dirname(dirname(__file__)))


def read_software_cfg(filename):
    """
    Return data from harmonisation software configuration file

    :type filename: str
    :param filename: path of software configuration file

    :return:
        :software: *str*

        Software implementation abbreviation, i.e.:
        FO - FastOpt
        EV - NPL Errors in Variables
        OM - NPL ODR+MC Approach

        :version: *str*

        Version number of software (format V.V)

        :tag: *str*

        VCS tag for software implementation (format TTTTTTT)
    """

    # First take string of whole configuration file
    with open(filename) as f:
        software_text = f.read()

    # Open file
    config = ConfigParser.RawConfigParser()
    config.read(filename)

    # Get naming info
    software = config.get('DEFAULT', 'software')
    software_version = config.get('DEFAULT', 'software_version')
    software_tag = config.get('DEFAULT', 'software_tag')

    return software, software_version, software_tag, software_text


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

    def main():
        return 0

    main()
