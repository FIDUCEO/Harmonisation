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
from os.path import isfile, abspath, split
from os.path import join as pjoin
import sys


def read_software_cfg(filename):
    """
    Return data from harmonisation software configuration file

    :param filename: str
        path of software configuration file

    :return:
        :software: str
            Software implementation abbreviation, i.e.:
            FO - FastOpt
            EV - NPL Errors in Variables
            OM - NPL ODR+MC Approach

        :version: str
            Version number of software (format V.V)

        :tag: str
            VCS tag for software implementation (format TTTTTTT)
    """

    # Open file
    config = ConfigParser.RawConfigParser()
    config.read(filename)

    # Get naming info
    software = config.get('DEFAULT', 'software')
    software_version = config.get('DEFAULT', 'software_version')
    software_tag = config.get('DEFAULT', 'software_tag')

    return software, software_version, software_tag


def read_job_cfg(filename):
    """
    Return data from harmonisation job configuration file

    :param filename: str
        path of job configuration file

    :return:
        :job_id: str
            ID of harmonisation job

        :matchup_dataset: str
            Abbreviated harmonisation dataset name, of the form:
            INSTR_TYPE_M_CORREL_MCN_YYYYMMDD_YYYYMMDD

            where,
            INSTR - instrument name
            TYPE - dataset type
            M - number of harmonised parameters in sensor model
            CORREL - correlation structures present in data
            MCN - Monte Carlo trial number (filled with ___ if not a MC dataset)

        :dataset_dir:  str
            Path of directory containing matchup data files only (as named in matchup_dataset)

        :parameter_path: str
            Path of parameters to be used as starting point for solver and comparison in diagnostics

        :output_dir: str
            Path of directory to store output data files in

        :sensor_functions_path: str
            Path of data reader (Default harm_data_reader.py, however, some test datasets require specific reader,
            for example to open averaging uncertainties

        :data_reader_path: str
            Path to sensor functions

    """

    # Open file
    config = ConfigParser.RawConfigParser()
    config.read(filename)

    # Get naming info
    job_id = config.get('DEFAULT', 'job_id')
    matchup_dataset = config.get('DEFAULT', 'matchup_dataset')

    # Get data directories and paths
    dataset_dir = abspath(config.get('DATA', 'dataset_dir'))
    parameter_path = abspath(config.get('DATA', 'parameter_path'))
    output_dir = abspath(config.get('DATA', 'output_dir'))

    # Get functions paths
    data_reader_path = abspath(config.get('FUNCTIONS', 'data_reader_path'))
    sensor_functions_path = abspath(config.get('FUNCTIONS', 'sensor_functions_path'))

    return job_id, matchup_dataset, dataset_dir, parameter_path, output_dir, sensor_functions_path, data_reader_path


def get_dataset_paths(dataset_dir):
    """
    Return list of matchup data files in

    :param dataset_dir: str
        path of directory containing only matchup data files

    :return:
        :dataset_paths: list: str
            Paths of matchup series files in matchup dataset directory
    """

    dataset_paths = []
    for f in listdir(dataset_dir):
        if (isfile(pjoin(dataset_dir, f))) and (f[-3:] == ".nc"):
            dataset_paths.append(pjoin(dataset_dir, f))

    return dataset_paths


def get_harm_paths(output_dir):
    """
    Return path of harmonisation output file and list of paths of harmonisation residuals files in directory

    :param output_dir: str
        path of directory containing only harmonisation output and residual data files

    :return:
        :dataset_paths: list: str
            Paths of matchup series files in matchup dataset directory
    """

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
