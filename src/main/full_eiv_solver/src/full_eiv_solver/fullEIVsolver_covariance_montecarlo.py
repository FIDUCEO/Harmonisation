"""
Compute covariance matrix for completed full EIV harmonisation run by Monte Carlo
"""

'''___Python Modules___'''
from os import makedirs
from os.path import basename
from sys import argv

'''___Third Party Modules___'''
from numpy import append, zeros, loadtxt
from numpy.random import normal
from netCDF4 import Dataset

'''___Harmonisation Modules___'''
from common import *
from HarmonisationOp import HarmonisationOp

main_directory = dirname(dirname(__file__))
sys.path.append(main_directory)
from nplcore import GNAlgo
from nplcore import HarmonisationResult

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def main(job_cfg_fname, trial_number=None):
    ################################################################################################################
    # Process configuration data
    ################################################################################################################

    # 1. Read configuration data
    conf = {}  # dictionary to store data

    #  a. Read software config file
    software_cfg_fname = "software.cfg"
    conf['software'], conf['version'], conf['tag'], conf['software_text'] = read_software_cfg(software_cfg_fname)

    # b. Read job config file
    conf['job_id'], conf['matchup_dataset'], dataset_dir, parameter_path, output_dir, \
    sensor_functions_path, data_reader_path, conf['job_text'] = read_job_cfg(job_cfg_fname)

    # 2. Get matchup data paths from directory
    dataset_paths = get_dataset_paths(dataset_dir)

    # 3. Import required specified functions
    sensor_functions = import_file(sensor_functions_path)
    sensor_model = sensor_functions.sensor_model
    try:
        adjustment_model = sensor_functions.adjustment_model
    except AttributeError:
        adjustment_model = None

    if basename(data_reader_path) == "DEFAULT":
        harm_data_reader = None
    else:
        harm_data_reader = import_file(data_reader_path).MatchUp

    # 4. Make output directory if it doesn't exist
    try:
        makedirs(output_dir)
    except OSError:
        pass

    ################################################################################################################
    # Run harmonisation
    ################################################################################################################

    # Initialise object
    HOp = HarmonisationOp(dataset_paths=dataset_paths,
                          output_dir=output_dir,
                          software_cfg=conf,
                          data_reader=harm_data_reader)

    # Run algorithm
    if trial_number is not None:
        HOp.run_mc_trial(trial_number)

    else:
        parameter_covaraince_matrix = HOp.combine_parameter_covariance_ij()

        dataset = Dataset(hout_path, "a")
        dataset.variables["parameter_covariance_matrix"][:] = parameter_covaraince_matrix
        dataset.close()

        # rmtree(HOp.temp_directory)

    return 0

if __name__ == "__main__":

    # Combine trial output
    if len(argv) == 2:
        main(os.path.abspath(argv[1]))

    # Run one trial
    elif len(argv) == 3:
        main(os.path.abspath(argv[1]), int(argv[2]))
