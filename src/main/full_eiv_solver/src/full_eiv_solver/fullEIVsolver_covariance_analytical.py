"""
Compute covariance matrix for completed full EIV harmonisation run
"""

'''___Python Modules___'''
from os import makedirs
from os.path import basename, dirname
from os.path import join as pjoin
import sys
from sys import argv
from glob import glob
from shutil import rmtree

'''___Third Party Modules___'''
from numpy import append, array, savetxt, zeros, loadtxt
from netCDF4 import Dataset

'''___Harmonisation Modules___'''
from config_functions import *
from harmonisation import HarmonisationOp

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def main(job_cfg_fname, i=None, j=None):

    ################################################################################################################
    # Process configuration data
    ################################################################################################################

    element = str(i)+", "+str(j)
    if (i is None) and (j is None):
        element = "Combination"

    print "Errors-in-Variables Sensor Harmonisation - Parameter Covariance Element "+element

    print "\nReading Configuration Data:"
    print "> " + job_cfg_fname + "\n"

    # 1. Read configuration data
    conf = {}  # dictionary to store data

    #  a. Read software config file
    software_cfg_fname = "software.cfg"
    conf['software'], conf['version'], conf['tag'], conf['software_text'] = read_software_cfg(software_cfg_fname)

    # b. Read job config file
    conf['job_id'], conf['matchup_dataset'], dataset_dir, \
    sensor_data_path, output_dir, data_reader_path, conf['job_text'] = read_job_cfg(job_cfg_fname)

    # 2. Get matchup data paths from directory
    dataset_paths = get_dataset_paths(dataset_dir)

    # 3. Import required specified functions
    if basename(data_reader_path) == "DEFAULT":
        harm_data_reader = None
    else:
        harm_data_reader = import_file(data_reader_path).MatchUp

    # 4. Get harmonisation output files paths if Monte Carlo run
    hout_path, hres_paths = get_harm_paths(output_dir)

    ####################################################################################################################
    # Compute Parameter Covariance Matrix
    ####################################################################################################################

    # Initialise object
    HOp = HarmonisationOp(dataset_paths=dataset_paths,
                          sensor_data_path=sensor_data_path,
                          output_dir=output_dir,
                          software_cfg=conf,
                          data_reader=harm_data_reader)

    # Run algorithm
    if (i is not None) and (j is not None):
        parameter_covariance_ij = HOp.calculate_parameter_covariance_ij(i, j)
        save_path = pjoin(HOp.temp_directory, "parameter_covariance_matrix[" + str(i) + "," + str(j) + "].dat")

        print "Writing result to:", save_path
        savetxt(save_path, array([parameter_covariance_ij]))

    else:
        parameter_covaraince_matrix = HOp.combine_parameter_covariance_ij()

        dataset = Dataset(hout_path, "a")
        dataset.variables["parameter_covariance_matrix"][:] = parameter_covaraince_matrix
        dataset.close()

    return 0

if __name__ == "__main__":

    # Combine elements
    if len(argv) == 2:
        main(os.path.abspath(argv[1]))

    # Calculate ijth element
    elif len(argv) == 4:
        main(os.path.abspath(argv[1]), int(argv[2]), int(argv[3]))
