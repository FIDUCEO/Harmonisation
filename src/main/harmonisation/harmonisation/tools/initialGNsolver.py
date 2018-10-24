"""
Main Script to Determine Preconditioner of NPL Harmonsation Implementation

Usage:
python fullEIVsolver.py /path/to/job.cfg
"""

'''___Built-In Modules___'''
from os import makedirs
from os.path import basename
from sys import argv

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from common import *
from HarmonisationOp import HarmonisationOp


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Arta Dillo", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


'''___Constants___'''

# Tolerances
TOLPC = 1e-6    # Preconditioner alogrithm convergence tolerance
TOL = 1e-6      # Gauss-Newton algorithm convergence tolerance
TOLA = 1e-8     # Gauss-Newton algorithm LSMR tolerance tolA
TOLB = 1e8      # Gauss-Newton algorithm LSMR tolerance tolB
TOLU = 1e-8     # Gauss-Newton algorithm Minres rtol


def main(job_cfg_fname):

    ################################################################################################################
    # Process configuration data
    ################################################################################################################

    print "Errors-in-Variables Sensor Harmonisation - Initial GN"

    print "\nReading Configuration Data:"
    print "> " + job_cfg_fname + "\n"

    # 1. Read configuration data
    conf = {}   # dictionary to store data

    #  a. Read software config file
    software_cfg_fname = "software.cfg"
    conf['software'], conf['version'], conf['tag'], conf['software_text'] = read_software_cfg(software_cfg_fname)

    # b. Read job config file
    conf['job_id'], conf['matchup_dataset'], dataset_dir,\
        sensor_data_path, output_dir, data_reader_path, conf['job_text'] = read_job_cfg(job_cfg_fname)

    # 2. Get matchup data paths from directory
    dataset_paths = get_dataset_paths(dataset_dir)

    # 3. Import required specified functions
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
    H = HarmonisationOp(dataset_paths=dataset_paths,
                        sensor_data_path=sensor_data_path,
                        output_dir=output_dir,
                        software_cfg=conf,
                        data_reader=harm_data_reader)

    # Run algorithm
    H.calculate_GNOp_initial()

    return 0

if __name__ == "__main__":
    main(os.path.abspath(argv[1]))

