"""
Main Operator of NPL Harmonsation Implementation

Usage:
python fullEIVsolver.py /path/to/job.cfg <0/1 - return covariance (optional)>
"""

'''___Python Modules___'''
from os import makedirs
from os.path import basename, dirname
from os.path import join as pjoin
from sys import argv

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''
from config_functions import *
from harmonisation import HarmonisationOp

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


def main(job_cfg_fname, return_covariance=True, PC_directory=None, GNOp_directory=None):

    ################################################################################################################
    # Process configuration data
    ################################################################################################################

    print "Errors-in-Variables Sensor Harmonisation"

    print "\nReading Configuration Data:"
    print "> " + job_cfg_fname + "\n"

    # 1. Read configuration data
    conf = {}   # dictionary to store data

    #  a. Read software config file
    software_cfg_fname = pjoin(dirname(__file__), "software.cfg")
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
                        PC_dir=PC_directory,
                        GN_dir=GNOp_directory,
                        software_cfg=conf,
                        data_reader=harm_data_reader)

    # Run algorithm
    H.run(tolPC=TOLPC, tol=TOL, tolA=TOLA, tolB=TOLB, tolU=TOLU, show=True, return_covariance=return_covariance)

    return 0

if __name__ == "__main__":

    # Parameter + Covariance
    if len(argv) == 2:
        main(os.path.abspath(argv[1]))

    # Parameter + Covariance switch
    elif len(argv) == 3:
        main(os.path.abspath(argv[1]), bool(int(argv[2])))

    # Parameter + Covariance + PC Solution Path
    elif len(argv) == 4:
        main(os.path.abspath(argv[1]), bool(int(argv[2])), os.path.abspath(argv[3]))

    # Parameter + Covariance + PC Solution Directory + Initial GNOp Directory
    elif len(argv) == 5:
        main(os.path.abspath(argv[1]), bool(int(argv[2])), os.path.abspath(argv[3]), os.path.abspath(argv[4]))

    # main("/home/seh2/src/DeployedProjects/full_EIV/src/test/test_configs/AVHRR_RSIM_3_test_newdata/AVHRR_RSIM_3_test_newdata.cfg")
