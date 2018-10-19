#!/opt/anaconda2/bin/python2.7

"""
Command line tool for NPL Harmonsation Implementation
"""

'''___Python Modules___'''
from os import makedirs
from os.path import basename, dirname
from os.path import join as pjoin
from sys import argv
import argparse
import logging

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''
from common import *
from harmonisation import HarmonisationOp

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Arta Dilo", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


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

parsed_cmdline = parse_cmdline()


def configure_logging(fname, verbose=False, quiet=False):
    """
    Configure logger

    :param path_to_log_directory:  path to directory to write log file in
    :return:
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

logger = configure_logging(parsed_cmdline.log, parsed_cmdline.verbose, parsed_cmdline.quiet)

'''___Constants___'''

# Tolerances
TOLPC = 1e-6    # Preconditioner alogrithm convergence tolerance
TOL = 1e-6      # Gauss-Newton algorithm convergence tolerance
TOLA = 1e-8     # Gauss-Newton algorithm LSMR tolerance tolA
TOLB = 1e8      # Gauss-Newton algorithm LSMR tolerance tolB
TOLU = 1e-8     # Gauss-Newton algorithm Minres rtol

preamble = "\nErrors-in-Variables Sensor Harmonisation"

def main(parsed_cmdline):

    job_cfg_fname = parsed_cmdline.config_file

    ####################################################################################################################
    # Process configuration data
    ####################################################################################################################

    logger.info(preamble)
    logger.info("\nReading Configuration Data: "+parsed_cmdline.config_file)

    # 1. Read configuration data
    conf = {}   # dictionary to store data

    #  a. Read software config file
    software_cfg_fname = pjoin(dirname(__file__), "software.cfg")
    conf['software'], conf['version'], conf['tag'], conf['software_text'] = read_software_cfg(software_cfg_fname)

    # b. Read job config file
    conf['job_id'], conf['matchup_dataset'], dataset_dir,\
        sensor_data_path, output_dir, conf['job_text'] = read_job_cfg(job_cfg_fname)

    # 4. Make output directory if it doesn't exist
    try:
        makedirs(output_dir)
    except OSError:
        pass

    # Logging options
    show = 1
    if parsed_cmdline.quiet:
        show = 0
    elif parsed_cmdline.verbose:
        show = 2

    ################################################################################################################
    # Run harmonisation
    ################################################################################################################

    print "Match-up Dataset Directory:", dataset_dir
    print "Sensor Data File:", sensor_data_path

    # Run algorithm
    H = HarmonisationOp()
    H.run(dataset_dir=dataset_dir,
          sensor_data_path=sensor_data_path,
          output_dir=output_dir,
          pc_input=parsed_cmdline.pc_input,
          save_pc=parsed_cmdline.save_pc,
          gn_input=parsed_cmdline.gn_input,
          save_gn=parsed_cmdline.save_gn,
          software_cfg=conf,
          tolPC=TOLPC, tol=TOL, tolA=TOLA, tolB=TOLB, tolU=TOLU,
          show=show,
          return_covariance=(not parsed_cmdline.no_covariance))

    return 0

if __name__ == "__main__":
    main(parsed_cmdline)
