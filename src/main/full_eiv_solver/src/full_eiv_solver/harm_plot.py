#!/opt/anaconda2/bin/python2.7

"""
Command line tool for generating harmonisation result diagnostic plots

For usage try:
python full_eiv_solver.py --help
"""

'''___Python Modules___'''
from os import makedirs

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''
from common import *
from HarmonisationPlottingOp import HarmonisationPlottingOp


def parse_cmdline():
    parser = argparse.ArgumentParser(
        description="Run harmonisation of match-up dataset",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("config_file", action="store",
                        help="Path of harmonisation configuration file")

    log_options = parser.add_mutually_exclusive_group()
    log_options.add_argument("--verbose", action="store_true",
                             help="Option for verbose output")

    log_options.add_argument("--quiet", action="store_true",
                             help="Option for quiet output")

    parser.add_argument("--version", action="version", version='v%s' % __version__)

    return parser.parse_args()


def try_makedirs(directory):
    try:
        makedirs(directory)
    except OSError:
        pass
    return 0


def main(p):

    ################################################################################################################
    # Process configuration data
    ################################################################################################################

    job_cfg_fname = p.config_file

    print "Harmonisation Output Plotting \n"

    print "Reading Job Config:", job_cfg_fname, "\n"

    # 1. Read configuration data
    conf = {}   # dictionary to store data
    conf['job_id'], conf['matchup_dataset'], dataset_dir, sensor_data_path,\
        output_dir, conf['job_text'] = read_job_cfg(job_cfg_fname)

    # 2. Get harmonisation result paths
    hout_path, hres_paths = get_harm_paths(output_dir)

    ################################################################################################################
    # Run harmonisation
    ################################################################################################################

    # Initialise object
    H = HarmonisationPlottingOp(dataset_dir=dataset_dir,
                                sensor_data_path=sensor_data_path,
                                output_dir=output_dir,
                                software_cfg=conf,
                                hout_path=hout_path,
                                hres_paths=hres_paths)

    # Run algorithm
    H.run()

    return 0

