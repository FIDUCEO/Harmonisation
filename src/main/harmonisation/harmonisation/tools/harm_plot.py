#!/opt/anaconda2/bin/python2.7

"""
Command line tool for generating harmonisation result diagnostic plots

For usage try:
python harm_plot.py --help
"""

'''___Built-In Modules___'''
from os import makedirs

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from harmonisation.version import __version__, __tag__
from common import *
from HarmonisationPlottingOp import HarmonisationPlottingOp


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "20/10/2018"
__credits__ = ["Jon Mittaz"]
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


parsed_cmdline = parse_cmdline(solver_options=False)


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

if __name__ == "__main__":
    main(parsed_cmdline)
