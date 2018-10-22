#!/opt/anaconda2/bin/python2.7

"""
Command line tool for generating harmonisation result diagnostic plots

For usage try:
python harm_plot.py --help
"""

'''___Python Modules___'''
from os import makedirs

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''
from common import *
from HarmonisationPlottingOp import HarmonisationPlottingOp


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
