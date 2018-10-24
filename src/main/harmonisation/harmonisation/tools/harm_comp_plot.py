#!/opt/anaconda2/bin/python2.7

"""
Command line tool for generating harmonisation result comparison plots

For usage try:
python harm_comp_plot.py --help
"""

'''___Python Modules___'''
from os import makedirs
from os.path import basename
from sys import argv

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from common import *
from HarmonisationCompPlottingOp import HarmonisationComparisonPlottingOp


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "12/12/2017"
__credits__ = ["Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def try_makedirs(directory):
    try:
        makedirs(directory)
    except OSError:
        pass
    return 0


def main(job_cfg_fname):

    ################################################################################################################
    # Process configuration data
    ################################################################################################################

    print "Harmonisation Output Plotting \n"

    print "Reading Job Config:", job_cfg_fname, "\n"

    # 1. Read configuration data
    conf = {}   # dictionary to store data

    # b. Read job config file
    conf['job_id'], conf['matchup_dataset'], dataset_dir, sensor_data_path,\
        output_dir, data_reader_path, conf['job_text'] = read_job_cfg(job_cfg_fname)

    # 2. Get matchup data paths from directory
    dataset_paths = get_dataset_paths(dataset_dir)

    # 3. Import required specified functions
    if basename(data_reader_path) == "DEFAULT":
        harm_data_reader = None
    else:
        harm_data_reader = import_file(data_reader_path).MatchUp

    # 4. Get harmonisation result paths
    hout_path, hres_paths = get_harm_paths(output_dir)

    ################################################################################################################
    # Run harmonisation
    ################################################################################################################

    # Initialise object
    H = HarmonisationComparisonPlottingOp(dataset_paths=dataset_paths,
                                          sensor_data_path=sensor_data_path,
                                          output_dir=output_dir,
                                          data_reader=harm_data_reader,
                                          hout_path=hout_path,
                                          hres_paths=hres_paths)

    # Run algorithm
    H.run()

    return 0

if __name__ == "__main__":
    main(os.path.abspath(argv[1]))
