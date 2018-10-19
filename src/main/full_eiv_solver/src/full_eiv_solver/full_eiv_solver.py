#!/opt/anaconda2/bin/python2.7

"""
Command line tool for Error-in-Variables Harmonisation Implementation

For usage try:
python full_eiv_solver.py --help
"""

'''___Python Modules___'''
from os import makedirs

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''
from version import version, tag
from common import *
from HarmonisationOp import HarmonisationOp

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Jon Mittaz"]
__version__ = version
__tag__ = tag
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

parsed_cmdline = parse_cmdline()
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
    conf = {}

    #  a. Read software config file
    conf['software'], conf['version'], conf['tag'] = software_short_name, __version__, __tag__

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
