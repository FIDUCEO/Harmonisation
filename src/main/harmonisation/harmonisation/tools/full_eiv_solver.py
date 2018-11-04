"""
Command line tool for Error-in-Variables Harmonisation Implementation

For usage try:
$ full_eiv_solver --help
"""

'''___Built-In Modules___'''

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from harmonisation.version import __version__, __tag__, _short_name_
from common import *
from HarmonisationOp import HarmonisationOp


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__version__ = __version__
__tag__ = __tag__
__credits__ = ["Jon Mittaz"]
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


parsed_cmdline = parse_cmdline()
logger = configure_logging(parsed_cmdline.log, parsed_cmdline.verbose, parsed_cmdline.quiet)

preamble = "\nErrors-in-Variables Sensor Harmonisation"


def main():

    logger.info(preamble)

    dataset_dir, sensor_data, output_dir, pc_input, save_pc, gn_input, save_gn, show = read_parsed_cmdline(parsed_cmdline)

    logger.info("Match-up Dataset Directory: " + dataset_dir)
    logger.info("Harmonisation Result Directory: " + output_dir)

    conf = {"software": _short_name_, "version": __version__, "tag": __tag__}

    # Initialise harmonisation operator
    harm_op = HarmonisationOp()

    # Run harmonisation operator
    harm_op.run(dataset_dir=dataset_dir,
                sensor_data=sensor_data,
                output_dir=output_dir,
                pc_input=pc_input,
                save_pc=save_pc,
                gn_input=gn_input,
                save_gn=save_gn,
                software_cfg=conf,
                tolPC=TOLPC, tol=TOL, tolA=TOLA, tolB=TOLB, tolU=TOLU,
                show=show,
                return_covariance=(not parsed_cmdline.no_covariance))

    logger.info("Complete")

    return 0

if __name__ == "__main__":
    main()
