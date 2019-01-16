"""
Command line tool for generating harmonisation result diagnostic plots

For usage try:
$ harm_plot --help
"""

'''___Built-In Modules___'''

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from harmonisation.version import __version__
from common import *
from HarmonisationPlottingOp import HarmonisationPlottingOp


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "20/10/2018"
__version__ = __version__
__credits__ = ["Jon Mittaz"]
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


parsed_cmdline = parse_cmdline(solver_options=False)
logger = configure_logging(parsed_cmdline.log, parsed_cmdline.verbose, parsed_cmdline.quiet)

preamble = "\nErrors-in-Variables Sensor Harmonisation Result Plotting"


def main():

    logger.info(preamble)

    dataset_dir, sensor_data, output_dir, show = read_parsed_cmdline(parsed_cmdline, solver_options=False)

    logger.info("Match-up Dataset Directory: "+dataset_dir)
    logger.info("Harmonisation Result Directory: "+output_dir)

    conf = {}

    # Initialise plotting operator
    harm_plot_op = HarmonisationPlottingOp(dataset_dir=dataset_dir,
                                           sensor_data=sensor_data,
                                           output_dir=output_dir,
                                           logger=logger)

    # Run plotting operator
    harm_plot_op.run()

    logger.info("Complete")

    return 0

if __name__ == "__main__":
    main()
