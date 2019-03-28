"""
Command line tool for generating harmonisation result comparison plots

For usage try:
$ harm_comp_plot --help
"""

'''___Built-In Modules___'''

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from harmonisation.version import __version__, __tag__
from common import *
from HarmonisationCompPlottingOp import HarmonisationComparisonPlottingOp


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "12/12/2017"
__version__ = __version__
__tag__ = __tag__
__credits__ = ["Jon Mittaz"]
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


parsed_cmdline = parse_cmdline(description="Generate plots to compare harmonisation outputs", solver_options=False)
logger = configure_logging(parsed_cmdline.log, parsed_cmdline.verbose, parsed_cmdline.quiet)

preamble = "\nErrors-in-Variables Sensor Harmonisation Result Comparison Plotting"


def main():

    logger.info(preamble)

    dataset_dir, sensor_data, output_dir, show = read_parsed_cmdline(parsed_cmdline, solver_options=False)

    logger.info("Match-up Dataset Directory: "+dataset_dir)
    logger.info("Harmonisation Result Directory: "+output_dir)

    conf = {}

    # Initialise object
    harm_comp_op = HarmonisationComparisonPlottingOp(dataset_paths=dataset_dir,
                                                     sensor_data_path=sensor_data,
                                                     output_dir=output_dir)

    # Run algorithm
    harm_comp_op.run()

    return 0

if __name__ == "__main__":
    main()
