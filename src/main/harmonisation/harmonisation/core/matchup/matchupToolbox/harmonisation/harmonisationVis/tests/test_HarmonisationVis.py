"""
Tests for HarmonisationVis Class
"""

'''___Built-In Modules___'''
import unittest
from os.path import abspath
from os.path import join as pjoin
from os import makedirs
from shutil import rmtree

'''___Third-Party Modules___'''

'''___harmonisation Modules___'''
from harmonisation import open_matchup
from harmonisation.core.matchup.matchupToolbox.harmonisation.harmonisationIO.HarmonisationResult import HarmonisationResult


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "16/1/2019"
__credits__ = ["Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


# Constant
DATASET_PATHS = [abspath("../../../../../../../data/simulated_matchup/0_1.nc"),
                 abspath("../../../../../../../data/simulated_matchup/1_2.nc"),
                 abspath("../../../../../../../data/simulated_matchup/2_3.nc")]
OUTPUT_DIRECTORY = abspath("../../../../../../../data/simulated_matchup_result")
PLOTS_DIRECTORY = pjoin(OUTPUT_DIRECTORY, "plots")


def setup_open():
    from harmonisation.core.matchup.matchupToolbox.harmonisation.harmonisationVis import HarmonisationVis

    matchup = open_matchup(DATASET_PATHS)
    harmonisation_result = HarmonisationResult(OUTPUT_DIRECTORY)

    harm_vis = HarmonisationVis(matchup, harmonisation_result)

    try:
        makedirs(PLOTS_DIRECTORY)
    except OSError:
        pass

    return harm_vis


def teardown_rm():
    rmtree(PLOTS_DIRECTORY)
    return 0


class TestHarmonisationVis(unittest.TestCase):

    def test_plot_L1_nom_v_L2_nom_scatter(self):
        harm_vis = setup_open()
        harm_vis.plot_L1_nom_v_L2_nom_scatter(PLOTS_DIRECTORY)

        teardown_rm()

    def test_plot_L1_harm_v_L2_harm_scatter(self):
        harm_vis = setup_open()
        harm_vis.plot_L1_harm_v_L2_harm_scatter(PLOTS_DIRECTORY)

        # teardown_rm()



if __name__ == '__main__':
    unittest.main()
