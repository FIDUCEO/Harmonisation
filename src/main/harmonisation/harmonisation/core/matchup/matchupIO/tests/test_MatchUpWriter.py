"""
Tests for MatchUpWriter Class
"""

'''___Built-In Modules___'''
import unittest
from os.path import abspath

'''___Third-Party Modules___'''

'''___harmonisation Modules___'''
from harmonisation import open_matchup


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "16/1/2019"
__credits__ = ["Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


# Constant
DATASET_PATHS = [abspath("../../../../../data/simulated_matchup/0_1.nc"),
                 abspath("../../../../../data/simulated_matchup/1_2.nc"),
                 abspath("../../../../../data/simulated_matchup/2_3.nc")]
OUTPUT_DIRECTORY = abspath("../../../../../data/test_output_matchup")


def setup_open_matchup():
    from harmonisation.core.matchup.matchupIO.MatchUpWriter import MatchUpWriter
    writer = MatchUpWriter()
    matchup = open_matchup(DATASET_PATHS)
    return writer, matchup


class TestMatchUpWriter(unittest.TestCase):

    def test_writeMatchUp(self):

        writer, matchup = setup_open_matchup()
        writer.writeMatchUp(matchup, OUTPUT_DIRECTORY)


if __name__ == '__main__':
    unittest.main()
