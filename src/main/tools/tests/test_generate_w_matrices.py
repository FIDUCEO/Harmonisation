"""
Test of MatchUpData Class
"""

'''___Built-In Modules___'''
import sys
from os.path import dirname
import unittest

'''___Third-Party Modules___'''
from numpy import array, arange
from scipy.sparse import csr_matrix

'''___NPL Modules___'''
tools_directory = dirname(dirname(__file__))
sys.path.append(tools_directory)
from generate_w_matrices import generate_rolling_average_w_matrix

'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "18/11/2017"
__credits__ = ["Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class TestGenerateWMatrices(unittest.TestCase):
    def test_generate_rolling_average_w_matrix_uniform_ufull(self):
        """
        Test for function which generates a w matrix for a rolling average operation
        - one observation per scanline
        """

        ################################################################################################################
        # 1. Define input data
        ################################################################################################################

        matchup_times = array([0.0, 1.0, 2.0, 3.0, 4.0])
        scanline_time = 1.0
        kernel_uncertainty = array([[1.0, 1.1, 1.2, 1.3, 1.4],
                                    [1.1, 1.2, 1.3, 1.4, 1.5],
                                    [1.2, 1.3, 1.4, 1.5, 1.6],
                                    [1.3, 1.4, 1.5, 1.6, 1.7],
                                    [1.4, 1.5, 1.6, 1.7, 1.8]])

        ################################################################################################################
        # 2. Define expected output
        ################################################################################################################

        w_matrix_expected = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                                   [0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        

        u_matrix_expected = array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])

        ################################################################################################################
        # 3. W matrix generation script
        ################################################################################################################

        w_matrix_test, u_matrix_test = generate_rolling_average_w_matrix(matchup_times, scanline_time,
                                                                         kernel_uncertainty)
        w_matrix_test = w_matrix_test.toarray()

        ################################################################################################################
        # 3. Compare output
        ################################################################################################################

        # a. w matrix
        for row_test, row_expected in zip(w_matrix_test, w_matrix_expected):
            self.assertEquals(row_test.tolist(), row_expected.tolist())

        # b. u matrix
        self.assertEquals(u_matrix_expected.tolist(), u_matrix_test.tolist())

    def test_generate_rolling_average_w_matrix_nonununiform_ufull(self):
        """
        Test for function which generates a w matrix for a rolling average operation
        - not equal time gap between match-ups
        """

        ################################################################################################################
        # 1. Define input data
        ################################################################################################################

        matchup_times = array([0.0, 0.0, 2.0, 4.0, 4.0])
        scanline_time = 1.0
        kernel_uncertainty = array([[1.0, 1.1, 1.2, 1.3, 1.4],
                                    [1.0, 1.1, 1.2, 1.3, 1.4],
                                    [1.2, 1.3, 1.4, 1.5, 1.6],
                                    [1.4, 1.5, 1.6, 1.7, 1.8],
                                    [1.4, 1.5, 1.6, 1.7, 1.8]])

        ################################################################################################################
        # 2. Define expected output
        ################################################################################################################

        w_matrix_expected = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                                   [0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2],
                                   [0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])

        u_matrix_expected = array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8])

        ################################################################################################################
        # 3. W matrix generation script
        ################################################################################################################

        w_matrix_test, u_matrix_test = generate_rolling_average_w_matrix(matchup_times, scanline_time,
                                                                         kernel_uncertainty)
        w_matrix_test = w_matrix_test.toarray()

        ################################################################################################################
        # 3. Compare output
        ################################################################################################################

        # a. w matrix
        for row_test, row_expected in zip(w_matrix_test, w_matrix_expected):
            self.assertEquals(row_test.tolist(), row_expected.tolist())

        # b. u matrix
        self.assertEquals(u_matrix_expected.tolist(), u_matrix_test.tolist())

    def test_generate_rolling_average_w_matrix_break_ufull(self):
        """
        Test for function which generates a w matrix for a rolling average operation
        - break where time between match-ups >> kernel time
        """

        ################################################################################################################
        # 1. Define input data
        ################################################################################################################

        matchup_times = array([0.0, 3.0, 20.0, 20.0, 21.0])
        scanline_time = 1.0
        kernel_uncertainty = array([[1.0, 1.1, 1.2, 1.3, 1.4],
                                    [1.3, 1.4, 1.5, 1.6, 1.7],
                                    [2.0, 2.1, 2.2, 2.3, 2.4],
                                    [2.0, 2.1, 2.2, 2.3, 2.4],
                                    [2.1, 2.2, 2.3, 2.4, 2.5]])

        ################################################################################################################
        # 2. Define expected output
        ################################################################################################################

        w_matrix_expected = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])

        u_matrix_expected = array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5])

        ################################################################################################################
        # 3. W matrix generation script
        ################################################################################################################

        w_matrix_test, u_matrix_test = generate_rolling_average_w_matrix(matchup_times, scanline_time,
                                                                         kernel_uncertainty)
        w_matrix_test = w_matrix_test.toarray()

        ################################################################################################################
        # 3. Compare output
        ################################################################################################################

        # a. w matrix
        for row_test, row_expected in zip(w_matrix_test, w_matrix_expected):
            self.assertEquals(row_test.tolist(), row_expected.tolist())

        # b. u matrix
        self.assertEquals(u_matrix_expected.tolist(), u_matrix_test.tolist())

    def test_generate_rolling_average_w_matrix_unordered_ufull(self):
        """
        Test for function which generates a w matrix for a rolling average operation
        - match-ups not in time order
        """

        ################################################################################################################
        # 1. Define input data
        ################################################################################################################

        matchup_times = array([102.0, 101.0, 100.0, 0.0, 104.0])
        scanline_time = 1.0
        kernel_uncertainty = array([[1.7, 1.8, 1.9, 2.0, 2.1],
                                    [1.6, 1.7, 1.8, 1.9, 2.0],
                                    [1.5, 1.6, 1.7, 1.8, 1.9],
                                    [1.0, 1.1, 1.2, 1.3, 1.4],
                                    [1.9, 2.0, 2.1, 2.2, 2.3]])

        ################################################################################################################
        # 2. Define expected output
        ################################################################################################################

        w_matrix_expected = array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                                   [0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])


        u_matrix_expected = array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3])

        ################################################################################################################
        # 3. W matrix generation script
        ################################################################################################################

        w_matrix_test, u_matrix_test = generate_rolling_average_w_matrix(matchup_times, scanline_time,
                                                                         kernel_uncertainty)
        w_matrix_test = w_matrix_test.toarray()

        ################################################################################################################
        # 3. Compare output
        ################################################################################################################

        # a. w matrix
        for row_test, row_expected in zip(w_matrix_test, w_matrix_expected):
            self.assertEquals(row_test.tolist(), row_expected.tolist())

        # b. u matrix
        self.assertEquals(u_matrix_expected.tolist(), u_matrix_test.tolist())

    def test_generate_rolling_average_w_matrix_uniform_uint(self):
        """
        Test for function which generates a w matrix for a rolling average operation
        - one observation per scanline
        """

        ################################################################################################################
        # 1. Define input data
        ################################################################################################################

        matchup_times = array([0.0, 1.0, 2.0, 3.0, 4.0])
        scanline_time = 1.0
        kernel_uncertainty = 5.0
        weights = array([0.2, 0.2, 0.2, 0.2, 0.2])

        ################################################################################################################
        # 2. Define expected output
        ################################################################################################################

        w_matrix_expected = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                                   [0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])

        u_matrix_expected = array([5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0])

        ################################################################################################################
        # 3. W matrix generation script
        ################################################################################################################

        w_matrix_test, u_matrix_test = generate_rolling_average_w_matrix(matchup_times, scanline_time,
                                                                         kernel_uncertainty, weights=weights)
        w_matrix_test = w_matrix_test.toarray()

        ################################################################################################################
        # 3. Compare output
        ################################################################################################################

        # a. w matrix
        for row_test, row_expected in zip(w_matrix_test, w_matrix_expected):
            self.assertEquals(row_test.tolist(), row_expected.tolist())

        # b. u matrix
        self.assertEquals(u_matrix_expected.tolist(), u_matrix_test.tolist())

    def test_generate_rolling_average_w_matrix_unordered_uint(self):
        """
        Test for function which generates a w matrix for a rolling average operation
        - match-ups not in time order
        """

        ################################################################################################################
        # 1. Define input data
        ################################################################################################################

        matchup_times = array([102.0, 101.0, 100.0, 0.0, 104.0])
        scanline_time = 1.0
        kernel_uncertainty = 1.0
        weights = array([0.2, 0.2, 0.2, 0.2, 0.2])

        ################################################################################################################
        # 2. Define expected output
        ################################################################################################################

        w_matrix_expected = array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                                   [0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])

        u_matrix_expected = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

        ################################################################################################################
        # 3. W matrix generation script
        ################################################################################################################

        w_matrix_test, u_matrix_test = generate_rolling_average_w_matrix(matchup_times, scanline_time,
                                                                         kernel_uncertainty, weights=weights)
        w_matrix_test = w_matrix_test.toarray()

        ################################################################################################################
        # 3. Compare output
        ################################################################################################################

        # a. w matrix
        for row_test, row_expected in zip(w_matrix_test, w_matrix_expected):
            self.assertEquals(row_test.tolist(), row_expected.tolist())

        # b. u matrix
        self.assertEquals(u_matrix_expected.tolist(), u_matrix_test.tolist())

if __name__ == '__main__':
    unittest.main()
