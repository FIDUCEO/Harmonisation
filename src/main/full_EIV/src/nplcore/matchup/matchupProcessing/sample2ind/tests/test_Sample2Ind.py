"""
Test of Sample2Ind
"""

'''___Built-In Modules___'''
import unittest
from os.path import dirname
from os.path import join as pjoin
import sys

'''___Third-Party Modules___'''
from numpy import array, ndarray
from scipy.sparse import csr_matrix

'''___NPL Modules___'''
sample2ind_directory = dirname(dirname(__file__))
sys.path.append(sample2ind_directory)
from Sample2Ind import Sample2Ind

matchupIO_directory = pjoin(dirname(dirname(dirname(dirname(__file__)))), "matchupIO")
sys.path.append(matchupIO_directory)
from Uncertainty import Uncertainty
from MatchUp import MatchUp

'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "23/11/2017"
__credits__ = ["Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def return_MatchUpTest_r__():
    """
    Return a MatchUp dataset object for testing

    :return:
        :MatchUpTest: *eopy.matchup.matchupIO.MatchUp*

        Test match-up dataset
    """

    ####################################################################################################################
    # 1. Initialise test data
    ####################################################################################################################

    values = array([470.5, 720.56, 450.9, 295.6, 315.23,
                    70.5, 70.6, 70.3, 70.7, 70.5,
                    71.5, 71.6, 71.3, 71.7,
                    80.5, 80.6, 80.3, 80.7,
                    150.5, 151.1, 149.8, 150.2, 151.4,
                    140.5, 141.1, 139.8, 140.2,
                    160.5, 161.1, 169.8, 160.2,
                    30.2, 20.4, 28.2, 50.7, 45.6,
                    29.2, 37.4, 28.2, 50.7,
                    28.2, 32.4, 22.2, 53.7, ])
    unc = [Uncertainty(1, array([1.6, 1.5, 1.5, 1.3, 1.5])),
           Uncertainty(1, array([3.1, 3.2, 3.2, 3.1, 3.0])),
           Uncertainty(1, array([3.3, 3.4, 3.1, 3.2])),
           Uncertainty(1, array([2.1, 2.2, 2.2, 2.1])),
           Uncertainty(1, array([5.0, 4.7, 5.1, 5.2, 5.3])),
           Uncertainty(1, array([4.2, 4.3, 4.4, 4.3])),
           Uncertainty(1, array([4.0, 3.7, 4.4, 4.7])),
           Uncertainty(1, array([2.2, 1.7, 2.0, 4.3, 2.6])),
           Uncertainty(1, array([2.3, 1.2, 2.3, 4.4])),
           Uncertainty(1, array([3.2, 2.7, 3.0, 5.3]))]
    ks = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
    unck = [Uncertainty(1, array([0.25, 0.25, 0.25, 0.25, 0.25])),
            Uncertainty(1, array([0.2644, 0.2644, 0.2644, 0.2644]))]
    idx = {"Nm": [5, 4],
           "cNm": [0, 5, 9],
           "Im": [[0, 1], [1, 2]],
           "sensors": [-1, 1, 2],
           "sensor_ms": [1, 3, 3],
           "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
           "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
           "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
           "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
           "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44],
           "Ia": [1, 1, 1, 2, 2, 2]}
    a = array([1., 1.3, 0.002, 0.5, 1.1, 0.0005])

    ####################################################################################################################
    # 3. Initialise MatchUp object
    ####################################################################################################################

    MatchUpTest = MatchUp()
    MatchUpTest.values = values
    MatchUpTest.unc = unc
    MatchUpTest.ks = ks
    MatchUpTest.unck = unck
    MatchUpTest.idx = idx
    MatchUpTest.a = a

    return MatchUpTest


def return_MatchUpTest_rsw():
    """
    Return a MatchUp dataset object for testing

    :return:
        :MatchUpTest: *eopy.matchup.matchupIO.MatchUp*

        Test match-up dataset
    """

    ####################################################################################################################
    # 1. Initialise test data
    ####################################################################################################################

    w2_matchup1 = array([[0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5]])
    w1_matchup2 = array([[0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5]])
    w2_matchup2 = array([[0.5, 0.5, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.5, 0.5, 0.0],
                         [0.0, 0.0, 0.5, 0.5, 0.0],
                         [0.0, 0.0, 0.0, 0.5, 0.5]])
    u2_matchup1 = array([1.0, 0.9, 0.78, 1.0, 0.9, 0.68, 1.0])
    u1_matchup2 = array([1.0, 0.9, 0.9, 0.5, 0.9, 0.58, 1.1])
    u2_matchup2 = array([1.0, 0.9, 0.9, 0.5, 0.8])

    values = array([470.5, 720.56, 450.9, 295.6, 315.23,
                    70.5, 70.6, 70.3, 70.7, 70.5,
                    71.5, 71.6, 71.3, 71.7,
                    80.5, 80.6, 80.3, 80.7,
                    150.5, 151.1, 149.8, 150.2, 151.4,
                    140.5, 141.1, 139.8, 140.2,
                    160.5, 161.1, 169.8, 160.2,
                    20.2, 21.0, 19.7, 19.7, 20.7,
                    13.1, 13.2, 11.8, 15.2,
                    12.6, 13.7, 13.7, 11.3])
    unc = [Uncertainty(1, array([1.6, 1.5, 1.5, 1.3, 1.5])),
           Uncertainty(1, array([3.1, 3.2, 3.2, 3.1, 3.0])),
           Uncertainty(1, array([3.3, 3.4, 3.1, 3.2])),
           Uncertainty(1, array([2.1, 2.2, 2.2, 2.1])),
           Uncertainty(1, array([5.0, 4.7, 5.1, 5.2, 5.3])),
           Uncertainty(1, array([4.2, 4.3, 4.4, 4.3])),
           Uncertainty(1, array([4.0, 3.7, 4.4, 4.7])),
           Uncertainty(3, (0, 0)),
           Uncertainty(3, (1, 1)),
           Uncertainty(3, (2, 2))]
    ks = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
    unck = [Uncertainty(1, array([0.25, 0.25, 0.25, 0.25, 0.25])),
            Uncertainty(1, array([0.2644, 0.2644, 0.2644, 0.2644]))]
    idx = {"Nm": [5, 4],
           "cNm": [0, 5, 9],
           "Im": [[0, 1], [1, 2]],
           "sensors": [-1, 1, 2],
           "sensor_ms": [1, 3, 3],
           "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
           "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
           "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
           "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
           "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44],
           "Ia": [1, 1, 1, 2, 2, 2]}
    a = array([1., 1.3, 0.002, 0.5, 1.1, 0.0005])
    w_matrices = [csr_matrix(w2_matchup1), csr_matrix(w1_matchup2), csr_matrix(w2_matchup2)]
    u_matrices = [u2_matchup1, u1_matchup2, u2_matchup2]

    ####################################################################################################################
    # 3. Initialise MatchUp object
    ####################################################################################################################

    MatchUpTest = MatchUp()
    MatchUpTest.values = values
    MatchUpTest.unc = unc
    MatchUpTest.ks = ks
    MatchUpTest.unck = unck
    MatchUpTest.idx = idx
    MatchUpTest.a = a
    MatchUpTest.w_matrices = w_matrices
    MatchUpTest.u_matrices = u_matrices

    return MatchUpTest


class TestSample2Ind(unittest.TestCase):
    def test_run_multi_r___sf1(self):
        """
        Test for Sample2Ind.run() for test ``eopy.matchup.matchupIO.MatchUp`` object with data for multiple matchup
        series
        - scale factor 1
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares sample to expected value
        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest_r__()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        values_expected = MatchUpTest.values
        unc_expected = MatchUpTest.unc
        w_matrices_expected = MatchUpTest.w_matrices
        u_matrices_expected = MatchUpTest.u_matrices
        ks_expected = MatchUpTest.ks
        unck_expected = MatchUpTest.unck
        idx_expected = MatchUpTest.idx

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        Sample2IndOp = Sample2Ind()
        MatchUpSample = Sample2IndOp.run(MatchUpTest, sf=1)

        values_test = MatchUpSample.values
        unc_test = MatchUpSample.unc
        w_matrices_test = MatchUpSample.w_matrices
        u_matrices_test = MatchUpSample.u_matrices
        ks_test = MatchUpSample.ks
        unck_test = MatchUpSample.unck
        idx_test = MatchUpSample.idx

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.typeID, block_unc_test.typeID)
            self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.typeID, block_unck_test.typeID)
            self.assertEqual(block_unck_expected.uR.tolist(), block_unck_test.uR.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

    def test_run_multi_rsw_sf1(self):
        """
        Test for Sample2Ind.run() for test ``eopy.matchup.matchupIO.MatchUp`` object with data for multiple matchup
        series
        - scale factor 1
        """

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest_rsw()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        values_expected = array([470.5, 450.9, 315.23,
                                 70.5, 70.3, 70.5,
                                 71.5, 71.3,
                                 80.5, 80.3,
                                 150.5, 149.8, 151.4,
                                 140.5, 139.8,
                                 160.5, 169.8,
                                 20.2, 19.7, 20.7,
                                 13.1, 11.8,
                                 12.6, 13.7])
        unc_expected = [Uncertainty(1, array([1.6, 1.5, 1.5])),
                        Uncertainty(1, array([3.1, 3.2, 3.0])),
                        Uncertainty(1, array([3.3, 3.1])),
                        Uncertainty(1, array([2.1, 2.2])),
                        Uncertainty(1, array([5.0, 5.1, 5.3])),
                        Uncertainty(1, array([4.2, 4.4])),
                        Uncertainty(1, array([4.0, 4.4])),
                        Uncertainty(1, array([0.951314593, 0.951314593, 0.855102109])),
                        Uncertainty(1, array([0.951314593, 0.728010979])),
                        Uncertainty(1, array([0.951314593, 0.728010979]))]
        ks_expected = array([1.2, 1.3, 1.3, 3.2, 3.3])
        unck_expected = [Uncertainty(1, array([0.25, 0.25, 0.25])),
                         Uncertainty(1, array([0.2644, 0.2644]))]
        idx_expected = {"Nm": [3, 2],
                        "cNm": [0, 3, 5],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [3, 3, 2, 2, 3, 2, 2, 3, 2, 2],
                        "idx": [0, 3, 6, 8, 10, 13, 15, 17, 20, 22, 24],
                        "Ia": [1, 1, 1, 2, 2, 2]}
        w_matrices_expected = []
        u_matrices_expected = []

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        Sample2IndOp = Sample2Ind()
        MatchUpSample = Sample2IndOp.run(MatchUpTest)

        values_test = MatchUpSample.values
        unc_test = MatchUpSample.unc
        w_matrices_test = MatchUpSample.w_matrices
        u_matrices_test = MatchUpSample.u_matrices
        ks_test = MatchUpSample.ks
        unck_test = MatchUpSample.unck
        idx_test = MatchUpSample.idx

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.typeID, block_unc_test.typeID)
            for uR_expected, uR_test in zip(block_unc_expected.uR, block_unc_test.uR):
                self.assertAlmostEqual(uR_expected, uR_test, places=5)

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.typeID, block_unck_test.typeID)
            for uR_expected, uR_test in zip(block_unck_expected.uR, block_unck_test.uR):
                self.assertAlmostEqual(uR_expected, uR_test, places=5)

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)


if __name__ == '__main__':
    unittest.main()
