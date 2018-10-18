"""
Tests for Transform2NormInd
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
transform2normind_directory = dirname(dirname(__file__))
sys.path.append(transform2normind_directory)
from Transform2NormInd import Transform2NormInd

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


def return_MatchUpTest___w():
    """
    Return a MatchUp dataset object for testing

    :return:
        :MatchUpTest: *eopy.matchup.matchupIO.MatchUp*

        Test match-up dataset
    """

    ####################################################################################################################
    # 1. Initialise test data
    ####################################################################################################################

    w1 = array([[0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                [0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                [0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25]])
    w2 = array([[0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                [0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00],
                [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25]])

    u1 = array([1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0])
    u2 = array([2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0])

    values = array([5.0, 3.0, 3.0, 2.5, 6.0, 3.0, 2.0, 4.0, 3.0, 4.0])
    unc = [Uncertainty(3, (0, 0)),
           Uncertainty(3, (1, 1))]
    ks = array([1.2, 1.7, 1.3, 1.4, 1.3])
    unck = [Uncertainty(1, array([0.25, 0.25, 0.25, 0.25, 0.25]))]
    idx = {"Nm": [5],
           "cNm": [0, 5, 10],
           "Im": [[1, 2]],
           "sensors": [1, 2],
           "sensor_ms": [1],
           "n_sensor": [1, 2],
           "n_mu": [1, 1],
           "n_cov": [1, 1],
           "N_var": [5, 5],
           "idx": [0, 5, 10],
           "Ia": [1, 1, 1, 2, 2, 2]}
    a = array([1., 1.3, 0.002, 0.5, 1.1, 0.0005])
    w_matrices = [csr_matrix(w1), csr_matrix(w2)]
    u_matrices = [u1, u2]

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


class TestTransform2NormInd(unittest.TestCase):
    def test_run_multi_r__(self):
        """
        Test for Transform2NormInd.run() for test ``eopy.matchup.matchupIO.MatchUp`` object with data for multiple
        matchup series
        - random uncertainty type only
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compare transformed dataset to expected value

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest_r__()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Original dataset values (should be unchanged)
        MatchUpOriginal_expected = return_MatchUpTest_r__()

        # Transformed dataset
        values_expected = array([294.0625, 480.3733333, 300.6, 227.3846154, 210.1533333,
                                 22.74193548, 22.0625, 21.96875, 22.80645161, 23.5,
                                 21.66666667, 21.05882353, 23, 22.40625,
                                 38.33333333, 36.63636364, 36.5, 38.42857143,
                                 30.1, 32.14893617, 29.37254902, 28.88461538, 28.56603774,
                                 33.45238095, 32.81395349, 31.77272727, 32.60465116,
                                 40.125, 43.54054054, 38.59090909, 34.08510638,
                                 13.72727273, 12, 14.1, 11.79069767, 17.53846154,
                                 12.69565217, 31.16666667, 12.26086957, 11.52272727,
                                 8.8125, 12, 7.4, 10.13207547])
        unc_expected = [Uncertainty(1, array([1.6, 1.5, 1.5, 1.3, 1.5])),
                        Uncertainty(1, array([3.1, 3.2, 3.2, 3.1, 3.0])),
                        Uncertainty(1, array([3.3, 3.4, 3.1, 3.2])),
                        Uncertainty(1, array([2.1, 2.2, 2.2, 2.1])),
                        Uncertainty(1, array([5.0, 4.7, 5.1, 5.2, 5.3])),
                        Uncertainty(1, array([4.2, 4.3, 4.4, 4.3])),
                        Uncertainty(1, array([4.0, 3.7, 4.4, 4.7])),
                        Uncertainty(1, array([2.2, 1.7, 2.0, 4.3, 2.6])),
                        Uncertainty(1, array([2.3, 1.2, 2.3, 4.4])),
                        Uncertainty(1, array([3.2, 2.7, 3.0, 5.3]))]
        ks_expected = array([4.8, 6.8, 5.2, 5.6, 5.2, 12.10287443, 13.99394856, 12.48108926, 12.85930408])
        unck_expected = [Uncertainty(1, array([0.25, 0.25, 0.25, 0.25, 0.25])),
                         Uncertainty(1, array([0.2644, 0.2644, 0.2644, 0.2644]))]
        idx_expected = {"Nm": [5, 4],
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
        a_expected = array([1., 1.3, 0.002, 0.5, 1.1, 0.0005])
        w_matrices_expected = []
        u_matrices_expected = []

        ################################################################################################################
        # 3. Run Transform2NormInd.run()
        ################################################################################################################

        Transform2NormIndOp = Transform2NormInd()
        MatchUpTransform = Transform2NormIndOp.run(MatchUpTest)

        values_test = MatchUpTransform.values
        unc_test = MatchUpTransform.unc
        w_matrices_test = MatchUpTransform.w_matrices
        u_matrices_test = MatchUpTransform.u_matrices
        ks_test = MatchUpTransform.ks
        unck_test = MatchUpTransform.unck
        idx_test = MatchUpTransform.idx

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test transformed data object attribute by attribute

        # a. values
        for i, (value_expected, value_test) in enumerate(zip(values_expected, values_test)):
            self.assertAlmostEqual(value_expected, value_test, places=5, msg=str(i))

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.typeID, block_unc_test.typeID)
            self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        for k_expected, k_test in zip(ks_expected, ks_test):
            self.assertAlmostEqual(k_test.tolist(), k_expected.tolist(), places=5)

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

        # Test original data object preserved attribute by attribute

        # a. values
        for i, (value_original_expected, value_original_test) in enumerate(zip(MatchUpOriginal_expected.values, MatchUpTest.values)):
            self.assertAlmostEqual(value_original_expected, value_original_test, places=5)

        # b. unc
        for block_unc_original_expected, block_unc_original_test in zip(MatchUpOriginal_expected.unc, MatchUpTest.unc):
            self.assertEqual(block_unc_original_expected.typeID, block_unc_original_test.typeID)
            self.assertEqual(block_unc_original_expected.uR.tolist(), block_unc_original_test.uR.tolist())

        # c. w_matrices
        self.assertEqual(MatchUpOriginal_expected.w_matrices, MatchUpTest.w_matrices)

        # d. u_matrices
        self.assertEqual(MatchUpOriginal_expected.u_matrices, MatchUpTest.u_matrices)

        # e. ks
        for k_original_expected, k_original_test in zip(MatchUpOriginal_expected.ks, MatchUpTest.ks):
            self.assertAlmostEqual(k_original_test.tolist(), k_original_expected.tolist(), places=5)

        # f. unck
        for block_unck_original_expected, block_unck_original_test in zip(MatchUpOriginal_expected.unck, MatchUpTest.unck):
            self.assertEqual(block_unck_original_expected.typeID, block_unck_original_test.typeID)
            self.assertEqual(block_unck_original_expected.uR.tolist(), block_unck_original_test.uR.tolist())

        # h. idx
        self.assertEqual(set(MatchUpOriginal_expected.idx), set(MatchUpTest.idx))
        for key in MatchUpOriginal_expected.idx.keys():
            idx_i_original_test = MatchUpTest.idx[key]
            idx_i_original_expected = MatchUpOriginal_expected.idx[key]
            if isinstance(idx_i_original_expected, ndarray):
                self.assertEqual(idx_i_original_test.tolist(), idx_i_original_expected.tolist())
            else:
                self.assertEqual(idx_i_original_test, idx_i_original_expected)

    def test_run_single___w(self):
        """
        Test for Transform2NormInd.run() for test ``eopy.matchup.matchupIO.MatchUp`` object with data for multiple
        matchup series
        - random uncertainty type only
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compare transformed dataset to expected value

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest___w()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        w1 = array([[0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                    [0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                    [0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                    [0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00],
                    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25]])
        w2 = array([[0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                    [0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                    [0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                    [0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00],
                    [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25]])

        u1 = array([1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0])
        u2 = array([2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0])

        # Original dataset values (should be unchanged)
        MatchUpOriginal_expected = return_MatchUpTest___w()

        # Transformed dataset
        values_expected = array([5.0, 2.5, 5.0, 2.5, 1.0, 0.5, 3.0, 1.5, 3.0, 3.0, 6.0, 3.0, 6.0,
                                 0.5, 1.0, 1.5, 3.0, 1.5, 3.0, 3.5, 1.0, 0.5, 4.0, 2.0, 4.0, 2.0])
        unc_expected = [Uncertainty(3, (0, 0)),
                        Uncertainty(3, (1, 1))]
        ks_expected = array([4.8, 6.8, 5.2, 5.6, 5.2])
        unck_expected = [Uncertainty(1, array([0.25, 0.25, 0.25, 0.25, 0.25]))]
        idx_expected = {"Nm": [5],
                        "cNm": [0, 5, 10],
                        "Im": [[1, 2]],
                        "sensors": [1, 2],
                        "sensor_ms": [1],
                        "n_sensor": [1, 2],
                        "n_mu": [1, 1],
                        "n_cov": [1, 1],
                        "N_var": [13, 13],
                        "idx": [0, 13, 26],
                        "Ia": [1, 1, 1, 2, 2, 2]}
        a_expected = array([1., 1.3, 0.002, 0.5, 1.1, 0.0005])

        ################################################################################################################
        # 3. Run Transform2NormInd.run()
        ################################################################################################################

        Transform2NormIndOp = Transform2NormInd()
        MatchUpTransform = Transform2NormIndOp.run(MatchUpTest)

        values_test = MatchUpTransform.values
        unc_test = MatchUpTransform.unc
        w_matrices_test = MatchUpTransform.w_matrices
        u_matrices_test = MatchUpTransform.u_matrices
        ks_test = MatchUpTransform.ks
        unck_test = MatchUpTransform.unck
        idx_test = MatchUpTransform.idx

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test transformed data object attribute by attribute

        # a. values
        for i, (value_expected, value_test) in enumerate(zip(values_expected, values_test)):
            self.assertAlmostEqual(value_expected, value_test, places=5, msg=str(i)+": Expected "+str(value_expected)+
                                                                             " != Test "+str(value_test))

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.typeID, block_unc_test.typeID)
            self.assertEqual(block_unc_expected.w_i, block_unc_test.w_i)
            self.assertEqual(block_unc_expected.u_i, block_unc_test.u_i)

        # e. ks
        for k_expected, k_test in zip(ks_expected, ks_test):
            self.assertAlmostEqual(k_test.tolist(), k_expected.tolist(), places=5)

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

        # Test original data object preserved attribute by attribute

        # a. values
        for i, (value_original_expected, value_original_test) in enumerate(
                zip(MatchUpOriginal_expected.values, MatchUpTest.values)):
            self.assertAlmostEqual(value_original_expected, value_original_test, places=5)

        # b. unc
        for block_unc_original_expected, block_unc_original_test in zip(MatchUpOriginal_expected.unc, MatchUpTest.unc):
            self.assertEqual(block_unc_original_expected.typeID, block_unc_original_test.typeID)
            self.assertEqual(block_unc_original_expected.w_i, block_unc_original_test.w_i)
            self.assertEqual(block_unc_original_expected.u_i, block_unc_original_test.u_i)

        # e. ks
        for k_original_expected, k_original_test in zip(MatchUpOriginal_expected.ks, MatchUpTest.ks):
            self.assertAlmostEqual(k_original_test.tolist(), k_original_expected.tolist(), places=5)

        # f. unck
        for block_unck_original_expected, block_unck_original_test in zip(MatchUpOriginal_expected.unck,
                                                                          MatchUpTest.unck):
            self.assertEqual(block_unck_original_expected.typeID, block_unck_original_test.typeID)
            self.assertEqual(block_unck_original_expected.uR.tolist(), block_unck_original_test.uR.tolist())

        # h. idx
        self.assertEqual(set(MatchUpOriginal_expected.idx), set(MatchUpTest.idx))
        for key in MatchUpOriginal_expected.idx.keys():
            idx_i_original_test = MatchUpTest.idx[key]
            idx_i_original_expected = MatchUpOriginal_expected.idx[key]
            if isinstance(idx_i_original_expected, ndarray):
                self.assertEqual(idx_i_original_test.tolist(), idx_i_original_expected.tolist())
            else:
                self.assertEqual(idx_i_original_test, idx_i_original_expected)

    def test_run_multi_rs_(self):

        # todo - write TestTransform2NormInd.test_run_multi_rs_
        pass

    def test_run_multi_rsw(self):

        # todo - write TestTransform2NormInd.test_run_multi_rs_
        pass

    def test_reverse_multi_r__(self):

        # todo - write TestTransform2NormInd.test_reverse_multi_r__
        pass

    def test_reverse_multi_rs_(self):

        # todo - write TestTransform2NormInd.test_reverse_multi_rs_
        pass

    def test_reverse_multi_rsw(self):

        # todo - write TestTransform2NormInd.test_reverse_multi_rsw
        pass


if __name__ == '__main__':
    unittest.main()
