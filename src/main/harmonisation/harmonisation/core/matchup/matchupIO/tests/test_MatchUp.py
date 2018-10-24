"""
Test of MatchUpData Class
"""

'''___Built-In Modules___'''
import unittest
from os.path import join as pjoin
from os import makedirs, getcwd
from os.path import dirname, exists
import sys
from shutil import rmtree
from datetime import datetime as dt
from random import random

'''___Third-Party Modules___'''
from numpy import array, ndarray, savetxt
from scipy.sparse import csr_matrix

'''___NPL Modules___'''
from test_functions.W_matrix_functions import write_input_file, return_w_matrix_variables, append_W_to_input_file
from harmonisation import MatchUp, Uncertainty


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "18/11/2017"
__credits__ = ["Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

# Constant
temp_data_directory = pjoin(getcwd(), "temp")


class TestMatchUp(unittest.TestCase):
    def test_openMatchUpData_multi_r__(self):
        """
        Test for MatchUpData.openMatchUp() method for case with multiple matchup series datasets
        - with only random uncertainty type
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, each with uncertainty types random
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(dirname(__file__), "temp_data_directory")
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")
        parameter_fname = pjoin(test_data_directory, "parameter.csv")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 1, 1])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 2.3],
                              [3.4, 4.3, 1.2],
                              [3.1, 4.4, 2.3],
                              [3.2, 4.3, 4.4]])

        Us1_matchup2 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 1, 1])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 3.2],
                              [2.2, 3.7, 2.7],
                              [2.2, 4.4, 3.0],
                              [2.1, 4.7, 5.3]])

        Us2_matchup2 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 1, 1])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # --------------------------------------------------------------------------------------------------------------

        # d. Test data to file -----------------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = [Uncertainty("r", array([1.6, 1.5, 1.5, 1.3, 1.5])),
                        Uncertainty("r", array([3.1, 3.2, 3.2, 3.1, 3.0])),
                        Uncertainty("r", array([3.3, 3.4, 3.1, 3.2])),
                        Uncertainty("r", array([2.1, 2.2, 2.2, 2.1])),
                        Uncertainty("r", array([5.0, 4.7, 5.1, 5.2, 5.3])),
                        Uncertainty("r", array([4.2, 4.3, 4.4, 4.3])),
                        Uncertainty("r", array([4.0, 3.7, 4.4, 4.7])),
                        Uncertainty("r", array([2.2, 1.7, 2.0, 4.3, 2.6])),
                        Uncertainty("r", array([2.3, 1.2, 2.3, 4.4])),
                        Uncertainty("r", array([3.2, 2.7, 3.0, 5.3]))]
        w_matrices_expected = []
        u_matrices_expected = []
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = [Uncertainty("r", array([0.25, 0.25, 0.25, 0.25, 0.25])),
                         Uncertainty("r", array([0.2644, 0.2644, 0.2644, 0.2644]))]
        time1_expected = array([dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4),
                                dt(1970, 1, 1, 1, 0, 5),
                                dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4)])
        time2_expected = array([dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000),
                                dt(1970, 1, 1, 1, 0, 5, 100000),
                                dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000)])
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2], parameter_fname)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.form, block_unc_test.form)
            self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.form, block_unck_test.form)
            self.assertEqual(block_unck_expected.uR.tolist(), block_unck_test.uR.tolist())

        # g. time
        self.assertEqual(time1_test.tolist(), time1_expected.tolist())
        self.assertEqual(time2_test.tolist(), time2_expected.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_rs_(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with multiple match-up series datasets
        - with random and random and systematic uncertainty type
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, with covariate 1 and 3 random uncertainty type and
        #    covariate 2 with systematic uncertainty type
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(temp_data_directory, str(int(random()*1000000)))
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 2, 1])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 2.3],
                              [3.4, 4.3, 1.2],
                              [3.1, 4.4, 2.3],
                              [3.2, 4.3, 4.4]])

        Us1_matchup2 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 2, 1])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 3.2],
                              [2.2, 3.7, 2.7],
                              [2.2, 4.4, 3.0],
                              [2.1, 4.7, 5.3]])

        Us2_matchup2 = array([[0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 2, 1])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # --------------------------------------------------------------------------------------------------------------

        # c. Test data to file -----------------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = [Uncertainty("r", array([1.6, 1.5, 1.5, 1.3, 1.5])),
                        Uncertainty("r", array([3.1, 3.2, 3.2, 3.1, 3.0])),
                        Uncertainty("r", array([3.3, 3.4, 3.1, 3.2])),
                        Uncertainty("r", array([2.1, 2.2, 2.2, 2.1])),
                        Uncertainty("rs", (array([5.0, 4.7, 5.1, 5.2, 5.3]), 1.0)),
                        Uncertainty("rs", (array([4.2, 4.3, 4.4, 4.3]), 1.0)),
                        Uncertainty("rs", (array([4.0, 3.7, 4.4, 4.7]), 2.0)),
                        Uncertainty("r", array([2.2, 1.7, 2.0, 4.3, 2.6])),
                        Uncertainty("r", array([2.3, 1.2, 2.3, 4.4])),
                        Uncertainty("r", array([3.2, 2.7, 3.0, 5.3]))]
        w_matrices_expected = []
        u_matrices_expected = []
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = [Uncertainty("r", array([0.25, 0.25, 0.25, 0.25, 0.25])),
                         Uncertainty("r", array([0.2644, 0.2644, 0.2644, 0.2644]))]
        time1_expected = array([dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4),
                                dt(1970, 1, 1, 1, 0, 5),
                                dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4)])
        time2_expected = array([dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000),
                                dt(1970, 1, 1, 1, 0, 5, 100000),
                                dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000)])
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2])

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.form, block_unc_test.form)
            self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())
            if block_unc_expected == "rs":
                self.assertEqual(block_unc_expected.uS.tolist(), block_unc_test.uS.tolist())

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.form, block_unck_test.form)
            self.assertEqual(block_unck_expected.uR.tolist(), block_unck_test.uR.tolist())

        # g. time
        self.assertEqual(time1_test.tolist(), time1_expected.tolist())
        self.assertEqual(time2_test.tolist(), time2_expected.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_rsw(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with mulitple match-up series datasets
        - with random, random and systematic and w matrix uncertainty typea
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, with covariate 1 w matrix uncertainty type, covariate
        #    2 systematic uncertainty type and covariate 3 w matrix uncertainty type
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(temp_data_directory, str(int(random()*1000000)))
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 2, 3])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # iv. w and u matrices
        w_matrix_use1_matchup1 = array([0])
        w_matrix_use2_matchup1 = array([0, 0, 1])
        u_matrix_use1_matchup1 = array([0])
        u_matrix_use2_matchup1 = array([0, 0, 1])

        w2_matchup1 = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        w2_matchup1 = csr_matrix(w2_matchup1)

        u_matrix2_matchup1 = array([0.2, 0.3, 0.2, 0.1, 0.3, 0.5, 0.3, 0.7, 0.8])

        w_matrix_val_matchup1, w_matrix_row_matchup1, \
            w_matrix_col_matchup1, w_matrix_nnz_matchup1, \
                u_matrix_row_count_matchup1, u_matrix_matchup1 \
                    = return_w_matrix_variables([w2_matchup1], [u_matrix2_matchup1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 0.0],
                              [3.4, 4.3, 0.0],
                              [3.1, 4.4, 0.0],
                              [3.2, 4.3, 0.0]])

        Us1_matchup2 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 2, 3])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 0.0],
                              [2.2, 3.7, 0.0],
                              [2.2, 4.4, 0.0],
                              [2.1, 4.7, 0.0]])

        Us2_matchup2 = array([[0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 2, 3])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # iv. w and u matrices
        w_matrix_use1_matchup2 = array([0, 0, 1])
        w_matrix_use2_matchup2 = array([0, 0, 1])
        u_matrix_use1_matchup2 = array([0, 0, 1])
        u_matrix_use2_matchup2 = array([0, 0, 2])

        w12_matchup2 = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0],
                              [0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                              [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                              [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        w12_matchup2 = csr_matrix(w12_matchup2)

        u_matrix1_matchup2 = array([0.2, 0.4, 0.3, 0.4, 0.3, 0.2, 0.5, 0.3])
        u_matrix2_matchup2 = array([0.2, 0.3, 0.2, 0.1, 0.3, 0.5, 0.3, 0.7, 0.8])

        w_matrix_val_matchup2, w_matrix_row_matchup2,\
            w_matrix_col_matchup2, w_matrix_nnz_matchup2, \
                u_matrix_row_count_matchup2, u_matrix_matchup2\
                    = return_w_matrix_variables([w12_matchup2],
                                                [u_matrix1_matchup2, u_matrix2_matchup2])

        # --------------------------------------------------------------------------------------------------------------

        # c. Write match-up data to file -------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)
        append_W_to_input_file(fname_matchup1,
                               w_matrix_val_matchup1, w_matrix_row_matchup1,
                               w_matrix_col_matchup1, w_matrix_nnz_matchup1,
                               u_matrix_row_count_matchup1, u_matrix_matchup1,
                               w_matrix_use1_matchup1, w_matrix_use2_matchup1,
                               u_matrix_use1_matchup1, u_matrix_use2_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)
        append_W_to_input_file(fname_matchup2,
                               w_matrix_val_matchup2, w_matrix_row_matchup2,
                               w_matrix_col_matchup2, w_matrix_nnz_matchup2,
                               u_matrix_row_count_matchup2, u_matrix_matchup2,
                               w_matrix_use1_matchup2, w_matrix_use2_matchup2,
                               u_matrix_use1_matchup2, u_matrix_use2_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = [Uncertainty("r", array([1.6, 1.5, 1.5, 1.3, 1.5])),
                        Uncertainty("r", array([3.1, 3.2, 3.2, 3.1, 3.0])),
                        Uncertainty("r", array([3.3, 3.4, 3.1, 3.2])),
                        Uncertainty("r", array([2.1, 2.2, 2.2, 2.1])),
                        Uncertainty("rs", (array([5.0, 4.7, 5.1, 5.2, 5.3]), 1.0)),
                        Uncertainty("rs", (array([4.2, 4.3, 4.4, 4.3]), 1.0)),
                        Uncertainty("rs", (array([4.0, 3.7, 4.4, 4.7]), 2.0)),
                        Uncertainty("ave", (0, 0)),
                        Uncertainty("ave", (1, 1)),
                        Uncertainty("ave", (2, 2))]
        w_matrices_expected = [w2_matchup1, w12_matchup2]
        u_matrices_expected = [u_matrix2_matchup1,
                                        u_matrix1_matchup2,
                                        u_matrix2_matchup2]
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = [Uncertainty("r", array([0.25, 0.25, 0.25, 0.25, 0.25])),
                         Uncertainty("r", array([0.2644, 0.2644, 0.2644, 0.2644]))]
        time1_expected = array([dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4),
                                dt(1970, 1, 1, 1, 0, 5),
                                dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4)])
        time2_expected = array([dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000),
                                dt(1970, 1, 1, 1, 0, 5, 100000),
                                dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000)])
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2])

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.form, block_unc_test.form)
            if block_unc_expected == "r":
                self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())
            if block_unc_expected == "rs":
                self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())
                self.assertEqual(block_unc_expected.uS.tolist(), block_unc_test.uS.tolist())
            if block_unc_expected == "ave":
                self.assertEqual(block_unc_expected.u_i, block_unc_test.u_i)
                self.assertEqual(block_unc_expected.w_i, block_unc_test.w_i)

        # c. w_matrices
        for w_matrix_expected, w_matrix_test in zip(w_matrices_expected, w_matrices_test):
            self.assertEqual(w_matrix_expected.data.tolist(), w_matrix_test.data.tolist())
            self.assertEqual(w_matrix_expected.indices.tolist(), w_matrix_test.indices.tolist())
            self.assertEqual(w_matrix_expected.indptr.tolist(), w_matrix_test.indptr.tolist())

        # d. u_matrices
        for u_matrix_expected, u_matrix_test in \
                zip(u_matrices_expected, u_matrices_test):
            self.assertEqual(u_matrix_expected.tolist(), u_matrix_test.tolist())

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.form, block_unck_test.form)
            self.assertEqual(block_unck_expected.uR.tolist(), block_unck_test.uR.tolist())

        # g. time
        self.assertEqual(time1_test.tolist(), time1_expected.tolist())
        self.assertEqual(time2_test.tolist(), time2_expected.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_r_w_reusedw(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with mulitple match-up series datasets
        - with random, random and systematic and w matrix uncertainty typea
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, with covariate 1 w matrix uncertainty type, covariate
        #    2 systematic uncertainty type and covariate 3 w matrix uncertainty type
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(temp_data_directory, str(int(random()*1000000)))
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 0.0, 0.0],
                              [3.2, 0.0, 0.0],
                              [3.2, 0.0, 0.0],
                              [3.1, 0.0, 0.0],
                              [3.0, 0.0, 0.0]])

        Us2_matchup1 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 3, 3])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # iv. w and u matrices
        w_matrix_use1_matchup1 = array([0])
        w_matrix_use2_matchup1 = array([0, 1, 1])
        u_matrix_use1_matchup1 = array([0])
        u_matrix_use2_matchup1 = array([0, 1, 1])

        w2_matchup1 = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        w2_matchup1 = csr_matrix(w2_matchup1)

        u_matrix2_matchup1 = array([0.2, 0.3, 0.2, 0.1, 0.3, 0.5, 0.3, 0.7, 0.8])

        w_matrix_val_matchup1, w_matrix_row_matchup1, \
        w_matrix_col_matchup1, w_matrix_nnz_matchup1, \
        u_matrix_row_count_matchup1, u_matrix_matchup1 \
            = return_w_matrix_variables([w2_matchup1], [u_matrix2_matchup1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 0.0, 0.0],
                              [3.4, 0.0, 0.0],
                              [3.1, 0.0, 0.0],
                              [3.2, 0.0, 0.0]])

        Us1_matchup2 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 3, 3])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 0.0, 0.0],
                              [2.2, 0.0, 0.0],
                              [2.2, 0.0, 0.0],
                              [2.1, 0.0, 0.0]])

        Us2_matchup2 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 3, 3])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # iv. w and u matrices
        w_matrix_use1_matchup2 = array([0, 1, 1])
        w_matrix_use2_matchup2 = array([0, 2, 2])
        u_matrix_use1_matchup2 = array([0, 1, 2])
        u_matrix_use2_matchup2 = array([0, 3, 4])

        w1_matchup2 = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0],
                             [0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                             [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                             [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        w1_matchup2 = csr_matrix(w1_matchup2)
        w2_matchup2 = array([[0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.0, 0.0, 0.0],
                             [0.0, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.0, 0.0],
                             [0.0, 0.0, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.0],
                             [0.0, 0.0, 0.0, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429, 0.1429]])
        w2_matchup2 = csr_matrix(w2_matchup2)

        u_matrix1_matchup2 = array([0.2, 0.4, 0.3, 0.4, 0.3, 0.2, 0.5, 0.3])
        u_matrix2_matchup2 = array([0.2, 0.4, 0.5, 0.4, 0.7, 0.3, 0.5, 0.3])
        u_matrix3_matchup2 = array([0.2, 0.3, 0.8, 0.2, 0.3, 0.5, 0.8, 0.7, 0.8, 0.6])
        u_matrix4_matchup2 = array([0.2, 0.3, 0.2, 0.1, 0.3, 0.5, 0.3, 0.7, 0.8, 0.3])

        w_matrix_val_matchup2, w_matrix_row_matchup2, \
        w_matrix_col_matchup2, w_matrix_nnz_matchup2, \
        u_matrix_row_count_matchup2, u_matrix_matchup2 \
            = return_w_matrix_variables([w1_matchup2, w2_matchup2],
                                        [u_matrix1_matchup2, u_matrix2_matchup2,
                                         u_matrix3_matchup2, u_matrix4_matchup2])

        # --------------------------------------------------------------------------------------------------------------

        # c. Write match-up data to file -------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)
        append_W_to_input_file(fname_matchup1,
                               w_matrix_val_matchup1, w_matrix_row_matchup1,
                               w_matrix_col_matchup1, w_matrix_nnz_matchup1,
                               u_matrix_row_count_matchup1, u_matrix_matchup1,
                               w_matrix_use1_matchup1, w_matrix_use2_matchup1,
                               u_matrix_use1_matchup1, u_matrix_use2_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)
        append_W_to_input_file(fname_matchup2,
                               w_matrix_val_matchup2, w_matrix_row_matchup2,
                               w_matrix_col_matchup2, w_matrix_nnz_matchup2,
                               u_matrix_row_count_matchup2, u_matrix_matchup2,
                               w_matrix_use1_matchup2, w_matrix_use2_matchup2,
                               u_matrix_use1_matchup2, u_matrix_use2_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = [Uncertainty("r", array([1.6, 1.5, 1.5, 1.3, 1.5])),
                        Uncertainty("r", array([3.1, 3.2, 3.2, 3.1, 3.0])),
                        Uncertainty("r", array([3.3, 3.4, 3.1, 3.2])),
                        Uncertainty("r", array([2.1, 2.2, 2.2, 2.1])),
                        Uncertainty("ave", (0, 0)),
                        Uncertainty("ave", (0, 1)),
                        Uncertainty("ave", (1, 2)),
                        Uncertainty("ave", (1, 3)),
                        Uncertainty("ave", (2, 4)),
                        Uncertainty("ave", (2, 5))]
        w_matrices_expected = [w2_matchup1, w1_matchup2, w2_matchup2]
        u_matrices_expected = [u_matrix2_matchup1,
                                        u_matrix1_matchup2,
                                        u_matrix2_matchup2,
                                        u_matrix3_matchup2,
                                        u_matrix4_matchup2]
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = [Uncertainty("r", array([0.25, 0.25, 0.25, 0.25, 0.25])),
                         Uncertainty("r", array([0.2644, 0.2644, 0.2644, 0.2644]))]
        time1_expected = array([dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4),
                                dt(1970, 1, 1, 1, 0, 5),
                                dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4)])
        time2_expected = array([dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000),
                                dt(1970, 1, 1, 1, 0, 5, 100000),
                                dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000)])
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test, \
        u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
            idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2])

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.form, block_unc_test.form)
            if block_unc_expected == "r":
                self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())
            if block_unc_expected == "ave":
                self.assertEqual(block_unc_expected.u_i, block_unc_test.u_i)
                self.assertEqual(block_unc_expected.w_i, block_unc_test.w_i)

        # c. w_matrices
        for w_matrix_expected, w_matrix_test in zip(w_matrices_expected, w_matrices_test):
            self.assertEqual(w_matrix_expected.data.tolist(), w_matrix_test.data.tolist())
            self.assertEqual(w_matrix_expected.indices.tolist(), w_matrix_test.indices.tolist())
            self.assertEqual(w_matrix_expected.indptr.tolist(), w_matrix_test.indptr.tolist())

        # d. u_matrices
        for u_matrix_expected, u_matrix_test in \
                zip(u_matrices_expected, u_matrices_test):
            self.assertEqual(u_matrix_expected.tolist(), u_matrix_test.tolist())

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.form, block_unck_test.form)
            self.assertEqual(block_unck_expected.uR.tolist(), block_unck_test.uR.tolist())

        # g. time
        self.assertEqual(time1_test.tolist(), time1_expected.tolist())
        self.assertEqual(time2_test.tolist(), time2_expected.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_r___notime(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with multiple matchup series datasets
        - with only random uncertainty type
        - option to open time switched off
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, each with uncertainty types random
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(dirname(__file__), "temp_data_directory")
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 1, 1])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 2.3],
                              [3.4, 4.3, 1.2],
                              [3.1, 4.4, 2.3],
                              [3.2, 4.3, 4.4]])

        Us1_matchup2 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 1, 1])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 3.2],
                              [2.2, 3.7, 2.7],
                              [2.2, 4.4, 3.0],
                              [2.1, 4.7, 5.3]])

        Us2_matchup2 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 1, 1])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # --------------------------------------------------------------------------------------------------------------

        # d. Test data to file -----------------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = [Uncertainty("r", array([1.6, 1.5, 1.5, 1.3, 1.5])),
                        Uncertainty("r", array([3.1, 3.2, 3.2, 3.1, 3.0])),
                        Uncertainty("r", array([3.3, 3.4, 3.1, 3.2])),
                        Uncertainty("r", array([2.1, 2.2, 2.2, 2.1])),
                        Uncertainty("r", array([5.0, 4.7, 5.1, 5.2, 5.3])),
                        Uncertainty("r", array([4.2, 4.3, 4.4, 4.3])),
                        Uncertainty("r", array([4.0, 3.7, 4.4, 4.7])),
                        Uncertainty("r", array([2.2, 1.7, 2.0, 4.3, 2.6])),
                        Uncertainty("r", array([2.3, 1.2, 2.3, 4.4])),
                        Uncertainty("r", array([3.2, 2.7, 3.0, 5.3]))]
        w_matrices_expected = []
        u_matrices_expected = []
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = [Uncertainty("r", array([0.25, 0.25, 0.25, 0.25, 0.25])),
                         Uncertainty("r", array([0.2644, 0.2644, 0.2644, 0.2644]))]
        time1_expected = None
        time2_expected = None
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2], open_time=False)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.form, block_unc_test.form)
            self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.form, block_unck_test.form)
            self.assertEqual(block_unck_expected.uR.tolist(), block_unck_test.uR.tolist())

        # g. time
        self.assertEqual(time1_test, time1_expected)
        self.assertEqual(time2_test, time2_expected)

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_rs__notime(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with multiple match-up series datasets
        - with random and random and systematic uncertainty type
        - option to open time switched off
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, with covariate 1 and 3 random uncertainty type and
        #    covariate 2 with systematic uncertainty type
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(temp_data_directory, str(int(random()*1000000)))
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 2, 1])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 2.3],
                              [3.4, 4.3, 1.2],
                              [3.1, 4.4, 2.3],
                              [3.2, 4.3, 4.4]])

        Us1_matchup2 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 2, 1])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 3.2],
                              [2.2, 3.7, 2.7],
                              [2.2, 4.4, 3.0],
                              [2.1, 4.7, 5.3]])

        Us2_matchup2 = array([[0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 2, 1])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # --------------------------------------------------------------------------------------------------------------

        # d. Test data to file -----------------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = [Uncertainty("r", array([1.6, 1.5, 1.5, 1.3, 1.5])),
                        Uncertainty("r", array([3.1, 3.2, 3.2, 3.1, 3.0])),
                        Uncertainty("r", array([3.3, 3.4, 3.1, 3.2])),
                        Uncertainty("r", array([2.1, 2.2, 2.2, 2.1])),
                        Uncertainty("rs", (array([5.0, 4.7, 5.1, 5.2, 5.3]), 1.0)),
                        Uncertainty("rs", (array([4.2, 4.3, 4.4, 4.3]), 1.0)),
                        Uncertainty("rs", (array([4.0, 3.7, 4.4, 4.7]), 2.0)),
                        Uncertainty("r", array([2.2, 1.7, 2.0, 4.3, 2.6])),
                        Uncertainty("r", array([2.3, 1.2, 2.3, 4.4])),
                        Uncertainty("r", array([3.2, 2.7, 3.0, 5.3]))]
        w_matrices_expected = []
        u_matrices_expected = []
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = [Uncertainty("r", array([0.25, 0.25, 0.25, 0.25, 0.25])),
                         Uncertainty("r", array([0.2644, 0.2644, 0.2644, 0.2644]))]
        time1_expected = None
        time2_expected = None
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2], open_time=None)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.form, block_unc_test.form)
            self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())
            if block_unc_expected == "rs":
                self.assertEqual(block_unc_expected.uS.tolist(), block_unc_test.uS.tolist())

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.form, block_unck_test.form)
            self.assertEqual(block_unck_expected.uR.tolist(), block_unck_test.uR.tolist())

        # g. time
        self.assertEqual(time1_test, time1_expected)
        self.assertEqual(time2_test, time2_expected)

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_rsw_notime(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with mulitple match-up series datasets
        - with random, random and systematic and w matrix uncertainty type
        - option to open time switched off
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, with covariate 1 w matrix uncertainty type, covariate
        #    2 systematic uncertainty type and covariate 3 w matrix uncertainty type
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(temp_data_directory, str(int(random()*1000000)))
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 2, 3])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # iv. w and u matrices
        w_matrix_use1_matchup1 = array([0])
        w_matrix_use2_matchup1 = array([0, 0, 1])
        u_matrix_use1_matchup1 = array([0])
        u_matrix_use2_matchup1 = array([0, 0, 1])

        w2_matchup1 = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        w2_matchup1 = csr_matrix(w2_matchup1)

        u_matrix2_matchup1 = array([0.2, 0.3, 0.2, 0.1, 0.3, 0.5, 0.3, 0.7, 0.8])

        w_matrix_val_matchup1, w_matrix_row_matchup1, \
            w_matrix_col_matchup1, w_matrix_nnz_matchup1, \
                u_matrix_row_count_matchup1, u_matrix_matchup1 \
                    = return_w_matrix_variables([w2_matchup1], [u_matrix2_matchup1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 0.0],
                              [3.4, 4.3, 0.0],
                              [3.1, 4.4, 0.0],
                              [3.2, 4.3, 0.0]])

        Us1_matchup2 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 2, 3])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 0.0],
                              [2.2, 3.7, 0.0],
                              [2.2, 4.4, 0.0],
                              [2.1, 4.7, 0.0]])

        Us2_matchup2 = array([[0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 2, 3])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # iv. w and u matrices
        w_matrix_use1_matchup2 = array([0, 0, 1])
        w_matrix_use2_matchup2 = array([0, 0, 1])
        u_matrix_use1_matchup2 = array([0, 0, 1])
        u_matrix_use2_matchup2 = array([0, 0, 2])

        w12_matchup2 = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0],
                              [0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                              [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                              [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        w12_matchup2 = csr_matrix(w12_matchup2)

        u_matrix1_matchup2 = array([0.2, 0.4, 0.3, 0.4, 0.3, 0.2, 0.5, 0.3])
        u_matrix2_matchup2 = array([0.2, 0.3, 0.2, 0.1, 0.3, 0.5, 0.3, 0.7, 0.8])

        w_matrix_val_matchup2, w_matrix_row_matchup2,\
            w_matrix_col_matchup2, w_matrix_nnz_matchup2, \
                u_matrix_row_count_matchup2, u_matrix_matchup2\
                    = return_w_matrix_variables([w12_matchup2],
                                                [u_matrix1_matchup2, u_matrix2_matchup2])

        # --------------------------------------------------------------------------------------------------------------

        # c. Write match-up data to file -------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)
        append_W_to_input_file(fname_matchup1,
                               w_matrix_val_matchup1, w_matrix_row_matchup1,
                               w_matrix_col_matchup1, w_matrix_nnz_matchup1,
                               u_matrix_row_count_matchup1, u_matrix_matchup1,
                               w_matrix_use1_matchup1, w_matrix_use2_matchup1,
                               u_matrix_use1_matchup1, u_matrix_use2_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)
        append_W_to_input_file(fname_matchup2,
                               w_matrix_val_matchup2, w_matrix_row_matchup2,
                               w_matrix_col_matchup2, w_matrix_nnz_matchup2,
                               u_matrix_row_count_matchup2, u_matrix_matchup2,
                               w_matrix_use1_matchup2, w_matrix_use2_matchup2,
                               u_matrix_use1_matchup2, u_matrix_use2_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = [Uncertainty("r", array([1.6, 1.5, 1.5, 1.3, 1.5])),
                        Uncertainty("r", array([3.1, 3.2, 3.2, 3.1, 3.0])),
                        Uncertainty("r", array([3.3, 3.4, 3.1, 3.2])),
                        Uncertainty("r", array([2.1, 2.2, 2.2, 2.1])),
                        Uncertainty("rs", (array([5.0, 4.7, 5.1, 5.2, 5.3]), 1.0)),
                        Uncertainty("rs", (array([4.2, 4.3, 4.4, 4.3]), 1.0)),
                        Uncertainty("rs", (array([4.0, 3.7, 4.4, 4.7]), 2.0)),
                        Uncertainty("ave", (0, 0)),
                        Uncertainty("ave", (1, 1)),
                        Uncertainty("ave", (2, 2))]
        w_matrices_expected = [w2_matchup1, w12_matchup2]
        u_matrices_expected = [u_matrix2_matchup1,
                                        u_matrix1_matchup2,
                                        u_matrix2_matchup2]
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = [Uncertainty("r", array([0.25, 0.25, 0.25, 0.25, 0.25])),
                         Uncertainty("r", array([0.2644, 0.2644, 0.2644, 0.2644]))]
        time1_expected = array([dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4),
                                dt(1970, 1, 1, 1, 0, 5),
                                dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4)])
        time2_expected = array([dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000),
                                dt(1970, 1, 1, 1, 0, 5, 100000),
                                dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000)])
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2])

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        for block_unc_test, block_unc_expected in zip(unc_test, unc_expected):
            self.assertEqual(block_unc_expected.form, block_unc_test.form)
            if block_unc_expected == "r":
                self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())
            if block_unc_expected == "rs":
                self.assertEqual(block_unc_expected.uR.tolist(), block_unc_test.uR.tolist())
                self.assertEqual(block_unc_expected.uS.tolist(), block_unc_test.uS.tolist())
            if block_unc_expected == "ave":
                self.assertEqual(block_unc_expected.u_i, block_unc_test.u_i)
                self.assertEqual(block_unc_expected.w_i, block_unc_test.w_i)

        # c. w_matrices
        for w_matrix_expected, w_matrix_test in zip(w_matrices_expected, w_matrices_test):
            self.assertEqual(w_matrix_expected.data.tolist(), w_matrix_test.data.tolist())
            self.assertEqual(w_matrix_expected.indices.tolist(), w_matrix_test.indices.tolist())
            self.assertEqual(w_matrix_expected.indptr.tolist(), w_matrix_test.indptr.tolist())

        # d. u_matrices
        for u_matrix_expected, u_matrix_test in \
                zip(u_matrices_expected, u_matrices_test):
            self.assertEqual(u_matrix_expected.tolist(), u_matrix_test.tolist())

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        for block_unck_test, block_unck_expected in zip(unck_test, unck_expected):
            self.assertEqual(block_unck_expected.form, block_unck_test.form)
            self.assertEqual(block_unck_expected.uR.tolist(), block_unck_test.uR.tolist())

        # g. time
        self.assertEqual(time1_test.tolist(), time1_expected.tolist())
        self.assertEqual(time2_test.tolist(), time2_expected.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_r___nounc(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with multiple matchup series datasets
        - with only random uncertainty type
        - option to open uncertainty switched off
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, each with uncertainty types random
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(dirname(__file__), "temp_data_directory")
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 1, 1])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 2.3],
                              [3.4, 4.3, 1.2],
                              [3.1, 4.4, 2.3],
                              [3.2, 4.3, 4.4]])

        Us1_matchup2 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 1, 1])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 3.2],
                              [2.2, 3.7, 2.7],
                              [2.2, 4.4, 3.0],
                              [2.1, 4.7, 5.3]])

        Us2_matchup2 = array([[0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0],
                              [0.0, 0.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 1, 1])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # --------------------------------------------------------------------------------------------------------------

        # c. Parameter data --------------------------------------------------------------------------------------------

        parameter = array([[1.0, 2.0],
                           [3.0, 4.0]])

        # --------------------------------------------------------------------------------------------------------------

        # d. Test data to file -----------------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = None
        w_matrices_expected = None
        u_matrices_expected = None
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = None
        time1_expected = array([dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4),
                                dt(1970, 1, 1, 1, 0, 5),
                                dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4)])
        time2_expected = array([dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000),
                                dt(1970, 1, 1, 1, 0, 5, 100000),
                                dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000)])
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2], open_uncertainty=False)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        self.assertEqual(unc_expected, unc_test)

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        self.assertEqual(unck_expected, unck_test)

        # g. time
        self.assertEqual(time1_test.tolist(), time1_expected.tolist())
        self.assertEqual(time2_test.tolist(), time2_expected.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_rs__nounc(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with multiple match-up series datasets
        - with random and random and systematic uncertainty type
        - with option to open uncertainties switched off
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, with covariate 1 and 3 random uncertainty type and
        #    covariate 2 with systematic uncertainty type
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(temp_data_directory, str(int(random()*1000000)))
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 2, 1])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 2.3],
                              [3.4, 4.3, 1.2],
                              [3.1, 4.4, 2.3],
                              [3.2, 4.3, 4.4]])

        Us1_matchup2 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 2, 1])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 3.2],
                              [2.2, 3.7, 2.7],
                              [2.2, 4.4, 3.0],
                              [2.1, 4.7, 5.3]])

        Us2_matchup2 = array([[0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 2, 1])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # --------------------------------------------------------------------------------------------------------------

        # d. Test data to file -----------------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = None
        w_matrices_expected = None
        u_matrices_expected = None
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = None
        time1_expected = array([dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4),
                                dt(1970, 1, 1, 1, 0, 5),
                                dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4)])
        time2_expected = array([dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000),
                                dt(1970, 1, 1, 1, 0, 5, 100000),
                                dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000)])
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2], open_uncertainty=False)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        self.assertEqual(unc_expected, unc_test)

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        self.assertEqual(unck_expected, unck_test)

        # g. time
        self.assertEqual(time1_test.tolist(), time1_expected.tolist())
        self.assertEqual(time2_test.tolist(), time2_expected.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)

    def test_openMatchUpData_multi_rsw_nounc(self):
        """
        Test for MatchUpData.openMatchUpData() method for case with mulitple match-up series datasets
        - with random, random and systematic and w matrix uncertainty type
        - with option to open uncertainties switched off
        """

        # Test Description
        # ================
        #
        # 1. This test writes two test match-up data files:
        #    + Reference - Sensor A
        #    + Sensor A - Sensor B
        #    Sensor A/B have three measurement function variables, with covariate 1 w matrix uncertainty type, covariate
        #    2 systematic uncertainty type and covariate 3 w matrix uncertainty type
        #
        # 2. The files are read by the match-up data reader
        #
        # 3. The MatchUpData object in memory is compared to the expected values

        ################################################################################################################
        # 1. Write test match-up data files
        ################################################################################################################

        # Define file paths
        test_data_directory = pjoin(temp_data_directory, str(int(random()*1000000)))
        if not exists(test_data_directory):
            makedirs(test_data_directory)

        fname_matchup1 = pjoin(test_data_directory, "matchup1.nc")
        fname_matchup2 = pjoin(test_data_directory, "matchup2.nc")

        # a. Reference - Sensor A match-up data ------------------------------------------------------------------------

        # i. Reference Sensor data
        sensor_1_name_matchup1 = -1

        X1_matchup1 = array([[16.2, 11.2, 15.1, 20.3, 18.1]]).T
        Ur1_matchup1 = array([1.6, 1.5, 1.5, 1.3, 1.5])
        Us1_matchup1 = array([0.0, 0.0, 0.0, 0.0, 0.0])
        uncertainty_type1_matchup1 = array([1])

        # ii. Sensor A data
        sensor_2_name_matchup1 = 1

        X2_matchup1 = array([[70.5, 150.5, 30.2],
                             [70.6, 151.1, 20.4],
                             [70.3, 149.8, 28.2],
                             [70.7, 150.2, 50.7],
                             [70.5, 151.4, 45.6]])

        Ur2_matchup1 = array([[3.1, 5.0, 2.2],
                              [3.2, 4.7, 1.7],
                              [3.2, 5.1, 2.0],
                              [3.1, 5.2, 4.3],
                              [3.0, 5.3, 2.6]])

        Us2_matchup1 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type2_matchup1 = array([1, 2, 3])

        # iii. Match-up data
        K_matchup1 = array([1.2, 1.7, 1.3, 1.4, 1.3])
        Kr_matchup1 = array([0.3, 0.3, 0.3, 0.3, 0.3])
        Ks_matchup1 = array([0.4, 0.4, 0.4, 0.4, 0.4])
        time1_matchup1 = array([1.0, 2.0, 3.0, 4.0, 5.0])
        time2_matchup1 = array([1.1, 2.1, 3.1, 4.1, 5.1])

        # iv. w and u matrices
        w_matrix_use1_matchup1 = array([0])
        w_matrix_use2_matchup1 = array([0, 0, 1])
        u_matrix_use1_matchup1 = array([0])
        u_matrix_use2_matchup1 = array([0, 0, 1])

        w2_matchup1 = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        w2_matchup1 = csr_matrix(w2_matchup1)

        u_matrix2_matchup1 = array([0.2, 0.3, 0.2, 0.1, 0.3, 0.5, 0.3, 0.7, 0.8])

        w_matrix_val_matchup1, w_matrix_row_matchup1, \
            w_matrix_col_matchup1, w_matrix_nnz_matchup1, \
                u_matrix_row_count_matchup1, u_matrix_matchup1 \
                    = return_w_matrix_variables([w2_matchup1], [u_matrix2_matchup1])

        # --------------------------------------------------------------------------------------------------------------

        # b. Sensor A - Sensor B match-up data -------------------------------------------------------------------------

        # i. Sensor A data
        sensor_1_name_matchup2 = 1

        X1_matchup2 = array([[71.5, 140.5, 29.2],
                             [71.6, 141.1, 37.4],
                             [71.3, 139.8, 28.2],
                             [71.7, 140.2, 50.7]])

        Ur1_matchup2 = array([[3.3, 4.2, 0.0],
                              [3.4, 4.3, 0.0],
                              [3.1, 4.4, 0.0],
                              [3.2, 4.3, 0.0]])

        Us1_matchup2 = array([[0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0],
                              [0.0, 1.0, 0.0]])

        uncertainty_type1_matchup2 = array([1, 2, 3])

        # ii. Sensor B data
        sensor_2_name_matchup2 = 2

        X2_matchup2 = array([[80.5, 160.5, 28.2],
                             [80.6, 161.1, 32.4],
                             [80.3, 169.8, 22.2],
                             [80.7, 160.2, 53.7]])

        Ur2_matchup2 = array([[2.1, 4.0, 0.0],
                              [2.2, 3.7, 0.0],
                              [2.2, 4.4, 0.0],
                              [2.1, 4.7, 0.0]])

        Us2_matchup2 = array([[0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0],
                              [0.0, 2.0, 0.0]])

        uncertainty_type2_matchup2 = array([1, 2, 3])

        # iii. Match-up data
        K_matchup2 = array([3.2, 3.7, 3.3, 3.4])
        Kr_matchup2 = array([0.5, 0.5, 0.5, 0.5])
        Ks_matchup2 = array([0.12, 0.12, 0.12, 0.12])
        time1_matchup2 = array([1.0, 2.0, 3.0, 4.0])
        time2_matchup2 = array([1.1, 2.1, 3.1, 4.1])

        # iv. w and u matrices
        w_matrix_use1_matchup2 = array([0, 0, 1])
        w_matrix_use2_matchup2 = array([0, 0, 1])
        u_matrix_use1_matchup2 = array([0, 0, 1])
        u_matrix_use2_matchup2 = array([0, 0, 2])

        w12_matchup2 = array([[0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0],
                              [0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0],
                              [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0],
                              [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2]])
        w12_matchup2 = csr_matrix(w12_matchup2)

        u_matrix1_matchup2 = array([0.2, 0.4, 0.3, 0.4, 0.3, 0.2, 0.5, 0.3])
        u_matrix2_matchup2 = array([0.2, 0.3, 0.2, 0.1, 0.3, 0.5, 0.3, 0.7, 0.8])

        w_matrix_val_matchup2, w_matrix_row_matchup2,\
            w_matrix_col_matchup2, w_matrix_nnz_matchup2, \
                u_matrix_row_count_matchup2, u_matrix_matchup2\
                    = return_w_matrix_variables([w12_matchup2],
                                                [u_matrix1_matchup2, u_matrix2_matchup2])

        # --------------------------------------------------------------------------------------------------------------

        # c. Write match-up data to file -------------------------------------------------------------------------------

        # i. Match-up 1
        write_input_file(fname_matchup1,
                         X1_matchup1, X2_matchup1,
                         Ur1_matchup1, Ur2_matchup1, Us1_matchup1, Us2_matchup1,
                         uncertainty_type1_matchup1, uncertainty_type2_matchup1,
                         K_matchup1, Kr_matchup1, Ks_matchup1,
                         time1_matchup1, time2_matchup1,
                         sensor_1_name_matchup1, sensor_2_name_matchup1)
        append_W_to_input_file(fname_matchup1,
                               w_matrix_val_matchup1, w_matrix_row_matchup1,
                               w_matrix_col_matchup1, w_matrix_nnz_matchup1,
                               u_matrix_row_count_matchup1, u_matrix_matchup1,
                               w_matrix_use1_matchup1, w_matrix_use2_matchup1,
                               u_matrix_use1_matchup1, u_matrix_use2_matchup1)

        # ii. Match-up 2
        write_input_file(fname_matchup2,
                         X1_matchup2, X2_matchup2,
                         Ur1_matchup2, Ur2_matchup2, Us1_matchup2, Us2_matchup2,
                         uncertainty_type1_matchup2, uncertainty_type2_matchup2,
                         K_matchup2, Kr_matchup2, Ks_matchup2,
                         time1_matchup2, time2_matchup2,
                         sensor_1_name_matchup2, sensor_2_name_matchup2)
        append_W_to_input_file(fname_matchup2,
                               w_matrix_val_matchup2, w_matrix_row_matchup2,
                               w_matrix_col_matchup2, w_matrix_nnz_matchup2,
                               u_matrix_row_count_matchup2, u_matrix_matchup2,
                               w_matrix_use1_matchup2, w_matrix_use2_matchup2,
                               u_matrix_use1_matchup2, u_matrix_use2_matchup2)

        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        # Expected value of HData attributes
        values_expected = array([16.2, 11.2, 15.1, 20.3, 18.1,
                                 70.5, 70.6, 70.3, 70.7, 70.5,
                                 71.5, 71.6, 71.3, 71.7,
                                 80.5, 80.6, 80.3, 80.7,
                                 150.5, 151.1, 149.8, 150.2, 151.4,
                                 140.5, 141.1, 139.8, 140.2,
                                 160.5, 161.1, 169.8, 160.2,
                                 30.2, 20.4, 28.2, 50.7, 45.6,
                                 29.2, 37.4, 28.2, 50.7,
                                 28.2, 32.4, 22.2, 53.7, ])
        unc_expected = None
        w_matrices_expected = None
        u_matrices_expected = None
        ks_expected = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
        unck_expected = None
        time1_expected = array([dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4),
                                dt(1970, 1, 1, 1, 0, 5),
                                dt(1970, 1, 1, 1, 0, 1),
                                dt(1970, 1, 1, 1, 0, 2),
                                dt(1970, 1, 1, 1, 0, 3),
                                dt(1970, 1, 1, 1, 0, 4)])
        time2_expected = array([dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000),
                                dt(1970, 1, 1, 1, 0, 5, 100000),
                                dt(1970, 1, 1, 1, 0, 1, 100000),
                                dt(1970, 1, 1, 1, 0, 2, 100000),
                                dt(1970, 1, 1, 1, 0, 3, 100000),
                                dt(1970, 1, 1, 1, 0, 4, 100000)])
        idx_expected = {"Nm": [5, 4],
                        "cNm": [0, 5, 9],
                        "Im": [[0, 1], [1, 2]],
                        "sensors": [-1, 1, 2],
                        "sensor_ms": [1, 3, 3],
                        "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
                        "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
                        "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
                        "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
                        "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44]}

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        MatchUpOp = MatchUp()

        values_test, unc_test, w_matrices_test,\
            u_matrices_test, ks_test, unck_test, time1_test, time2_test, \
                idx_test = MatchUpOp.openMatchUpData([fname_matchup1, fname_matchup2], open_uncertainty=False)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        # Test HData object attribute by attribute

        # a. values
        self.assertEqual(values_test.tolist(), values_expected.tolist())

        # b. unc
        self.assertEqual(unc_expected, unc_test)

        # c. w_matrices
        self.assertEqual(w_matrices_test, w_matrices_expected)

        # d. u_matrices
        self.assertEqual(u_matrices_test, u_matrices_expected)

        # e. ks
        self.assertEqual(ks_test.tolist(), ks_expected.tolist())

        # f. unck
        self.assertEqual(unck_expected, unck_test)

        # g. time
        self.assertEqual(time1_test.tolist(), time1_expected.tolist())
        self.assertEqual(time2_test.tolist(), time2_expected.tolist())

        # h. idx
        self.assertEqual(set(idx_expected.keys()), set(idx_test.keys()))
        for key in idx_expected.keys():
            idx_i_test = idx_test[key]
            idx_i_expected = idx_expected[key]
            if isinstance(idx_i_expected, ndarray):
                self.assertEqual(idx_i_test.tolist(), idx_i_expected.tolist())
            else:
                self.assertEqual(idx_i_test, idx_i_expected)

        ################################################################################################################
        # 5. Remove Test Data
        ################################################################################################################

        rmtree(test_data_directory)


if __name__ == '__main__':
    unittest.main()
