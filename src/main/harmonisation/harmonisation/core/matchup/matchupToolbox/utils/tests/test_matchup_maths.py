"""
Test of matchup_maths module functions
"""

'''___Built-In Modules___'''
import unittest

'''___Third-Party Modules___'''
from numpy import array, vstack, ones, dot, column_stack, eye

'''___NPL Modules___'''
from harmonisation import MatchUp, Uncertainty, evaluate_measurand, evaluate_adjusted_measurand, evaluate_K


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "23/11/2017"
__credits__ = ["Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def return_MatchUpTest():
    """
    Return a MatchUp dataset object for testing

    :return:
        :MatchUpTest: *eopy.matchup.matchupIO.MatchUp*

        Test match-up dataset
    """

    ####################################################################################################################
    # 1. Define test sensor function
    ####################################################################################################################

    def test_reference_function(X, a, sensor_c, sensor_xt_i, sensor_at_i, sensor_t):
        return X[:, 0], None

    def test_sensor_function(X, a, sensor_c, sensor_xt_i, sensor_at_i, sensor_t):
        """
        Arbitary sensor function

        :type a: numpy.ndarray
        :param a: calibration parameters

        :type X: list
        :param X: List of arrays of sensor observed parameters

        :return:
            :measurand: *numpy.ndarray*

            Computed sensor measurand
        """

        # Extract observed parameters
        X1 = X[:, 0]
        X2 = X[:, 1]
        X3 = X[:, 2]
        M = len(X1)

        # Evaluate measurand
        parameter_derivatives = vstack((ones(M), X1 * X2 / X3, X1 ** 2)).T
        parameter = [a[0], a[1], a[2]]
        measurand = dot(parameter_derivatives, parameter)

        # Evaluate derivatives
        # In the following order:
        # > d(measurand)/dX1
        # > d(measurand)/dX2
        # > d(measurand)/dX3
        # > d(measurand)/da0
        # > d(measurand)/da1
        # > d(measurand)/da2

        derivatives = column_stack((parameter[1] * X1 / X2 + 2 * parameter[2] * X1,
                                    parameter[1] * X1 / X3,
                                    -parameter[1] * X2 * X1 / X3 ** 2,
                                    parameter_derivatives))

        return measurand, derivatives

    def test_adjustment_model(measurand):
        """
        Arbitary sensor adjustment function to sample sensor 1

        :type measurand: numpy.ndarray
        :param measurand: measurand data

        :return:
            :adjusted_measurand: *numpy.ndarray*

            Adjusted sensor measurand

            :adjusted_measurand_derivatives: *numpy.ndarray*

            Adjusted sensor measurand derivatives
        """

        adjusted_measurand = 2 * measurand
        adjusted_measurand_derivatives = ones(len(adjusted_measurand))
        return adjusted_measurand, adjusted_measurand_derivatives

    ####################################################################################################################
    # 2. Initialise test data
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
    unc = [Uncertainty("r", array([1.6, 1.5, 1.5, 1.3, 1.5])),
           Uncertainty("r", array([3.1, 3.2, 3.2, 3.1, 3.0])),
           Uncertainty("r", array([3.3, 3.4, 3.1, 3.2])),
           Uncertainty("r", array([2.1, 2.2, 2.2, 2.1])),
           Uncertainty("r", array([5.0, 4.7, 5.1, 5.2, 5.3])),
           Uncertainty("r", array([4.2, 4.3, 4.4, 4.3])),
           Uncertainty("r", array([4.0, 3.7, 4.4, 4.7])),
           Uncertainty("r", array([2.2, 1.7, 2.0, 4.3, 2.6])),
           Uncertainty("r", array([2.3, 1.2, 2.3, 4.4])),
           Uncertainty("r", array([3.2, 2.7, 3.0, 5.3]))]
    ks = array([1.2, 1.7, 1.3, 1.4, 1.3, 3.2, 3.7, 3.3, 3.4])
    unck = [Uncertainty("r", array([0.25, 0.25, 0.25, 0.25, 0.25])),
            Uncertainty("r", array([0.2644, 0.2644, 0.2644, 0.2644]))]
    idx = {"Nm": [5, 4],
           "cNm": [0, 5, 9],
           "Im": [[0, 1], [1, 2]],
           "sensors": [-1, 1, 2],
           "sensor_m": [1, 3, 3],
           "n_sensor": [0, 1, 1, 2, 1, 1, 2, 1, 1, 2],
           "n_mu": [1, 1, 2, 2, 1, 2, 2, 1, 2, 2],
           "n_cov": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
           "N_var": [5, 5, 4, 4, 5, 4, 4, 5, 4, 4],
           "idx": [0, 5, 10, 14, 18, 23, 27, 31, 36, 40, 44],
           "parameter_sensor": [1, 1, 1, 2, 2, 2],
           "sensor_model_constant_sensor": [],
           "sensor_model_contant": None}
    a = array([1., 1.3, 0.002, 0.5, 1.1, 0.0005])

    ####################################################################################################################
    # 2. Initialise MatchUp object
    ####################################################################################################################

    MatchUpTest = MatchUp()
    MatchUpTest.values = values
    MatchUpTest.unc = unc
    MatchUpTest.ks = ks
    MatchUpTest.unck = unck
    MatchUpTest.idx = idx
    MatchUpTest.a = a
    MatchUpTest.sensor_model = [test_reference_function, test_sensor_function, test_sensor_function]
    MatchUpTest.adjustment_model = [test_adjustment_model, test_adjustment_model, test_adjustment_model]

    return MatchUpTest


class TestMatchupMaths(unittest.TestCase):
    def test_evaluate_measurand_multi_mu0(self):
        """
        Test for match_up.evaluate_measurand() function for an example match-up dataset with multiple match-up series
        - computing measurand for all match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed measurand

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        measurand_expected = array([[470.5, 467.6731159],
                                    [720.56, 690.7705827],
                                    [450.9, 496.3530452],
                                    [295.6, 283.282621],
                                    [315.23, 315.2345789],
                                    [458.4667945, 507.721508],
                                    [362.4185745, 444.5854022],
                                    [470.6731247, 679.3291801],
                                    [269.0340877, 268.5784796]])

        ################################################################################################################
        # 3. Run MatchUpData.read_daUncertaintyta()
        ################################################################################################################

        measurand_test = evaluate_measurand(MatchUpTest)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for row_expected, row_test in zip(measurand_expected, measurand_test):
            for measurand_expected_i, measurand_test_i in zip(row_expected, row_test):
                self.assertAlmostEqual(measurand_expected_i, measurand_test_i, places=3)

    def test_evaluate_measurand_multi_mu0_unc(self):
        """
        Test for match_up.evaluate_measurand() function for an example match-up dataset with multiple match-up series
        - computing measurand for all match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed measurand

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()
        parameter_covariance_matrix = eye(len(MatchUpTest.a))

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        measurand_expected = array([[470.5, 467.6731159],
                                    [720.56, 690.7705827],
                                    [450.9, 496.3530452],
                                    [295.6, 283.282621],
                                    [315.23, 315.2345789],
                                    [458.4667945, 507.721508],
                                    [362.4185745, 444.5854022],
                                    [470.6731247, 679.3291801],
                                    [269.0340877, 268.5784796]])
        uncertainty_expected = [[0., 4982.65185547],
                                [0., 5011.71582031],
                                [0., 4956.17919922],
                                [0., 5002.87597656],
                                [0., 4975.75878906],
                                [5123.81298828, 6496.42626953],
                                [5133.67138672, 6508.70947266],
                                [5095.96386719, 6477.2754]]

        ################################################################################################################
        # 3. Run MatchUpData.read_daUncertaintyta()
        ################################################################################################################

        measurand_test, uncertainty_test = evaluate_measurand(MatchUpTest, parameter_covariance_matrix=parameter_covariance_matrix)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for row_expected, row_test in zip(measurand_expected, measurand_test):
            for measurand_expected_i, measurand_test_i in zip(row_expected, row_test):
                self.assertAlmostEqual(measurand_expected_i, measurand_test_i, places=3)

        for row_expected, row_test in zip(uncertainty_expected, uncertainty_test):
            for measurand_expected_i, measurand_test_i in zip(row_expected, row_test):
                self.assertAlmostEqual(measurand_expected_i, measurand_test_i, places=3)


    def test_evaluate_measurand_multi_muint(self):
        """
        Test for match_up.evaluate_measurand() function for an example match-up dataset with multiple match-up series
        - computing measurand for 1 selected match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed measurand

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        measurand_expected = array([[458.4667945, 507.721508],
                                    [362.4185745, 444.5854022],
                                    [470.6731247, 679.3291801],
                                    [269.0340877, 268.5784796]])

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        measurand_test = evaluate_measurand(MatchUpTest, matchup_series=2)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for row_expected, row_test in zip(measurand_expected, measurand_test):
            for measurand_expected_i, measurand_test_i in zip(row_expected, row_test):
                self.assertAlmostEqual(measurand_expected_i, measurand_test_i, places=3)

    def test_evaluate_measurand_multi_mulst(self):
        """
        Test for match_up.evaluate_measurand() function for an example match-up dataset with multiple match-up series
        - computing measurand for a selection of a list of match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed measurand

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        measurand_expected = array([[470.5, 467.6731159],
                                    [720.56, 690.7705827],
                                    [450.9, 496.3530452],
                                    [295.6, 283.282621],
                                    [315.23, 315.2345789],
                                    [458.4667945, 507.721508],
                                    [362.4185745, 444.5854022],
                                    [470.6731247, 679.3291801],
                                    [269.0340877, 268.5784796]])

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        measurand_test = evaluate_measurand(MatchUpTest, matchup_series=[1, 2])

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for row_expected, row_test in zip(measurand_expected, measurand_test):
            for measurand_expected_i, measurand_test_i in zip(row_expected, row_test):
                self.assertAlmostEqual(measurand_expected_i, measurand_test_i, places=3)

    def test_evaluate_adjusted_measurand_multi_mu0(self):
        """
        Test for match_up.evaluate_adjusted_measurand() function for an example match-up dataset with multiple match-up
        series
        - computing adjusted measurand for all match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed adjusted measurand

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        adjusted_measurand_expected = array([[941, 935.3462318],
                                             [1441.12, 1381.541165],
                                             [901.8, 992.7060905],
                                             [591.2, 566.5652421],
                                             [630.46, 630.4691579],
                                             [916.933589, 1015.443016],
                                             [724.8371491, 889.1708044],
                                             [941.3462494, 1358.65836],
                                             [538.0681754, 537.1569593]])

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        adjusted_measurand_test = evaluate_adjusted_measurand(MatchUpTest)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for row_expected, row_test in zip(adjusted_measurand_expected, adjusted_measurand_test):
            for adjusted_measurand_expected_i, adjusted_measurand_test_i in zip(row_expected, row_test):
                self.assertAlmostEqual(adjusted_measurand_expected_i, adjusted_measurand_test_i, places=3)

    def test_evaluate_adjusted_measurand_multi_muint(self):
        """
        Test for match_up.evaluate_adjusted_measurand() function for an example match-up dataset with multiple match-up
        series
        - computing measurand for 1 selected match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed measurand

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        adjusted_measurand_expected = array([[916.933589, 1015.443016],
                                             [724.8371491, 889.1708044],
                                             [941.3462494, 1358.65836],
                                             [538.0681754, 537.1569593]])

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        adjusted_measurand_test = evaluate_adjusted_measurand(MatchUpTest, matchup_series=2)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for row_expected, row_test in zip(adjusted_measurand_expected, adjusted_measurand_test):
            for adjusted_measurand_expected_i, adjusted_measurand_test_i in zip(row_expected, row_test):
                self.assertAlmostEqual(adjusted_measurand_expected_i, adjusted_measurand_test_i, places=3)

    def test_evaluate_adjusted_measurand_multi_mulst(self):
        """
        Test for match_up.evaluate_adjusted_measurand() function for an example match-up dataset with multiple match-up
        series
        - computing measurand for a selection of a list of match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed measurand

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        adjusted_measurand_expected = array([[941, 935.3462318],
                                             [1441.12, 1381.541165],
                                             [901.8, 992.7060905],
                                             [591.2, 566.5652421],
                                             [630.46, 630.4691579],
                                             [916.933589, 1015.443016],
                                             [724.8371491, 889.1708044],
                                             [941.3462494, 1358.65836],
                                             [538.0681754, 537.1569593]])

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        adjusted_measurand_test = evaluate_adjusted_measurand(MatchUpTest, matchup_series=[1,2])

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for row_expected, row_test in zip(adjusted_measurand_expected, adjusted_measurand_test):
            for adjusted_measurand_expected_i, adjusted_measurand_test_i in zip(row_expected, row_test):
                self.assertAlmostEqual(adjusted_measurand_expected_i, adjusted_measurand_test_i, places=3)

    def test_evaluate_K_multi_mu0(self):
        """
        Test for match_up.evaluate_K() function for an example match-up dataset with multiple match-up
        series
        - computing K for all match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed K

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        K_expected = array([-5.653768212, -59.57883451, 90.9060905, -24.63475795, 0.009157895,
                            98.50942692, 164.3336554, 417.3121109, -0.911216111])

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        K_test = evaluate_K(MatchUpTest)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for K_expected_i, K_test_i in zip(K_expected, K_test):
            self.assertAlmostEqual(K_expected_i, K_test_i, places=3)

    def test_evaluate_K_multi_muint(self):
        """
        Test for match_up.evaluate_K() function for an example match-up dataset with multiple match-up
        series
        - computing K for 1 selected match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed K

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        K_expected = array([98.50942692, 164.3336554, 417.3121109, -0.911216111])

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        K_test = evaluate_K(MatchUpTest, matchup_series=2)

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for K_expected_i, K_test_i in zip(K_expected, K_test):
            self.assertAlmostEqual(K_expected_i, K_test_i, places=3)

    def test_evaluate_K_multi_mulst(self):
        """
        Test for match_up.evaluate_K() function for an example match-up dataset with multiple match-up
        series
        - computing K for a selection of a list of match-up series
        - one sensor function
        """

        # Test Description
        # ================
        #
        # 1. This test intialises an example *eopy.matchup.matchupIO.MatchUp* object
        #
        # 2. Compares computed K

        ################################################################################################################
        # 1. Initialise Test Data Object
        ################################################################################################################

        MatchUpTest = return_MatchUpTest()

        ################################################################################################################
        # 2. Define expected values
        ################################################################################################################

        K_expected = array([-5.653768212, -59.57883451, 90.9060905, -24.63475795, 0.009157895,
                            98.50942692, 164.3336554, 417.3121109, -0.911216111])

        ################################################################################################################
        # 3. Run MatchUpData.read_data()
        ################################################################################################################

        K_test = evaluate_K(MatchUpTest, matchup_series=[1, 2])

        ################################################################################################################
        # 4. Compare retrieve values to expect values
        ################################################################################################################

        for K_expected_i, K_test_i in zip(K_expected, K_test):
            self.assertAlmostEqual(K_expected_i, K_test_i, places=3)

if __name__ == '__main__':
    unittest.main()
