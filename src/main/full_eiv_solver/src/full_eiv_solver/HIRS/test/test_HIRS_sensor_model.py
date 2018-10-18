"""
HIRS sensor model test
"""

'''___Built-In Modules___'''
import unittest
import sys
from os.path import dirname

'''___Third-Party Modules___'''
from numpy import array

'''___NPL Modules___'''
sys.path.append(dirname(dirname(__file__)))

'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "20/02/2018"
__credits__ = [""]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class TestHIRSSensorModels(unittest.TestCase):
    def test_HIRSSensorModels3(self):
        from HIRS__RCh2_3_sensor_data import sensor_data

        # Open test
        sensor_A_data = sensor_data['ma']
        sensor_model = sensor_A_data['sensor_model']
        parameter = array([1.00e-15, 0.005, 2.00e-18])
        constant = sensor_A_data['sensor_model_constant']

        # Define input test vector
        X = array([[-2.1e1, -9.89e2, -4.42e2, 2.8338e2, -1.797204e-14],
                   [-2.1e1, -9.89e2, -4.42e2, 2.8338e2, -1.797204e-14]])

        # Define expected output
        L_expected = array([1.3087887633843708e-12, 1.3087887633843708e-12])
        J_expected = array([[1.25544266e-15,   2.65025477e-15,  -3.90569742e-15,   4.40301219e-14,
                             -1.0,   1.0,   1.77704631e-12,  -2.30287000e5],
                            [1.25544266e-15, 2.65025477e-15, -3.90569742e-15, 4.40301219e-14,
                             -1.0, 1.0, 1.77704631e-12, -2.30287000e5]])

        # Evaluate test vector
        L, J = sensor_model(X, parameter, constant, None, None, None)

        # Test results
        for element, element_expected in zip(L, L_expected):
            self.assertAlmostEquals(element, element_expected, places=5)

        for row, row_expected in zip(J, J_expected):
            for element, element_expected in zip(row, row_expected):
                self.assertAlmostEquals(element, element_expected, places=5)


if __name__ == "__main__":
    unittest.main()