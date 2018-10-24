"""
MW sensor model test
"""

'''___Built-In Modules___'''
import unittest
import sys
from os.path import dirname

'''___Third-Party Modules___'''
from numpy import array

'''___harmonisation Modules___'''
sys.path.append(dirname(dirname(__file__)))

'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "18/02/2018"
__credits__ = [""]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class TestMWSensorModels(unittest.TestCase):
    def test_MWSensorModels3(self):
        from MW____RCh3_3_sensor_data import sensor_data

        # Open test
        sensor_A_data = sensor_data['A']
        sensor_model = sensor_A_data['sensor_model']
        parameter = array([-15.5, 0.001, 0.71])
        constant = array(sensor_A_data['sensor_model_constant'])
        constant[:90] = 0.98
        index = array([1])

        # Define input test vector
        X = array([[10787.171875, 48380.703125, 283.285375, 283.2660625, 283.2778125,
                    283.28325, 283.38725, 41701, 331.5454321062],
                   [10787.171875, 48380.703125, 283.285375, 283.2660625, 283.2778125,
                    283.28325, 283.38725, 41701, 331.5454321062]])

        # Define expected output
        BT_expected = array([2.34186e2, 2.34186e2])
        J_expected = array([[1.20639e-2, 6.20918e-3, 5.59355e-1, 2.79677e-1, 2.79677e-1, 2.79677e-1,
                             2.79677e-1, 1.82988e2, 0, -1.19117e-13,  17.3055, 1.75557e16],
                            [1.20639e-2, 6.20918e-3, 5.59355e-1, 2.79677e-1, 2.79677e-1, 2.79677e-1,
                             2.79677e-1, 1.82988e2, 0, -1.19117e-13, 17.3055, 1.75557e16]])

        # Evaluate test vector
        BT, J = sensor_model(X, parameter, constant, index, None, None)

        # Test results
        for element, element_expected in zip(BT, BT_expected):
            self.assertAlmostEquals(element, element_expected, places=3)

        J[:,-1] /= 1e16
        J_expected[:, -1] /= 1e16
        for row, row_expected in zip(J, J_expected):
            for element, element_expected in zip(row, row_expected):
                self.assertAlmostEquals(element, element_expected, places=3)

        self.assertEquals(True, True)


if __name__ == "__main__":
    unittest.main()