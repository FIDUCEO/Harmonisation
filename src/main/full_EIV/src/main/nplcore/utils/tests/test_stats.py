"""
Test of stats module functions
"""

'''___Built-In Modules___'''
import unittest
from os.path import dirname
from os.path import join as pjoin
import sys
from datetime import datetime as dt

'''___Third-Party Modules___'''
from numpy import array, nan, isnan

'''___NPL Modules___'''
util_directory = dirname(dirname(__file__))
sys.path.append(util_directory)
from stats import time_average

'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "13/12/2017"
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class TestStats(unittest.TestCase):
    def test_time_average_monthly(self):

        # 1. Define input
        data = array([1.0, 2.0, 3.0, 2.0, 4.0, 1.0])
        times = array([dt(2000, 1, 3),
                       dt(2000, 1, 10),
                       dt(2000, 1, 15),
                       dt(2000, 2, 14),
                       dt(2000, 2, 16),
                       dt(2000, 4, 19)]).astype(dt)

        # 2. Define expected output
        averages_expected = array([2.0, 3.0, nan, 1.0])
        average_times_expected = array([dt(2000, 1, 1),
                                        dt(2000, 2, 1),
                                        dt(2000, 3, 1)]).astype(dt)

        # 3. Run averaging
        averages_test, average_times_test = time_average(data, times, time_period="month")

        # 4. Test output
        self.assertTrue(((averages_expected == averages_test) | (isnan(averages_expected) & isnan(averages_test))).all())

        for time_test, time_expected in zip(average_times_test, average_times_expected):
            self.assertEquals(time_test, time_expected)

    def test_time_average_monthly_unordered(self):

        # 1. Define input
        data = array([1.0, 2.0, 3.0, 2.0, 4.0, 1.0])
        times = array([dt(2000, 4, 3),
                       dt(2000, 2, 10),
                       dt(2000, 1, 15),
                       dt(2000, 1, 14),
                       dt(2000, 2, 16),
                       dt(2000, 1, 19)]).astype(dt)

        # 2. Define expected output
        averages_expected = array([2.0, 3.0, nan, 1.0])
        average_times_expected = array([dt(2000, 1, 1),
                                        dt(2000, 2, 1),
                                        dt(2000, 3, 1)]).astype(dt)

        # 3. Run averaging
        averages_test, average_times_test = time_average(data, times, time_period="month")

        # 4. Test output
        self.assertTrue(((averages_expected == averages_test) | (isnan(averages_expected) & isnan(averages_test))).all())

        for time_test, time_expected in zip(average_times_test, average_times_expected):
            self.assertEquals(time_test, time_expected)


if __name__ == '__main__':
    unittest.main()
