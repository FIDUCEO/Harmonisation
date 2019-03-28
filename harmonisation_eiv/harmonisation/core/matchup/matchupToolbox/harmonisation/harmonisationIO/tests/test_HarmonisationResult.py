"""
Test of MatchUpData Class
"""

'''___Built-In Modules___'''
import unittest
from os.path import join as pjoin
from os import makedirs, getcwd
from shutil import rmtree

'''___Third-Party Modules___'''
from netCDF4 import Dataset

'''___harmonisation Modules___'''
from harmonisation import HarmonisationResult
from harmonisation.version import __version__


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "5/11/2018"
__version__ = __version__
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

# Constant
temp_data_directory = pjoin(getcwd(), "temp")


def write_ncfile(path):
    rootgrp = Dataset(path, 'w')
    rootgrp.attr = "test"
    rootgrp.close()
    return 0


def setup_read(output_dir):
    makedirs(pjoin(output_dir, "plots"))

    path_result = pjoin(output_dir, "harm_test.nc")
    path_res1 = pjoin(output_dir, "harm_test_res_0_1.nc")
    path_res2 = pjoin(output_dir, "harm_test_res_1_2.nc")
    path_txt = pjoin(output_dir, "test.cfg")

    write_ncfile(path_result)
    write_ncfile(path_res1)
    write_ncfile(path_res2)

    with open(path_txt, "w") as f:
        f.write("test")

    harm_res = HarmonisationResult()

    return harm_res, path_result, [path_res1, path_res2]


def teardown_read(output_dir):
    rmtree(output_dir)
    return 0


class TestHarmonisationResult(unittest.TestCase):
    def test__get_harm_output_path(self):
        harm_res, path_result, paths_res = setup_read(temp_data_directory)

        test_path_result = harm_res._get_harm_output_path(temp_data_directory)

        self.assertEqual(test_path_result, path_result)

        teardown_read(temp_data_directory)

    def test__get_harm_res_paths(self):
        harm_res, path_result, paths_res = setup_read(temp_data_directory)

        test_paths_res = harm_res._get_harm_res_paths(temp_data_directory)

        self.assertItemsEqual(test_paths_res, paths_res)

        teardown_read(temp_data_directory)


if __name__ == '__main__':
    unittest.main()
