"""
Generate a new dataset with evaluated K

Usage:
python simulated_AVHRR.py <input dataset> <sensor data> <output dataset>
"""

'''___Python Modules___'''
import sys
from sys import argv
from distutils.dir_util import copy_tree

'''___Third Party Modules___'''
from netCDF4 import Dataset
from os import listdir
from os.path import join as pjoin

'''___Harmonisation Modules___'''
core_directory = "/home/seh2/src/FIDUCEO/WP2/Harmonisation/src/main/full_EIV/src/main"
sys.path.append(core_directory)
from nplcore import MatchUp
from nplcore import evaluate_K, evaluate_measurand

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "14/03/2018"
__credits__ = ["Ralf Quast", "Jon Mittaz", "Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def main(dataset_in_directory, sensor_model_path, dataset_out_directory):
    dataset_in_paths = [pjoin(dataset_in_directory, fname) for fname in listdir(dataset_in_directory)]
    print "Files:"
    for dataset_in_path in dataset_in_paths:
        print ">", dataset_in_path

    print "\nOpening data..."
    dataset_in = MatchUp(dataset_in_directory, sensor_model_path, open_uncertainty=False)
    print "Done"

    print "\nEvaluating 'True' K..."
    k_true = evaluate_K(dataset_in)
    print "Done"

    print "\nGenerating 'True' K dataset..."
    copy_tree(dataset_in_directory, dataset_out_directory)
    dataset_out_paths = [pjoin(dataset_in_directory, fname) for fname in listdir(dataset_in_directory)]
    for i_mu, fname in enumerate(dataset_out_paths):
        istart = dataset_in.idx['cNm'][i_mu]
        iend = dataset_in.idx['cNm'][i_mu+1]

        dataset_out_i = Dataset(fname, 'a')
        dataset_out_i.variables['K'][:] = k_true[istart:iend]
        dataset_out_i.close()
    print "Done"

    return 0

if __name__ == "__main__":
    main(argv[1], argv[2], argv[3])
