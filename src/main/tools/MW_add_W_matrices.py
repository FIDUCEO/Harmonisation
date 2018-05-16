"""
Generate the required W matrices for a given MW match-up file and add them to the file

Usage:
python MW_add_W_matrices.py path/to/matchup/file_in.nc path/to/matchup/file_out.nc
"""

'''___Python Modules___'''
from sys import argv
from shutil import copyfile
from os import listdir
from os.path import join as pjoin

'''___Third Party Modules___'''
from netCDF4 import Dataset
from numpy import array, int8, ones

'''___NPL Modules___'''
from generate_w_matrices import generate_rolling_average_w_matrix
from W_matrix_functions import return_w_matrix_variables, append_W_to_input_file
from harmonisation_input_checker import main as test_input_file

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "04/12/2017"
__credits__ = ["Martin Burgdorf"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

AMSUB_sensors = ["15"]
MHS_sensors = ["19", "A", "B"]
ref_sensors = ["18"]


def read_matchup_data(file_path):
    """
    Return required data from MW match-up file

    :type file_path: str
    :param file_path: path of match-up file

    :return:
        :sensor_1_name: *int*

        sensor 1 ID

        :sensor_2_name: *int*

        sensor 2 ID

        :time1: *numpy.ndarray*

        Match-up times for sensor 1

        :time2: *numpy.ndarray*

        Match-up times for sensor 2

        :uR1: *float*

        Random uncertainty for X1, first row

        :uR2: *float*

        Random uncertainty for X2, first row

        :scanline_time: *float*

        Time interval between consecutive scanlines

        :weights: *numpy.ndarray*

        Averaging kernel weights for rolling average operator
    """

    dataset = Dataset(file_path, "r")

    sensor_1_name = dataset.sensor_1_name
    sensor_2_name = dataset.sensor_2_name

    time1 = dataset.variables['time1'][:]
    time2 = dataset.variables['time1'][:]

    uR1 = dataset.variables['Ur1'][0, :]
    uR2 = dataset.variables['Ur2'][0, :]

    scanline_time = 8.0/3.0

    dataset.close()

    weights = array([0.0625, 0.125, 0.1875, 0.25, 0.1875, 0.125, 0.0625])

    return sensor_1_name, sensor_2_name, time1, time2, uR1, uR2, scanline_time, weights


def return_w_matrix_variables_MW(time1, time2, scanline_time, uR1, uR2, weights, sensor_1_instrument, sensor_2_instrument):
    """

    :type time1: numpy.ndarray
    :param time1: Matchup time for sensor 1

    :type time2: numpy.ndarray
    :param time2: Matchup time for sensor 2

    :type scanline_time: float
    :param scanline_time: time interval between consecutive scanlines

    :type uR1: float
    :param uR1: Random uncertainty for X1, first row

    :type uR2: float
    :param uR2: Random uncertainty for X1, first row

    :type sensor_1_name: str
    :param sensor_1_name: Sensor ID for match-up sensor 1

    :type weights: numpy.ndarray
    :param weights:

    :return:
        :w_matrix_val: *numpy.ndarray*

        Concatenated array of non-zero elements of W matrices

        :w_matrix_row: *numpy.ndarray*

        Concatenated array of row indices of W matrices

        :w_matrix_col: *numpy.ndarray*

        Concatenated array of row indices of W matrices

        :w_matrix_nnz: *numpy.ndarray*

        Number of non-zeros elements per W matrix

        :u_matrix_row_count: *numpy.ndarray*

        Number of non-zero elements per u matrix

        :u_matrix_val: *numpy.ndarray*

        Concatenated array of u matrix non-zero values

        :w_matrix_use1: *numpy.ndarray*

        Mapping from X1 array column index to W

        :w_matrix_use2: *numpy.ndarray*

        Mapping from X2 array column index to W

        :u_matrix_use1: *numpy.ndarray*

        Mapping from X1 array column index to U

        :u_matrix_use2: *numpy.ndarray*

        Mapping from X2 array column index to U
    """

    # 1. Generate lists of W matrices and uncertainty vectors

    # a. Build w and u matrix vectors for Sensor 2
    w_matrix2, u_matrix20 = generate_rolling_average_w_matrix(time2, scanline_time, float(uR2[0]), weights=weights)
    w_matrices = [w_matrix2]
    u_matrix21 = ones(len(u_matrix20)) * uR2[1]
    u_matrix22 = ones(len(u_matrix20)) * uR2[2]

    u_matrices = [u_matrix20, u_matrix21, u_matrix22]

    # b. Build w matrix and uncertainty vectors for Sensor 1 (if not reference)
    if len(uR1) != 1:
        # a. Compute W data
        w_matrix1, u_matrix10 = generate_rolling_average_w_matrix(time1, scanline_time, float(uR1[0]), weights=weights)
        w_matrices = [w_matrix1, w_matrix2]
        u_matrix11 = ones(len(u_matrix20)) * uR1[1]
        u_matrix12 = ones(len(u_matrix20)) * uR1[2]

        u_matrices = [u_matrix10, u_matrix11, u_matrix12, u_matrix20, u_matrix21, u_matrix22]

    # 2. Convert W matrices and uncertainty vector to input file variables
    w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz, \
    u_matrix_row_count, u_matrix_val = return_w_matrix_variables(w_matrices, u_matrices)

    # 3. Assign type/use matrices
    if (sensor_1_instrument == "ref") and (sensor_2_instrument == "AMSUB"):
        w_matrix_use1 = array([0], dtype=int8)
        w_matrix_use2 = array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0], dtype=int8)
        u_matrix_use1 = array([0], dtype=int8)
        u_matrix_use2 = array([1, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0], dtype=int8)

    if (sensor_1_instrument == "ref") and (sensor_2_instrument == "MHS"):
        w_matrix_use1 = array([0], dtype=int8)
        w_matrix_use2 = array([1, 1, 1, 1, 1, 1, 1, 0, 0], dtype=int8)
        u_matrix_use1 = array([0], dtype=int8)
        u_matrix_use2 = array([1, 2, 3, 3, 3, 3, 3, 0, 0], dtype=int8)

    if (sensor_1_instrument == "MHS") and (sensor_2_instrument == "MHS"):
        w_matrix_use1 = array([1, 1, 1, 1, 1, 1, 1, 0, 0], dtype=int8)
        w_matrix_use2 = array([2, 2, 2, 2, 2, 2, 2, 0, 0], dtype=int8)
        u_matrix_use1 = array([1, 2, 3, 3, 3, 3, 3, 0, 0], dtype=int8)
        u_matrix_use2 = array([4, 5, 6, 6, 6, 6, 6, 0, 0], dtype=int8)

    if (sensor_1_instrument == "AMSUB") and (sensor_2_instrument == "MHS"):
        w_matrix_use1 = array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0], dtype=int8)
        w_matrix_use2 = array([2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0], dtype=int8)
        u_matrix_use1 = array([1, 2, 3, 3, 3, 3, 3, 0, 0], dtype=int8)
        u_matrix_use2 = array([4, 5, 6, 6, 6, 6, 6, 0, 0], dtype=int8)

    if (sensor_1_instrument == "MHS") and (sensor_2_instrument == "AMSUB"):
        w_matrix_use1 = array([1, 1, 1, 1, 1, 1, 1, 0, 0], dtype=int8)
        w_matrix_use2 = array([2, 2, 2, 2, 2, 2, 2, 0, 0], dtype=int8)
        u_matrix_use1 = array([1, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0], dtype=int8)
        u_matrix_use2 = array([4, 5, 6, 6, 6, 6, 6, 6, 6, 0, 0], dtype=int8)

    if (sensor_1_instrument == "AMSUB") and (sensor_2_instrument == "AMSUB"):
        w_matrix_use1 = array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0], dtype=int8)
        w_matrix_use2 = array([2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0], dtype=int8)
        u_matrix_use1 = array([1, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0], dtype=int8)
        u_matrix_use2 = array([4, 5, 6, 6, 6, 6, 6, 6, 6, 0, 0], dtype=int8)

    return w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz,\
            u_matrix_row_count, u_matrix_val, w_matrix_use1, w_matrix_use2, u_matrix_use1, u_matrix_use2


def run(input_file_path, output_file_path):
    """
    Routine to update MW matchup file to include w variables and match newest spec

    :type input_file_path: str
    :param input_file_path: path of input harmonisation input file

    :type output_file_path: str
    :param output_file_path: path of output harmonisation input file
    """

    # 1. Read required data to generate W matrices
    print "Reading file:", input_file_path
    sensor_1_name, sensor_2_name, time1, time2, uR1, uR2, scanline_time, weights = read_matchup_data(input_file_path)

    if sensor_2_name in MHS_sensors:
        sensor_2_instrument = "MHS"
    elif sensor_2_name in AMSUB_sensors:
        sensor_2_instrument = "AMSUB"
    else:
        raise ValueError("Sensor 2 Unknown - "+sensor_2_name)

    if sensor_1_name in MHS_sensors:
        sensor_1_instrument = "MHS"
    elif sensor_1_name in AMSUB_sensors:
        sensor_1_instrument = "AMSUB"
    elif sensor_1_name in ref_sensors:
        sensor_1_instrument = "ref"
    else:
        raise ValueError("Sensor 1 Unknown - "+sensor_1_name)

    # 2. Generate required W matrix variables
    print "Generating W Matrix Variables from data..."

    w_matrix_val, w_matrix_row, w_matrix_col,\
        w_matrix_nnz, u_matrix_row_count, u_matrix_val, \
            w_matrix_use1, w_matrix_use2, u_matrix_use1, \
                u_matrix_use2 = return_w_matrix_variables_MW(time1, time2, scanline_time, uR1, uR2, weights,
                                                             sensor_1_instrument, sensor_2_instrument)

    # 3. Update file to include W matrix variables
    print "Writing to new file:", output_file_path

    # a. Copy old file to new location
    copyfile(input_file_path, output_file_path)

    # b.  Amend existing variables
    dataset = Dataset(output_file_path, 'a')

    # + uncertainty_type1/2
    if sensor_1_instrument == "ref":
        uncertainty_type1 = array([1], dtype=int8)
    elif sensor_1_instrument == "AMSUB":
        uncertainty_type1 = array([3, 3, 4, 4, 4, 4, 4, 4, 4, 1, 2], dtype=int8)
    elif sensor_1_instrument == "MHS":
        uncertainty_type1 = array([3, 3, 4, 4, 4, 4, 4, 1, 2], dtype=int8)
    if sensor_2_instrument == "AMSUB":
        uncertainty_type2 = array([3, 3, 4, 4, 4, 4, 4, 4, 4, 1, 2], dtype=int8)
    elif sensor_2_instrument == "MHS":
        uncertainty_type2 = array([3, 3, 4, 4, 4, 4, 4, 1, 2], dtype=int8)

    dataset.variables['uncertainty_type1'][:] = uncertainty_type1
    dataset.variables['uncertainty_type2'][:] = uncertainty_type2

    # + Ur1/2
    Ur2 = dataset.variables['Ur2'][:]
    if sensor_2_instrument == "AMSUB":
        Ur2[:, 0:9] = 0
    elif sensor_2_instrument == "MHS":
        Ur2[:, 0:7] = 0
    dataset.variables['Ur2'][:] = Ur2

    Ur1 = dataset.variables['Ur1'][:]
    if sensor_1_instrument == "AMSUB":
        Ur1[:, :9] = 0
    elif sensor_1_instrument == "MHS":
        Ur1[:, :7] = 0
    dataset.variables['Ur1'][:] = Ur1

    dataset.close()

    # c. Add w matrix variables to updated file
    append_W_to_input_file(output_file_path,
                           w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz,
                           u_matrix_row_count, u_matrix_val,
                           w_matrix_use1, w_matrix_use2, u_matrix_use1, u_matrix_use2)

    # 4. Test harmonisation
    test_input_file(output_file_path)

    return 0


def main(input_directory, output_directory):

    for fname in listdir(input_directory):
        run(pjoin(input_directory, fname), pjoin(output_directory, fname))

    return 0


if __name__ == "__main__":
    # main(argv[1], argv[2])
    main("/home/seh2/data/FIDUCEO/matchups/MW/MW____RCh3_1_RS______",
         "/home/seh2/data/FIDUCEO/matchups/MW/MW____RCh3_1_RSAS____")

