"""
Generate the required W matrices for a given MW match-up file and add them to the file

Usage:
python MW_add_W_matrices.py path/to/matchup/file_in.nc path/to/matchup/file_out.nc
"""

'''___Python Modules___'''
from sys import argv
from shutil import copyfile

'''___Third Party Modules___'''
from netCDF4 import Dataset
from numpy import array, int8

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

        :uR_PRT1: *float*

        Random uncertainty for single PRT measurements for sensor 1

        :uR_PRT2: *float*

        Random uncertainty for single PRT measurements for sensor 2

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

    uR_PRT1 = None
    if sensor_1_name != -1:
        uR_PRT1 = float(dataset.variables['Ur1'][1, 4]*2)
    uR_PRT2 = float(dataset.variables['Ur2'][1, 4]*2)

    scanline_time = 8.0/3.0

    dataset.close()

    weights = array([0.0625, 0.125, 0.1875, 0.25, 0.1875, 0.125, 0.0625])

    return sensor_1_name, sensor_2_name, time1, time2, uR_PRT1, uR_PRT2, scanline_time, weights


def return_w_matrix_variables_MW(time1, time2, scanline_time, uR_PRT1, uR_PRT2, sensor_1_name, weights):
    """

    :type time1: numpy.ndarray
    :param time1: Matchup time for sensor 1

    :type time2: numpy.ndarray
    :param time2: Matchup time for sensor 2

    :type scanline_time: float
    :param scanline_time: time interval between consecutive scanlines

    :type uR_PRT1: float
    :param uR_PRT1: Random uncertainty for single PRT measurements for sensor 1

    :type uR_PRT1: float
    :param uR_PRT1: Random uncertainty for single PRT measurements for sensor 1

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
    w_matrix2, u_matrix2 = generate_rolling_average_w_matrix(time2, scanline_time, uR_PRT2, weights=weights)
    w_matrices = [w_matrix2]
    u_matrices = [u_matrix2]

    # b. Build w matrix and uncertainty vectors for Sensor 1 (if not reference)
    if sensor_1_name != -1:
        # a. Compute W data
        w_matrix1, u_matrix1 = generate_rolling_average_w_matrix(time1, scanline_time, uR_PRT1, weights=weights)
        w_matrices = [w_matrix1, w_matrix2]
        u_matrices = [u_matrix1, u_matrix2]

    # 2. Convert W matrices and uncertainty vector to input file variables
    w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz, \
    u_matrix_row_count, u_matrix_val = return_w_matrix_variables(w_matrices, u_matrices)

    # 3. Assign type/use matrices
    w_matrix_use1 = array([0, 0, 1, 1, 1, 1, 1, 0, 0], dtype=int8)
    w_matrix_use2 = array([0, 0, 2, 2, 2, 2, 2, 0, 0], dtype=int8)
    u_matrix_use1 = array([0, 0, 1, 1, 1, 1, 1, 0, 0], dtype=int8)
    u_matrix_use2 = array([0, 0, 2, 2, 2, 2, 2, 0, 0], dtype=int8)

    if sensor_1_name == -1:
        w_matrix_use1 = array([0], dtype=int8)
        w_matrix_use2 = array([0, 0, 2, 2, 2, 2, 2, 0, 0], dtype=int8)
        u_matrix_use1 = array([0], dtype=int8)
        u_matrix_use2 = array([0, 0, 2, 2, 2, 2, 2, 0, 0], dtype=int8)

    return w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz,\
            u_matrix_row_count, u_matrix_val, w_matrix_use1, w_matrix_use2, u_matrix_use1, u_matrix_use2


def main(input_file_path, output_file_path):
    """
    Routine to update MW matchup file to include w variables and match newest spec

    :type input_file_path: str
    :param input_file_path: path of input harmonisation input file

    :type output_file_path: str
    :param output_file_path: path of output harmonisation input file
    """

    # 1. Read required data to generate W matrices
    print "Reading file:", input_file_path
    sensor_1_name, sensor_2_name, time1, time2, uR_PRT1, uR_PRT2, scanline_time, weights = read_matchup_data(input_file_path)

    # 2. Generate required W matrix variables
    print "Generating W Matrix Variables from data..."

    w_matrix_val, w_matrix_row, w_matrix_col,\
        w_matrix_nnz, u_matrix_row_count, u_matrix_val, \
            w_matrix_use1, w_matrix_use2, u_matrix_use1, \
                u_matrix_use2 = return_w_matrix_variables_MW(time1, time2, scanline_time,
                                                             uR_PRT1, uR_PRT2, sensor_1_name, weights)

    # 3. Update file to include W matrix variables
    print "Writing to new file:", output_file_path

    # a. Copy old file to new location
    copyfile(input_file_path, output_file_path)

    # b.  Amend existing variables
    dataset = Dataset(output_file_path, 'a')

    # + uncertainty_type1/2
    uncertainty_type1 = array([1, 1, 4, 4, 4, 4, 4, 1, 2], dtype=int8)
    uncertainty_type2 = array([1, 1, 4, 4, 4, 4, 4, 1, 2], dtype=int8)
    if sensor_1_name == -1:
        uncertainty_type1 = array([1], dtype=int8)
        uncertainty_type2 = array([1, 1, 4, 4, 4, 4, 4, 1, 2], dtype=int8)

    dataset.variables['uncertainty_type1'][:] = uncertainty_type1
    dataset.variables['uncertainty_type2'][:] = uncertainty_type2

    # + Ur1/2
    Ur2 = dataset.variables['Ur2'][:]
    Ur2[:, 2:7] = 0
    dataset.variables['Ur2'][:] = Ur2

    if sensor_1_name != -1:
        Ur1 = dataset.variables['Ur1'][:]
        Ur1[:, 2:7] = 0
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


if __name__ == "__main__":
    main(argv[1], argv[2])

