"""
Generate the required W matrices for a given harmonisation match-up file and add them to the file

Usage:
python add_W_to_matchup_file.py path/to/matchup/file.nc
"""

'''___Python Modules___'''
from sys import argv

'''___Third Party Modules___'''
from numpy import append, array, vstack, int8, trim_zeros, zeros, ones, arange, count_nonzero, sum, asarray, where
from scipy.sparse import csr_matrix
from netCDF4 import Dataset

'''___Harmonisation Modules___'''
from W_matrix_functions import return_w_matrix_variables, write_input_file, append_W_to_input_file
from harmonisation_input_checker import main as test_input_file

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "23/08/2017"
__credits__ = ["Ralf Quast", "Jon Mittaz", "Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

def read_matchup_data(file_path):
    """
    Return time data from harmonisation match-up file

    :type file_path: str
    :param file_path: path of match-up file

    :return:
        :lm: *numpy.ndarray*

        Stores satellite pairs with number of entries

        :H: *numpy.ndarray*

        Radiances and counts per matchup

        :Us: *numpy.ndarray*

        Systematic uncertainties for H array

        :Ur: *numpy.ndarray*

        Random uncertainties for H array

        :K: *numpy.ndarray*

        K (sensor-to-sensor differences) for zero shift case

        :Kr: *numpy.ndarray*

        K (sensor-to-sensor differences) random uncertainties (matchup uncertainty)

        :Ks: *numpy.ndarray*

        K (sensor-to-sensor differences) systematic uncertainties for zero shift case

        :sensor_1_name: *int*

        sensor 1 ID

        :sensor_2_name: *int*

        sensor 2 ID

        :time_sensor_1: *numpy.ndarray*

        Match-up times for sensor 1

        :time_sensor_2: *numpy.ndarray*

        Match-up times for sensor 2

        :width_sensor_1: *float*

        Time per scanline of sensor 1

        :width_sensor_2: *float*

        Time per scanline of sensor 2

        :u_C_S_sensor_1: *numpy.ndarray*

        Scanline uncertainties for space counts for the reference sensor

        :u_C_ICT_sensor_1: *numpy.ndarray*

        Scanline uncertainties for internal calibration target counts for the reference sensor

        u_C_S_sensor_2: *numpy.ndarray*

        Scanline uncertainties for space counts for the sensor

        :u_C_ICT_sensor_2: *numpy.ndarray*

        Scanline uncertainties for internal calibration target counts for the sensor
    """

    rootgrp = Dataset(file_path)
    H = rootgrp.variables["H"][:]
    Us = rootgrp.variables["Us"][:]
    Ur = rootgrp.variables["Ur"][:]
    K = rootgrp.variables["K"][:]
    Kr = rootgrp.variables["Kr"][:]
    Ks = rootgrp.variables["Ks"][:]
    sensor_1_name = rootgrp.variables["lm"][0, 0]
    sensor_2_name = rootgrp.variables["lm"][0, 1]
    time_sensor_1 = rootgrp.variables["ref_time_matchup"][:]
    time_sensor_2 = rootgrp.variables["time_matchup"][:]
    width_sensor_1 = rootgrp.variables["corrData"][:]
    width_sensor_2 = rootgrp.variables["corrData"][:]
    u_C_S_sensor_1 = rootgrp.variables['ref_cal_Sp_Ur'][:, :]
    u_C_ICT_sensor_1 = rootgrp.variables['ref_cal_BB_Ur'][:, :]
    u_C_S_sensor_2 = rootgrp.variables['cal_Sp_Ur'][:, :]
    u_C_ICT_sensor_2 = rootgrp.variables['cal_BB_Ur'][:, :]
    rootgrp.close()

    return H, Us, Ur, K, Kr, Ks, sensor_1_name, sensor_2_name, time_sensor_1, time_sensor_2, \
           width_sensor_1, width_sensor_2, u_C_S_sensor_1, u_C_ICT_sensor_1, u_C_S_sensor_2, u_C_ICT_sensor_2


def return_valid_averages_mask(u11, u12, u21, u22, sensor_1_name):

    valid_averages = ones(u11.shape[0], dtype=bool)
    if sensor_1_name == -1:
        for i, (row21, row22) in enumerate(zip(u21, u22)):
            if (0 in row21) or (0 in row22):
                valid_averages[i] = False

    else:
        for i, (row11, row12, row21, row22) in enumerate(zip(u11, u12, u21, u22)):
            if (0 in row11) or (0 in row12) or (0 in row21) or (0 in row22):
                valid_averages[i] = False

    return valid_averages


def return_W_matrix(times, width, uncertainty):
    """
    Return weighting matrix in sparse representation and uncertainty vector.

    :type times: numpy.ndarray
    :param times: match-up times for match-ups in match-up series

    :type width: numpy.ndarray
    :param width: match-up time data

    :return:
        :W: *scipy.sparse.csr_matrix*
        weighting matrix
        :uncertainty_vector: *numpy.ndarray*
        scanline uncertainty vector
    """
    n_var = len(times)  # number of match_ups
    N_W = len(uncertainty[0])  # length of full averaging kernel

    # 1. Determine first column index of each row
    rel_col = []
    col = 0
    for i in xrange(n_var):
        col_step = 0
        if i > 0:
            diff = times[i] - times[i - 1]
            s = 1
            if diff < 0:
                s = -1
            if abs(diff) > width:
                corr_val = 0.
            else:
                corr_val = 1. - (abs(diff) / width)
            col_step = s * int(round(N_W * (1 - corr_val)))
        col += col_step
        rel_col.append(col)
    rel_col = asarray(rel_col)
    col = rel_col - min(rel_col)

    # 2. Generate sparse matrix indices
    # initialise sparse matrix index and values arrays (of maximum size, i.e. if all windows are N_W)
    ir = zeros(n_var * N_W)
    jc = zeros(n_var * N_W)
    ws = zeros(n_var * N_W)
    uncertainty_vector = zeros(n_var * N_W)

    # Populate arrays
    iend = 0
    for i, col_i in enumerate(col):
        ui = uncertainty[i]  # scanline non-zero uncertainties
        # Fill sparse matrix index and value arrays
        istart = iend
        iend = istart + N_W
        ir[istart:iend] = ones(N_W) * i
        jc[istart:iend] = arange(N_W) + col_i
        ws[istart:iend] = ones(N_W) / N_W
        uncertainty_vector[col_i:col_i + N_W] = ui
    # Trim off trailing zeros if maximum size not required, i.e. if all windows are not N_W in length
    ir = trim_zeros(ir, trim='b')
    jc = trim_zeros(jc, trim='b')
    ws = trim_zeros(ws, trim='b')
    uncertainty_vector = trim_zeros(uncertainty_vector, trim='b')

    # 3. Build sparse matrix
    W = csr_matrix((ws, (ir, jc)))

    return W, uncertainty_vector


def return_w_matrices(u_C_S_sensor_1, u_C_ICT_sensor_1, u_C_S_sensor_2, u_C_ICT_sensor_2,
                      time_sensor_1, width_sensor_1, time_sensor_2, width_sensor_2, sensor_i_name):

    """
    Return w matrices and uncertainty vectors for a match-up

    :type u_C_S_sensor_1: numpy.ndarray
    :param u_C_S_sensor_1: Scanline uncertainties for space counts for sensor 1

    :type u_C_ICT_sensor_1: numpy.ndarray
    :param u_C_ICT_sensor_1: Scanline uncertainties for ICT counts for sensor 1

    :type u_C_S_sensor_2: numpy.ndarray
    :param u_C_S_sensor_2: Scanline uncertainties for space counts for sensor 2

    :type u_C_ICT_sensor_2: numpy.ndarray
    :param u_C_ICT_sensor_2: Scanline uncertainties for ICT counts for sensor 2

    :type time_sensor_1: numpy.ndarray
    :param time_sensor_1: match-up time for sensor 1

    :type width_sensor_1: int
    :param width_sensor_1: scanline width for sensor 1

    :type time_sensor_2: numpy.ndarray
    :param time_sensor_2: match-up time for sensor 2

    :type width_sensor_2: int
    :param width_sensor_2: scanline width for sensor 2

    :type sensor_i_name: int
    :param sensor_i_name: sensor 1 ID

    :return:
        :W: *scipy.sparse.csr_matrix*
        weighting matrix
        :uncertainty_vector: *numpy.ndarray*
        scanline uncertainty vector
    """

    # 1. Build w matrix and uncertainty vectors for Sensor 2
    w_C_S_sensor_2, uncertainty_vector_S_sensor_2 = return_W_matrix(time_sensor_2, width_sensor_2, u_C_S_sensor_2)
    _, uncertainty_vector_ICT_sensor_2 = return_W_matrix(time_sensor_2, width_sensor_2, u_C_ICT_sensor_2)

    w_matrices = [w_C_S_sensor_2]
    uncertainty_vectors = [uncertainty_vector_S_sensor_2, uncertainty_vector_ICT_sensor_2]

    # 2. Build w matrix and uncertainty vectors for Sensor 1 (if not reference)
    if sensor_i_name != -1:
        # a. Compute W data
        w_C_S_sensor_1, uncertainty_vector_S_sensor_1 = return_W_matrix(time_sensor_1, width_sensor_1, u_C_S_sensor_1)
        _, uncertainty_vector_ICT_sensor_1 = return_W_matrix(time_sensor_1, width_sensor_1, u_C_ICT_sensor_1)
        w_matrices = [w_C_S_sensor_1, w_C_S_sensor_2]
        uncertainty_vectors = [uncertainty_vector_S_sensor_1, uncertainty_vector_ICT_sensor_1,
                               uncertainty_vector_S_sensor_2, uncertainty_vector_ICT_sensor_2]

    return w_matrices, uncertainty_vectors


def main(input_file_path, output_file_path):
    """
    Routine to update Jon's AVHRR harmonisation matchup file to include w variables and match newest spec

    :type input_file_path: str
    :param input_file_path: path of input harmonisation input file

    :type output_file_path: str
    :param output_file_path: path of output harmonisation input file
    """

    # 1. Read required data to generate W matrices
    print "Reading file:", input_file_path
    H, Us, Ur, K, Kr, Ks, sensor_1_name, sensor_2_name, time_sensor_1, time_sensor_2, \
    width_sensor_1, width_sensor_2, u_C_S_sensor_1, u_C_ICT_sensor_1, \
    u_C_S_sensor_2, u_C_ICT_sensor_2 = read_matchup_data(input_file_path)

    # Find bad data to remove:
    # > Determine valid scanlines (i.e. full averaging kernal available)
    valid_averages = return_valid_averages_mask(u_C_S_sensor_1, u_C_ICT_sensor_1, u_C_S_sensor_2, u_C_ICT_sensor_2,
                                                sensor_1_name)
    # > Remove incorrect assignment in Ur array
    if sensor_1_name != -1:
        Ur[:, [0, 1, 5, 6]] = 0
    else:
        Ur[:, [5, 6]] = 0
        Us[:, 0] = 0
    # > Define X1, X2, Ur1, Ur2, Us1 and Us2 arrays
    m1 = 5
    m2 = 5
    if sensor_1_name == -1:
        m1 = 1
    elif (H[:, 9] == 0).all():
        m1 = 4
    if (H[:, 9] == 0).all():
        m2 = 4

    X1 = H[valid_averages, :m1]
    Ur1 = Ur[valid_averages, :m1]
    Us1 = Us[valid_averages, :m1]
    X2 = H[valid_averages, 5:5+m2]
    Ur2 = Ur[valid_averages, 5:5+m2]
    Us2 = Us[valid_averages, 5:5+m2]

    del H, Ur, Us

    # 2. Generate required W matrix variables
    print "Generating W Matrix Variables from data..."

    # a. Generate lists of W matrices and uncertainty vectors
    w_matrices, uncertainty_vectors = return_w_matrices(u_C_S_sensor_1[valid_averages, :],
                                                        u_C_ICT_sensor_1[valid_averages, :],
                                                        u_C_S_sensor_2[valid_averages, :],
                                                        u_C_ICT_sensor_2[valid_averages, :],
                                                        time_sensor_1[valid_averages], width_sensor_1,
                                                        time_sensor_2[valid_averages], width_sensor_2,
                                                        sensor_1_name)
    del u_C_S_sensor_1, u_C_ICT_sensor_1, u_C_S_sensor_2, u_C_ICT_sensor_2

    # b. Convert W matrices and uncertainty vector to input file variables
    w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz, \
        uncertainty_vector_row_count, uncertainty_vector = return_w_matrix_variables(w_matrices, uncertainty_vectors)
    del w_matrices, uncertainty_vectors

    # c. Assign type/use matrices
    w_matrix_use1 = array([1, 1, 0, 0, 0], dtype=int8)
    w_matrix_use2 = array([2, 2, 0, 0, 0], dtype=int8)
    uncertainty_vector_use1 = array([1, 2, 0, 0, 0], dtype=int8)
    uncertainty_vector_use2 = array([3, 4, 0, 0, 0], dtype=int8)
    uncertainty_type1 = array([3, 3, 1, 2, 2], dtype=int8)
    uncertainty_type2 = array([3, 3, 1, 2, 2], dtype=int8)
    if sensor_1_name == -1:
        w_matrix_use1 = array([0], dtype=int8)
        w_matrix_use2 = array([1, 1, 0, 0, 0], dtype=int8)
        uncertainty_vector_use1 = array([0], dtype=int8)
        uncertainty_vector_use2 = array([1, 2, 0, 0, 0], dtype=int8)
        uncertainty_type1 = array([1], dtype=int8)
        uncertainty_type2 = array([3, 3, 1, 2, 2], dtype=int8)

    if X2.shape[1] == 4:
        w_matrix_use1 = array([1, 1, 0, 0], dtype=int8)
        w_matrix_use2 = array([2, 2, 0, 0], dtype=int8)
        uncertainty_vector_use1 = array([1, 2, 0, 0], dtype=int8)
        uncertainty_vector_use2 = array([3, 4, 0, 0], dtype=int8)
        uncertainty_type1 = array([3, 3, 1, 1], dtype=int8)
        uncertainty_type2 = array([3, 3, 1, 1], dtype=int8)
        if sensor_1_name == -1:
            w_matrix_use1 = array([0], dtype=int8)
            w_matrix_use2 = array([1, 1, 0, 0], dtype=int8)
            uncertainty_vector_use1 = array([0], dtype=int8)
            uncertainty_vector_use2 = array([1, 2, 0, 0], dtype=int8)
            uncertainty_type1 = array([1], dtype=int8)
            uncertainty_type2 = array([3, 3, 1, 1], dtype=int8)

    # 3. Update file to include W matrix variables
    print "Writing to new file:", output_file_path

    # a. Update existing file to comply with new specification
    write_input_file(output_file_path,
                     X1, X2,
                     Ur1, Ur2, Us1, Us2, uncertainty_type1, uncertainty_type2,
                     K[valid_averages], Kr[valid_averages], Ks[valid_averages],
                     time_sensor_1[valid_averages], time_sensor_2[valid_averages],
                     sensor_1_name, sensor_2_name)

    # b. Add w matrix variables to updated file
    append_W_to_input_file(output_file_path,
                           w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz,
                           uncertainty_vector_row_count, uncertainty_vector,
                           w_matrix_use1, w_matrix_use2, uncertainty_vector_use1, uncertainty_vector_use2)

    # 4. Test harmonisation
    test_input_file(output_file_path)

    print "Done"
    return 0

if __name__ == "__main__":
    main("/home/seh2/downloads/n09_n08.nc", "/home/seh2/downloads/n09_n08_test.nc")
    #main(argv[1], argv[2])



