"""Script to convert PHM match-up files to FIDUCEO format, try 'python phm2fiduco.py --help' for usage"""

'''___Built-In Modules___'''
import sys
import argparse
import logging
import os

'''___Third-Party Modules___'''
from numpy import loadtxt, savetxt, zeros, full, array, arange, int8, ones, trim_zeros, asarray, arange, unique, int32
from scipy.sparse import csr_matrix
from netCDF4 import Dataset

'''___NPL Modules___'''


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "18/06/2018"
__credits__ = []
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

PMH_FNAMES = ["data_pams.txt", "data_vals.txt", "data_uncs.txt"]


def parse_cmdline():
    parser = argparse.ArgumentParser(
        description="Tool to convert PHM's simulated sensor match-up files to FIDUCEO file format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--input_directory", action="store",
                        help="Path of directory containing PHM format files")

    parser.add_argument("-o", "--output_directory", action="store", default=".",
                        help="Directory to write output file to")

    logging_options = parser.add_mutually_exclusive_group()

    logging_options.add_argument("--verbose", action="store_true",
                                 help="Option for verbose output")

    logging_options.add_argument("--quiet", action="store_true",
                                 help="Option for quiet output")

    parser.add_argument("--log", action="store", type=str,
                        help="Log file to write to. Leave out for stdout.")

    parser.add_argument("--version", action="version", version='v%s' % __version__)

    return parser.parse_args()

parsed_cmdline = parse_cmdline()


def configure_logging(fname, verbose=False, quiet=False):
    """
    Configure logger

    :param path_to_log_directory:  path to directory to write log file in
    :return:
    """

    logger = logging.getLogger(__name__)
    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)

    if fname is not None:
        file_formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(message)s')
        file_handler = logging.FileHandler(fname)
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    stream_formatter = logging.Formatter('%(message)s')
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(stream_formatter)
    logger.addHandler(stream_handler)

    return logger

logger = configure_logging(parsed_cmdline.log, parsed_cmdline.verbose, parsed_cmdline.quiet)


def prepare_directories(input_directory, output_directory):

    logger.info("Input directory: "+input_directory)
    logger.info("Output directory: "+output_directory)

    # Check directories
    logger.debug("Preparing directories...")

    # Check input directory exists
    input_directory_exists = os.path.isdir(input_directory)
    if not input_directory_exists:
        logger.debug("Input directory does not exist")
        raise IOError("Input directory does not exist '"+input_directory+"'")
    else:
        logger.debug("Input directory already exists")

    # Check required files in input directory
    input_directory_contents = os.listdir(input_directory)
    pmh_files_present = [True if f in input_directory_contents else False for f in PMH_FNAMES]

    if all(f is True for f in pmh_files_present):
        logger.debug("Input directory contains required PMH files")
    else:
        missing_files = [f for f, p in zip(PMH_FNAMES, pmh_files_present) if p is False]
        logger.debug("Input directory missing PMH files '"+str(missing_files)+"'")
        raise IOError("Input directory missing PMH files '"+str(missing_files)+"'")

    # If required create output directory
    output_directory_exists = os.path.isdir(output_directory)
    if output_directory_exists:
        logger.debug("Output directory already exists")
    else:
        os.makedirs(output_directory)
        logger.debug("Output directory created")

    return 0


def open_data_vals(input_directory):
    path = os.path.join(input_directory, "data_vals.txt")
    logger.debug("Opening: "+path)
    return loadtxt(path)


def open_data_pams(input_directory):
    path = os.path.join(input_directory,"data_pams.txt")
    logger.debug("Opening: " + path)
    with open(path) as data_pams:
        Im = []
        Nm = []
        data = []
        for line in data_pams:
            if "Im" in line:
                pass
            elif "Nm" in line:
                Im = data
                data = []
            elif "atrue" in line:
                Nm = data
                data = []
            else:
                data.append(line.rstrip().split())
        atrue_list = data

    Nm = [int(n[0]) for n in Nm]

    atrue = zeros((len(atrue_list), len(atrue_list[0])))
    for i, row in enumerate(atrue_list):
        for j, elem in enumerate(row):
            atrue[i,j] = elem
    return Im, Nm, atrue


def open_data_uncs(input_directory):
    path = os.path.join(input_directory, "data_uncs.txt")
    logger.debug("Opening: " + path)
    with open(path) as data_uncs:
        for i_tot, l in enumerate(data_uncs):
            pass

    uR = []
    uS = []
    nW = []
    uncertainty_type = []
    with open(path) as data_uncs:
        for i, line in enumerate(data_uncs):
            if (i == 0) or (i == i_tot-1):
                pass
            elif i == 1:
                uR_ref = float(line)
            elif i == i_tot:
                Kr = float(line)
            elif i % 2 == 0:
                unc_type = line.split()[2][:-1]
            else:
                if unc_type == "R":
                    uR.append(float(line.rstrip()))
                    uS.append(0.0)
                    nW.append(0)
                    uncertainty_type.append(1)
                elif unc_type == "RS":
                    uR.append(float(line.split()[0]))
                    uS.append(float(line.split()[1]))
                    nW.append(0)
                    uncertainty_type.append(2)
                elif unc_type == "MA":
                    uR.append(float(line.split()[0]))
                    uS.append(0.0)
                    nW.append(int(line.split()[1]))
                    uncertainty_type.append(3)
                else:
                    raise TypeError("Unrecognised uncertainty type in input data")

    return uR_ref, uR, uS, nW, Kr, uncertainty_type


def generate_rolling_average_w_matrix(matchup_times, scanline_time, kernel_uncertainty, weights=None, kernel_size=None):
    """
    Return w matrix in sparse representation and diagonal of u matrix for a rolling average operation, where the
    covariance matrix is w_matrix^T u_matrix^T u_matrix w_matrix

    :type matchup_times: numpy.ndarray
    :param matchup_times: match-up times for match-ups in match-up series

    :type scanline_time: numpy.ndarray
    :param scanline_time: time interval between consecutive scanlines

    :type kernel_uncertainty: float/numpy.ndarray
    :param kernel_uncertainty: uncertainty per unaveraged value - either full array for each value or one value for all

    :type weights: str
    :param weights: kernel weighting of averages

    :type kernel_size: int
    :param kernel_size: length of averaging kernel, only required if weights is None and kernel_uncertainty is float

    :return:
        :w_matrix: *scipy.sparse.csr_matrix*

        W matrix

        :u_matrix: *numpy.ndarray*

        Uncertainty matrix
    """

    # Useful parameters
    num_matchups = int(len(matchup_times))          # number of match_ups

    if weights is not None:
        kernel_size = len(weights)

    elif (kernel_size is not None) and (weights is None):
        weights = ones(kernel_size) / kernel_size

    elif (kernel_size is None) and (weights is None):
        kernel_size = int(len(kernel_uncertainty[0]))
        weights = ones(kernel_size) / kernel_size

    else:
        raise RuntimeError("Error: Averaging kernel not appropriately specified")

    single_uncertainty = False
    if type(kernel_uncertainty) == float:
        single_uncertainty = True

    ####################################################################################################################
    # 1. Determine first column index of each row
    ####################################################################################################################

    # column index of first non-zero element of each w row relative to first match-up
    relative_column_postions = []

    for i in range(num_matchups):
        scanline_pos = 0
        if i > 0:
            rel_time = round(scanline_time*round(float(matchup_times[i]-matchup_times[0])/scanline_time), 2)
            scanline_pos = int(round(rel_time/scanline_time))

        relative_column_postions.append(int(scanline_pos))

    # Alter relative columns positions so that zero is the lowest value to give absolute column position
    relative_column_postions = asarray(relative_column_postions)
    column_positions = relative_column_postions - min(relative_column_postions)

    ####################################################################################################################
    # 2. Generate matrices
    ####################################################################################################################

    # 1. W Matrix

    # Populate sparse W matrix index/data vectors
    JA = zeros(num_matchups * kernel_size, dtype=int32)
    WA = zeros(num_matchups * kernel_size)
    UA = zeros(num_matchups * kernel_size)
    for i, column_position in enumerate(column_positions):
        WA[i*kernel_size:i*kernel_size+kernel_size] = weights
        JA[i*kernel_size:i*kernel_size+kernel_size] = arange(column_position, column_position+kernel_size)

    IA = zeros(num_matchups+1, dtype=int32)
    for i in range(num_matchups+1):
        IA[i] = i*kernel_size

    w_shape = (num_matchups, max(JA)+1)
    w_matrix = csr_matrix((WA, JA, IA), shape=w_shape)

    # Remove zero columns
    w_zero_col = unique(w_matrix.nonzero()[1])
    w_matrix = w_matrix[:, w_zero_col]

    # 2. U Matrix

    # First generate as sparse
    if single_uncertainty is True:
        u_matrix = kernel_uncertainty*ones(w_matrix.shape[1])

    else:
        UA = zeros(num_matchups * kernel_size)
        for i, column_position in enumerate(column_positions):
            UA[i * kernel_size:i * kernel_size + kernel_size] = kernel_uncertainty[i]
        u_matrix = csr_matrix((UA, JA, IA), shape=w_shape)
        u_matrix = u_matrix[:, w_zero_col]
        u_matrix = u_matrix.max(axis=0).toarray()[0]

    return w_matrix.sorted_indices(), u_matrix


def write_input_file(file_path, X1, X2, Ur1, Ur2, Us1, Us2, uncertainty_type1, uncertainty_type2, K, Kr, Ks,
                     time1, time2, sensor_1_name, sensor_2_name, additional_variables=None):
    """
    Write harmonisation input file from input data arrays (no w matrix variables, see func append_W_to_input_file(...)
    for this functionality)

    :type file_path: str
    :param file_path: match-up file path

    :type X1: numpy.ndarray
    :param X1: Radiances and counts per matchup for sensor 1

    :type X2: numpy.ndarray
    :param X2: Radiances and counts per matchup for sensor 2

    :type Ur1: numpy.ndarray
    :param Ur1: Random uncertainties for X1 array

    :type Ur2: numpy.ndarray
    :param Ur2: Random uncertainties for X2 array

    :type Us1: numpy.ndarray
    :param Us1: Systematic uncertainties for X1 array

    :type Us2: numpy.ndarray
    :param Us2: Systematic uncertainties for X2 array

    :type uncertainty_type1: numpy.ndarray
    :param uncertainty_type1: Uncertainty correlation type per X1 column

    :type uncertainty_type2: numpy.ndarray
    :param uncertainty_type2: Uncertainty correlation type per X2 column

    :type K: numpy.ndarray
    :param K: K (sensor-to-sensor differences) for zero shift case

    :type Kr: numpy.ndarray
    :param Kr: K (sensor-to-sensor differences) random uncertainties (matchup uncertainty)

    :type Ks: numpy.ndarray
    :param Ks: K (sensor-to-sensor differences) systematic uncertainties for zero shift case

    :type time1: numpy.ndarray
    :param time1: "Match-up time sensor 1, seconds since 1970-01-01"

    :type time2: numpy.ndarray
    :param time2: "Match-up time sensor 2, seconds since 1970-01-01"

    :type sensor_1_name: int
    :param sensor_1_name: sensor i ID

    :type sensor_2_name: int
    :param sensor_2_name: sensor j ID

    :type additional_variables: dict
    :param additional_variables: dictionary of additional, non-required variable to add to harmonisation input files. To
    be for e.g. testing, diagnostics etc.

    One dictionary entry per additional variable, with each entry as:

    "variable_name": {"data": data_array, "dtype": dtype_str, "dim": dim_tuple, "Description", desc_str}

    where:
    * data_array(*numpy.ndarray*) - array variable data
    * dtype_str(*str*) - netCDF variable data type (e.g. "i4", "f8", ...)
    * dim_tuple(*tuple:str*) - tuple of the variable dimension names (e.g. ('M',) )
    * desc_str(*str*) - description of the variable
    """

    # 1. Open file
    rootgrp = Dataset(file_path, mode='w')

    # 2. Create attributes
    rootgrp.sensor_1_name = sensor_1_name
    rootgrp.sensor_2_name = sensor_2_name

    # 2. Create dimensions
    # > M - number of matchups
    M_dim = rootgrp.createDimension("M", X1.shape[0])

    # > m - number of columns in X1 and X2 arrays
    m1_dim = rootgrp.createDimension("m1", X1.shape[1])
    m2_dim = rootgrp.createDimension("m2", X2.shape[1])

    # 3. Create new variables

    # > X1 - Radiances and counts per matchup for sensor 1
    X1_var = rootgrp.createVariable('X1', 'f8', ('M', 'm1'), zlib=True, complevel=9)
    X1_var.description = "Radiances and counts per matchup for sensor 1"
    X1_var[:] = X1[:]

    # > X2 - Radiances and counts per matchup for sensor 2
    X2_var = rootgrp.createVariable('X2', 'f8', ('M', 'm2'), zlib=True, complevel=9)
    X2_var.description = "Radiances and counts per matchup for sensor 2"
    X2_var[:] = X2[:]

    # > Ur1 - Random uncertainties for X1 array
    Ur1_var = rootgrp.createVariable('Ur1', 'f8', ('M', 'm1'), zlib=True, complevel=9)
    Ur1_var.description = "Random uncertainties for X1 array"
    Ur1_var[:] = Ur1[:]

    # > Ur2 - Random uncertainties for X2 array
    Ur2_var = rootgrp.createVariable('Ur2', 'f8', ('M', 'm2'), zlib=True, complevel=9)
    Ur2_var.description = "Random uncertainties for X2 array"
    Ur2_var[:] = Ur2[:]

    # > Us1 - Systematic uncertainties for X1 array
    Us1_var = rootgrp.createVariable('Us1', 'f8', ('M', 'm1'), zlib=True, complevel=9)
    Us1_var.description = "Systematic uncertainties for X1 array"
    Us1_var[:] = Us1[:]

    # > Us2 - Systematic uncertainties for X2 array
    Us2_var = rootgrp.createVariable('Us2', 'f8', ('M', 'm2'), zlib=True, complevel=9)
    Us2_var.description = "Systematic uncertainties for X2 array"
    Us2_var[:] = Us2[:]

    # > uncertainty_type1 - Uncertainty correlation type per X1 column
    uncertainty_type1_var = rootgrp.createVariable('uncertainty_type1', 'i4', ('m1',), zlib=True, complevel=9)
    uncertainty_type1_var.description = "Uncertainty correlation type per X1 column, labelled as, " + \
                                        "(1) Independent Error Correlation, " + \
                                        "(2) Independent + Systematic Error Correlation, " + \
                                        "(3) Structured Error Correlation, or" + \
                                        "(4) Structured + Systematic Error Correlation"
    uncertainty_type1_var[:] = uncertainty_type1[:]

    # > uncertainty_type2 - Uncertainty correlation type per X2 column
    uncertainty_type2_var = rootgrp.createVariable('uncertainty_type2', 'i4', ('m2',), zlib=True, complevel=9)
    uncertainty_type2_var.description = "Uncertainty correlation type per X2 column, labelled as, " + \
                                        "(1) Independent Error Correlation, " + \
                                        "(2) Independent + Systematic Error Correlation, " + \
                                        "(3) Structured Error Correlation, or" + \
                                        "(4) Structured + Systematic Error Correlation"
    uncertainty_type2_var[:] = uncertainty_type2[:]

    # > K - K (sensor-to-sensor differences) for zero shift case
    K_var = rootgrp.createVariable('K', 'f8', ('M',), zlib=True, complevel=9)
    K_var.description = "K (sensor-to-sensor differences) for zero shift case"
    K_var[:] = K[:]

    # > Kr - K (sensor-to-sensor differences) random uncertainties (matchup uncertainty)
    Kr_var = rootgrp.createVariable('Kr', 'f8', ('M',), zlib=True, complevel=9)
    Kr_var.description = "K (sensor-to-sensor differences) random uncertainties (matchup uncertainty)"
    Kr_var[:] = Kr[:]

    # > Ks - K (sensor-to-sensor differences) systematic uncertainties for zero shift case
    Ks_var = rootgrp.createVariable('Ks', 'f8', ('M',), zlib=True, complevel=9)
    Ks_var.description = "K (sensor-to-sensor differences) systematic uncertainties for zero shift case"
    Ks_var[:] = Ks[:]

    # > time1 - Sensor 1 time of match-up
    time1_var = rootgrp.createVariable('time1', 'f8', ('M',), zlib=True, complevel=9)
    time1_var.description = "Match-up time sensor 1, seconds since 1970-01-01"
    time1_var[:] = time1[:]

    # > time 2 - Sensor 2 time of match-up
    time2_var = rootgrp.createVariable('time2', 'f8', ('M',), zlib=True, complevel=9)
    time2_var.description = "Match-up time sensor 2, seconds since 1970-01-01"
    time2_var[:] = time2[:]

    # > additional variables - non-required variable to add to harmonisation input files
    if additional_variables is not None:
        for variable in additional_variables.keys():
            additional_var = rootgrp.createVariable(variable, additional_variables[variable]['dtype'],
                                                    additional_variables[variable]['dim'], zlib=True, complevel=9)
            additional_var.Description = additional_variables[variable]['Description']
            additional_var[:] = additional_variables[variable]['data'][:]

    # 5. Close file
    rootgrp.close()

    return 0


def return_w_matrix_variables(w_matrices, u_matrices):
    """
    Produce harmonisation input file W matrix variable arrays from lists of W and U matrices.
    NB: Also additionally required for the file would be w_matrix_use1/2 and u_matrix_use1/2 variable arrays.

    :type w_matrices: list:scipy.sparse.csr_matrix
    :param w_matrices: List of scipy.sparse.csr_matrix object W matrices
                       (In the order they will be indexed in harmonisation input file w_matrix_use1/2 variable)

    :type u_matrices: list:numpy.ndarray
    :param u_matrices: List of U matrices
                         (In the order they will be indexed in harmonisation input file u_matrices_use1/2 variable)

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

        Number of elements per uncertainty vector

        :u_matrix_val: *numpy.ndarray*

        Concatenated array of u matrix diagonals
    """

    # 1. Generate W matrix variable arrays

    # > w_matrix_nnz
    w_matrix_nnz = zeros(len(w_matrices), dtype=int)
    for i, w in enumerate(w_matrices):
        w_matrix_nnz[i] = w.nnz

    # > w_matrix_val & w_matrix_col
    total_nnz = sum(w_matrix_nnz)
    w_matrix_val = zeros(total_nnz)
    w_matrix_col = zeros(total_nnz)

    istart = 0
    iend = 0
    for w, w_nnz in zip(w_matrices, w_matrix_nnz):
        iend += w_nnz
        w_matrix_val[istart:iend] = w.data
        w_matrix_col[istart:iend] = w.indices
        istart = iend

    # > w_matrix_row
    w_matrix_row = zeros((len(w_matrices), len(w_matrices[0].indptr)))
    for i, w in enumerate(w_matrices):
        w_matrix_row[i, :] = w.indptr

    # 2. Generate u matrix variable arrays

    # > u_matrix_row_count
    u_matrix_row_count = zeros(len(u_matrices), dtype=int)
    for i, u in enumerate(u_matrices):
        u_matrix_row_count[i] = len(u)

    # > u_matrix_val
    total_row = sum(u_matrix_row_count)
    u_matrix_val = zeros(total_row)

    istart = 0
    iend = 0
    for u, u_row in zip(u_matrices, u_matrix_row_count):
        iend += u_row
        u_matrix_val[istart:iend] = u
        istart = iend

    return w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz, u_matrix_row_count, u_matrix_val


def append_W_to_input_file(filepath,
                           w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz,
                           u_matrix_row_count, u_matrix_val,
                           w_matrix_use1, w_matrix_use2, u_matrix_use1, u_matrix_use2):
    """
    Append set of W matrix variable arrays to a given harmonisation input file

    :type w_matrix_val: numpy.ndarray
    :param w_matrix_val: Concatenated array of non-zero elements of W matrices

    :type w_matrix_row: numpy.ndarray
    :param w_matrix_row: Concatenated array of row indices of W matrices

    :type w_matrix_col: numpy.ndarray
    :param w_matrix_col: Concatenated array of row indices of W matrices

    :type w_matrix_nnz: numpy.ndarray
    :param w_matrix_nnz: Number of non-zeros elements per W matrix

    :type u_matrix_row_count: numpy.ndarray
    :param u_matrix_row_count: Number of non-zero elements per u matrix

    :type u_matrix_val: numpy.ndarray
    :param u_matrix_val: Concatenated array of u matrix non-zero values

    :type w_matrix_use1: numpy.ndarray
    :param w_matrix_use1: mapping from X1 array column index to W

    :type w_matrix_use2: numpy.ndarray
    :param w_matrix_use2: mapping from X2 array column index to W

    :type u_matrix_use1: numpy.ndarray
    :param u_matrix_use1: a mapping from X1 array column index to U

    :type u_matrix_use2: numpy.ndarray
    :param u_matrix_use2: a mapping from X2 array column index to U
    """

    # 1. Open file
    rootgrp = Dataset(filepath, mode='a')

    # 2. Create dimensions
    # > w_matrix_count - number of W matrices
    w_matrix_count_dim = rootgrp.createDimension("w_matrix_count", len(w_matrix_nnz))

    # > w_matrix_row - number of rows in each w matrix
    if len(w_matrix_row.shape) > 1:
        num_row = w_matrix_row.shape[1]
    else:
        num_row = len(w_matrix_row)
    w_matrix_row_count_dim = rootgrp.createDimension("w_matrix_row_count", num_row)

    # > w_matrix_sum_nnz - sum of non-zero elements in all W matrices
    w_matrix_nnz_sum_dim = rootgrp.createDimension("w_matrix_nnz_sum", sum(w_matrix_nnz))

    # > u_matrix_count - number of u matrices
    u_matrix_count_dim = rootgrp.createDimension("u_matrix_count", len(u_matrix_row_count))

    # > u_matrix_row_count_sum - sum of rows in u matrices
    u_matrix_row_count_sum_dim = rootgrp.createDimension("u_matrix_row_count_sum", sum(u_matrix_row_count))

    # 3. Create new variables
    # > w_matrix_nnz - number of non-zero elements for each W matrix
    w_matrix_nnz_var = rootgrp.createVariable('w_matrix_nnz', 'i4', ('w_matrix_count',), zlib=True, complevel=9)
    w_matrix_nnz_var.description = "number of non-zero elements for each W matrix"

    # > w_matrix_row - CSR row numbers for each W matrix
    row_dims = ('w_matrix_row_count',)
    if len(w_matrix_row.shape) > 1:
        row_dims = ('w_matrix_count', 'w_matrix_row_count')
    w_matrix_row_var = rootgrp.createVariable('w_matrix_row', 'i4', row_dims, zlib=True, complevel=9)
    w_matrix_row_var.description = "CSR row numbers for each W matrix"

    # > w_matrix_col - CSR column numbers for all W matrices
    w_matrix_col_var = rootgrp.createVariable('w_matrix_col', 'i4', ('w_matrix_nnz_sum',), zlib=True, complevel=9)
    w_matrix_col_var.description = "CSR column numbers for all W matrices"

    # > w_matrix_val - CSR values for all W matrices
    w_matrix_val_var = rootgrp.createVariable('w_matrix_val', 'f4', ('w_matrix_nnz_sum',), zlib=True, complevel=9)
    w_matrix_val_var.description = "CSR values for all W matrices"

    # > w_matrix_use1 - a mapping from X2 array column index to W
    w_matrix_use1_var = rootgrp.createVariable('w_matrix_use1', 'i4', ('m1',), zlib=True, complevel=9)
    w_matrix_use1_var.description = "mapping from X1 array column index to W"

    # > w_matrix_use2 - a mapping from X2 array column index to W
    w_matrix_use2_var = rootgrp.createVariable('w_matrix_use2', 'i4', ('m2',), zlib=True, complevel=9)
    w_matrix_use2_var.description = "mapping from X2 array column index to W"

    # > u_matrix_row_count - number of rows of each u matrix
    u_matrix_row_count_var = rootgrp.createVariable('u_matrix_row_count', 'i4',
                                                              ('u_matrix_count',), zlib=True, complevel=9)
    u_matrix_row_count_var.description = "number of rows of each u matrix"

    # > u matrix val - uncertainty of each scanline value
    u_matrix_val_var = rootgrp.createVariable('u_matrix_val', 'f4', ('u_matrix_row_count_sum',),
                                                    zlib=True, complevel=9)
    u_matrix_val_var.description = "u matrix non-zero diagonal elements"

    # > u_matrix_use1 - a mapping from X1 array column index to U
    u_matrix_use1_var = rootgrp.createVariable('u_matrix_use1', 'i4', ('m1',), zlib=True, complevel=9)
    u_matrix_use1_var.description = "mapping from X1 array column index to U"

    # > u_matrix_use2 - a mapping from X2 array column index to U
    u_matrix_use2_var = rootgrp.createVariable('u_matrix_use2', 'i4', ('m2',), zlib=True, complevel=9)
    u_matrix_use2_var.description = "mapping from X2 array column index to U"

    # 4. Add data
    w_matrix_nnz_var[:] = w_matrix_nnz[:]
    w_matrix_row_var[:] = w_matrix_row[:]
    w_matrix_col_var[:] = w_matrix_col[:]
    w_matrix_val_var[:] = w_matrix_val[:]
    w_matrix_use1_var[:] = w_matrix_use1[:]
    w_matrix_use2_var[:] = w_matrix_use2[:]

    u_matrix_row_count_var[:] = u_matrix_row_count[:]
    u_matrix_val_var[:] = u_matrix_val[:]
    u_matrix_use1_var[:] = u_matrix_use1[:]
    u_matrix_use2_var[:] = u_matrix_use2[:]

    # 5. Close file
    rootgrp.close()

    return 0


def run_conversion(input_directory, output_directory):

    # Open data
    logger.info("Opening Input Data...")

    H = open_data_vals(input_directory)
    logger.debug("H.shape: " + str(H.shape))

    Im, Nm, atrue = open_data_pams(input_directory)
    logger.debug("Im: " + str(Im))
    logger.debug("Nm: " + str(Nm))
    logger.debug("atrue: " + str(atrue))

    uR_ref, uR, uS, nW, Kr, uncertainty_type = open_data_uncs(input_directory)
    logger.debug("uR_ref: "+str(uR_ref))
    logger.debug("uR: "+str(uR))
    logger.debug("uS: "+str(uS))
    logger.debug("nW: "+str(nW))
    logger.debug("Kr: "+str(Kr))
    logger.debug("uncertainty_type: "+str(uncertainty_type))

    # Write data to file
    logger.info("Writing Output Data...")
    i_start = 0
    i_end = 0
    nX = len(uR)

    output_file_paths = [os.path.join(output_directory, "_".join(s)+".nc") for s in Im]

    for i, (sensor_pair, n_mu, file_path) in enumerate(zip(Im, Nm, output_file_paths)):
        sensor_1_name = sensor_pair[0]
        sensor_2_name = sensor_pair[1]
        logger.debug("Writing: "+file_path)

        i_start = i_end
        i_end = i_start + n_mu

        if sensor_pair[0] == '0':
            X1 = array([H[i_start:i_end, 0]]).T
            Ur1 = full(n_mu, uR_ref)
            Us1 = zeros(n_mu)
            uncertainty_type1 = array([1])
        else:
            X1 = H[i_start:i_end, :nX]
            Ur1 = zeros((n_mu, nX))
            Us1 = zeros((n_mu, nX))
            for j, u in enumerate(uR):
                if (uncertainty_type[j] == 1) or (uncertainty_type[j] == 2):
                    Ur1[:,j] = u
            for j, u in enumerate(uS):
                Us1[:,j] = u
            uncertainty_type1 = array(uncertainty_type)

        X2 = H[i_start:i_end, nX:2*nX]
        Ur2 = zeros((n_mu, nX))
        Us2 = zeros((n_mu, nX))
        for j, u in enumerate(uR):
            if (uncertainty_type[j] == 1) or (uncertainty_type[j] == 2):
                Ur2[:, j] = u
        for j, u in enumerate(uS):
            Us2[:, j] = u
        uncertainty_type2 = array(uncertainty_type)

        K = H[i_start:i_end, -1]
        Kr = full(n_mu, Kr)
        Ks = zeros(n_mu)

        time = arange(i_start, i_end)

        write_input_file(file_path, X1, X2, Ur1, Ur2, Us1, Us2, uncertainty_type1, uncertainty_type2, K, Kr, Ks,
                         time, time, sensor_1_name, sensor_2_name, additional_variables=None)

    logger.info("Adding W matrices...")
    uR_w = [u for u,t in zip(uR, uncertainty_type) if t == 3]
    nW_w = [n for n,t in zip(nW, uncertainty_type) if t == 3]

    # Check all W matrices the same
    if (nW_w.count(nW_w[0]) != len(nW_w)) and (uR_w.count(uR_w[0]) != len(uR_w)):
        raise ValueError("Expected all W matrices to be the same, Sam needs to do more coding")
    else:
        w_matrix, uncertainty_matrix = generate_rolling_average_w_matrix(time, 1.0, uR_w[0], kernel_size=nW_w[0])
        w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz, \
        uncertainty_vector_row_count, uncertainty_vector = return_w_matrix_variables([w_matrix], [uncertainty_matrix])

    w_matrix_use2 = array([1 if t == 3 else 0 for t in uncertainty_type], int8)
    uncertainty_vector_use2 = array([1 if t == 3 else 0 for t in uncertainty_type], int8)

    for o, sensor_pair in zip(output_file_paths, Im):
        if sensor_pair[0] == '0':
            w_matrix_use1 = array([0], int8)
            uncertainty_vector_use1 = array([0], int8)
        else:
            w_matrix_use1 = array([1 if t == 3 else 0 for t in uncertainty_type], int8)
            uncertainty_vector_use1 = array([1 if t == 3 else 0 for t in uncertainty_type], int8)
        append_W_to_input_file(o,
                               w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz,
                               uncertainty_vector_row_count, uncertainty_vector,
                               w_matrix_use1, w_matrix_use2, uncertainty_vector_use1, uncertainty_vector_use2)

    logger.info("Saving True Parameters...")
    savetxt(os.path.join(output_directory, "atrue.txt"), atrue)

    logger.info("Done")
    return 0


def main(parsed_cmdline):

    logger.info("Converting Files")

    # Prepare directories
    input_directory = parsed_cmdline.input_directory
    output_directory = parsed_cmdline.output_directory
    prepare_directories(input_directory, output_directory)

    # Run conversion
    run_conversion(input_directory, output_directory)

    return 0

if __name__ == "__main__":
    #main(parsed_cmdline)
    main()
