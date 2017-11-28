"""
Functions convert N scipy.sparse.csr_matrix object W matrices into corresponding variables for harmonisation input files
and attach them to files
"""

'''___Python Modules___'''
from numpy import zeros

'''___Third Party Modules___'''
from netCDF4 import Dataset

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "28/09/2017"
__credits__ = ["Ralf Quast", "Jon Mittaz", "Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


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
                                       "(2) Independent + Systematic Error Correlation, or " + \
                                       "(3) Structured Error Correlation"
    uncertainty_type1_var[:] = uncertainty_type1[:]

    # > uncertainty_type2 - Uncertainty correlation type per X2 column
    uncertainty_type2_var = rootgrp.createVariable('uncertainty_type2', 'i4', ('m2',), zlib=True, complevel=9)
    uncertainty_type2_var.description = "Uncertainty correlation type per X2 column, labelled as, " + \
                                        "(1) Independent Error Correlation, " + \
                                        "(2) Independent + Systematic Error Correlation, or " + \
                                        "(3) Structured Error Correlation"
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
    w_matrix_val_var = rootgrp.createVariable('w_matrix_val', 'f8', ('w_matrix_nnz_sum',), zlib=True, complevel=9)
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
    u_matrix_val_var = rootgrp.createVariable('u_matrix_val', 'f8', ('u_matrix_row_count_sum',),
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

if __name__ == "__main__":
    pass
