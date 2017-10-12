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


def write_input_file(file_path, H, Us, Ur, K, Kr, Ks, sensor_i_name, sensor_j_name, additional_variables):
    """
    Write harmonisation input file from input data arrays (no w matrix variables, see func append_W_to_input_file(...)
    for this functionality)

    :type file_path: str
    :param file_path: match-up file path

    :type H: numpy.ndarray
    :param H: Radiances and counts per matchup

    :type Us: numpy.ndarray
    :param Us: Systematic uncertainties for H array

    :type Ur: numpy.ndarray
    :param Ur: Random uncertainties for H array

    :type K: numpy.ndarray
    :param K: K (sensor-to-sensor differences) for zero shift case

    :type Kr: numpy.ndarray
    :param Kr: K (sensor-to-sensor differences) random uncertainties (matchup uncertainty)

    :type Ks: numpy.ndarray
    :param Ks: K (sensor-to-sensor differences) systematic uncertainties for zero shift case

    :type sensor_i_name: int
    :param sensor_i_name: sensor i ID

    :type sensor_j_name: int
    :param sensor_j_name: sensor j ID

    :type additional_variables: dict
    :param additional_variables: dictionary of additional, non-required variable to add to harmonisation input files.
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
    rootgrp.sensor_i_name = sensor_i_name
    rootgrp.sensor_j_name = sensor_j_name

    # 2. Create dimensions
    # > M - number of matchups
    M_dim = rootgrp.createDimension("M", H.shape[0])

    # > m - number of columns in H
    m_dim = rootgrp.createDimension("m", H.shape[1])

    # 3. Create new variables

    # > H - Radiances and counts per matchup
    H_var = rootgrp.createVariable('H', 'f8', ('M', 'm'), zlib=True, complevel=9)
    H_var.description = "H array (M,m). Radiances and counts per matchup"
    H_var[:] = H[:]

    # > Us - Systematic uncertainties for H array
    Us_var = rootgrp.createVariable('Us', 'f8', ('M', 'm'), zlib=True, complevel=9)
    Us_var.description = "Systematic uncertainties for H array"
    Us_var[:] = Us[:]

    # > Ur - Random uncertainties for H array
    Ur_var = rootgrp.createVariable('Ur', 'f8', ('M', 'm'), zlib=True, complevel=9)
    Ur_var.description = "Random uncertainties for H array"
    Ur_var[:] = Ur[:]

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

    for variable in additional_variables.keys():
        additional_var = rootgrp.createVariable(variable, additional_variables[variable]['dtype'],
                                                additional_variables[variable]['dim'], zlib=True, complevel=9)
        additional_var.Description = additional_variables[variable]['Description']
        additional_var[:] = additional_variables[variable]['data'][:]

    # 5. Close file
    rootgrp.close()

    return 0

def return_w_matrix_variables(w_matrices, uncertainty_vectors):
    """
    Produce harmonisation input file W matrix variable arrays from lists of Ws matrices and uncertainty vectors.
    NB: Also additionally required for the file would be w_matrix_use and uncertainty_vector_use variable arrays.

    :type w_matrices: list:scipy.sparse.csr_matrix
    :param w_matrices: List of scipy.sparse.csr_matrix object W matrices
                       (In the order they will be indexed in harmonisation input file w_matrix_use variable)

    :type uncertainty_vectors: list:numpy.ndarray
    :param uncertainty_vectors: List of uncertainty vectors
                         (In the order they will be indexed in harmonisation input file uncertainty_vector_use variable)

    :return:
        :w_matrix_val: *numpy.ndarray*

        Concatenated array of non-zero elements of W matrices

        :w_matrix_row: *numpy.ndarray*

        Concatenated array of row indices of W matrices

        :w_matrix_col: *numpy.ndarray*

        Concatenated array of row indices of W matrices

        :w_matrix_nnz: *numpy.ndarray*

        Number of non-zeros elements per W matrix

        :uncertainty_vector_row_count: *numpy.ndarray*

        Number of elements per uncertainty vector

        :uncertainty_vector: *numpy.ndarray*

        Concatenated array of uncertainty vectors
    """

    # 1. Generate W matrix variable arrays

    # > w_matrix_nnz
    w_matrix_nnz = zeros(len(w_matrices))
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

    # 2. Generate uncertainty vector variable arrays

    # > uncertainty_vector_row_count
    uncertainty_vector_row_count = zeros(len(uncertainty_vectors))
    for i, u in enumerate(uncertainty_vectors):
        uncertainty_vector_row_count[i] = len(u)

    # > uncertainty_vector
    total_row = sum(uncertainty_vector_row_count)
    uncertainty_vector = zeros(total_row)

    istart = 0
    iend = 0
    for u, u_row in zip(uncertainty_vectors, uncertainty_vector_row_count):
        iend += u_row
        uncertainty_vector[istart:iend] = u
        istart = iend

    return w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz, uncertainty_vector_row_count, uncertainty_vector


def append_W_to_input_file(filepath,
                           w_matrix_val, w_matrix_row, w_matrix_col, w_matrix_nnz,
                           uncertainty_vector_row_count, uncertainty_vector,
                           w_matrix_use, uncertainty_vector_use):
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

    :type uncertainty_vector_row_count: numpy.ndarray
    :param uncertainty_vector_row_count: Number of elements per uncertainty vector

    :type uncertainty_vector: numpy.ndarray
    :param uncertainty_vector: Concatenated array of uncertainty vectors
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
    w_matrix_num_row_dim = rootgrp.createDimension("w_matrix_num_row", num_row)

    # > w_matrix_sum_nnz - sum of non-zero elements in all W matrices
    w_matrix_sum_nnz_dim = rootgrp.createDimension("w_matrix_sum_nnz", sum(w_matrix_nnz))

    # > uncertainty_vector_count - number of uncertainty vectors
    uncertainty_vector_count_dim = rootgrp.createDimension("uncertainty_vector_count", len(uncertainty_vector_row_count))

    # > uncertainty_vector_sum_row - sum of rows in uncertainty vector
    uncertainty_vector_sum_row_dim = rootgrp.createDimension("uncertainty_vector_sum_row", sum(uncertainty_vector_row_count))

    # 3. Create new variables
    # > w_matrix_nnz - number of non-zero elements for each W matrix
    w_matrix_nnz_var = rootgrp.createVariable('w_matrix_nnz', 'i4', ('w_matrix_count',), zlib=True, complevel=9)
    w_matrix_nnz_var.description = "number of non-zero elements for each W matrix"

    # > w_matrix_row - CSR row numbers for each W matrix
    row_dims = ('w_matrix_num_row',)
    if len(w_matrix_row.shape) > 1:
        row_dims = ('w_matrix_count', 'w_matrix_num_row')
    w_matrix_row_var = rootgrp.createVariable('w_matrix_row', 'i4', row_dims, zlib=True, complevel=9)
    w_matrix_row_var.description = "CSR row numbers for each W matrix"

    # > w_matrix_col - CSR column numbers for all W matrices
    w_matrix_col_var = rootgrp.createVariable('w_matrix_col', 'i4', ('w_matrix_sum_nnz',), zlib=True, complevel=9)
    w_matrix_col_var.description = "CSR column numbers for all W matrices"

    # > w_matrix_val - CSR values for all W matrices
    w_matrix_val_var = rootgrp.createVariable('w_matrix_val', 'f8', ('w_matrix_sum_nnz',), zlib=True, complevel=9)
    w_matrix_val_var.description = "CSR values for all W matrices"

    # > w_matrix_use - a mapping from H array column index to W
    w_matrix_use_var = rootgrp.createVariable('w_matrix_use', 'i4', ('m',), zlib=True, complevel=9)
    w_matrix_use_var.description = "mapping from H array column index to W"

    # > uncertainty_vector_row_count - number of rows of each uncertainty vector
    uncertainty_vector_row_count_var = rootgrp.createVariable('uncertainty_vector_row_count', 'i4', ('uncertainty_vector_count',), zlib=True, complevel=9)
    uncertainty_vector_row_count_var.description = "number of rows of each uncertainty vector"

    # > uncertainty_vector - uncertainty of each scanline value
    uncertainty_vector_var = rootgrp.createVariable('uncertainty_vector', 'f8', ('uncertainty_vector_sum_row',), zlib=True, complevel=9)
    uncertainty_vector_var.description = "uncertainty of each scanline value"

    # > uncertainty_vector_use - a mapping from H array column index to uncertainty vector
    uncertainty_vector_use_var = rootgrp.createVariable('uncertainty_vector_use', 'i4', ('m',), zlib=True, complevel=9)
    uncertainty_vector_use_var.description = "mapping from H array column index to uncertainty vector"

    # 4. Add data
    w_matrix_nnz_var[:] = w_matrix_nnz[:]
    w_matrix_row_var[:] = w_matrix_row[:]
    w_matrix_col_var[:] = w_matrix_col[:]
    w_matrix_val_var[:] = w_matrix_val[:]
    w_matrix_use_var[:] = w_matrix_use[:]

    uncertainty_vector_row_count_var[:] = uncertainty_vector_row_count[:]
    uncertainty_vector_var[:] = uncertainty_vector[:]
    uncertainty_vector_use_var[:] = uncertainty_vector_use[:]

    # 5. Close file
    rootgrp.close()

    return 0

if __name__ == "__main__":
    pass