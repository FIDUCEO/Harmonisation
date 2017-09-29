"""
Check harmonisation NPL files for internal consistancy

Usage:
python add_W_to_matchup_file.py path/to/matchup/file.nc
"""

'''___Python Modules___'''
from numpy import zeros

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "28/09/2017"
__credits__ = ["Ralf Quast", "Jon Mittaz", "Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

def write_matchup_file(file_path, Ws, H_matrix_column_W_ID, Us, H_matrix_U_ID):
    """
    Update given file with W matrix variables

    :type file_path: str
    :param file_path: match-up file path

    :type w_matrix_val: numpy.ndarray
    :param w_matrix_val: Concatenated array of non-zero elements of W matrices

    :type w_matrix_row: numpy.ndarray
    :param w_matrix_row: Concatenated array of row indices of W matrices

    :type w_matrix_col: numpy.ndarray
    :param w_matrix_col: Concatenated array of row indices of W matrices

    :type w_matrix_nnz: numpy.ndarray
    :param w_matrix_nnz: Number of non-zeros elements per W matrix

    :type w_matrix_use: numpy.ndarray
    :param w_matrix_use: Array mapping H array columns to w matrices

    :type uncertainty_vector_row_count: numpy.ndarray
    :param uncertainty_vector_row_count: Number of elements per uncertainty vector

    :type uncertainty_vector: numpy.ndarray
    :param uncertainty_vector: Concatenated array of uncertainty vectors

    :type uncertainty_vector_use: numpy.ndarray
    :param uncertainty_vector_use: Array mapping H array columns to uncertainty vectors
    """

    ####################################################################################################################
    # 1. Generate W matrix input variable
    ####################################################################################################################

    w_matrix_nnz = zeros(len(Ws))
    for i, W in enumerate(Ws):
        w_matrix_nnz[i] = W.nnz

    total_nnz = sum(w_matrix_nnz)

    w_matrix_values = zeros(total_nnz)
    w_matrix_col_

    istart = 0
    iend = 0
    for W, W_nnz in zip(Ws, w_matrix_nnz):
        iend += W_nnz
        w_matrix_nnz[istart:iend] = W.values
        w_matrix_col[istart:iend] = W.indices
        istart = iend

    # 1. Open file
    rootgrp = Dataset(file_path, mode='w')

    # 2. Create dimensions
    # > L - number of matchup series
    L_dim = rootgrp.createDimension("L", 1)

    # > nl - width of lm
    nl_dim = rootgrp.createDimension("nl", 3)

    # > M - number of matchups
    M_dim = rootgrp.createDimension("M", H.shape[0])

    # > m - number of columns in H
    m_dim = rootgrp.createDimension("m", H.shape[1])

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
    # > lm - Stores satellite pairs with number of entries
    lm_var = rootgrp.createVariable('lm', 'i4', ('L', 'nl'), zlib=True, complevel=9)
    lm_var.description = "lm variable (L,nl). Stores satellite pairs with number of entries"

    # > H - Radiances and counts per matchup
    H_var = rootgrp.createVariable('H', 'f8', ('M', 'm'), zlib=True, complevel=9)
    H_var.description = "H array (M,m). Radiances and counts per matchup"

    # > Us - Systematic uncertainties for H array
    Us_var = rootgrp.createVariable('Us', 'f8', ('M', 'm'), zlib=True, complevel=9)
    Us_var.description = "Systematic uncertainties for H array"

    # > Ur - Random uncertainties for H array
    Ur_var = rootgrp.createVariable('Ur', 'f8', ('M', 'm'), zlib=True, complevel=9)
    Ur_var.description = "Random uncertainties for H array"

    # > K - K (sensor-to-sensor differences) for zero shift case
    K_var = rootgrp.createVariable('K', 'f8', ('M',), zlib=True, complevel=9)
    K_var.description = "K (sensor-to-sensor differences) for zero shift case"

    # > Kr - K (sensor-to-sensor differences) random uncertainties (matchup uncertainty)
    Kr_var = rootgrp.createVariable('Kr', 'f8', ('M',), zlib=True, complevel=9)
    Kr_var.description = "K (sensor-to-sensor differences) random uncertainties (matchup uncertainty)"

    # > Ks - K (sensor-to-sensor differences) systematic uncertainties for zero shift case
    Ks_var = rootgrp.createVariable('Ks', 'f8', ('M',), zlib=True, complevel=9)
    Ks_var.description = "K (sensor-to-sensor differences) systematic uncertainties for zero shift case"

    # > time_matchup - Matchup time
    time_matchup_var = rootgrp.createVariable('time_matchup', 'f8', ('M',), zlib=True, complevel=9)
    time_matchup_var.description = "Matchup time"
    time_matchup_var.unit = "seconds since 1970-01-01"

    # > ref_time_matchup - Reference matchup time
    ref_time_matchup_var = rootgrp.createVariable('ref_time_matchup', 'f8', ('M',), zlib=True, complevel=9)
    ref_time_matchup_var.description = "Reference matchup time"
    ref_time_matchup_var.unit = "seconds since 1970-01-01"

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
    lm_var[:] = lm[:]
    H_var[:] = H[:]
    Ur_var[:] = Ur[:]
    Us_var[:] = Us[:]
    K_var[:] = K[:]
    Ks_var[:] = Ks[:]
    Kr_var[:] = Kr[:]
    time_matchup_var[:] = time_matchup[:]
    ref_time_matchup[:] = ref_time_matchup[:]

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
