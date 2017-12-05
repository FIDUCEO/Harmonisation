"""
Module with functions to build a W matrix
"""

'''___Python Modules___'''


'''___Third Party Modules___'''
from numpy import zeros, ones, trim_zeros, asarray, arange, unique, int32
from scipy.sparse import csr_matrix

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "27/11/2017"
__credits__ = ["Ralf Quast", "Jon Mittaz", "Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


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

    for i in xrange(num_matchups):
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

    return w_matrix, u_matrix


if __name__ == "__main__":
    pass
