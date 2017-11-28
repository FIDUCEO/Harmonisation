"""
Module with functions to build a W matrix
"""

'''___Python Modules___'''


'''___Third Party Modules___'''
from numpy import zeros, ones, trim_zeros, asarray, arange
from scipy.sparse import csr_matrix

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "27/11/2017"
__credits__ = ["Ralf Quast", "Jon Mittaz", "Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def generate_rolling_average_w_matrix(matchup_times, scanline_time, kernel_uncertainty):
    """
    Return w matrix in sparse representation and diagonal of u matrix for a rolling average operation, where the
    covariance matrix is w_matrix^T u_matrix^T u_matrix w_matrix

    :type matchup_times: numpy.ndarray
    :param matchup_times: match-up times for match-ups in match-up series

    :type scanline_time: numpy.ndarray
    :param scanline_time: time interval between consecutive scanlines

    :type kernel_uncertainty: numpy.ndarray
    :param kernel_uncertainty: uncertainty per

    :return:
        :w_matrix: *scipy.sparse.csr_matrix*

        W matrix

        :u_matrix: *numpy.ndarray*

        Uncertainty matrix
    """

    # Useful parameters
    num_matchups = int(len(matchup_times))        # number of match_ups
    kernel_size = int(len(kernel_uncertainty[0]))   # length of full averaging kernel

    ####################################################################################################################
    # 1. Determine first column index of each row
    ####################################################################################################################

    # column index of first non-zero element of each w row relative to first match-up
    relative_column_postions = []

    current_column = 0
    for i in xrange(num_matchups):
        scanline_step = 0
        if i > 0:
            time_difference = matchup_times[i] - matchup_times[i-1]
            scanline_step = time_difference/scanline_time
            if scanline_step > kernel_size:
                scanline_step = kernel_size
        current_column += scanline_step
        relative_column_postions.append(int(current_column))

    # Alter relative columns positions so that zero is the lowest value to give absolute column position
    relative_column_postions = asarray(relative_column_postions)
    column_positions = relative_column_postions - min(relative_column_postions)

    ####################################################################################################################
    # 2. Generate sparse matrix indices
    ####################################################################################################################

    # Initialise sparse matrix index and values arrays (maximum size, i.e. no shared scanline values between averages)
    ir = zeros(num_matchups * kernel_size)
    jc = zeros(num_matchups * kernel_size)
    ws = zeros(num_matchups * kernel_size)
    u_matrix = zeros(num_matchups * kernel_size)

    # Populate arrays
    iend = 0
    for i, col_i in enumerate(column_positions):
        ui = kernel_uncertainty[i]  # scanline non-zero uncertainties

        # Fill sparse matrix index and value arrays
        istart = iend
        iend = istart + kernel_size
        ir[istart:iend] = ones(kernel_size) * i
        jc[istart:iend] = arange(kernel_size) + col_i
        ws[istart:iend] = ones(kernel_size) / kernel_size
        u_matrix[int(col_i):int(col_i + kernel_size)] = ui

    # Trim off trailing zeros if maximum size not required, i.e. if sharing of scanline values between averages
    ir = trim_zeros(ir, trim='b')
    jc = trim_zeros(jc, trim='b')
    ws = trim_zeros(ws, trim='b')
    u_matrix = trim_zeros(u_matrix, trim='b')

    # 3. Build sparse matrix
    w_matrix = csr_matrix((ws, (ir, jc)))

    return w_matrix, u_matrix


if __name__ == "__main__":
    pass
