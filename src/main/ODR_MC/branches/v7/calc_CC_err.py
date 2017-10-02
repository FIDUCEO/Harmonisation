"""
Created 2017/02/22

@author: Sam Hunt, NPL\ENV
"""

"""___Python Modules___"""
from scipy.sparse import csr_matrix
from numpy import zeros, ones, trim_zeros, arange
from netCDF4 import Dataset
import numpy.random as random


def calc_CC_err(u, times, corrData):
    """
    Return weighting matrix: sparse representation.

    :param u: float
        standard uncertainties
    :param times: numpy.ndarray
        match-up times for match-ups in match-up series
    :param corrData: numpy.ndarray
        match-up time data

    :return:
        :C_ICT_raw_err: numpy.ndarray
            error for raw scanline level C_ICT counts, uncertainty standardised to 1
        :C_ICT_err: numpy.ndarray
            error for averaged C_ICT counts
    """

    n_var = len(times)  # number of match_ups
    N_W = len(u[0])     # length of maximum averaging kernel

    # initialise sparse matrix index and values arrays (of maximum size, i.e. if all windows are n_w)
    ir = zeros(n_var * N_W)
    jc = zeros(n_var * N_W)
    ws = zeros(n_var * N_W)

    col = 0  # column number
    iend = 0
    for i in xrange(n_var):
        ui = u[i][u[i] != 0]  # scanline non-zero uncertainties
        n_w = len(ui)  # width of match-up specific averaging window

        # find col_step of match-up compared to last match-up (if not first match-up)
        col_step = 0
        if i > 0:
            corr_val = return_correlation(times, corrData, i, i - 1)
            col_step = int(round(n_w * (1 - corr_val)))
        col += col_step

        # fill sparse matrix index and value arrays
        istart = iend
        iend = istart + n_w

        ir[istart:iend] = ones(n_w) * i
        jc[istart:iend] = arange(n_w) + col
        ws[istart:iend] = ui * ones(n_w) / n_w

    # trim off trailing zeros if maximum size not required, i.e. if all windows are not n_w in length
    ir = trim_zeros(ir, trim='b')
    jc = trim_zeros(jc, trim='b')
    ws = trim_zeros(ws, trim='b')

    # build sparse matrix
    W = csr_matrix((ws, (ir, jc)))

    # generate raw scanline errors (uncertainy normalised to 1)
    C_ICT_raw_err = random.normal(loc=zeros(W.indices[-1]+1))

    # average raw errors to generate C_ICT_err (scaling by raw C_ICT uncertainty)
    C_ICT_err = W.dot(C_ICT_raw_err)

    return C_ICT_err, C_ICT_raw_err


def return_correlation(index, corr_array, cent_pos, req_pos):
    """
    Function to give error correlation between scan lines
    Operated on a case by case basis (does not assume inputs are arrays)

    :param index: numpy.ndarray
        CorrIndexArray data from file (this name)
    :param corr_array: numpy.ndarray
        corrData auxiliary data from file (this name)
    :param cent_pos: int
        central scanline of averaging
    :param req_pos: int
        outer scanline of interest

    :return
        :corr_val:
            correlation between scanlines

    """

    diff = abs(index[cent_pos] - index[req_pos])
    if diff > corr_array[0]:
        return 0.
    else:
        return 1. - (diff / corr_array[0])


if __name__ == "__main__":

    def code_sample():
        """
        Sample code for using calc_C_ICT_err to generate averaged C_ICT errors
        """

        # open data file
        dirMU = "D:\Projects\FIDUCEO\Data\Simulated\m02_n15.nc"
        rootgrp = Dataset(dirMU, 'r')

        # access required data (first 100 match-ups for this example)
        u_raw = rootgrp.variables['cal_BB_Ur'][:100, :]
        times = rootgrp.variables['CorrIndexArray'][:100]
        corrData = rootgrp.variables['corrData'][:100]

        rootgrp.close()

        # run function
        C_ICT_err, C_ICT_raw_err = calc_CC_err(u_raw, times, corrData)

        return 0

    code_sample()

