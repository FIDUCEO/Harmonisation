"""
Created on Day Month Date  Year Hour:Minute:Second

@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules____'''
from numpy import zeros, ones, arange, trim_zeros
from scipy.sparse import csr_matrix


class CorrelForm:
    """
    Allow users to create an object containing information of a covariates uncertainty values and form

    :Attributes:
        :form: str
            specifics correlation form of uncertainty, can have values:
            > "r" - random form
            > "rs" - random + systematic form
            > "ave" - averaging form

        + additional attribute/s containing uncertainty data depending on correlation form:

        > "r" form
            :uR: float
                where uR is the standard uncertainty per match-up

        > "rs" form
            :uR: float
                where uR is the standard uncertainty per match-up
            :uS: float
                where uS is single valued standard uncertainty per match-up caused by systematic effects

        > "ave" form
            :uR: float
                where uR is the array of n_w random uncertainties per match-up, where n_w is the width of the widest
                averaging kernel for that match-up series. Scanlines in this window not used in the a particular
                match-up average should be valued 0
            :W: scipy.sparse.csr_matrix
                averaging operator matrix

    """

    def __init__(self, form, data_tuple):
        """
        Take user input covariate uncertainty correlation form and data and apply as attributes of the class

        :param form: str
            user specifies correlation form of covariate uncertainty (specified in class description)
        :param data_tuple:
            user specifies uncertainty values, structure depends on correlation form, as follows,
            > "r" - uR
                :uR: float
                    where uR is the standard uncertainty per match-up
            > "rs" - (uR, uS)
                :uR: float
                    where uR is the standard uncertainty per match-up
                :uS: float
                    where uS is single valued standard uncertainty per match-up caused by systematic effects
            > "ave" - (uR, times, corrData)
                :uR: float
                    where uR is the array of n_w random uncertainties per match-up, where n_w is the width of the widest
                    averaging kernel for that match-up series. Scanlines in this window not used in the a particular
                    match-up average should be valued 0
                :times: numpy.ndarray
                    match-up times for match-ups in match-up series
                :corrData: numpy.ndarray
                    match-up time data

        """

        # set form attribute
        self.form = form

        # set data value attributes from input data_tuple
        if form == "r":
            self.uR = data_tuple

        elif form == 'rs':
            self.uR = data_tuple[0]
            self.uS = data_tuple[1]

        elif form == 'ave':
            self.uR = data_tuple[0]
            self.W = self.calc_W(data_tuple[0], data_tuple[1], data_tuple[2])

    def calc_W(self, u, times, corrData):
        """
        Return weighting matrix in sparse representation.

        :type u: numpy.ndarray
        :param u: standard uncertainties

        :type times: numpy.ndarray
        :param times: match-up times for match-ups in match-up series

        :type corrData: numpy.ndarray
        :param corrData: match-up time data

        :return:
            :W: *scipy.sparse.csr_matrix*

            weighting matrix
        """

        n_var = len(times)    # number of match_ups
        N_W = len(u[0])       # length of full averaging kernel

        # initialise sparse matrix index and values arrays (of maximum size, i.e. if all windows are N_W)
        ir = zeros(n_var * N_W)
        jc = zeros(n_var * N_W)
        ws = zeros(n_var * N_W)

        col = 0  # column number
        iend = 0
        for i in xrange(n_var):
            ui = u[i][u[i] != 0]      # scanline non-zero uncertainties
            n_w = len(ui)             # width of match-up specific averaging window

            # find col_step of match-up compared to last match-up (if not first match-up)
            col_step = 0
            if i > 0:
                corr_val = self.return_correlation(times, corrData, i, i - 1)
                col_step = int(round(n_w * (1 - corr_val)))
            col += col_step

            # fill sparse matrix index and value arrays
            istart = iend
            iend = istart + n_w

            ir[istart:iend] = ones(n_w) * i
            jc[istart:iend] = arange(n_w) + col
            ws[istart:iend] = ui * ones(n_w)/n_w

        # trim off trailing zeros if maximum size not required, i.e. if all windows are not N_W in length
        ir = trim_zeros(ir, trim='b')
        jc = trim_zeros(jc, trim='b')
        ws = trim_zeros(ws, trim='b')

        # build sparse matrix
        W = csr_matrix((ws, (ir, jc)))

        return W

    def return_correlation(self, times, width, i_1, i_2):
        """
        Function to give error correlation between scan lines

        :type times: numpy.ndarray
        :param times: CorrIndexArray data from file (this name)

        :type width: numpy.ndarray
        :param width: time width of averaging window corrData auxiliary data from file (this name)

        :type i_1: int
        :param i_1: index in times of central scanline of averaging

        :type i_1: int
        :param i_2: index in times of outer scanline of interest

        :return
            :corr_val: *float*

            correlation between scanlines
        """

        diff = abs(times[i_1] - times[i_2])
        if diff > width:
            return 0.
        else:
            return 1. - (diff / width)

if __name__ == "__main__":

    def main():
        return 0

    main()
