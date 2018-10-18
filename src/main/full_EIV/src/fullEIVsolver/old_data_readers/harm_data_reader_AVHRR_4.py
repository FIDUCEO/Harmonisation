"""
Created on Mon Jan 9  2017 17:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
from netCDF4 import Dataset
from numpy import array, zeros, ones, loadtxt, append, vstack, asarray, arange, trim_zeros
from scipy.sparse import csr_matrix

'''___Harmonisation Modules___'''
main_directory = dirname(dirname(dirname(__file__)))
sys.path.append(main_directory)
from nplcore import MatchUp as MatchUp_template
from nplcore import CorrelForm


class MatchUp(MatchUp_template):
    """
    Class to open and store harmonisation match-up data from file

    :Attributes:
        :values: numpy.ndarray, dtype=float
            harmonisation variables in 1D array (see open_PH method for description of structure)
        :unc: list:CorrelForm
            uncertainty information for harmonisation data
        :ks: numpy.ndarray, dtype=float
            adjustment factors between sensor match-ups orders by match-up series
        :unck: list:CorrelForm
            uncertainty on adjustment factors
        :a: numpy.ndarray, dtype=float
            calibration parameters
        :idx: dict:list
            dictionary of data indices listed by data block, following are provided: n_sensor, n_mu, n_cov, N_var
            and idx (see open_PH method for description of structure)
    """

    def openMatchUp(self, paths, path_parameters, open_uncertainty=True, open_time=True):
        """
        Function to open shortened version simulated AVHRR data produced by Jon Mittaz for harmonisation

        :param paths: list
            list containing the paths of the harmonisation match-up data in netCDF file
        :param path_parameters: str
            path of file containing

        :return:
            :values: numpy.ndarray, dtype=float
                harmonisation variables (excluding calibration parameters) - structure described below
            :unc: list:CorrelForms
                uncertainties associated with values
            :ks: numpy.ndarray
                adjustment factors between sensors in match-ups
            :uncks: list:CorrelForms
                uncertainty on adjustment factors
            :a: numpy.ndarray, dtype=float
                calibration parameters
            :idx: dict:lists
                dictionary describing structure of harmonisation data
            :idx_orig: dict:list
                reference copy of original state of idx

        """

        # Harmonisation requires input match-up data to be restructured into a 1D array structured in blocks
        # of data in the following way:
        #
        # values = [ Rref | X_1(1,2) |...| X_1(L,S-1) | X_1(L ,S) |   ...   | X_m(1,2) |...| X_m(L,S-1)| X_m(L,S) ]
        #
        # where:
        # - Rref - reference sensor measured radiance
        # - X_m - model covariate of number m of a total m covariates
        # - (L, s) - covariate indices:
        #            > L - match-up series index of L of a total L match-up series
        #            > s - sensor number of S sensors
        #
        # Adjustment Factor, K, data is considered separately and stored in a separate 1D array in blocks by
        # match-up series in the following way:
        #
        # ks = [ K(1) |...| K(L)]
        #
        # Calibration parameters are given a0 to ap (p+1 total parameters) per sensor (for S sensors)
        # and are stored in 1D arrays in the following way:
        #
        # a = [ a_0 (s=0),..., a_p (s=0), ... , a_0 (s=S),...,a_p(s=S) ]
        #
        # Each covariate (+ks) in each match-up has its own uncertainty values and type, create a list to contain this
        # information for the consecuative blocks in the values data array.

        ################################################################################################################
        # 1. Build idx dictionary which describes required data structure
        ################################################################################################################

        # a. open match-up series description data

        # initialise match-up info lists
        lm = zeros((len(paths), 3))  # list to contain one array per match-ups series, which each contain data:
        # [sensor_1_name, sensor_2_name, number_of_matchups]
        ms = []  # list to contain number of parameters in sensor model per match-up series

        for i, dir in enumerate(paths):
            rootgrp = Dataset(dir, 'r')
            lm[i, :] = rootgrp.variables['lm'][:][0]
            lm[i,2] = int(lm[i,2])
            ms.append(rootgrp.variables['H'][:].shape[1] / 2)
            rootgrp.close()

        lm =lm.astype(int)

        # check all sensor models require same number of covariates
        if len(set(ms)) != 1:
            exit('Sensor Model Mismatch - sensor_model per match-up series much take the same number of variables')
        m = ms[0]  # set number of covariates in each sensor model

        # b. Obtaining lists of indices to describe data structure

        # i. per match-up series indices
        Nm = []  # number of match-up per match-up series
        cNm = [0]  # cumulative number of match-ups by match-up series
        Im = []  # sensors per match up by ID (starting from 0 and increasing)
        sensors = []  # sensor name corresponding to ID used in Im

        # loop through match-up series description data to populate index lists
        numbers = []  # sensor numbers used
        total = 0  # number of match-ups per match up series cumulative total
        for info in lm:

            # add number of match-ups in match-up series to list
            Nm.append(info[2])

            # calculate cumulative total of match-ups
            total += info[2]
            cNm.append(total)

            # initialise pair for sensor IDs in match-up series
            pair = [0, 0]

            # determine sensor IDs
            for i, sensor in enumerate(info[:2]):
                if sensor in sensors:
                    pair[i] = numbers[sensors.index(sensor)]
                else:
                    sensors.append(sensor)

                    # check if new sensor needing new ID
                    if numbers == []:
                        number = 0
                    else:
                        number = numbers[-1] + 1
                    numbers.append(number)

                    pair[i] = number

            Im.append(pair)

        # ii. per data block indices

        sensor_ms = [1 if s == -1 else m for s in sensors]

        # generate required data
        rs_list = [num for pair in Im for num in pair]  # full list of sensors by match-up series
        n_ref = rs_list.count(0)  # number of reference-sensor match-up series
        s_list = [value for value in rs_list if value != 0]  # list of sensors

        # Create New Indices:
        # > n_sensor - list of sensor number of consecuative data blocks
        n_sensor = [0] * n_ref + s_list * m

        # > n_mu - list of match-up series number of consecuative data blocks
        ref_mus = []
        sensor_mus = []
        for i, pair in enumerate(Im):
            if 0 in pair:
                ref_mus.append(i + 1)
                sensor_mus.append(i + 1)
            else:
                sensor_mus.append(i + 1)
                sensor_mus.append(i + 1)

        n_mu = ref_mus + sensor_mus * m

        # > n_cov - list of covariate number of consecuative blocks of data
        n_cov = [1] * n_ref
        for cov in range(1, m + 1):
            n_cov += [cov] * (2 * len(Nm) - n_ref)

        # - N_var - list of total number of variables in consecuative data blocks
        #          (initially number of match-ups before modification)
        N_var = [Nm[n - 1] for n in n_mu]

        # - idxs - index in 1D harmonisation data array of first element of each data block (and last element index)
        idxs = [0]
        total = 0
        for N in N_var:
            total += N
            idxs.append(int(total))

        # iii. compile indices lists into dictionary
        idx = {"Nm": Nm,
               "cNm": cNm,
               "Im": Im,
               'lm': lm,
               "sensors": sensors,
                "sensor_ms": sensor_ms,
               "n_sensor": n_sensor,
               "n_mu": n_mu,
               "n_cov": n_cov,
               "N_var": N_var,
               "idx": idxs}

        ################################################################################################################
        # 2. Open match-up data files and parameter estimates
        ################################################################################################################

        # initialise arrays
        ks = zeros(cNm[-1])
        uKarray = zeros(cNm[-1])
        times = zeros(cNm[-1])
        times_ref = zeros(cNm[-1])
        corr_data = array([])
        Darray = zeros((cNm[-1], 2 * m))
        uRarray = zeros((cNm[-1], 2 * m))
        uSarray = zeros((cNm[-1], 2 * m))

        for i, dir in enumerate(paths):
            rootgrp = Dataset(dir, 'r')

            istart = cNm[i]
            iend = cNm[i + 1]

            # initialise data arrays for first match-up series file
            ks[istart:iend] = rootgrp.variables['K'][:]
            uKarray[istart:iend] = (rootgrp.variables['Kr'][:]**2 + rootgrp.variables['Ks'][:]**2)**0.5
            times[istart:iend] = rootgrp.variables['time_matchup'][:]
            times_ref[istart:iend] = rootgrp.variables['ref_time_matchup'][:]
            corr_data = append(corr_data, rootgrp.variables['corrData'][:])

            # per covariate per match-up data
            Darray[istart:iend, :] = rootgrp.variables['H'][:, :]
            uRarray[istart:iend, :] = rootgrp.variables['Ur'][:, :]
            uSarray[istart:iend, :] = rootgrp.variables['Us'][:, :]

            rootgrp.close()

        # open parameter best estimates
        a = loadtxt(path_parameters, delimiter=',')

        # one sensor
        if a.ndim == 1:
            idx['Ia'] = asarray([idx['sensors'][int(j) + 1] for j in zeros(a.shape[0])])

        # multiple sensors
        else:
            idx['Ia'] = asarray([idx['sensors'][j + 1] for j in
                                 vstack([arange(a.shape[0]) for i in range(a.shape[1])]).flatten('F')])

        a = a.flatten()

        ######################### AHVRR SPECIFIC CODE ##################################################################
        n_ws = []

        # Open data match-up series file by match-up series file
        for i, dir in enumerate(paths):
            rootgrp = Dataset(dir, 'r')
            n_ws.append(rootgrp.getncattr('Calibraton_Average_No_Scanline'))
            rootgrp.close()
        n_w = max(n_ws)

        cal_uRarray = zeros((cNm[-1], n_w))
        ref_cal_uRarray = zeros((cNm[-1], n_w))
        space_uRarray = zeros((cNm[-1], n_w))
        ref_space_uRarray = zeros((cNm[-1], n_w))

        for i, dir in enumerate(paths):
            rootgrp = Dataset(dir, 'r')
            istart = cNm[i]
            iend = cNm[i + 1]

            cal_uRarray[istart:iend] = rootgrp.variables['cal_BB_Ur'][:, :]
            ref_cal_uRarray[istart:iend] = rootgrp.variables['ref_cal_BB_Ur'][:, :]
            space_uRarray[istart:iend] = rootgrp.variables['cal_Sp_Ur'][:, :]
            ref_space_uRarray[istart:iend] = rootgrp.variables['ref_cal_Sp_Ur'][:, :]

            rootgrp.close()

        ################################################################################################################

        ################################################################################################################
        # 3. Reformat uncertainty data
        ################################################################################################################

        # uncertainties for values array
        unc = [0] * len(idx['n_cov'])
        w_matrices = []
        uncertainty_vectors = []

        for i, (cov, mu) in enumerate(zip(idx['n_cov'], idx['n_mu'])):

            # indices defining first and last positions in data matrix
            istartm = idx['cNm'][mu - 1]
            iendm = idx['cNm'][mu]

            # index defining column in data matrix
            if idx['Im'][mu - 1][0] == idx['n_sensor'][i]:  # if the sensor is the first sensor in the match-up series
                col_idx = cov - 1
            if idx['Im'][mu - 1][1] == idx['n_sensor'][i]:  # if the sensor is the second sensor in the match-up series
                col_idx = cov + m - 1

            # assign CorrelForm object per block of data depend on type of correlation

            # > if reference sensor - correlation form random
            if idx['n_sensor'][i] == 0:
                unc[i] = CorrelForm("r", uRarray[istartm:iendm, col_idx])

            ######################### AHVRR SPECIFIC CODE ##############################################################

            ######################### NO A3 SPECIFIC CODE #########################
            elif (col_idx == 4) or (col_idx == 9):
                unc[i] = CorrelForm("r", uRarray[istartm:iendm, col_idx])
            #######################################################################

            # > if C_Space/C_ICT - averaging correlation (location hard coded)
            elif (col_idx == 0) or (col_idx == 5) or (col_idx == 1) or (col_idx == 6):

                if col_idx == 0:
                    uR_w = ref_space_uRarray[istartm:iendm]

                elif col_idx == 5:
                    uR_w = space_uRarray[istartm:iendm]

                elif col_idx == 1:
                    uR_w = ref_cal_uRarray[istartm:iendm]

                elif col_idx == 6:
                    uR_w = cal_uRarray[istartm:iendm]

                corr_data = 25.0
                W = self.calc_W(uR_w, times[istartm:iendm], corr_data)
                w_matrices.append(W)
                uncertainty_vectors.append(ones(W.shape[1]))

                unc[i] = CorrelForm("ave", (len(w_matrices)-1, len(uncertainty_vectors)-1))
            ############################################################################################################

            # > random effect - if systematic component is zero
            elif uSarray[istartm, col_idx] == 0:
                unc[i] = CorrelForm("r", uRarray[istartm:iendm, col_idx])

            # > systematic effect
            elif uSarray[istartm, col_idx] != 0:
                unc[i] = CorrelForm("rs", (uRarray[istartm:iendm, col_idx],
                                           sum(uSarray[istartm:iendm, col_idx] ** 2) ** 0.5 / (iendm - istartm)**0.5))

        # uncertainties for ks array
        unck = [0] * len(idx['Nm'])
        for i in range(len(idx['Nm'])):
            istart = idx['cNm'][i]
            iend = idx['cNm'][i + 1]

            unck[i] = CorrelForm('r', uKarray[istart:iend])

        # Flatten values into required 1d form
        values = self.flatten_values(Darray, idx)

        return values, unc, w_matrices, uncertainty_vectors, ks, unck, a, idx

    def flatten_values(self, values, idx):
        """
        Return 1d form of 2d match-up data array and descriptive dictionary of indices

        :param values: numpy.ndarray
            2d harmonisation match-up data array
        :param idx: dict
            dict of lists describing the structure of values

        :return:
            :values_flat: numpy.ndarray
                1d harmonisation match-up data array
        """

        m = values.shape[1]/2

        # initialise 1D numpy array of size of all variables + reference radiance data
        values_flat = zeros(idx['idx'][-1])

        # fill values array with from Darray according to structure specified in idx dict
        for i in range(len(idx['n_mu'])):

            # indices defining first and last positions in data matrix

            istartm = idx['cNm'][idx['n_mu'][i] - 1]
            iendm = idx['cNm'][idx['n_mu'][i]]

            # indices defining first and last positions in data vector
            istartv = idx['idx'][i]
            iendv = idx['idx'][i + 1]

            # if the sensor is the first sensor in the match-up series
            if idx['Im'][idx['n_mu'][i] - 1][0] == idx['n_sensor'][i]:
                values_flat[istartv:iendv] = values[istartm:iendm, idx['n_cov'][i] - 1]

            # if the sensor is the second sensor in the match-up series:
            if idx['Im'][idx['n_mu'][i] - 1][1] == idx['n_sensor'][i]:
                values_flat[istartv:iendv] = values[istartm:iendm, idx['n_cov'][i] + m - 1]

        return values_flat

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
