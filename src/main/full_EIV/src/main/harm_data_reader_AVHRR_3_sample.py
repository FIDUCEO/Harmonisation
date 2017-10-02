"""
Created on Mon Jan 9  2017 17:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
from numpy import array, zeros, ones, loadtxt, append, mean, delete, vstack, asarray, arange
from netCDF4 import Dataset
from copy import deepcopy

'''___Harmonisation Modules___'''
from harm_data_reader import HarmData as HarmData_template
from correl_forms import CorrelForm


class HarmData(HarmData_template):
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

    def open_data(self, paths, path_parameters):
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

        max_len = 1000  # i.e. all match-ups
        sel_cov = [0, 1, 2, 3, 5, 6, 7, 8] # i.e. ignore temperature columns

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
            lm[i, 2] = max_len
            ms.append(rootgrp.variables['H'][:].shape[1] / 2)
            rootgrp.close()

        lm =lm.astype(int)

        # check all sensor models require same number of covariates
        if len(set(ms)) != 1:
            exit('Sensor Model Mismatch - sensor_model per match-up series much take the same number of variables')
        m = len(sel_cov)/2  # set number of covariates in each sensor model

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
            ks[istart:iend] = rootgrp.variables['K'][:max_len]
            uKarray[istart:iend] = (rootgrp.variables['Kr'][:max_len]**2 + rootgrp.variables['Ks'][:max_len]**2)**0.5
            times[istart:iend] = rootgrp.variables['time_matchup'][:max_len]
            times_ref[istart:iend] = rootgrp.variables['ref_time_matchup'][:max_len]
            corr_data = append(corr_data, rootgrp.variables['corrData'][:])

            # per covariate per match-up data
            Darray[istart:iend, :] = rootgrp.variables['H'][:max_len, sel_cov]
            uRarray[istart:iend, :] = rootgrp.variables['Ur'][:max_len, sel_cov]
            uSarray[istart:iend, :] = rootgrp.variables['Us'][:max_len, sel_cov]

            rootgrp.close()

        # pre-processing checks
        bad_mus = self.find_bad_mus(uSarray, uRarray)
        Darray, uRarray, uSarray, ks, uKarray, idx = self.remove_bad_mus(bad_mus, Darray, uRarray,
                                                                         uSarray, ks, uKarray, idx)

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

            cal_uRarray[istart:iend] = rootgrp.variables['cal_BB_Ur'][:max_len, :]
            ref_cal_uRarray[istart:iend] = rootgrp.variables['ref_cal_BB_Ur'][:max_len, :]
            space_uRarray[istart:iend] = rootgrp.variables['cal_Sp_Ur'][:max_len, :]
            ref_space_uRarray[istart:iend] = rootgrp.variables['ref_cal_Sp_Ur'][:max_len, :]

            rootgrp.close()

        cal_uRarray = delete(cal_uRarray, bad_mus, axis=0)
        ref_cal_uRarray = delete(ref_cal_uRarray, bad_mus, axis=0)
        space_uRarray = delete(space_uRarray, bad_mus, axis=0)
        ref_space_uRarray = delete(ref_space_uRarray, bad_mus, axis=0)
        ################################################################################################################

        ################################################################################################################
        # 3. Reformat uncertainty data
        ################################################################################################################

        # uncertainties for values array
        unc = [0] * len(idx['n_cov'])

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
            elif (col_idx == 3) or (col_idx == 7):
                unc[i] = CorrelForm("r", uRarray[istartm:iendm, col_idx])
            #######################################################################

            # > if C_Space/C_ICT - averaging correlation (location hard coded)
            elif (col_idx == 0) or (col_idx == 4) or (col_idx == 1) or (col_idx == 5):

                if col_idx == 0:
                    uR_w = ref_space_uRarray[istartm:iendm]

                elif col_idx == 4:
                    uR_w = space_uRarray[istartm:iendm]

                elif col_idx == 1:
                    uR_w = ref_cal_uRarray[istartm:iendm]

                elif col_idx == 5:
                    uR_w = cal_uRarray[istartm:iendm]

                corr_data = 25.0
                unc[i] = CorrelForm("ave", (uR_w, times[istartm: iendm], corr_data))
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

        return Darray, unc, ks, unck, a, idx, times


if __name__ == "__main__":

    def main():
        return 0

    main()
