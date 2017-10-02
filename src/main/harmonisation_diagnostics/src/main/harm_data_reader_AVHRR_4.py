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
        times = zeros(cNm[-1])
        times_ref = zeros(cNm[-1])
        corr_data = array([])
        Darray = zeros((cNm[-1], 2 * m))

        for i, dir in enumerate(paths):
            rootgrp = Dataset(dir, 'r')

            istart = cNm[i]
            iend = cNm[i + 1]

            # initialise data arrays for first match-up series file
            ks[istart:iend] = rootgrp.variables['K'][:]
            times[istart:iend] = rootgrp.variables['time_matchup'][:]
            times_ref[istart:iend] = rootgrp.variables['ref_time_matchup'][:]
            corr_data = append(corr_data, rootgrp.variables['corrData'][:])

            # per covariate per match-up data
            Darray[istart:iend, :] = rootgrp.variables['H'][:, :]

            rootgrp.close()

        # open parameter best estimates
        a = loadtxt(path_parameters, delimiter=',')

        # one sensor
        if a.ndim == 1:
            idx['Ia'] = asarray([idx['sensors'][int(j) + 1] for j in zeros(a.shape[0])])

        # multiple sensors
        else:
            idx['Ia'] = asarray([idx['sensors'][j + 1] for j in
                                 vstack([arange(a.shape[0]) for i in range(a.shape[0])]).flatten('F')])

        a = a.flatten()

        unc = None
        unck = None

        return Darray, unc, ks, unck, a, idx, times


if __name__ == "__main__":

    def main():
        return 0

    main()
