"""
Created on Mon Jan 9  2017 17:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
from numpy import array, zeros, loadtxt, append, delete, isnan, vstack, arange, asarray
from datetime import datetime
from netCDF4 import Dataset
from copy import deepcopy
from sys import exit
from os.path import join as pjoin

'''___Harmonisation Modules___'''
from correl_forms import CorrelForm


class HarmData:
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

    def __init__(self, path=None, path_parameters=None, sensor_model=None, adjustment_model=None, flatten=True):
        """
        Initialise harmonisation data object, opening data from directory if specified

        :param directory: str
            directory of harmonisation data to be opened
        """

        # initialise attributes
        self.values = array([])
        self.unc = array([])
        self.ks = array([])
        self.unck = array([])
        self.a = array([])
        self.idx = {}
        self.sensor_model = None
        self.adjustment_model = None

        # open data
        if path is not None:
            if path_parameters is not None:
                if sensor_model is not None:
                    if adjustment_model is not None:
                        self.sensor_model = sensor_model
                        self.adjustment_model = adjustment_model
                        self.values, self.unc, self.ks, self.unck,\
                            self.a, self.idx, times = self.open_data(path, path_parameters)

                        self.times = self.seconds2date(times)

                        # save separate copy of original indices for future reference
                        self.idx_orig = deepcopy(self.idx)

                    else:
                        exit('Missing Parameter - adjustment_model missing from HarmData')
                else:
                    exit('Missing Parameter - sensor_model missing from HarmData')
            else:
                exit('Missing Parameter - parameter_path missing from HarmData')

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

        # open match-up info from match-up data files

        for i, dir in enumerate(paths):
            rootgrp = Dataset(dir, 'r')
            lm[i, :] = rootgrp.variables['lm'][:][0]
            lm[i,2] = int(lm[i,2])
            ms.append(rootgrp.variables['H'][:].shape[1] / 2)
            rootgrp.close()

        lm = lm.astype(int)

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

            uKarray[istart:iend] = (rootgrp.variables['Kr'][:]**2+rootgrp.variables['Ks'][:]**2)**0.5

            # initialise data arrays for first match-up series file
            ks[istart:iend] = rootgrp.variables['K'][:]
            times[istart:iend] = rootgrp.variables['time_matchup'][:]
            times_ref[istart:iend] = rootgrp.variables['ref_time_matchup'][:]
            corr_data = append(corr_data, rootgrp.variables['corrData'][:])

            # per covariate per match-up data
            Darray[istart:iend, :] = rootgrp.variables['H'][:]
            uRarray[istart:iend, :] = rootgrp.variables['Ur'][:]
            uSarray[istart:iend, :] = rootgrp.variables['Us'][:]

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
                                 vstack([arange(a.shape[0]) for i in range(a.shape[0])]).flatten('F')])

        a = a.flatten()

        ################################################################################################################
        # 4. Reformat uncertainty data
        ################################################################################################################

        # uncertainties for values array
        unc = [0] * len(idx['n_cov'])

        for i, (cov, mu) in enumerate(zip(idx['n_cov'], idx['n_mu'])):

            # indices defining first and last positions in data matrix
            istartm = idx['cNm'][mu - 1]
            iendm = idx['cNm'][mu]
            ib = idx['idx'][i]  # start idx of data block
            ie = idx['idx'][i + 1]  # end idx of data block

            # index defining column in data matrix
            if idx['Im'][mu - 1][0] == idx['n_sensor'][i]:  # if the sensor is the first sensor in the match-up series
                col_idx = cov - 1
            if idx['Im'][mu - 1][1] == idx['n_sensor'][i]:  # if the sensor is the second sensor in the match-up series
                col_idx = cov + m - 1

            # assign CorrelForm object per block of data depend on type of correlation

            # > if reference sensor - correlation form random
            if idx['n_sensor'][i] == 0:
                unc[i] = CorrelForm("r", uRarray[istartm:iendm, col_idx])

            # > random effect - if systematic component is zero
            elif uSarray[istartm, col_idx] == 0:
                unc[i] = CorrelForm("r", uRarray[istartm:iendm, col_idx])

            # > systematic effect
            elif uSarray[istartm, col_idx] != 0:
                unc[i] = CorrelForm("rs", (uRarray[istartm:iendm, col_idx],
                                           sum(uSarray[istartm:iendm, col_idx] ** 2)**0.5 / (iendm - istartm)**0.5))

        # uncertainties for ks array
        unck = [0] * len(idx['Nm'])
        for i in range(len(idx['Nm'])):
            istart = idx['cNm'][i]
            iend = idx['cNm'][i + 1]

            unck[i] = CorrelForm('r', uKarray[istart:iend])

        return Darray, unc, ks, unck, a, idx, times

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

    def unflatten_values(self, values_flat, idx):
        """
        Return 1d form of 2d match-up data array and descriptive dictionary of indices

        :param values_flat: numpy.ndarray
            1d harmonisation match-up data array (unconverted)
        :param idx_flat: dict
            dict of lists describing the structure of values (original)

        :return:
            :values: numpy.ndarray
                2d harmonisation match-up data array
        """

        # required parameters
        N_cov = idx['n_cov'][-1]  # total number of covariates
        N_mu = idx['cNm'][-1]

        # Reformat into data arrays
        values = zeros((N_mu, 2 * N_cov))

        for k in range(len(idx['n_mu'])):

            # indices defining first and last positions in data matrix
            istartm = idx['cNm'][idx['n_mu'][k] - 1]
            iendm = idx['cNm'][idx['n_mu'][k]]

            # indices defining first and last positions in data vector
            istartv = idx['idx'][k]
            iendv = idx['idx'][k + 1]

            # if the sensor is the first sensor in the match-up series
            if idx['Im'][idx['n_mu'][k] - 1][0] == idx['n_sensor'][k]:
                values[istartm:iendm, idx['n_cov'][k] - 1] = values_flat[istartv:iendv]

            # if the sensor is the second sensor in the match-up series:
            if idx['Im'][idx['n_mu'][k] - 1][1] == idx['n_sensor'][k]:
                values[istartm:iendm, idx['n_cov'][k] + N_cov - 1] = values_flat[istartv:iendv]

        return values

    def find_bad_mus(self, uSarray, uRarray):
        """
        Find indices of match-ups in combined array of match-up series data that contains invalid values

        :param uRarray: numpy.ndarray
            combined match-up series covariate random uncertainty data
        :param uSarray: numpy.ndarray
            combined match-up series covariate systematic uncertainty data

        :return:
            :bad_mus: list
                list of indices of invalid match-ups in combined match-up series data
        """

        bad_mus = []

        #check for negative uncertainties
        for i, (uR, uS) in enumerate(zip(uRarray, uSarray)):
            #filter nans
            uR = uR[~isnan(uR)]
            uS = uS[~isnan(uS)]
            if any(uR < 0) or any(uS < 0):
                bad_mus.append(i)

        return bad_mus

    def remove_bad_mus(self, bad_mus, Darray, uRarray, uSarray, ks, uKarray, idx):
        """
        Remove invalid match-ups from data

        :param bad_mus: numpy.ndarray
            list of indices of invalid match-ups in combined match-up series data
        :param Darray: numpy.ndarray
            combined match-up series covariate data
        :param uRarray: numpy.ndarray
            combined match-up series covariate random uncertainty data
        :param uSarray: numpy.ndarray
            combined match-up series covariate systematic uncertainty data
        :param ks: numpy.ndarray
            combined match-up series sensor adjustment factor data
        :param uKarray: numpy.ndarray
            combined match-up series sensor adjustment factor random uncertainty data
        :param idx: dict:lists
            dictionary describing structure of harmonisation data

        :return:
            :Darray: numpy.ndarray
                combined match-up series covariate data with invalid match-ups removed
            :uRarray: numpy.ndarray
                combined match-up series covariate random uncertainty data with invalid match-ups removed
            :uSarray: numpy.ndarray
                combined match-up series covariate systematic uncertainty data with invalid match-ups removed
            :ks: numpy.ndarray
                combined match-up series sensor adjustment factor data with invalid match-ups removed
            :uKarray: numpy.ndarray
                combined match-up series sensor adjustment factor random uncertainty data with invalid match-ups removed
            :idx: dict:lists
                dictionary describing structure of harmonisation data adjusted to reflect new structure
        """

        if bad_mus != []:

            # a. remove bad data from arrays
            Darray = delete(Darray, bad_mus, axis=0)
            uRarray = delete(uRarray, bad_mus, axis=0)
            uSarray = delete(uSarray, bad_mus, axis=0)
            ks = delete(ks, bad_mus)
            uKarray = delete(uKarray, bad_mus)

            # b. update idx

            # find which match-up series bad match-ups belong to
            n_bad_mus = zeros(len(idx['Nm']))

            for i_bad in bad_mus:
                for j in xrange(len(idx['Nm'])):
                    if i_bad <= idx['cNm'][j+1]:
                        n_bad_mus[j] += 1
                        break

            # update Nm
            for i in xrange(len(idx['Nm'])):
                idx['Nm'][i] = int(idx['Nm'][i] - n_bad_mus[i])

            # update cNm
            total = 0
            cNm = [0]
            for N in idx['Nm']:
                total += N
                cNm.append(total)
            idx['cNm'] = cNm

            # update N_var
            N_var = [idx['Nm'][n - 1] for n in idx['n_mu']]
            idx['N_var'] = N_var

            # update idx
            idxs = [0]
            total = 0
            for N in idx['N_var']:
                total += N
                idxs.append(int(total))
            idx['idx'] = idxs

        return Darray, uRarray, uSarray, ks, uKarray, idx

    def seconds2date(self, times):
        """
        Return matchup times in datetime format

        :param times: numpy.ndarray: float
            Times in seconds since 1/1/1970

        :return:
            :dates: numpy.ndarray: datetime.datetime
                Times in datetime formate
        """

        dates = zeros(times.shape[0]).astype(datetime)
        for i, time in enumerate(times):
            dates[i] = datetime.fromtimestamp(time)

        return dates

    def save_data(self, path, a, Va, f, values_est):
        """
        Function to save harmonisation data to netCDF file at path

        :param path: str
            path of file to store data to
        :param a: numpy.ndarray
            parameter estimates for GNAlgo
        :param Va: numpy.ndarray
            covariance matrices for parameter estimates
        :param f:
            final value of f from Gauss Newton Algorithm
        :param values_est:
            variable estimates from Gauss Newton Algorithm
        """

        N_cov = self.idx['n_cov'][-1]
        N_mu = self.idx['cNm'][-1]
        N_L = len(self.idx['Nm'])
        N_var = self.idx['idx'][-1]
        N_a = len(self.a)

        nl = 3

        # Reformat into data arrays
        values_est_array = zeros((N_mu, 2*N_cov))

        for k in range(len(self.idx['n_mu'])):

            # indices defining first and last positions in data matrix
            istartm = self.idx['cNm'][self.idx['n_mu'][k]-1]
            iendm = self.idx['cNm'][self.idx['n_mu'][k]]

            # indices defining first and last positions in data vector
            istartv = self.idx_orig['idx'][k]
            iendv = self.idx_orig['idx'][k+1]

            # if the sensor is the first sensor in the match-up series
            if self.idx['Im'][self.idx['n_mu'][k]-1][0] == self.idx['n_sensor'][k]:
                values_est_array[istartm:iendm, self.idx['n_cov'][k] - 1] = values_est[istartv:iendv]

            # if the sensor is the second sensor in the match-up series:
            if self.idx['Im'][self.idx['n_mu'][k]-1][1] == self.idx['n_sensor'][k]:
                values_est_array[istartm:iendm, self.idx['n_cov'][k]+N_cov-1] = values_est[istartv:iendv]

        # make lm
        lm = zeros((N_L, nl))
        for i in range(len(lm)):
            lm[i, 0] = self.idx['Nm'][i]
            lm[i, 1] = self.idx['sensors'][self.idx['Im'][i][0]]
            lm[i, 2] = self.idx['sensors'][self.idx['Im'][i][1]]

        # save results
        rootgrp = Dataset(path, 'w')

        # make dimensions
        N_p_dim = rootgrp.createDimension('N_p_dim', a.shape[0])
        N_mu_dim = rootgrp.createDimension('N_mu_dim', N_mu)
        xyza_dim = rootgrp.createDimension('xyza_dim', len(f))
        N_sensors_dim = rootgrp.createDimension('N_sensors_dim', a.shape[1])
        col_dim = rootgrp.createDimension('col_dim', 2*N_cov)
        L_dim = rootgrp.createDimension('L_dim', N_L)
        nl_dim = rootgrp.createDimension('nl_dim', nl)

        lm_var = rootgrp.createVariable('lm', 'f4', ('L_dim', 'nl_dim',))
        lm_var.Description = 'lm variable (L,nl). Stores satellite pairs with number of entries'
        a_var = rootgrp.createVariable('a', 'f4', ('N_p_dim', 'N_sensors_dim',))
        a_var.Description = 'Parameter estimates'
        Va_var = rootgrp.createVariable('Va', 'f4', ('N_p_dim', 'N_p_dim', 'N_sensors_dim',))
        Va_var.Description = 'Covariance matrices for parameter estimates'
        f_var = rootgrp.createVariable('f', 'f4', ('xyza_dim',))
        f_var.Description = 'Final value of f from Gauss Newton Algorithm'
        H_est_var = rootgrp.createVariable('H_est', 'f4', ('N_mu_dim', 'col_dim',))
        H_est_var.Description = 'Residuals of the variable estimates'

        lm_var[:, :] = lm[:, :]
        a_var[:, :] = a[:, :]
        Va_var[:, :, :] = Va[:, :, :]
        f_var[:] = f[:]
        H_est_var[:, :] = values_est_array[:, :]

        rootgrp.close()

        return 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        attrs = vars(self).keys()
        for attr in attrs:
            self.__delattr__(attr)

        self.closed = True


if __name__ == "__main__":

    def main():
        return 0

    main()
