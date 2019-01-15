"""
MatchUp Data Class
"""

'''___Built-In Modules___'''
from copy import deepcopy
from datetime import datetime
from os import listdir
from os.path import join as pjoin
from os.path import split, splitext, basename, isdir, isfile
import sys
from importlib import import_module

'''___Third-Party Modules___'''
from numpy import array, zeros, append, asarray, cumsum, int32, ones, nan, float32
from scipy.sparse import csr_matrix
from netCDF4 import Dataset

'''___harmonisation Modules___'''
from Uncertainty import Uncertainty

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "18/11/2017"
__credits__ = ["Jonathan Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


STANDARD_VARIABLE_NAMES = ["time1", "time2", "K", "Kr", "Ks", "across_track_index1", "across_track_index2",
                           "along_track_index1", "along_track_index2"]


class MatchUp(object):
    """
    Class to handle sensor series match-up data from file (format specified within FIDUCEO project)

    :Attributes:
        .. py:attribute:: values

            *numpy.ndarray*

            Match-up data in 1D data structure

        .. py:attribute:: unc

            *list:CorrelForm*

            Uncertainty information for match-up data

        .. py:attribute:: w_matrices

            *list*

            List of dataset w matrices (as ``scipy.sparse.csr_matrix``)

        .. py:attribute:: u_matrices

            *list*

            List of dataset u matrices (as ``numpy.ndarray`` of matrix diagonals)

        .. py:attribute:: ks

            *numpy.ndarray*

            Adjustment factors between sensors in match-up in 1D data structure

        .. py:attribute:: unck

            *list:CorrelForm*

            Uncertainty on adjustment factors


        .. py:attribute:: additional_values

            *numpy.ndarray*

            Per match-up data additional to sensor measurement function data

        .. py:attribute:: a

            *list:CorrelForm*

            Calibration parameters for sensor_function

        .. py:attribute:: idx

            *dict:list*

            Dictionary of indexing lists which define the 1D structure of the match-up data

        .. py:attribute:: _original_idx

            *dict:list*

            Copy of self.idx on opening

    :Methods:
        .. py:method:: openMatchUp(...):

            Reads the data product specified by the given file path.

        .. py:method:: save_to_netcdf(...):

            Writes entire data product to netcdf file at specified path

        .. py:method:: write_matchup_file(...):

            Write match-up file from variable arrays

        .. py:method:: _open_w_matrix(...):

            Open a particular numbered W matrix from match-up file in sparse CSR representation

        .. py:method:: _open_u_matrix(...):

            Open a particular numbered uncertainty from match-up file in sparse CSR representation

        .. py:method:: _seconds2date(self, times):

            Return times formatted as seconds since 1/1/1970 as ``datetime.datetime`` objects

        .. py:method:: _date2seconds(self, times):

            Return times formatted as ``datetime.datetime`` objects as seconds since 1/1/1970

        .. py:method:: _return_w_matrix_variables(...):

            Produce match-up file W matrix variable arrays from lists of W and U matrices data.
    """

    def __init__(self, paths_matchup_data=None, sensor_data=None, open_uncertainty=True, open_additional_values=True):
        """
        Initialise match-up data object, opening match-up data from file if specified

        :type paths_matchup_data: list:str
        :param paths_matchup: paths of match-up file to open. If one file string of path, if multiple files provide list
        of strings of path

        :type sensor_data: str
        :param sensor_data: sensor data

        :type open_uncertainty: bool
        :param open_uncertainty: Switch to determine whether or not to open uncertainty data

        :type open_additional_values: bool
        :param open_additional_values: Switch to determine whether to open data additional to measurement function data
        """

        # 1. Initialise attributes

        # Match-up data
        self.values = array([])
        self.unc = []
        self.w_matrices = []
        self.u_matrices = []
        self.ks = array([])
        self.unck = []
        self.time1 = None
        self.time2 = None
        self.across_track_index1 = None
        self.across_track_index2 = None
        self.along_track_index1 = None
        self.along_track_index2 = None
        self.additional_values = None

        # Sensor data
        self.a = array([], dtype=float32)
        self.sensor_data = None
        self.sensor_model = None
        self.sensor_model_constant = array([])
        self.adjustment_model = None

        # Data indexing
        self.idx = {}
        self._original_idx = {}

        # 2. Populate attributes if appropriate data provided
        if paths_matchup_data is not None:
            self.values, self.unc,\
                self.w_matrices, self.u_matrices,\
                    self.ks, self.unck, \
                        self.time1, self.time2, \
                            self.across_track_index1, self.across_track_index2, \
                                self.along_track_index1, self.along_track_index2, \
                                    self.additional_values, self.idx \
                                        = self.openMatchUpData(paths_matchup_data,
                                                               open_uncertainty=open_uncertainty,
                                                               open_additional_values=open_additional_values)

        if sensor_data is not None:
            self.sensor_data = sensor_data

        if (paths_matchup_data is not None) and (sensor_data is not None):
            self.a, sensor_model_parameter_sensor, \
                self.sensor_model_constant, sensor_model_constant_sensor, \
                    self.sensor_model, self.adjustment_model = self._extract_matchup_sensor_data()

            self.idx["parameter_sensor"] = sensor_model_parameter_sensor
            self.idx["sensor_model_constant_sensor"] = sensor_model_constant_sensor

        # save separate copy of original indices for future reference
        self._original_idx = deepcopy(self.idx)

    def openMatchUpData(self, paths_matchup, open_uncertainty=True, open_additional_values=True):
        """
        Return match-up data from file

        :type paths_matchup: str
        :param paths_matchup: paths of match-up file to open. If one file string of path, if multiple files provide list
        of strings of path

        :type open_uncertainty: bool
        :param open_uncertainty: Switch to determine whether or not to open uncertainty data

        :type open_additional_values: bool
        :param open_additional_values: Switch to determine whether to open data additional to measurement function data

        :return:
            :values: *numpy.ndarray*

            Match-up data in 1D block data structure

            :unc: *list

            Uncertainties associated with values by data block (as ``harmonisation.core.matchup.matchupIO``).
            None if open_uncertainty is False.

            :w_matrices: *list*

            List of dataset w matrices (as ``scipy.sparse.csr_matrix``)

            :u_matrices: *list*

            List of dataset u matrices (as ``numpy.ndarray`` of matrix diagonals)

            :ks: *numpy.ndarray*

            Adjustment factors between sensors in match-ups

            :unck: *list*

            Uncertainty on adjustment factors. None if open_uncertainties is False.

            :time1: *numpy.ndarray*

            Sensor 1 match-up times

            :time2: *numpy.ndarray*

            Sensor 2 match-up times

            :a: *numpy.ndarray*

            Calibration parameters for sensor_function

            :idx: *dict*

            Dictionary of indexing lists which define the 1D structure of the match-up data

            :_original_idx: *dict*

            Reference copy of original state of idx
        """

        # Input match-up data is restructured into a 1D array structured in blocks of data in the following way:
        #
        # values = [ X_1(1,1) |...| X_1(N,S-1) | X_1(L ,S) |   ...   | X_m(1,2) |...| X_m(N,S-1)| X_m(L,S) ]
        #
        # where:
        # - X_m - model covariate of number m of a total m covariates
        # - (N, s) - covariate indices:
        #            > N - match-up series index of N of a total N match-up series
        #            > s - sensor number of S sensors
        #
        # Adjustment Factor, K, data is considered separately and stored in a separate 1D array in blocks by
        # match-up series in the following way:
        #
        # ks = [ K(1) |...| K(N)]
        #
        # Calibration parameters are given a0 to ap (p+1 total parameters) per sensor (for S sensors)
        # and are stored in 1D arrays in the following way:
        #
        # a = [ a_0 (s=0),..., a_p (s=0), ... , a_0 (s=S),...,a_p(s=S) ]
        #
        # Each covariate (+ks) in each match-up has its own uncertainty values and type, create a list to contain this
        # information for the consecuative blocks in the values data array.

        # Open datasets
        datasets = self.open_nc_datasets(paths_matchup)

        # Define 1D Data Structure
        idx = self.get_idx_from_nc_datasets(datasets)

        # Open match-up data
        values = self.open_values_from_nc_datasets(datasets, idx)
        ks = self.open_ks_from_nc_datasets(datasets, idx)

        time1, time2 = self.open_time_from_nc_datasets(datasets, idx)
        across_track_index1 = self.open_optional_index_from_nc_datasets("across_track_index1", datasets, idx)
        across_track_index2 = self.open_optional_index_from_nc_datasets("across_track_index2", datasets, idx)
        along_track_index1 = self.open_optional_index_from_nc_datasets("along_track_index1", datasets, idx)
        along_track_index2 = self.open_optional_index_from_nc_datasets("along_track_index2", datasets, idx)

        unc, unck, w_matrices, u_matrices = None, None, None, None
        if open_uncertainty is True:
            unc, w_matrices, u_matrices = self.open_unc_from_nc_datasets(datasets, idx)
            unck = self.open_unck_from_nc_datasets(datasets, idx)

        additional_values = None
        if open_additional_values:
            additional_values, idx = self.open_additional_values_from_nc_datasets(datasets, idx)

        # Close datasets
        for dataset in datasets:
            dataset.close()

        return values, unc, w_matrices, u_matrices, ks, unck, time1, time2, across_track_index1, across_track_index2,\
               along_track_index1, along_track_index2, additional_values, idx

    def open_nc_datasets(self, paths_matchup):
        """

        :param paths_matchup: str/list
        :param paths_matchup: matchup file paths, may be:

        * string of path of single file
        * string of path of directory containing files
        * list of file paths
        * list of paths of directories containing files

        :return:
        """
        # Determine file paths of match-up files
        if type(paths_matchup) == str:
            if isdir(paths_matchup):
                paths_matchup = sorted([pjoin(paths_matchup, fname)
                                        for fname in listdir(paths_matchup) if fname[-3:] == ".nc"])
            elif isfile(paths_matchup):
                paths_matchup = [paths_matchup]
            else:
                return None

        elif type(paths_matchup) == list:
            if all(isfile(path) for path in paths_matchup):
                pass
            elif all(isdir(path) for path in paths_matchup):
                dirs_matchup = paths_matchup
                paths_matchup = []
                for dir_matchup in dirs_matchup:
                    paths_matchup += [pjoin(dir_matchup, fname) for fname in listdir(dir_matchup)]
            else:
                return None

        return [Dataset(path, 'r') for path in paths_matchup]

    def get_idx_from_nc_datasets(self, datasets):

        Nm = self.get_idx_Nm_from_nc_datasets(datasets)
        cNm = self.get_idx_cNm_from_nc_datasets(datasets)
        sensors, sensor_ms = self.get_idx_sensors_from_nc_datasets(datasets)
        Im = self.get_idx_Im_from_nc_datasets(datasets, sensors)
        n_sensors = self.get_idx_n_sensor_from_nc_datasets(Im, sensor_ms)
        n_mu = self.get_idx_n_mu_from_nc_datasets(Im, sensor_ms)
        n_cov = self.get_idx_n_cov_from_nc_datasets(Im, sensor_ms)
        N_var = self.get_idx_N_var_from_nc_datasets(Nm, n_mu)
        idxs = self.get_idx_idxs_from_nc_datasets(N_var)

        idx = {"Nm": Nm,
               "cNm": cNm,
               "Im": Im,
               "sensors": sensors,
               "sensor_ms": sensor_ms,
               "n_sensor": n_sensors,
               "n_mu": n_mu,
               "n_cov": n_cov,
               "N_var": N_var,
               "idx": idxs}
        return idx

    def get_block_idxs(self, idx):
        return [(i, j, k) for i, j, k in zip(idx['n_sensor'], idx['n_mu'], idx['n_cov'])]

    def get_idx_Nm_from_nc_datasets(self, datasets):
        return [dataset.dimensions['M'].size for dataset in datasets]

    def get_idx_cNm_from_nc_datasets(self, datasets):
        sm = 0
        cNm = [sm]
        for dataset in datasets:
            sm += dataset.dimensions['M'].size
            cNm.append(sm)
        return cNm

    def get_idx_sensors_from_nc_datasets(self, datasets):
        sensors = []
        sensor_ms = []
        for i, dataset in enumerate(datasets):
            for j, (sensor, sensor_m) in enumerate(zip([dataset.sensor_1_name, dataset.sensor_2_name],
                                                       [dataset.dimensions['m1'].size, dataset.dimensions['m2'].size])):
                if sensor in sensors:
                    pass
                else:
                    sensors.append(sensor)
                    sensor_ms.append(sensor_m)
        return sensors, sensor_ms

    def get_idx_Im_from_nc_datasets(self, datasets, sensors):
        Im = [0] * len(datasets)
        for i, dataset in enumerate(datasets):
            Im[i] = [sensors.index(sensor) for sensor in [dataset.sensor_1_name, dataset.sensor_2_name]]
        return Im

    def get_idx_n_sensor_from_nc_datasets(self, Im, sensor_ms):
        sensor_list = [num for pair in Im for num in pair]
        n_sensor = [int(s) for m_i in range(1, max(sensor_ms) + 1)
                    for s in sensor_list if sensor_ms[s] >= m_i]
        return n_sensor

    def get_idx_n_mu_from_nc_datasets(self, Im, sensor_ms):
        n_mu = [int(i + 1) for m_i in range(1, max(sensor_ms) + 1)
                for i, pair in enumerate(Im)
                for s in pair if sensor_ms[s] >= m_i]
        return n_mu

    def get_idx_n_cov_from_nc_datasets(self, Im, sensor_ms):
        sensor_list = [num for pair in Im for num in pair]
        n_cov = [m_i for m_i in range(1, max(sensor_ms) + 1)
                 for s in sensor_list if sensor_ms[s] >= m_i]
        return n_cov

    def get_idx_N_var_from_nc_datasets(self, Nm, n_mu):
        return [Nm[n - 1] for n in n_mu]

    def get_idx_idxs_from_nc_datasets(self, N_var):
        idxs = [0]
        total = 0
        for N in N_var:
            total += N
            idxs.append(int(total))
        return idxs

    def get_idx_additional_values_name_from_nc_datasets(self, datasets):
        additional_variables = set()
        for dataset in datasets:
            matchup_additional_variables = set([str(variable) for variable in dataset.variables.keys()
                                                if (dataset.variables[variable].dimensions == ("M",)) and
                                                   (variable not in STANDARD_VARIABLE_NAMES)])
            additional_variables = additional_variables.union(matchup_additional_variables)

        return list(additional_variables)

    def open_values_from_nc_datasets(self, datasets, idx):

        values = zeros(idx['idx'][-1], dtype=float32)  # 1D organised match-up measurement data

        # sensor number, match up number and covariate number per data block as a tuple for indexing
        block_idxs = self.get_block_idxs(idx)

        # Open data by match-up series
        for i_matchup, dataset in enumerate(datasets):

            # a. Open data by sensor -----------------------------------------------------------------------------------
            for i_sensor_pair, n_sensor in enumerate(idx["Im"][i_matchup]):

                # Open data and uncertainties by covariate
                for i_cov in range(idx["sensor_ms"][n_sensor]):
                    # Find location of block for this covariate in data structure
                    block_idx = block_idxs.index((n_sensor, i_matchup + 1, i_cov + 1))

                    # i. X - add data to values array ------------------------------------------------------------------
                    i_start = idx['idx'][block_idx]
                    i_end = idx['idx'][block_idx + 1]
                    values[i_start:i_end] = dataset.variables['X' + str(i_sensor_pair + 1)][:, i_cov]

        return values

    def open_unc_from_nc_datasets(self, datasets, idx):

        unc = [0] * len(idx['N_var'])            # 1D organised match-up measurement data uncertainty

        w_matrix_idxs, u_matrix_idxs = self.open_matrix_idxs_from_nc_datasets(datasets)

        w_matrices = [0] * len(w_matrix_idxs)
        u_matrices = [0] * len(u_matrix_idxs)

        # sensor number, match up number and covariate number per data block as a tuple for indexing
        block_idxs = self.get_block_idxs(idx)

        # Open data by match-up series
        for i_matchup, dataset in enumerate(datasets):

            # a. Open data by sensor -----------------------------------------------------------------------------------
            for i_sensor_pair, n_sensor in enumerate(idx["Im"][i_matchup]):

                # Open sensor data labels
                uncertainty_type = dataset.variables['uncertainty_type' + str(i_sensor_pair + 1)][:]

                w_matrix_use = zeros(uncertainty_type.shape[0])
                u_matrix_use = zeros(uncertainty_type.shape[0])
                if {"w_matrix_use1", "w_matrix_use2", "u_matrix_use1", "u_matrix_use2"}.issubset(
                        dataset.variables.keys()):
                    w_matrix_use = dataset.variables['w_matrix_use' + str(i_sensor_pair + 1)][:]
                    u_matrix_use = dataset.variables['u_matrix_use' + str(i_sensor_pair + 1)][:]

                # Open data and uncertainties by covariate
                for i_cov, (uncertainty_type, w_j, u_j) in enumerate(
                        zip(uncertainty_type, w_matrix_use, u_matrix_use)):

                    # Find location of block for this covariate in data structure
                    block_idx = block_idxs.index((n_sensor, i_matchup + 1, i_cov + 1))

                    # > independent errors
                    if uncertainty_type == 1:
                        unc[block_idx] = Uncertainty(1, dataset.variables['Ur' + str(i_sensor_pair + 1)][:, i_cov])

                    # > independent errors + systematic error
                    if uncertainty_type == 2:
                        unc[block_idx] = Uncertainty(2, (
                        dataset.variables['Ur' + str(i_sensor_pair + 1)][:, i_cov],
                        dataset.variables['Us' + str(i_sensor_pair + 1)][:, i_cov]))

                    # > w matrix defined error correlation
                    if uncertainty_type == 3:
                        # Find location to store w and u matrix in w_matrices and u_matrices storage lists
                        w_matrix_idx = w_matrix_idxs.index([i_matchup + 1, w_j])
                        u_matrix_idx = u_matrix_idxs.index([i_matchup + 1, u_j])

                        # Open w and u matrix
                        w_matrices[w_matrix_idx] = self._open_w_matrix(dataset, w_j)
                        u_matrices[u_matrix_idx] = self._open_u_matrix(dataset, u_j)

                        unc[block_idx] = Uncertainty(3, (w_matrix_idx, u_matrix_idx))

                    # > w matrix + systematic error correlation
                    if uncertainty_type == 4:
                        # Find location to store w and u matrix in w_matrices and u_matrices storage lists
                        w_matrix_idx = w_matrix_idxs.index([i_matchup + 1, w_j])
                        u_matrix_idx = u_matrix_idxs.index([i_matchup + 1, u_j])

                        # Open w and u matrix
                        w_matrices[w_matrix_idx] = self._open_w_matrix(dataset, w_j)
                        u_matrices[u_matrix_idx] = self._open_u_matrix(dataset, u_j)

                        unc[block_idx] = Uncertainty(4, (w_matrix_idx, u_matrix_idx,
                                                         dataset.variables['Us' + str(i_sensor_pair + 1)][1, i_cov]))

        # Determine separate instances of systematic errors
        sys_indices = []
        for i, block_index in enumerate(zip(idx['n_sensor'], idx['n_cov'])):
            if (unc[i].typeID == 2) or (unc[i].typeID == 4):
                if block_index not in sys_indices:
                    sys_indices.append(block_index)
                    unc[i].uS_i = len(sys_indices)
                else:
                    unc[i].uS_i = sys_indices.index(block_index) + 1

        return unc, w_matrices, u_matrices

    def _open_w_matrix(self, dataset, w_i):
        """
        Return a particular numbered W matrix from match-up data file in sparse CSR representation

        :type dataset: netCDF4.Dataset
        :param dataset: in-memory representation of harmonisation input file

        :type w_i: int
        :param w_i: number of W matrix to open from file

        :return:
            :W: *scipy.sparse.csr_matrix*

            W matrix in sparse csr representation
        """

        # W matrix data indices
        i_nnz = append(zeros(1), cumsum(dataset.variables["w_matrix_nnz"][:])).astype(dtype=int32)

        # Read appropriate W matrix data
        data = dataset.variables["w_matrix_val"][i_nnz[w_i - 1]:i_nnz[w_i]]
        indices = dataset.variables["w_matrix_col"][i_nnz[w_i - 1]:i_nnz[w_i]]
        indptr = dataset.variables["w_matrix_row"][w_i - 1, :]

        # Determine W matrix shape
        shape = (len(indptr) - 1, max(indices) + 1)

        return csr_matrix((data, indices, indptr), shape=shape)

    def _open_u_matrix(self, dataset, u_i):
        """
        Return diagonal of a particular numbered u matrix from match-up data file

        :type dataset: netCDF4.Dataset
        :param dataset: in-memory representation of harmonisation input file

        :type u_i: int
        :param u_i: number of uncertainty to open from file

        :return:
            :u_matrix: *numpy.ndarray*

            diagonal of u matrix
        """

        # W matrix data indices
        i_vec = append(zeros(1), cumsum(dataset.variables["u_matrix_row_count"][:])).astype(dtype=int32)

        return dataset.variables["u_matrix_val"][i_vec[u_i - 1]:i_vec[u_i]]

    def open_matrix_idxs_from_nc_datasets(self, datasets):
        w_matrix_idxs = []
        u_matrix_idxs = []

        for i_matchup, dataset in enumerate(datasets):

            if {"w_matrix_use1", "w_matrix_use2", "u_matrix_use1", "u_matrix_use2"}.issubset(dataset.variables.keys()):
                w_matrix_use = append(dataset["w_matrix_use1"][:], dataset["w_matrix_use2"][:])
                for w_matrix_i in set(w_matrix_use[w_matrix_use != 0]):
                    w_matrix_idxs.append([i_matchup + 1, int(w_matrix_i)])

                u_matrix_use = append(dataset["u_matrix_use1"][:], dataset["u_matrix_use2"][:])
                for u_matrix_i in set(u_matrix_use[u_matrix_use != 0]):
                    u_matrix_idxs.append([i_matchup + 1, int(u_matrix_i)])

        return w_matrix_idxs, u_matrix_idxs

    def open_ks_from_nc_datasets(self, datasets, idx):

        ks = zeros(idx['cNm'][-1])  # 1D organised match-up k data

        for i_matchup, dataset in enumerate(datasets):

            i_start = idx['cNm'][i_matchup]
            i_end = idx['cNm'][i_matchup + 1]
            ks[i_start:i_end] = dataset.variables['K'][:]

        return ks

    def open_unck_from_nc_datasets(self, datasets, idx):
        return [Uncertainty(1, (d.variables['Kr'][:] ** 2 + d.variables['Ks'][:] ** 2) ** 0.5) for d in datasets]

    def open_time_from_nc_datasets(self, datasets, idx):

        time1 = zeros(idx['cNm'][-1]).astype(datetime)
        time2 = zeros(idx['cNm'][-1]).astype(datetime)

        for i_matchup, dataset in enumerate(datasets):
            i_start = idx['cNm'][i_matchup]
            i_end = idx['cNm'][i_matchup + 1]

            time1[i_start:i_end] = self._seconds2date(dataset.variables['time1'][:])
            time2[i_start:i_end] = self._seconds2date(dataset.variables['time2'][:])

        return time1, time2

    def _seconds2date(self, times):
        """
        Return match-up times in ``datetime.datetime`` format

        :type times: *numpy.ndarray*
        :param times: Times in seconds since 1/1/1970

        :return:
            :dates: *numpy.ndarray*

            Times in ``datetime.datetime`` format
        """

        dates = zeros(times.shape[0], dtype=datetime)
        for i, time in enumerate(times):
            dates[i] = datetime.fromtimestamp(int(time))

        return dates

    def open_optional_index_from_nc_datasets(self, index_name, datasets, idx):

        index = None

        index_req_mu = set([True if index_name in dataset.variables.keys() else False
                                          for dataset in datasets])
        if True in index_req_mu:
            index = zeros(idx['cNm'][-1])

            for i_matchup, dataset in enumerate(datasets):

                i_start = idx['cNm'][i_matchup]
                i_end = idx['cNm'][i_matchup + 1]

                if index_name in dataset.variables.keys():
                    index[i_start:i_end] = dataset.variables[index_name][:].astype(int32)

                else:
                    index[i_start:i_end] = nan

        return index

    def open_additional_values_from_nc_datasets(self, datasets, idx):

        idx["additional_values_name"] = self.get_idx_additional_values_name_from_nc_datasets(datasets)

        additional_values = zeros((idx['cNm'][-1], len(idx["additional_values_name"])))
        for i_matchup, dataset in enumerate(datasets):
            i_start = idx['cNm'][i_matchup]
            i_end = idx['cNm'][i_matchup + 1]

            for i_variable, variable in enumerate(idx["additional_values_name"]):

                if variable in dataset.variables.keys():
                    additional_values[i_start:i_end, i_variable] = dataset.variables[variable][:]

                    # Set fill data (value -1e30) to nan
                    if variable == "nominal_BT2":
                        pass
                    additional_values[i_start:i_end, i_variable]\
                        [additional_values[i_start:i_end, i_variable] == -1e30] = nan
                else:
                    additional_values[i_start:i_end, i_variable] = nan

        return additional_values, idx

    def _open_sensor_data(self, path_sensor_data):
        """
        Return match-up sensor data from file

        :type path_sensor_data: str
        :param path_sensor_data: path of sensor data file

        :return:
            :sensor_data: *dict*

            Sensor data dictionary
        """

        directory_sensor_data = split(path_sensor_data)[0]
        fname_sensor_data = splitext(basename(path_sensor_data))[0]

        sys.path.insert(0, directory_sensor_data)

        sensor_data = import_module(fname_sensor_data).sensor_data

        return sensor_data

    def _extract_matchup_sensor_data(self):
        """
        Return matchup relative sensor data from opened sensor file

        :return:
        :sensor_data: *dict *

        Sensor
        data
        dictionary
        """

        sensor_model_parameter = []
        sensor_model_parameter_sensor = []
        sensor_model_constant = []
        sensor_model_constant_sensor = []
        sensor_model = []
        adjustment_model = []

        for sensor in self.sensor_data.keys():
            if sensor in self.idx['sensors']:
                sensor_model_parameter += self.sensor_data[sensor]["sensor_model_parameter"]
                sensor_model_parameter_sensor += [sensor for a in self.sensor_data[sensor]["sensor_model_parameter"]]

                sensor_model_constant += self.sensor_data[sensor]["sensor_model_constant"]
                sensor_model_constant_sensor += [sensor for a in self.sensor_data[sensor]["sensor_model_constant"]]

        for sensor in self.idx['sensors']:
            sensor_model.append(self.sensor_data[sensor]["sensor_model"])

            if "adjustment_model" in self.sensor_data[sensor].keys():
                adjustment_model.append(self.sensor_data[sensor]["sensor_adjustment_model"])
            else:
                adjustment_model.append(self._default_adjustment_model)

        return asarray(sensor_model_parameter, dtype=float32), sensor_model_parameter_sensor,\
               sensor_model_constant, sensor_model_constant_sensor, sensor_model, adjustment_model

    def _default_adjustment_model(self, measurand):
        """
        Return default adjusted measurand array and derivatives from measurand values for default case,
        adjusted_measurand = measurand

        :type measurand: numpy.ndarray
        :param measurand: Measurand values

        :return:
            :adjusted_measurand: *numpy.ndarray*

            Adjusted measurand

            :adjusted_measurand_derivatives: *numpy.ndarray*

            Adjusted measurand derivatives
        """

        return measurand, ones(len(measurand))

    def save_to_netcdf(self, path):
        """
        Writes match-up data to netcdf file at specified path

        :type path: str
        :param path: Path of file to store data to
        """

        # todo - write save_to_netcdf()

        return 0

    def _return_w_matrix_variables(self, w_matrices, u_matrices):
        # todo - write this part!

        return 0

    def _date2seconds(self, times):
        """
        Return times in ``datetime.datetime`` format as seconds since 1/1/1970

        :type times: *numpy.ndarray*
        :param times: Times in ``datetime.datetime``

        :return:
            :dates: *numpy.ndarray*

            Times as seconds since 1/1/1970
        """

        # todo - write _date2seconds()

        return 0

    def getVariableData(self, variable_number, sensor_name, matchup_number):

        matchup_numbers = []
        if matchup_number is None:
            matchup_numbers = list(range(1, len(self.idx['Nm']) + 1))
        elif type(matchup_number) == float:
            matchup_numbers = [matchup_number]
        elif type(matchup_number) == list:
            matchup_numbers = matchup_number

        block_idxs = [(i, j, k) for i, j, k in zip([self.idx['sensors'][s] for s in self.idx['n_sensor']],
                                                   self.idx['n_mu'], self.idx['n_cov'])]

        idx = block_idxs.index((sensor_name, matchup_number, variable_number))
        istart = self.idx['idx'][idx]
        iend = self.idx['idx'][idx+1]
        return self.values[istart:iend]

    def getVariableDataMax(self, variable_number, sensor_name=None, matchup_number=None):

        try:
            return max(self.getVariableData(variable_number, sensor_name, matchup_number))
        except ValueError:
            return None

    def getVariableDataMin(self, variable_number, sensor_name,  matchup_number=None):

        try:
            return min(self.getVariableData(variable_number, sensor_name, matchup_number))
        except ValueError:
            return None

    def getAdditionalData(self, additional_variable_name, matchup_number=None):
        istart = 0
        iend = None
        if matchup_number is not None:
            istart = self.idx['cNm'][matchup_number-1]
            iend = self.idx['cNm'][matchup_number]
        additional_variable_col = self.idx['additional_values_name'].index(additional_variable_name)
        return self.additional_values[istart:iend, additional_variable_col]

    def getSensorAlongTrackIndex(self, n_mu, sensor_name):
        """
        :type n_mu: int
        :param n_mu: match-up number

        :type n_mu: str
        :param n_mu: name of sensor

        :return:
            :idx:
        """

        istart = self.idx['cNm'][n_mu-1]
        iend = self.idx['cNm'][n_mu]

        if sensor_name == self.idx['sensors'][self.idx['Im'][n_mu-1][0]]:
            if self.along_track_index1 is not None:
                return self.along_track_index1[istart:iend]
            else:
                return None
        elif sensor_name == self.idx['sensors'][self.idx['Im'][n_mu-1][1]]:
            if self.along_track_index2 is not None:
                return self.along_track_index2[istart:iend]
            else:
                return None
        else:
            raise IndexError

    def getSensorAcrossTrackIndex(self, n_mu, sensor_name):
        """
        :type n_mu: int
        :param n_mu: match-up number

        :type n_mu: str
        :param n_mu: name of sensor

        :return:
            :idx:
        """

        istart = self.idx['cNm'][n_mu-1]
        iend = self.idx['cNm'][n_mu]

        if sensor_name == self.idx['sensors'][self.idx['Im'][n_mu-1][0]]:
            if self.across_track_index1 is not None:
                return self.across_track_index1[istart:iend]
            else:
                return None
        elif sensor_name == self.idx['sensors'][self.idx['Im'][n_mu-1][1]]:
            if self.across_track_index2 is not None:
                return self.across_track_index2[istart:iend]
            else:
                return None
        else:
            raise IndexError

    def getSensorTime(self, n_mu, sensor_name):
        """
        :type n_mu: int
        :param n_mu: match-up number

        :type n_mu: str
        :param n_mu: name of sensor

        :return:
            :idx:
        """

        istart = self.idx['cNm'][n_mu-1]
        iend = self.idx['cNm'][n_mu]

        if sensor_name == self.idx['sensors'][self.idx['Im'][n_mu-1][0]]:
            if self.time1 is not None:
                return self.time1[istart:iend]
            else:
                return None
        elif sensor_name == self.idx['sensors'][self.idx['Im'][n_mu-1][1]]:
            if self.time2 is not None:
                return self.time2[istart:iend]
            else:
                return None
        else:
            raise IndexError

    def getAdditionalDataNames(self):
        return self.idx['additional_values_name']

    def getAdditionalDataMax(self, additional_variable_name, matchup_number=None):
        return max(self.getAdditionalData(additional_variable_name, matchup_number))

    def getAdditionalDataMin(self, additional_variable_name, matchup_number=None):
        return min(self.getAdditionalData(additional_variable_name, matchup_number))

    def write_input_file(self, file_path, X1, X2, Ur1, Ur2, Us1, Us2, uncertainty_type1, uncertainty_type2,
                               K, Kr, Ks, time1, time2, sensor_1_name, sensor_2_name,
                               w_matrix_variables=None, additional_variables=None):
        """
        Write match-up file from variable arrays

        :type file_path: str
        :param file_path: match-up file path

        :type X1: numpy.ndarray
        :param X1: Radiances and counts per matchup for sensor 1

        :type X2: numpy.ndarray
        :param X2: Radiances and counts per matchup for sensor 2

        :type Ur1: numpy.ndarray
        :param Ur1: Random uncertainties for X1 array

        :type Ur2: numpy.ndarray
        :param Ur2: Random uncertainties for X2 array

        :type Us1: numpy.ndarray
        :param Us1: Systematic uncertainties for X1 array

        :type Us2: numpy.ndarray
        :param Us2: Systematic uncertainties for X2 array

        :type uncertainty_type1: numpy.ndarray
        :param uncertainty_type1: Uncertainty correlation type per X1 column

        :type uncertainty_type2: numpy.ndarray
        :param uncertainty_type2: Uncertainty correlation type per X2 column

        :type K: numpy.ndarray
        :param K: K (sensor-to-sensor differences) for zero shift case

        :type Kr: numpy.ndarray
        :param Kr: K (sensor-to-sensor differences) random uncertainties (matchup uncertainty)

        :type Ks: numpy.ndarray
        :param Ks: K (sensor-to-sensor differences) systematic uncertainties for zero shift case

        :type time1: numpy.ndarray
        :param time1: "Match-up time sensor 1, seconds since 1970-01-01"

        :type time2: numpy.ndarray
        :param time2: "Match-up time sensor 2, seconds since 1970-01-01"

        :type sensor_1_name: int
        :param sensor_1_name: sensor i ID

        :type sensor_2_name: int
        :param sensor_2_name: sensor j ID

        :type w_matrix_variables: tuple
        :param w_matrix_variables: 10 variables required to build the file w matrix variables, arranged as:

        * w_matrix_nnz = w_matrix_variables[0]
        * w_matrix_row = w_matrix_variables[1]
        * w_matrix_col = w_matrix_variables[2]
        * w_matrix_val = w_matrix_variables[3]
        * w_matrix_use1 = w_matrix_variables[4]
        * w_matrix_use2 = w_matrix_variables[5]
        * uncertainty_vector_row_count = w_matrix_variables[6]
        * uncertainty_vector = w_matrix_variables[7]
        * uncertainty_vector_use1 = w_matrix_variables[8]
        * uncertainty_vector_use2 = w_matrix_variables[9]

        :type additional_variables: dict
        :param additional_variables: dictionary of additional, non-required variable to add to match up input files
        To be for e.g. testing, diagnostics etc.

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
        rootgrp.sensor_1_name = sensor_1_name
        rootgrp.sensor_2_name = sensor_2_name

        # 2. Create dimensions
        # > M - number of matchups
        M_dim = rootgrp.createDimension("M", X1.shape[0])

        # > m - number of columns in X1 and X2 arrays
        m1_dim = rootgrp.createDimension("m1", X1.shape[1])
        m2_dim = rootgrp.createDimension("m2", X2.shape[1])

        # 3. Create new variables

        # > X1 - Radiances and counts per matchup for sensor 1
        X1_var = rootgrp.createVariable('X1', 'f8', ('M', 'm1'), zlib=True, complevel=9)
        X1_var.description = "Radiances and counts per matchup for sensor 1"
        X1_var[:] = X1[:]

        # > X2 - Radiances and counts per matchup for sensor 2
        X2_var = rootgrp.createVariable('X2', 'f8', ('M', 'm2'), zlib=True, complevel=9)
        X2_var.description = "Radiances and counts per matchup for sensor 2"
        X2_var[:] = X2[:]

        # > Ur1 - Random uncertainties for X1 array
        Ur1_var = rootgrp.createVariable('Ur1', 'f8', ('M', 'm1'), zlib=True, complevel=9)
        Ur1_var.description = "Random uncertainties for X1 array"
        Ur1_var[:] = Ur1[:]

        # > Ur2 - Random uncertainties for X2 array
        Ur2_var = rootgrp.createVariable('Ur2', 'f8', ('M', 'm2'), zlib=True, complevel=9)
        Ur2_var.description = "Random uncertainties for X2 array"
        Ur2_var[:] = Ur2[:]

        # > Us1 - Systematic uncertainties for X1 array
        Us1_var = rootgrp.createVariable('Us1', 'f8', ('M', 'm1'), zlib=True, complevel=9)
        Us1_var.description = "Systematic uncertainties for X1 array"
        Us1_var[:] = Us1[:]

        # > Us2 - Systematic uncertainties for X2 array
        Us2_var = rootgrp.createVariable('Us2', 'f8', ('M', 'm2'), zlib=True, complevel=9)
        Us2_var.description = "Systematic uncertainties for X2 array"
        Us2_var[:] = Us2[:]

        # > uncertainty_type1 - Uncertainty correlation type per X1 column
        uncertainty_type1_var = rootgrp.createVariable('uncertainty_type1', 'i4', ('m1',), zlib=True, complevel=9)
        uncertainty_type1_var.description = "Uncertainty correlation type per X1 column, labelled as, " + \
                                            "(1) Independent Error Correlation, " + \
                                            "(2) Independent + Systematic Error Correlation, or " + \
                                            "(3) Structured Error Correlation"
        uncertainty_type1_var[:] = uncertainty_type1[:]

        # > uncertainty_type2 - Uncertainty correlation type per X2 column
        uncertainty_type2_var = rootgrp.createVariable('uncertainty_type2', 'i4', ('m2',), zlib=True, complevel=9)
        uncertainty_type2_var.description = "Uncertainty correlation type per X2 column, labelled as, " + \
                                            "(1) Independent Error Correlation, " + \
                                            "(2) Independent + Systematic Error Correlation, or " + \
                                            "(3) Structured Error Correlation"
        uncertainty_type2_var[:] = uncertainty_type2[:]

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

        # > time1 - Sensor 1 time of match-up
        time1_var = rootgrp.createVariable('time1', 'f8', ('M',), zlib=True, complevel=9)
        time1_var.description = "Match-up time sensor 1, seconds since 1970-01-01"
        time1_var[:] = time1[:]

        # > time 2 - Sensor 2 time of match-up
        time2_var = rootgrp.createVariable('time2', 'f8', ('M',), zlib=True, complevel=9)
        time2_var.description = "Match-up time sensor 2, seconds since 1970-01-01"
        time2_var[:] = time2[:]

        # > w_matrix variables
        if w_matrix_variables is not None:
            pass
            # todo - update this part!!!

        # > additional variables - non-required variable to add to harmonisation input files
        if additional_variables is not None:
            for variable in additional_variables.keys():
                additional_var = rootgrp.createVariable(variable, additional_variables[variable]['dtype'],
                                                        additional_variables[variable]['dim'], zlib=True, complevel=9)
                additional_var.Description = additional_variables[variable]['Description']
                additional_var[:] = additional_variables[variable]['data'][:]

        # 5. Close file
        rootgrp.close()

        return 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __getitem__(self, i, j, k):
        return

    def __repr__(self):

        string = "Match-Up Dataset\n"
        string += "".join(["="]*len(string))+"\n"

        string += "Sensors - "+str([str(sensor) for sensor in self.idx['sensors']])
        string += "\nTotal Match-Ups - " + str(self.idx['cNm'][-1])
        string += "\nMatch-Up Series - " + str(len(self.idx['cNm']))
        string +="\nTotal Match-Ups - "+str(self.idx['cNm'][-1])
        string +="\nTotal Sensor State Data Values - "+str(self.idx['idx'][-1])+"\n"

        for n_mu, Im in enumerate(self.idx['Im']):

            Nm_mu = self.idx['Nm'][n_mu]
            title = "\nMatch-up Series " + str(n_mu+1)
            string += title+"\n"
            underline = ["-"]*len(title)
            string += "".join(underline)+"\n"

            string += "Series Sensors: "+str(self.idx['sensors'][Im[0]])+", "+str(self.idx['sensors'][Im[1]])+"\n"
            string += "Series Match-ups: " + str(Nm_mu)+"\n"

            if Nm_mu < 10:
                shortened = False
                n_mu_sel = list(range(self.idx['Nm'][n_mu]))
            else:
                shortened = True
                n_mu_sel = [0,1,2,3,-1]

            for i_pair, n_sensor in enumerate(Im):
                string += "\n- " + str(self.idx['sensors'][n_sensor]) + " Sensor State Data \n"

                for i_n_mu_i, n_mu_i in enumerate(n_mu_sel):

                    if shortened:
                        if i_n_mu_i == len(n_mu_sel)-2:
                            string += "...\n"

                        else:
                            string += "["
                            for n_cov in range(self.idx['sensor_ms'][n_sensor]):
                                string += "{:.2f}".format(self.getVariableData(n_cov+1,
                                                                        self.idx['sensors'][n_sensor], n_mu+1)[n_mu_i])+", "
                            string += "]\n"
                    else:
                        string += "["
                        for n_cov in range(self.idx['sensor_ms'][n_sensor]):
                            string += "{:.2f}".format(self.getVariableData(n_cov + 1,
                                                                           self.idx['sensors'][n_sensor], n_mu + 1)[
                                                          n_mu_i]) + ", "
                        string += "]\n"

                if self.unc is not None:
                    string += "\n"+str(self.idx['sensors'][n_sensor])+" Sensor State Data Uncertainty Types: "
                    string += "["
                    indices = [(i1, i2, i3) for i1, i2, i3 in zip(self.idx['n_sensor'], self.idx['n_mu'], self.idx['n_cov'])]

                    for n_cov in range(self.idx['sensor_ms'][n_sensor]):
                        im = indices.index((n_sensor, n_mu + 1, n_cov + 1))
                        string += str(self.unc[im].typeID)+", "
                    string += "]\n"

            string += "\n- Series K\n"
            i_start = self.idx['cNm'][n_mu]
            i_end = self.idx['cNm'][n_mu+1]
            for i_n_mu_i, n_mu_i in enumerate(n_mu_sel):

                if shortened:
                    if i_n_mu_i == len(n_mu_sel) - 2:
                        string += "...\n"

                    else:
                        string += "["
                        string += "{:.2f}".format(self.ks[i_start+n_mu_i]) + ", "
                        string += "]\n"
                else:
                    string += "["
                    string += "{:.2f}".format(self.ks[i_start + n_mu_i]) + ", "
                    string += "]\n"

            if self.unck is not None:
                string+="\nSeries K Uncertainty Type: "+str(self.unck[n_mu].typeID)

            string+="\n"

        return string

    def close(self):
        attrs = vars(self).keys()
        for attr in attrs:
            self.__delattr__(attr)

        self.closed = True


if __name__ == "__main__":
    from numpy import where, isnan
    MatchUp("/home/data/satellite/AVHRR/L0/ch4/matchups/AVH11_REAL_4_RSA/m02_n18.nc")
    pass
