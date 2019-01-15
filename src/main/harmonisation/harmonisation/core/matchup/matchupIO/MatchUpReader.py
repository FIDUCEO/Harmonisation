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
from harmonisation import MatchUp
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


def open_matchup(paths_matchup, path_sensor_data=None, open_uncertainty=True, open_additional_values=True):
    MatchUpReaderOp = MatchUpReader()
    return MatchUpReaderOp.openMatchUp(paths_matchup, path_sensor_data,
                                       open_uncertainty=open_uncertainty, open_additional_values=open_additional_values)


class MatchUpReader:
    """
    Class to handle match-up data from file (format specified within FIDUCEO project)

    :Methods:
        .. py:method:: openMatchUp(...):

            Reads the data product specified by the given file path.

        .. py:method:: _open_w_matrix(...):

            Open a particular numbered W matrix from match-up file in sparse CSR representation

        .. py:method:: _open_u_matrix(...):

            Open a particular numbered uncertainty from match-up file in sparse CSR representation

        .. py:method:: _seconds2date(self, times):

            Return times formatted as seconds since 1/1/1970 as ``datetime.datetime`` objects

    """

    def openMatchUp(self, paths_matchup, path_sensor_data=None, open_uncertainty=True, open_additional_values=True):
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
            :matchup: *harmonisation.MatchUp*

            match up data object
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

        # Open match-up datasets
        datasets = self.open_nc_datasets(paths_matchup)

        # Initialise empty matchup object
        matchup = MatchUp()

        # Define 1D Data Structure
        matchup.idx = self.get_idx_from_nc_datasets(datasets)

        # Open match-up data
        matchup.values = self.open_values_from_nc_datasets(datasets, matchup.idx)
        matchup.ks = self.open_ks_from_nc_datasets(datasets, matchup.idx)

        matchup.time1, matchup.time2 = self.open_time_from_nc_datasets(datasets, matchup.idx)
        matchup.across_track_index1 = self.open_optional_index_from_nc_datasets("across_track_index1",
                                                                                datasets, matchup.idx)
        matchup.across_track_index2 = self.open_optional_index_from_nc_datasets("across_track_index2",
                                                                                datasets, matchup.idx)
        matchup.along_track_index1 = self.open_optional_index_from_nc_datasets("along_track_index1",
                                                                               datasets, matchup.idx)
        matchup.along_track_index2 = self.open_optional_index_from_nc_datasets("along_track_index2",
                                                                               datasets, matchup.idx)

        matchup.unc, matchup.unck, matchup.w_matrices, matchup.u_matrices = None, None, None, None
        if open_uncertainty is True:
            matchup.unc, matchup.w_matrices, matchup.u_matrices = self.open_unc_from_nc_datasets(datasets, matchup.idx)
            matchup.unck = self.open_unck_from_nc_datasets(datasets, matchup.idx)

        matchup.additional_values = None
        if open_additional_values:
            matchup.additional_values, matchup.idx = self.open_additional_values_from_nc_datasets(datasets, matchup.idx)

        # Close datasets
        for dataset in datasets:
            dataset.close()

        # Open sensor data
        if path_sensor_data is not None:
            sensor_data = self.open_sensor_data(path_sensor_data)
            matchup.a, matchup.sensor_model_constant, matchup.sensor_model,\
                matchup.adjustment_model, matchup.idx = self.extract_sensor_data(sensor_data, matchup.idx)

        matchup._original_idx = deepcopy(matchup.idx)

        return matchup

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

    def open_sensor_data(self, path_sensor_data):
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

        sensor_data = import_module(fname_sensor_data).SENSOR_DATA

        return sensor_data

    def extract_sensor_data(self, sensor_data, idx):
        """
        Return matchup relative sensor data from opened sensor file

        :return:
            :sensor_data: *dict*

            Sensor data dictionary
        """

        sensor_model_parameter = []
        sensor_model_parameter_sensor = []
        sensor_model_constant = []
        sensor_model_constant_sensor = []
        sensor_model = []
        adjustment_model = []

        for sensor in sensor_data.keys():
            if sensor in idx['sensors']:
                sensor_model_parameter += sensor_data[sensor]["sensor_model_parameter"]
                sensor_model_parameter_sensor += [sensor for a in sensor_data[sensor]["sensor_model_parameter"]]

                sensor_model_constant += sensor_data[sensor]["sensor_model_constant"]
                sensor_model_constant_sensor += [sensor for a in sensor_data[sensor]["sensor_model_constant"]]

        for sensor in idx['sensors']:
            sensor_model.append(sensor_data[sensor]["sensor_model"])

            if "adjustment_model" in sensor_data[sensor].keys():
                adjustment_model.append(sensor_data[sensor]["sensor_adjustment_model"])
            else:
                adjustment_model.append(self._default_adjustment_model)

        idx["parameter_sensor"] = sensor_model_parameter_sensor
        idx["sensor_model_constant_sensor"] = sensor_model_constant_sensor

        return asarray(sensor_model_parameter, dtype=float32), sensor_model_constant, sensor_model, \
               adjustment_model, idx

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


if __name__ == "__main__":
    pass
