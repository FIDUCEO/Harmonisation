"""
MatchUp Data Class
"""

'''___Built-In Modules___'''
from copy import deepcopy

'''___Third-Party Modules___'''
from numpy import array, float32
from netCDF4 import Dataset

'''___harmonisation Modules___'''
from harmonisation.core.matchup.matchupIO.MatchUpReader import MatchUpReader


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

        if sensor_data is not None:
            self.sensor_data = sensor_data

        if (paths_matchup_data is not None) and (sensor_data is not None):
            self.a, sensor_model_parameter_sensor, \
                self.sensor_model_constant, sensor_model_constant_sensor, \
                    self.sensor_model, self.adjustment_model = self.extract_sensor_data()

            self.idx["parameter_sensor"] = sensor_model_parameter_sensor
            self.idx["sensor_model_constant_sensor"] = sensor_model_constant_sensor

        # save separate copy of original indices for future reference
        self._original_idx = deepcopy(self.idx)

    def setSensorData(self, sensor_data):
        MatchUpReaderOp = MatchUpReader()
        self.a, self.sensor_model_constant, self.sensor_model, \
            self.adjustment_model, self.idx = MatchUpReaderOp.extract_sensor_data(sensor_data)

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
