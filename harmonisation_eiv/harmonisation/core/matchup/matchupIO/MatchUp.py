"""
MatchUp Data Class
"""

'''___Built-In Modules___'''

'''___Third-Party Modules___'''
from numpy import array, float32

'''___harmonisation Modules___'''
import MatchUpReader
import MatchUpWriter


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

    def __init__(self):
        """
        Initialise match-up data object, opening match-up data from file if specified
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

    def setSensorData(self, sensor_data):
        self.sensor_data = sensor_data

        MatchUpReaderOp = MatchUpReader.MatchUpReader()
        self.a, self.sensor_model_constant, self.sensor_model, \
            self.adjustment_model, self.idx = MatchUpReaderOp.extract_sensor_data(sensor_data, self.idx)

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

    def save_to_netcdf(self, path):
        """
        Writes match-up data to netcdf file at specified path

        :type path: str
        :param path: Path of file to store data to
        """

        matchupwriter_op = MatchUpWriter.MatchUpWriter()
        matchupwriter_op.writeMatchUp(self, path)
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
    pass
