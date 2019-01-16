"""
MatchUpWriter Class
"""

'''___Built-In Modules___'''

'''___Third-Party Modules___'''
from netCDF4 import Dataset

'''___harmonisation Modules___'''


'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "15/1/2019"
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class MatchUpWriter:
    """
    Class to handle sensor series match-up files from matchup data object (format specified within FIDUCEO project)

    :Methods:
        .. py:method:: writeMatchUp(...):

            Writes the matchup dataset to given file path.

        .. py:method:: write_matchup_file(...):

            Write match-up file from variable arrays

        .. py:method:: _date2seconds(self, times):

            Return times formatted as ``datetime.datetime`` objects as seconds since 1/1/1970

        .. py:method:: _return_w_matrix_variables(...):

            Produce match-up file W matrix variable arrays from lists of W and U matrices data.
    """

    def writeMatchUp(self, matchup, path):

        # todo - write this...
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


if __name__ == "__main__":
    pass
