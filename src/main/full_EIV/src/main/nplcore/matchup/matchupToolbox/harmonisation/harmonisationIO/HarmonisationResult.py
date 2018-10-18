"""
HarmonisationResult Class - In memory representation of the results of the harmonisation of a MatchUp data object
"""

'''___Built-In Modules___'''
from os.path import split, splitext
from os.path import join as pjoin

'''___Third-Party Modules___'''
from netCDF4 import Dataset, stringtochar, chartostring
from numpy import zeros, ndarray, array_str, asarray

'''___NPL Modules___'''

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "19/11/2017"
__credits__ = ["Jonathan Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class HarmonisationResult:
    """
    In-memory representation of the results of the harmonisation of a MatchUp data object

    Sample Code:

    .. code-block:: python

        # Read data
        HarmonisationOutput = HarmonisationResult("some/path/harm.nc", ["harm_res_S1_S2.nc", "harm_res_S2_S3.nc", ...])

        # Write data
        path = /some/directory/to/write/to
        HarmOutput.save(path)

    :Attributes:
        .. py:attribute:: parameter

            *numpy.ndarray*

            Harmonisation parameters

        .. py:attribute:: parameter_sensors

            *list*

            Mapping of sensor ID to associated parameters in parameter array

        .. py:attribute:: parameter_covariance_matrix

            *numpy.ndarray*

            Harmonisation parameter covariance matrix

        .. py:attribute:: cost

            *float*

            Sum of squares of the final value of the objective/cost function

        .. py:attribute:: cost_p_value

            *float*

            Chi-squared probability

        .. py:attribute:: cost_dof

            *float*

            Cost/objective function degrees of freedom


        .. py:attribute:: idx

            *dict:lists*

            Dictionary of indexing lists which define the 1D structure of the residual data

        .. py:attribute:: values_res

            *numpy.ndarray*

            Data residuals if available (not provided by all fitting routines) in 1D data structure

        .. py:attribute:: k_res

            *numpy.ndarray*

            k residuals in 1D data structure

        .. py:attribute:: additional_attributes

            *dict*

            Dictionary of additional attributes to give to write to harmonisation product file for labelling or
            diagnostic purposes

    :Methods:
        .. py:method:: open_harmonisation_output_file(...):

            Return data from harmonisation output file

        .. py:method:: open_harmonisation_residual_files(...):

            Return data from harmonisation residual file[s]

        .. py:method:: save(...):

            Write harmonisation output to file using functionality of self.write_harmonisation_output_file and
            self.write_harmonisation_residual file

        .. py:method:: write_harmonisation_output_file(...):

            Write data to harmonisation output file

        .. py:method:: write_harmonisation_residual_files(...):

            Write data to harmonisation residual file[s]
    """

    def __init__(self, harmonisation_output_file_path=None, harmonisation_residual_file_paths=None):
        """
        Initialise HarmonisationResult object to store harmonisation output data

        :type harmonisation_output_file_path: str
        :param harmonisation_output_file_path: path of harmonisation output file

        :type harmonisation_residual_file_paths: list:str
        :param harmonisation_residual_file_paths: list of path of harmonisation residual file paths
        """

        # 1. Initialise attributes -------------------------------------------------------------------------------------

        # a. Harmonisation results
        self.parameter = None
        self.parameter_sensors = None
        self.parameter_covariance_matrix = None

        # b. Statistics
        self.cost = None
        self.cost_dof = None
        self.cost_p_value = None

        # c. Residuals
        self.idx = None
        self.values_res = None
        self.ks_res = None

        # d. Additional labelling/test attributes
        self.additional_attributes = {}

        # e. Additional variables
        self.additional_variables = {}

        # --------------------------------------------------------------------------------------------------------------

        # 2. Open harmonisation residual files if path provided --------------------------------------------------------
        if harmonisation_residual_file_paths is not None:
            self.values_res, self.ks_res, \
                self.idx, self.additional_attributes \
                    = self.open_harmonisation_residual_files(harmonisation_residual_file_paths)
        # --------------------------------------------------------------------------------------------------------------

        # 3. Open harmonisation output file if path provided -----------------------------------------------------------
        if harmonisation_output_file_path is not None:
            self.parameter, self.parameter_sensors, \
                self.parameter_covariance_matrix, self.cost, \
                    self.cost_dof, self.cost_p_value, \
                        self.additional_attributes, self.additional_variables \
                            = self.open_harmonisation_output_file(harmonisation_output_file_path)
        # --------------------------------------------------------------------------------------------------------------

    def open_harmonisation_output_file(self, harmonisation_output_file_path):
        """
        Return data from harmonisation output file

        :type harmonisation_output_file_path: str
        :param harmonisation_output_file_path: harmonisation output file path

        :return:
            :parameter:

            *numpy.ndarray*

            Harmonisation parameters

            :parameter_sensors:

            *list*

            Sensors associated with harmonisation parameters

            :parameter_covariance_matrix:

            *numpy.ndarray*

            Harmonisation parameter covariance matrix

            :cost: *float*

            Sum of squares of the final value of the objective/cost function

            :cost_p_value: *float*

            Chi-squared probability

            :cost_dof: *float*

            Cost/objective function degrees of freedom

            :additional_attributes: *dict*

            Dictionary of additional attributes to give to write to harmonisation product file for labelling or
            diagnostic purposes
        """

        # 1. Open file -------------------------------------------------------------------------------------------------
        dataset = Dataset(harmonisation_output_file_path, "r")
        # --------------------------------------------------------------------------------------------------------------

        # 2. Get attributes --------------------------------------------------------------------------------------------
        # a. Open required attributes
        cost = dataset.cost
        cost_dof = dataset.cost_dof
        cost_p_value = dataset.cost_p_value

        # b. Open additinal attributes
        additional_attributes_keys = [str(attr) for attr in dataset.ncattrs()
                                                    if (attr not in ["cost", "cost_dof", "cost_p_value"])
                                                        and (attr[0] != "_")]
        additional_attributes = {}
        for additional_attributes_key in additional_attributes_keys:
            additional_attributes[additional_attributes_key] = dataset.getncattr(additional_attributes_key)
        # --------------------------------------------------------------------------------------------------------------

        # 3. Get variable data -----------------------------------------------------------------------------------------
        # > Parameter
        parameter = dataset.variables['parameter'][:]

        # > Parameter sensor
        parameter_sensors = [s.strip() for s in list(chartostring(dataset.variables['parameter_sensors'][:]))]

        # > Parameter covariance matrix
        parameter_covariance_matrix = dataset.variables['parameter_covariance_matrix'][:]
        # --------------------------------------------------------------------------------------------------------------

        # > 4. Open additional variables -------------------------------------------------------------------------------

        additional_variables = {}
        additional_variable_names = [variable for variable in dataset.variables.keys()
                                                  if variable not in ["parameter", "parameter_sensors",
                                                                      "parameter_covariance_matrix"]]

        for additional_variable_name in additional_variable_names:
            additional_variables[additional_variable_name] = dataset.variables[additional_variable_name][:]

        # --------------------------------------------------------------------------------------------------------------

        # 4. Close netCDF file -----------------------------------------------------------------------------------------
        dataset.close()
        # --------------------------------------------------------------------------------------------------------------

        return parameter, parameter_sensors, parameter_covariance_matrix,\
                    cost, cost_dof, cost_p_value, additional_attributes, additional_variables

    def open_harmonisation_residual_files(self, harmonisation_residual_file_paths):
        """
        Return data from harmonisation residual file

        :param harmonisation_residual_file_paths: str
            list of harmonisation residual file paths to open

        :return:
            :values_res: *numpy.ndarray*

            Data residuals if available (not provided by all fitting routines) in 1D data structure

            :ks_res: *numpy.ndarray*

            k residuals

            :idx: *dict*

            Dictionary of indexing lists which define the 1D structure of the residual data

            :additional_attributes: *dict*

            Dictionary of additional attributes to give to write to harmonisation product file for labelling or
            diagnostic purposes
        """

        # Input residual data is restructured into a 1D array structured in blocks of data in the following way:
        #
        # values_res = [ Xref_res | X_1_res(1,2) |...| X_1_res(N,S-1) |
        #                            | X_1_res(L,S) |...| X_m_res(1,2) |...| X_m_res(N,S-1)| X_m_res(L,S) ]
        #
        # where:
        # - Xref_res - reference sensor measurand fitting residual (if dataset contains reference sensor)
        # - X_m_res - model covariate fitting residual of number m of a total m covariates
        # - (N, s) - covariate indices:
        #            > N - match-up series index of N of a total N match-up series
        #            > s - sensor number of S sensors
        #
        # Adjustment Factor residuals, K_res, data is considered separately and stored in a separate 1D array in
        # blocks by match-up series in the following way:
        #
        # k_res = [ k_res(1) |...| k_res(N)]
        #
        # NB: Not all fitting routines return values_res

        ################################################################################################################
        # 1. Open additional_attributes
        ################################################################################################################

        # a. Open datasets
        datasets = [Dataset(path, 'r') for path in harmonisation_residual_file_paths]

        # b. Read attribute values from first harmonisation residual file
        additional_attributes_keys = [str(attr) for attr in datasets[0].ncattrs() if attr[0] != "_"]
        additional_attributes = {}
        for additional_attributes_key in additional_attributes_keys:
            additional_attributes[additional_attributes_key] = datasets[0].getncattr(additional_attributes_key)

        ################################################################################################################
        # 2. Define 1D Data Structure
        ################################################################################################################

        # Check if X1 and X2 residuals contained in harmonisation residual datasets
        values_res_available = False
        values_res_per_dataset = [True if ('X1_res' in d.variables) and ('X2_res' in d.variables) else False
                                  for d in datasets]
        if all(values_res_per_dataset):
            values_res_available = True

        # Obtaining lists of indices to describe data structure
        # a. per match-up series labels --------------------------------------------------------------------------------
        # Initialise idx lists
        num_matchups = len(datasets)
        Nm = [0] * num_matchups              # number of match-up per match-up series
        cNm = [0] * (num_matchups + 1)       # cumulative number of match-ups by match-up series
        Im = [0] * num_matchups              # sensor pairs per match up by ID (starting from 0 and increasing)
        sensors = []                         # sensor name corresponding to ID used in Im (fix reference sensor to be 0)
        sensor_ms = []

        # Populate idx lists by match-up series
        for i, dataset in enumerate(datasets):
            try:
                Nm[i] = dataset.dimensions['M'].size      # add number of match-ups in match-up series to list
            except KeyError:
                Nm[i] = dataset.dimensions['m'].size

            cNm[i+1] = cNm[i] + Nm[i]                 # calculate cumulative total of match-ups

            # Determine sensor IDs
            Im[i] = [0, 0]

            try:
                sensor_names_i = [dataset.sensor_1_name, dataset.sensor_2_name]
            except AttributeError:
                sensor_names_i = [dataset.sensor_i, dataset.sensor_j]

            if values_res_available:
                sensor_m_i = [dataset.dimensions['m1'].size, dataset.dimensions['m2'].size]
            else:
                sensor_m_i = [0, 0]

            for j, (sensor, sensor_m) in enumerate(zip(sensor_names_i, sensor_m_i)):
                if sensor in sensors:
                    Im[i][j] = sensors.index(sensor)
                else:
                    sensors.append(sensor)
                    sensor_ms.append(sensor_m)
                    Im[i][j] = sensors.index(sensor)

        # Form idx dictionary
        idx = {"Nm": Nm,
               "cNm": cNm,
               "Im": Im,
               "sensors": sensors,
               "sensor_ms": sensor_ms}
        # --------------------------------------------------------------------------------------------------------------

        # b. determine 1D data block structure definition indices ------------------------------------------------------

        if values_res_available:
            # generate required data
            full_sensor_list = [num for pair in Im for num in pair]  # full list of sensors (ordered by match-up series)
            num_ref_matchup = full_sensor_list.count(0)  # number of reference-sensor match-up series
            sensor_list = [sensor for sensor in full_sensor_list if
                           sensor != 0]  # list of sensors (excluding reference)

            # Create New Indices:
            # > n_sensor - list of sensor number of consecuative data blocks
            n_sensor = [int(0)] * num_ref_matchup + [int(sensor) for m_i in range(1, max(sensor_ms) + 1)
                                                     for sensor in sensor_list
                                                     if sensor_ms[sensor] >= m_i]

            # > n_mu - list of match-up series number of consecuative data blocks
            # For reference sensors
            n_mu_refs = [int(i + 1) for i, pair in enumerate(Im) if 0 in pair]

            # For series sensors
            n_mu_sensors = [int(i + 1) for m_i in range(1, max(sensor_ms) + 1)
                            for i, pair in enumerate(Im)
                            for sensor in pair
                            if (sensor_ms[sensor] >= m_i) and (sensor != 0)]

            # Combine
            n_mu = n_mu_refs + n_mu_sensors

            # > n_cov - list of covariate number of consecuative blocks of data
            # For reference sensors
            n_cov_refs = [int(1)] * num_ref_matchup

            # For series sensors
            n_cov_sensors = [m_i for m_i in range(1, max(sensor_ms) + 1)
                             for i, pair in enumerate(Im)
                             for sensor in pair
                             if (sensor_ms[sensor] >= m_i) and (sensor != 0)]

            # Combine
            n_cov = n_cov_refs + n_cov_sensors

            # - N_var - list of total number of variables in consecuative data blocks
            #          (initially number of match-ups before modification)
            N_var = [Nm[n - 1] for n in n_mu]

            # - idxs - index in 1D harmonisation data array of first element of each data block (and last element index)
            idxs = [0]
            total = 0
            for N in N_var:
                total += N
                idxs.append(int(total))

            # Add to idx dictionary
            idx["n_sensor"] = n_sensor
            idx["n_mu"] = n_mu
            idx["n_cov"] = n_cov
            idx["N_var"] = N_var
            idx["idx"] = idxs
        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 3. Open residual data
        ################################################################################################################

        # a. Read data residuals (if available) ------------------------------------------------------------------------
        values_res = None
        if values_res_available:
            values_res = zeros(idx['idx'][-1])  # initialise array to store 1D organised match-up residual data

            # sensor number, match up number and covariate number per data block as a tuple for indexing
            block_idxs = [(i, j, k) for i, j, k in zip(idx['n_sensor'], idx['n_mu'], idx['n_cov'])]

            # Read data by match-up series
            for i_matchup, dataset in enumerate(datasets):

                # Read match-up series data by sensor
                for i_sensor_pair, n_sensor in enumerate(idx["Im"][i_matchup]):

                    # Read sensor data by covariate
                    m = dataset.dimensions['m' + str(i_sensor_pair + 1)].size
                    for i_cov in range(m):

                        # Find location of block for this covariate in 1D data structure
                        block_idx = block_idxs.index((n_sensor, i_matchup + 1, i_cov + 1))

                        #  Add data to values_res array
                        i_start = idx['idx'][block_idx]
                        i_end = idx['idx'][block_idx + 1]
                        values_res[i_start:i_end] = dataset.variables['X' + str(i_sensor_pair + 1) + '_res'][:, i_cov]
        # --------------------------------------------------------------------------------------------------------------

        # b. Read adjustment factor residuals --------------------------------------------------------------------------
        ks_res = zeros(idx['cNm'][-1])  # 1D organised match-up k residual data

        # Read data by match-up series
        for i_matchup, dataset in enumerate(datasets):

                # add K_res data to ks_res array
                i_start = idx['cNm'][i_matchup]
                i_end = idx['cNm'][i_matchup + 1]
                ks_res[i_start:i_end] = dataset.variables['k_res'][:]
        # -------------------------------------------------- -----------------------------------------------------------

        ################################################################################################################
        # 4. Close all datasets
        ################################################################################################################

        for dataset in datasets:
            dataset.close()

        return values_res, ks_res, idx, additional_attributes

    def save(self, path, save_residuals=True):
        """
        Write harmonisation output to file using functionality of self.write_harmonisation_output_file and
        self.write_harmonisation_residual file

        :type directory: str
        :param directory: directory to write files to

        :type save_residuals: bool
        :param save_residuals: Switch to turn off writing harmonisation output residual file

        """

        self.write_harmonisation_output_file(path,
                                             self.parameter, self.parameter_sensors, self.parameter_covariance_matrix,
                                             self.cost, self.cost_dof, self.cost_p_value,
                                             self.additional_attributes, self.additional_variables)
        if save_residuals:
            self.write_harmonisation_residual_files(path, self.values_res, self.ks_res, self.idx,
                                                    self.additional_attributes)

    def write_harmonisation_output_file(self, path, parameter, parameter_sensors, parameter_covariance_matrix,
                                        cost, cost_dof, cost_p_value, additional_attributes, additional_variables):
        """
        Write harmonisation output file

        :type path: str
        :param path: path to write file to

        :type parameter: numpy.ndarray
        :param parameter: Harmonisation parameters

        :type parameter_sensors: list
        :param parameter_sensors: Sensors associated with harmonisation parameters

        :type parameter_covariance_matrix: numpy.ndarray
        :param parameter_covariance_matrix: Harmonisation parameter covariance matrix

        :type cost: float
        :param cost: Sum of squares of the final value of the objective/cost function

        :type cost_p_value: type
        :param cost_p_value: Chi-squared probability

        :type cost_dof: float
        :param cost_dof: Cost/objective function degrees of freedom

        :type additional_attributes: dict
        :param additional_attributes: Dictionary of additional attributes for labelling/diagnostic purposes
        """

        # 1. Open netCDF file
        rootgrp = Dataset(path, 'w')

        # 2. Set attributes
        rootgrp.cost = cost
        rootgrp.cost_dof = cost_dof
        rootgrp.cost_p_value = cost_p_value
        rootgrp.setncatts(additional_attributes)

        # 3. Create dimensions
        n = rootgrp.createDimension('n', len(parameter))
        max_l_name = 80
        l_name = rootgrp.createDimension('l_name', max_l_name)

        # 4. Create variables

        # > Harmonised parameter variable
        parameter_var = rootgrp.createVariable('parameter', 'f8', ('n',), zlib=True, complevel=9)
        parameter_var.description = "Harmonisation parameters"
        parameter_var[:] = parameter[:]

        # > Harmonised parameter covariance matrix
        parameter_covariance_matrix_var = rootgrp.createVariable('parameter_covariance_matrix', 'f8', ('n', 'n',),
                                                                 zlib=True, complevel=9)
        parameter_covariance_matrix_var.description = 'Harmonisation parameter covariance matrix'
        parameter_covariance_matrix_var[:] = parameter_covariance_matrix[:]

        # > Harmonisation parameter sensor name variable
        fmt_str = '{0: <'+str(max_l_name)+'}'
        parameter_sensors = asarray([fmt_str.format(p) for p in parameter_sensors])
        parameter_sensors_var = rootgrp.createVariable('parameter_sensors', 'S1', ('n', 'l_name'), zlib=True,
                                                       complevel=9)
        parameter_sensors_var.description = "Sensors associated with harmonisation parameters"
        parameter_sensors_var[:] = stringtochar(parameter_sensors[:])

        # > additional variables - non-required variable to add to harmonisation input files
        if additional_variables is not None:
            for variable in additional_variables.keys():
                for i, dim in enumerate(additional_variables[variable]['dim']):
                    if dim not in rootgrp.dimensions.keys():
                        new_dim = rootgrp.createDimension(dim, len(additional_variables[variable]['data']))

                if additional_variables[variable]['dtype'] == 'S1':
                    fmt_str = '{0: <' + str(max_l_name) + '}'
                    additional_variables[variable]['data'] = asarray([fmt_str.format(s) for s in
                                                                      additional_variables[variable]['data']])
                    additional_var = rootgrp.createVariable(variable, additional_variables[variable]['dtype'],
                                                            (additional_variables[variable]['dim'], 'l_name'),
                                                            zlib=True, complevel=9)
                    additional_var[:] = stringtochar(additional_variables[variable]['data'][:])
                else:
                    additional_var = rootgrp.createVariable(variable, additional_variables[variable]['dtype'],
                                                            additional_variables[variable]['dim'],
                                                            zlib=True, complevel=9)
                    additional_var[:] = additional_variables[variable]['data'][:]

                if "Description" in additional_variables[variable].keys():
                    additional_var.Description = additional_variables[variable]['Description']

        # 6. Close netCDF file
        rootgrp.close()

        return 0

    def write_harmonisation_residual_files(self, path, values_res, ks_res, idx, additional_attributes):
        """
        Write harmonisation set of residual files as defined in "Harmonisation Output File Format Definition" document

        :type path: str
        :param path: path of associated harmonisation output file, to be suffixed with residual name

        :type values_res: numpy.ndarray
        :param values_res: Data residuals if available (not provided by all fitting routines) in 1D data structure. Set
        as None if not available

        :type ks_res: numpy.ndarray
        :param ks_res: k residuals

        :type idx: dict
        :param idx: Dictionary of indexing lists which define the 1D structure of the residual data

        :type additional_attributes: dict
        :param additional_attributes: Dictionary of additional attributes to give to write to harmonisation product file
        for labelling or diagnostic purposes
        """

        # Residual file path components
        directory, output_fname = split(path)
        output_fname_noext, output_fname_ext = splitext(output_fname)

        # sensor number, match up number and covariate number per data block as a tuple for indexing
        block_idxs = [(i, j, k) for i, j, k in zip(idx['n_sensor'], idx['n_mu'], idx['n_cov'])]

        # Write file per match-up series
        for i_matchup, n_sensor_pair in enumerate(idx["Im"]):

            # 1. Initialise file ---------------------------------------------------------------------------------------
            sensor_name_pair = [idx['sensors'][n_sensor_pair[0]], idx['sensors'][n_sensor_pair[1]]]

            # Define filename
            fname_res = "_".join((output_fname_noext, "res", str(sensor_name_pair[0]), str(sensor_name_pair[1])))\
                        + output_fname_ext

            path_res = pjoin(directory, fname_res)

            # Open netCDF file
            dataset = Dataset(path_res, 'w')

            # Set attributes
            dataset.sensor_1_name = sensor_name_pair[0]
            dataset.sensor_2_name = sensor_name_pair[1]

            # Create dimensions
            M_matchup = idx['Nm'][i_matchup]
            M_matchup_dim = dataset.createDimension('M', M_matchup)
            # ----------------------------------------------------------------------------------------------------------

            # 2. Write values_res data if available --------------------------------------------------------------------
            if values_res is not None:

                # Write data per sensor
                for i_sensor_pair, (n_sensor, sensor_name) in enumerate(zip(n_sensor_pair, sensor_name_pair)):

                    # a. Build data array
                    # Initialise X_res array
                    m_sensor = idx['sensor_ms'][n_sensor]
                    X_res = zeros((M_matchup, m_sensor))

                    # Populate array by covariate
                    for i_cov in range(m_sensor):

                        # Find data location within 1D structure
                        block_idx = block_idxs.index((n_sensor, i_matchup + 1, i_cov + 1))
                        i_start = idx['idx'][block_idx]
                        i_end = idx['idx'][block_idx + 1]

                        # Populate array
                        X_res[:, i_cov] = values_res[i_start:i_end]

                    # ii. Add to file
                    # Create extra dimension
                    m_sensor_dim = dataset.createDimension("m"+str(i_sensor_pair+1), m_sensor)

                    # Create variable
                    X_res_var = dataset.createVariable('X'+str(i_sensor_pair+1)+'_res', 'f8',
                                                       ('M', "m"+str(i_sensor_pair+1)), zlib=True, complevel=9)
                    X_res_var.description = "X"+str(i_sensor_pair+1)+" harmonisation fitting residuals"

                    # Store data
                    X_res_var = X_res[:]

            # 3. Write ks_res data if available ------------------------------------------------------------------------

            # a. Create k_residual variable
            k_res_var = dataset.createVariable('k_res', 'f8', ('M',), zlib=True, complevel=9)
            k_res_var.description = "k harmonisation fitting residuals"

            # b. Store data
            # Locate data
            i_start = idx['cNm'][i_matchup]
            i_end = idx['cNm'][i_matchup + 1]

            # Add to file
            k_res_var[:] = ks_res[i_start:i_end]
            # ----------------------------------------------------------------------------------------------------------

            # Close file
            dataset.close()

        return 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        attrs = vars(self).keys()
        for attr in attrs:
            delattr(self, attr)

        self.closed = True


if __name__ == "__main__":
    pass
