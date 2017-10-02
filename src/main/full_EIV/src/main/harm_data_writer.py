"""
Created on Wed April 12  2017 17:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
from netCDF4 import Dataset
from os.path import join as pjoin
from numpy import zeros

class HarmOutput:
    """
    Class to generate object to store and write output harmonisation data

    Sample Code:

    > Opening data

    HOut = HarmOutput("some/path/harm.nc", ["harm_res_S1_S2.nc", "harm_res_S2_S3.nc", ...])

    > Writing data

    # Open empty HarmOutput object
    HOut = HarmOutput

    # Populate empty object with data
    HOut.parameter = a
    HOut.parameter_sensors = Ia
    HOut.parameter_covariance = V
    HOut.cost = F
    HOut.degrees_of_freedom = v
    HOut.cost_p_value = p
    HOut.lm = lm
    HOut.k_res = k_res
    HOut.H_res = H_res
    HOut.software = "MM"
    HOut.software_version = "V.V"
    HOut.software_tag = "TTTTTTT"
    HOut.job_id = "CC"
    HOut.matchup_dataset = "INSTR_TYPE_M_CORREL_MCN_YYYYMMDD_YYYYMMDD"

    # Write with save method
    HOut.save("some/directory")

    """

    def __init__(self, harmonisation_output_file_path=None, harmonisation_residual_file_paths=None):
        """
        Initialise HarmOut object to store harmonisation output data

        :type harmonisation_output_file_path: str
        :param harmonisation_output_file_path: path of harmonisation output file

        :type harmonisation_residual_file_paths: list:str
        :param harmonisation_residual_file_paths: list of path of harmonisation residual file paths

        :Attributes:
            :parameter: numpy.ndarray
                Harmonisation parameters

            :parameter_sensors: numpy.ndarray
                Sensors associated with harmonisation parameters

            :parameter_covariance_matrix: numpy.ndarray
                Harmonisation parameter covariance matrix

            :cost: float
                Sum of squares of the final value of the objective/cost function

            :cost_p_value: float
                Chi-squared probability

            :cost_dof:
                Cost/objective function degrees of freedom

            :software: str
                Software implementation abbreviation, i.e.:
                FO - FastOpt
                EV - NPL Errors in Variables
                OM - NPL ODR+MC Approach

            :software_version: str
                Version number of software (format V.V)

            :software_tag: str
                VCS tag for software implementation (format TTTTTTT)

            :job_id: str
                ID of harmonisation job

            :matchup_dataset: str
                Abbreviated harmonisation dataset name, of the form:
                INSTR_TYPE_M_CORREL_MCN_YYYYMMDD_YYYYMMDD

                where,
                INSTR - instrument name
                TYPE - dataset type
                M - number of harmonised parameters in sensor model
                CORREL - correlation structures present in data
                MCN - Monte Carlo trial number (filled with ___ if not a MC dataset)

            :lm: numpy.ndarray
                match-up series description

            :k_res: numpy.ndarray
                k residuals

            :H_res: numpy.ndarray
                data residual if available

            :additional_attributes: *dict*

            *(optional)* Dictionary of additional attributes to give to harmonisation product file for diagnostic
            purposes

        """

        # harmonisation results
        self.parameter = None
        self.parameter_sensors = None
        self.parameter_covariance_matrix = None

        # statistics
        self.cost = None
        self.cost_dof = None
        self.cost_p_value = None

        # residuals
        self.lm = None
        self.k_res = None
        self.H_res = None

        # software naming
        self.software = None
        self.software_version = None
        self.software_tag = None
        self.job_id = None

        # dataset naming
        self.matchup_dataset = None

        # addition test attributes
        self.additional_attributes = None

        if harmonisation_output_file_path is not None:
            self.parameter, self.parameter_sensors, self.parameter_covariance_matrix, \
            self.cost, self.cost_dof, self.cost_p_value, \
            self.software, self.software_version, self.software_tag, self.job_id,\
            self.matchup_dataset = self.open_harmonisation_output_file(harmonisation_output_file_path)

        if harmonisation_residual_file_paths is not None:
            self.lm, self.k_res, self.H_res, \
            self.software, self.software_version, self.software_tag, self.job_id, \
            self.matchup_dataset = self.open_harmonisation_residual_files(harmonisation_residual_file_paths)

    def open_harmonisation_output_file(self, harmonisation_output_file_path):
        """
        Return data from harmonisation output file

        :param harmonisation_output_file_path: str
            harmonisation output file path

        :return:
            :parameter: numpy.ndarray
                Harmonisation parameters

            :parameter_sensors: numpy.ndarray
                Sensors associated with harmonisation parameters

            :parameter_covariance_matrix: numpy.ndarray
                Harmonisation parameter covariance matrix

            :cost: float
                Sum of squares of the final value of the objective/cost function

            :cost_p_value: float
                Chi-squared probability

            :cost_dof:
                Cost/objective function degrees of freedom

            :software: str
                Software implementation abbreviation, i.e.:
                FO - FastOpt
                EV - NPL Errors in Variables
                OM - NPL ODR+MC Approach

            :software_version: str
                Version number of software (format V.V)

            :software_tag: str
                VCS tag for software implementation (format TTTTTTT)

            :job_id: str
                ID of harmonisation job

            :matchup_dataset: str
                Abbreviated harmonisation dataset name, of the form:
                INSTR_TYPE_M_CORREL_MCN_YYYYMMDD_YYYYMMDD

                where,
                INSTR - instrument name
                TYPE - dataset type
                M - number of harmonised parameters in sensor model
                CORREL - correlation structures present in data
                MCN - Monte Carlo trial number (filled with ___ if not a MC dataset)
        """

        # open netCDF file
        rootgrp = Dataset(harmonisation_output_file_path, "r")

        # get attribute data
        cost = rootgrp.cost
        cost_dof = rootgrp.cost_dof
        cost_p_value = rootgrp.cost_p_value
        software = rootgrp.software
        software_version = rootgrp.software_version
        software_tag = rootgrp.software_tag
        job_id = rootgrp.job_id
        matchup_dataset = rootgrp.matchup_dataset

        # get variable data
        parameter = rootgrp.variables['parameter'][:]
        parameter_sensors = rootgrp.variables['parameter_sensors'][:].astype(int)
        parameter_covariance_matrix = rootgrp.variables['parameter_covariance_matrix'][:, :]

        # Hacky solution to turning parameter_sensors into numbers from FastOpt
        if software == "FO":
            parameter_sensors_string = parameter_sensors[:]
            parameter_sensors = zeros(parameter_sensors_string.shape[0])
            for i, string in enumerate(parameter_sensors_string):
                parameter_sensors[i] = int("".join(string[1:3]))

        # close netCDF file
        rootgrp.close()

        return parameter, parameter_sensors, parameter_covariance_matrix, cost, cost_dof, cost_p_value,\
               software, software_version, software_tag, job_id, matchup_dataset

    def open_harmonisation_residual_files(self, harmonisation_residual_file_paths):
        """
        Return data from harmonisation residual file

        :param harmonisation_residual_file_paths: str
            list of harmonisation residual file paths to open

        :return:
            :lm: numpy.ndarray
                match-up series description

            :k_res: numpy.ndarray
                k residuals

            :H_res: numpy.ndarray
                data residual if available

             :software: str
                Software implementation abbreviation, i.e.:
                FO - FastOpt
                EV - NPL Errors in Variables
                OM - NPL ODR+MC Approach

            :software_version: str
                Version number of software (format V.V)

            :software_tag: str
                VCS tag for software implementation (format TTTTTTT)

            :job_id: str
                ID of harmonisation job

            :matchup_dataset: str
                Abbreviated harmonisation dataset name, of the form:
                INSTR_TYPE_M_CORREL_MCN_YYYYMMDD_YYYYMMDD

                where,
                INSTR - instrument name
                TYPE - dataset type
                M - number of harmonised parameters in sensor model
                CORREL - correlation structures present in data
                MCN - Monte Carlo trial number (filled with ___ if not a MC dataset)
        """

        # first get required descriptive data
        lm = zeros((len(harmonisation_residual_file_paths), 3))
        for i, harmonisation_residual_file_path in enumerate(harmonisation_residual_file_paths):

            # open file
            rootgrp = Dataset(harmonisation_residual_file_path, "r")

            # get attributes from first file
            if i == 0:
                software = rootgrp.software
                software_version = rootgrp.software_version
                software_tag = rootgrp.software_tag
                job_id = rootgrp.job_id
                matchup_dataset = rootgrp.matchup_dataset

                # check if H_res data available and get dimensions
                data_res = False
                if 'H_res' in rootgrp.variables.keys():
                    data_res = True
                    n_col = rootgrp.variables['H_res'].shape[1]

            # build lm array from data file
            sensor_i_name = rootgrp.sensor_i_name
            sensor_j_name = rootgrp.sensor_j_name

            # Hacky solution to turning parameter_sensors into numbers from FastOpt
            # This really will need changing
            if software == "FO":
                if matchup_dataset[:5] == "AVHRR":

                    if sensor_i_name == "m02":
                        sensor_i_name = -1
                    else:
                        sensor_i_name = int(sensor_i_name[1:])

                    sensor_j_name = sensor_j_name[1:]
                else:
                    print "Sensor name problem"

            lm[i, 0] = sensor_i_name
            lm[i, 1] = sensor_j_name
            lm[i, 2] = rootgrp.variables['k_res'].shape[0]

            rootgrp.close()

        # initialise data arrays
        total = 0
        idx = [0]
        for n in lm[:, 2]:
            total += int(n)
            idx.append(total)
        n_mu = idx[-1]

        k_res = zeros(n_mu)
        if data_res:
            H_res = zeros((n_mu, n_col))

        for i, harmonisation_residual_file_path in enumerate(harmonisation_residual_file_paths):

            istart = idx[i]
            iend = idx[i+1]

            # open netCDF file
            rootgrp = Dataset(harmonisation_residual_file_path, "r")

            # get variable data
            k_res[istart:iend] = rootgrp.variables['k_res'][:]
            if data_res:
                H_res[istart:iend, :] = rootgrp.variables['H_res'][:, :]
            else:
                H_res = None

            # close netCDF file
            rootgrp.close()

        return lm, k_res, H_res, software, software_version, software_tag, job_id, matchup_dataset

    def save(self, directory, res=True):
        """
        Write harmonisation output to file

        :type directory: str
        :param directory: directory to write files to

        :type res: bool
        :param res: Switch to turn off writing harmonisation output residual file

        """

        self.write_harmonisation_output_file(directory,
                                             self.parameter, self.parameter_sensors, self.parameter_covariance_matrix,
                                             self.cost, self.cost_dof, self.cost_p_value,
                                             self.software, self.software_version, self.software_tag, self.job_id,
                                             self.matchup_dataset, self.additional_attributes)
        if res:
            self.write_harmonisation_residual_files(directory,
                                                    self.software, self.software_version, self.software_tag, self.job_id,
                                                    self.matchup_dataset, self.lm, self.k_res, self.H_res,)

    def write_harmonisation_output_file(self, directory,
                                              parameter, parameter_sensors, parameter_covariance_matrix,
                                              cost, cost_dof, cost_p_value,
                                              software, software_version, software_tag, job_id,
                                              matchup_dataset, additional_attributes=None):
        """
        Write harmonisation output file as defined in "Harmonisation Output File Format Definition" document

        :type directory: str
        :param directory: directory to write file to

        :type parameter: numpy.ndarray
        :param parameter: Harmonisation parameters

        :param parameter_sensors: numpy.ndarray
            Sensors associated with harmonisation parameters

        :type parameter_covariance_matrix: numpy.ndarray
        :param parameter_covariance_matrix: Harmonisation parameter covariance matrix

        :type cost: float
        :param cost: Sum of squares of the final value of the objective/cost function

        :type cost_p_value: type
        :param cost_p_value: Chi-squared probability

        :type cost_dof: float
        :param cost_dof: Cost/objective function degrees of freedom

        :type software: str
        :param software: Software implementation abbreviation, i.e.:
        * FO - FastOpt
        * EV - NPL Errors in Variables
        * OM - NPL ODR+MC Approach

        :type software_version: str
        :param software_version: Version number of software (format V.V)

        :type software_tag: str
        :param software_tag: VCS tag for software implementation (format TTTTTTT)

        :type job_id: str
        :param job_id: ID of harmonisation job

        :type matchup_dataset: str
        :param matchup_dataset: Abbreviated harmonisation dataset name, of the form INSTR_TYPE_M_CORREL_MCN_YYYYMMDD_YYYYMMDD
        where:
        * INSTR - instrument name
        * TYPE - dataset type
        * M - number of harmonised parameters in sensor model
        * CORREL - correlation structures present in data
        * MCN - Monte Carlo trial number (filled with ___ if not a MC dataset)

        :type additional_attributes: dict
        :param additional_attributes: Dictionary of additional attributes to give to harmonisation product file for diagnostic purposes
        """

        # define file path
        fname = "_".join(("harm", software, software_version, software_tag, job_id, matchup_dataset)) + ".nc"
        path = pjoin(directory, fname)

        # open netCDF file
        rootgrp = Dataset(path, 'w')

        # set attributes
        rootgrp.cost = cost
        rootgrp.cost_dof = cost_dof
        rootgrp.cost_p_value = cost_p_value
        rootgrp.software = software
        rootgrp.software_version = software_version
        rootgrp.software_tag = software_tag
        rootgrp.job_id = job_id
        rootgrp.matchup_dataset = matchup_dataset

        if additional_attributes is not None:
            rootgrp.setncatts(additional_attributes)

        # create dimensions
        n = rootgrp.createDimension('n', len(parameter))

        # create variables

        # > harmonised parameter variable, a
        parameter_var = rootgrp.createVariable('parameter', 'f8', ('n',), zlib=True, complevel=9)
        parameter_var.description = "Harmonisation parameters"

        # > harmonised parameter covariance matrix, V
        parameter_covariance_matrix_var = rootgrp.createVariable('parameter_covariance_matrix', 'f8', ('n', 'n',), zlib=True, complevel=9)
        parameter_covariance_matrix_var.description = 'Harmonisation parameter covariance matrix'

        # harmonisation parameter sensor name variable, Ia
        parameter_sensors_var = rootgrp.createVariable('parameter_sensors', 'f8', ('n',), zlib=True, complevel=9)
        parameter_sensors_var.description = "Sensors associated with harmonisation parameters"

        # store data
        parameter_var[:] = parameter[:]
        parameter_covariance_matrix_var[:] = parameter_covariance_matrix[:]
        parameter_sensors_var[:] = parameter_sensors[:]

        # close netCDF file
        rootgrp.close()

        return 0

    def write_harmonisation_residual_files(self, directory,
                                           software, software_version, software_tag, job_id, matchup_dataset,
                                           lm, k_res, H_res=None):
        """
        Write harmonisation set of residual files as defined in "Harmonisation Output File Format Definition" document

        :param lm: numpy.ndarray
            match-up series description

        :param k_res: numpy.ndarray
            k residuals

        :param H_res: numpy.ndarray
            data residual if available

         :param software: str
            Software implementation abbreviation, i.e.:
            FO - FastOpt
            EV - NPL Errors in Variables
            OM - NPL ODR+MC Approach

        :param software_version: str
            Version number of software (format V.V)

        :param software_tag: str
            VCS tag for software implementation (format TTTTTTT)

        :param job_id: str
            ID of harmonisation job

        :param matchup_dataset: str
            Abbreviated harmonisation dataset name, of the form:
            INSTR_TYPE_M_CORREL_MCN_YYYYMMDD_YYYYMMDD

            where,
            INSTR - instrument name
            TYPE - dataset type
            M - number of harmonised parameters in sensor model
            CORREL - correlation structures present in data
            MCN - Monte Carlo trial number (filled with ___ if not a MC dataset)

        """

        total = 0
        idx = [0]
        for n in lm[:, 2]:
            total += n
            idx.append(total)
        n_mu = idx[-1]

        # Write file for each match-up series
        for i, pair in enumerate(lm):

            # Get required data from lm variable
            sensor_i = pair[0]
            sensor_j = pair[1]
            n_mu = pair[2]

            istart = idx[i]
            iend = idx[i+1]

            # define file path
            fname = "_".join(("harm", software, software_version, software_tag, job_id,
                              matchup_dataset, "res", str(sensor_i), str(sensor_j))) + ".nc"
            path = pjoin(directory, fname)

            # open netCDF file
            rootgrp = Dataset(path, 'w')

            # set attributes
            rootgrp.software = software
            rootgrp.software_version = software_version
            rootgrp.software_tag = software_tag
            rootgrp.job_id = job_id
            rootgrp.matchup_dataset = matchup_dataset
            rootgrp.sensor_i_name = sensor_i
            rootgrp.sensor_j_name = sensor_j

            # create dimensions
            m = rootgrp.createDimension('m', n_mu)
            if H_res is not None:
                n_col = rootgrp.createDimension('n_col', H_res.shape[1])

            # create variables

            # > Harmonisation match-up adjustment factor residuals
            k_res_var = rootgrp.createVariable('k_res', 'f8', ('m',), zlib=True, complevel=9)
            k_res_var.description = "k residuals"

            if H_res is not None:
                # > Harmonisation match-up date residuals
                H_res_var = rootgrp.createVariable('H_res', 'f8', ('m', 'n_col',), zlib=True, complevel=9)
                H_res_var.description = "Data residuals"

            # store data

            k_res_var[:] = k_res[istart:iend]
            if H_res is not None:
                H_res_var[:, :] = H_res[istart:iend, :]

            # close netCDF file
            rootgrp.close()

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

    def main():
        return 0

    main()