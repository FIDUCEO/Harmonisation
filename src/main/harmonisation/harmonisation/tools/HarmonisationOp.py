"""
Main Operator of FIDUCEO Harmonisation Implementation
"""

'''___Built-In Modules___'''
from glob import glob
from copy import deepcopy
from os import makedirs

'''___Third Party Modules___'''
from numpy import append, zeros, loadtxt, mean, std, cov, savetxt
from numpy.random import normal
from numpy import all as nall

'''___harmonisation Modules___'''
from common import *
from harmonisation import MatchUp, HarmonisationEIV, HarmonisationResult, Sample2Ind, Transform2NormInd
from harmonisation.core.matchup.matchupToolbox.harmonisation.\
    harmonisationProcessing.harmonisation_eiv.GN_algo import GNAlgo
from harmonisation.core.matchup.matchupToolbox.harmonisation.\
    harmonisationProcessing.harmonisation_eiv_pc.pc_algo import PCAlgo


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Arta Dillo", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


'''___Constants___'''

# Tolerances
TOLPC = 1e-6    # Preconditioner alogrithm convergence tolerance
TOL = 1e-6      # Gauss-Newton algorithm convergence tolerance
TOLA = 1e-8     # Gauss-Newton algorithm LSMR tolerance tolA
TOLB = 1e8      # Gauss-Newton algorithm LSMR tolerance tolB
TOLU = 1e-8     # Gauss-Newton algorithm Minres rtol


class HarmonisationOp:
    """
    This class runs the harmonisation process of Peter Harris for satellite sensor match-up data contained in the
    input directory

    Sample Code:

    .. code-block:: python

        H = HarmOp()
        H.run(...)

    :Methods:
        .. py:method:: run(...):

            Run harmonisation for a given input satellite sensor match-up data set configuration

        . py:method:: calculate_parameter_covariance_ij(...):

            Return element of harmonisation output parameter covariance matrix

        .. py:method:: combine_parameter_covariance_ij(...):

            Combine parameter covariance elements to form full matrix and save to file

    """

    def __init__(self):
        pass

    def run(self, dataset_dir, sensor_data, output_dir,
                  pc_input=None, save_pc=None,
                  gn_input=None, save_gn=None,
                  software_cfg=None,
                  tolPC=TOLPC, tol=TOL, tolA=TOLA, tolB=TOLB, tolU=TOLU,
                  show=1,
                  return_covariance=True):
        """
        Run harmonisation for a given input satellite sensor match-up data set configuration
        """

        # Default to save residual data
        res = True

        ################################################################################################################
        # 1.	Read Harmonisation Matchup Data
        ################################################################################################################

        print "\nOpening data..."
        HData = MatchUp(dataset_dir, sensor_data)
        print "Complete"

        print "\nData Info"
        print "=========="
        print "Reference Sensors - ", [str(sensor) for sensor in HData.idx['sensors']
                                       if sensor not in HData.idx["parameter_sensor"]]
        print "Harmonising Sensors - ", [str(sensor) for sensor in HData.idx['sensors']
                                         if sensor in HData.idx["parameter_sensor"]]
        print "Total Match-Ups - ", HData.idx['cNm'][-1]
        print "Total Sensor State Data Values - ", HData.idx['idx'][-1]
        print "Total Harmonisation Paramaters - ", len(HData.idx['parameter_sensor'])

        ################################################################################################################
        # 2.	Perform harmonisation
        ################################################################################################################

        print "\nBeginning Harmonisation..."
        Harmonisation = HarmonisationEIV()
        HarmonisationOutput = Harmonisation.run(HData,
                                                pc_input=pc_input, save_pc=save_pc,
                                                gn_input=gn_input, save_gn=save_gn,
                                                show=show)

        print "Final Solution:"
        print HarmonisationOutput.parameter
        print HarmonisationOutput.parameter_covariance_matrix

        ################################################################################################################
        # 3.	Write data to file
        ################################################################################################################

        print 'Writing data to file...'
        # Add metadata
        startDate = ""
        endDate = ""
        HarmonisationOutput.additional_attributes = {"software": software_cfg['software'],
                                                     "software_version": software_cfg['version'],
                                                     "software_tag": software_cfg['tag'],
                                                     "job_id": software_cfg["job_id"],
                                                     "matchup_dataset": "_".join((software_cfg['matchup_dataset'],
                                                                                  startDate, endDate))}

        # Write data
        fname_output = "_".join(("harm", software_cfg['software'], software_cfg['version'],
                                 software_cfg['tag'], software_cfg["job_id"],
                                 software_cfg['matchup_dataset'])) + ".nc"
        HarmonisationOutput.save(pjoin(output_dir, fname_output), save_residuals=res)

        print "\nOutput Data Written To:"
        print ">", output_dir

    def runPC(self, show=False):
        """
        This function runs the preconditioner only of the harmonisation of satellite instrument calibration parameters
        for group of sensors with a reference sensor from the match-up data located in the input directory.

        :globals:
            :self.dataDir: *str*

            directory of match-up data

        :return:
            :a: *numpy.ndarray*

            harmonised sensor calibration parameters

            :Va: *numpy.ndarray*

            covariance matrices for harmonised sensor calibration parameters
        """

        # Initialise

        # 1. Directories
        dataset_paths = self.dataset_paths
        sensor_data_path = self.sensor_data_path
        output_dir = self.output_dir

        # 2. Software Info
        software = self.software
        software_version = self.software_version
        software_tag = self.software_tag
        job_id = self.job_id
        matchup_dataset = self.matchup_dataset

        # Default to save residual data
        res = True

        ################################################################################################################
        # 1.	Read Harmonisation Matchup Data
        ################################################################################################################

        print "Match-up Dataset:"
        for path in dataset_paths:
            print ">", path

        print "\nSensor Data File:"
        print "> "+sensor_data_path

        print "\nOpening data..."
        HData = self.data_reader(dataset_paths, sensor_data_path, open_uncertainty=True)

        ################################################################################################################
        # 2.	Perform harmonisation
        ################################################################################################################

        print "Complete"

        print "\nData Info"
        print "=========="
        print "Reference Sensors - ", [str(sensor) for sensor in HData.idx['sensors']
                                       if sensor not in HData.idx["parameter_sensor"]]
        print "Harmonising Sensors - ", [str(sensor) for sensor in HData.idx['sensors']
                                         if sensor in HData.idx["parameter_sensor"]]
        print "Total Match-Ups - ", HData.idx['cNm'][-1]
        print "Total Sensor State Data Values - ", HData.idx['idx'][-1]
        print "Total Harmonisation Paramaters - ", len(HData.idx['parameter_sensor'])

        print "\nComputing Preconditioner Solution..."
        # a. sample data for preconditioning
        Sample2IndOp = Sample2Ind()
        HData_sample = Sample2IndOp.run(HData, sf=1, show=show)

        Transform2NormIndOp = Transform2NormInd()
        HData_sample = Transform2NormIndOp.run(HData_sample)

        # b. determine preconditioner solution
        PC = PCAlgo(HData_sample)
        a_PC, S_PC = PC.runPC(tol=1e-6)
        del HData_sample

        ################################################################################################################
        # 3.	Write data to file
        ################################################################################################################

        print '\nWriting data to file...'

        # Write data
        save_directory_PC = pjoin(output_dir, "temp", "PC")
        try:
            os.makedirs(save_directory_PC)
        except OSError:
            pass

        S_PC_path = pjoin(save_directory_PC, "S_PC.txt")
        a_PC_path = pjoin(save_directory_PC, "a_PC.txt")

        savetxt(S_PC_path, S_PC)
        savetxt(a_PC_path, a_PC)

        print "Output Data Written To:"
        print ">", save_directory_PC

    def calculate_GNOp_initial(self):
        """
        Writes the initial value of GNOp by transforming the input data to independent

        :globals:
            :self.dataDir: *str*

            directory of match-up data

        :return:
            :a: *numpy.ndarray*

            harmonised sensor calibration parameters

            :Va: *numpy.ndarray*

            covariance matrices for harmonised sensor calibration parameters
        """

        # Initialise
        Transform2NormIndOp = Transform2NormInd()

        # Directories
        dataset_paths = self.dataset_paths
        sensor_data_path = self.sensor_data_path
        output_dir = self.output_dir

        ################################################################################################################
        # 1.	Read Harmonisation Matchup Data
        ################################################################################################################

        print "Match-up Dataset:"
        for path in dataset_paths:
            print ">", path

        print "\nSensor Data File:"
        print "> "+sensor_data_path

        print "\nOpening data..."
        HData = self.data_reader(dataset_paths, sensor_data_path, open_uncertainty=True)

        ################################################################################################################
        # 2.	Perform harmonisation
        ################################################################################################################

        print "Complete"

        print "\nData Info"
        print "=========="
        print "Reference Sensors - ", [str(sensor) for sensor in HData.idx['sensors']
                                       if sensor not in HData.idx["parameter_sensor"]]
        print "Harmonising Sensors - ", [str(sensor) for sensor in HData.idx['sensors']
                                         if sensor in HData.idx["parameter_sensor"]]
        print "Total Match-Ups - ", HData.idx['cNm'][-1]
        print "Total Sensor State Data Values - ", HData.idx['idx'][-1]
        print "Total Harmonisation Paramaters - ", len(HData.idx['parameter_sensor'])

        print "\nComputing Initial GN Value..."
        # a. reparameterise input data such that output data are independent quantities
        HData = Transform2NormIndOp.run(HData)

        # b. run GN algorithm on modified data
        GNOp = GNAlgo(HData)

        ################################################################################################################
        # 3.	Write data to file
        ################################################################################################################

        print '\nWriting data to file...'

        # Write data
        save_directory_GNOp = pjoin(output_dir, "temp", "GNOp_initial")
        try:
            os.makedirs(save_directory_GNOp)
        except OSError:
            pass

        GNOp.save(save_directory_GNOp)

        print "Output Data Written To:"
        print ">", save_directory_GNOp

    def calculate_parameter_covariance_ij(self, i, j):
        """
        Return one element of harmonisation output parameter covariance matrix

        :type i: int
        :param i: covariance matrix element row

        :type j: int
        :param j: covariance matrix element column
        """

        # Initialise

        # 1. Input data
        dataset_paths = self.dataset_paths
        sensor_data_path = self.sensor_data_path

        # 2. Output Data
        temp_directory = self.temp_directory

        ################################################################################################################
        # 1.	Read Matchup Data and harmonisation output data
        ################################################################################################################

        # Input data
        HData = self.data_reader(dataset_paths, sensor_data_path,
                                 open_uncertainty=True, open_additional_values=False)

        # Re-open final solver
        GNOp = GNAlgo(HData)
        GNOp.open(temp_directory)

        ################################################################################################################
        # 2.	Perform harmonisation parameter covariance matrix element
        ################################################################################################################

        parameter_covariance_ij = GNOp.calculate_parameter_covariance_matrix_ij([i, j])

        return parameter_covariance_ij

    def combine_parameter_covariance_ij(self):
        """
        Combine parameter covariance elements to form full matrix and save to file
        """

        filelist = [file for file in glob(pjoin(self.temp_directory, "parameter_covariance_matrix*.dat"))]
        N_parameters = int(len(filelist)**0.5)

        parameter_covariance_matrix = zeros((N_parameters, N_parameters))

        for file in filelist:
            parameter_covariance_matrix[int(file[-8]), int(file[-6])] = loadtxt(file)

        parameter_covariance_matrix = (parameter_covariance_matrix + parameter_covariance_matrix.T) / 2

        return parameter_covariance_matrix

    def run_mc_trial(self, trial_number):
        """
        Run one Monte Carlo trial of harmonisation process

        :type trial_number: int
        :param trial_number: covariance matrix element row
        """

        # Initialise

        # 1. Input data
        dataset_paths = self.dataset_paths
        sensor_data_path = self.sensor_data_path

        # 2. Output Data
        temp_directory = self.temp_directory
        mc_directory = self.temp_directory

        # 3. Software Info
        software = self.software
        software_version = self.software_version
        software_tag = self.software_tag
        job_id = self.job_id
        matchup_dataset = self.matchup_dataset

        res = True

        ################################################################################################################
        # 1.	Read Matchup Data and harmonisation output data
        ################################################################################################################

        print "Opening Match-up data:"
        for path in dataset_paths:
            print ">", path

        # Input data
        HData = self.data_reader(dataset_paths, sensor_data_path,
                                 open_uncertainty=True, open_additional_values=False)

        print "Re-opening GN solver..."

        # Re-open final solver
        GNOp = GNAlgo(HData)
        GNOp.open(temp_directory)

        ################################################################################################################
        # 2.	Perform harmonisation parameter covariance matrix element
        ################################################################################################################

        # Set up MC trial
        N_var = HData.idx['idx'][-1]
        GNOp.HData.values = normal(GNOp.xyza[N_var], scale=1.0)
        GNOp.HData.ks = normal(GNOp.HData.ks, scale=1.0)
        GNOp.xyza = append(GNOp.HData.values, GNOp.HData.a)

        # Run trial
        HarmonisationOutput = HarmonisationResult()

        HarmonisationOutput.parameter_sensors = HData.idx["Ia"]
        HarmonisationOutput.idx = deepcopy(HData._original_idx)

        HarmonisationOutput.parameter, HarmonisationOutput.parameter_covariance_matrix, \
        HarmonisationOutput.cost, HarmonisationOutput.cost_dof, \
        HarmonisationOutput.cost_p_value, HarmonisationOutput.values_res, \
        HarmonisationOutput.ks_res \
            = GNOp.run(show=True, return_covariance=False)

        # Add metadata
        trial_number_str = str(trial_number).zfill(3)
        startDate = ""
        endDate = ""
        matchup_dataset[-3:] = trial_number_str
        HarmonisationOutput.additional_attributes = {"software": software,
                                                     "software_version": software_version,
                                                     "software_tag": software_tag,
                                                     "job_id": job_id,
                                                     "matchup_dataset": "_".join((matchup_dataset, startDate, endDate))}

        # Write data
        trial_directory = pjoin(mc_directory, trial_number_str)
        try:
            os.makedirs(trial_directory)
        except OSError:
            pass

        makedirs(trial_directory)
        fname_output = "_".join(("harm", software, software_version, software_tag, job_id, matchup_dataset)) + ".nc"
        HarmonisationOutput.save(pjoin(trial_directory, fname_output), save_residuals=res)

        print "Done"

        return 0

    def combine_mc_trials(self):
        """
        Run routine to combine and save MC trial output data
        """

        # Initialise
        output_dir = self.output_dir
        harm_out_path_best_estimate, harm_res_paths_best_estimate = get_harm_paths(output_dir)  # Get path
        mc_directory = self.mc_directory
        output_combine_dir = pjoin(mc_directory, 'combine')

        try:
            os.makedirs(output_combine_dir)
        except OSError:
            pass

        HarmComb = HarmonisationResult()

        ################################################################################################################
        # 1. Compute MC statistics
        ################################################################################################################

        # a. Get individual mc trial output directories
        mc_dirs = [mc_dir for mc_dir in glob(pjoin(mc_directory, "*/")) if split(split(mc_dir)[0])[1] != "combine"]

        # Compute parameter statistics
        a = self.get_parameters_mc(mc_dirs)  # Get parameters from individual mc trials
        parameter_mean = mean(a, axis=0)
        parameter_std = std(a, axis=0)
        parameter_covariance_matrix = cov(a, rowvar=0)

        # b. Compute normalised residuals
        # i. compute residual variance

        # Check if monte carlo trials have residual data
        mc_res = False
        if get_harm_paths(mc_dirs)[1] is not None:
            mc_res = True

        # Incrementally calculate variance
        H_res_00 = []
        if mc_res:
            for i, mc_dir in enumerate(mc_dirs):
                harm_res_paths_trial = get_harm_paths(mc_dir)[1]
                n_trial = i+1

                if n_trial == 1:
                    _, k_res_e, H_res_e, _, _, _, _, _ = HarmComb.open_harmonisation_residual_files(harm_res_paths_trial)
                    H_res_00.append(H_res_e[0,0])
                    k_res_variance = zeros(k_res_e.shape)
                    H_res_variance = zeros(H_res_e.shape)
                else:
                    _, k_res_i, H_res_i, _, _, _, _, _ = HarmComb.open_harmonisation_residual_files(harm_res_paths_trial)
                    k_res_mu = k_res_e + (k_res_i-k_res_e)/n_trial
                    H_res_mu = H_res_e + (H_res_i-H_res_e)/n_trial
                    H_res_00.append(H_res_i[0, 0])
                    k_res_v = (n_trial-2)*k_res_variance + (n_trial-1)*(k_res_e-k_res_mu)**2 + (k_res_i-k_res_mu)**2
                    H_res_v = (n_trial-2)*H_res_variance + (n_trial-1)*(H_res_e-H_res_mu)**2 + (H_res_i-H_res_mu)**2

                    k_res_e = k_res_mu
                    H_res_e = H_res_mu

                    k_res_variance = k_res_v/(n_trial-1)
                    H_res_variance = H_res_v/(n_trial-1)

            # ii. normalise best estimate residuals
            lm, k_res, H_res, _, _, _, _, _ = self.open_harmonisation_residual_files(harm_res_paths_best_estimate)
            k_res_normalised = k_res / k_res_variance**0.5
            H_res_normalised = H_res / H_res_variance**0.5

        ################################################################################################################
        # 2. Get initial data from first trial output
        ################################################################################################################

        # Open harmonisation output file
        self.parameter, self.parameter_sensors, _, self.cost, self.cost_dof, self.cost_p_value, _, \
        self.software_version, self.software_tag, self.job_id, self.matchup_dataset\
            = HarmComb.open_harmonisation_output_file(harm_out_path_best_estimate)

        # Update variables to MC specific values
        self.parameter_covariance_matrix = parameter_covariance_matrix
        self.software = "EM"
        self.lm = lm.astype(int)
        self.additional_attributes = {"parameter_trial_mean": parameter_mean,
                                      "parameter_trial_std": parameter_std}

        if mc_res:
            self.k_res = k_res_normalised
            self.H_res = H_res_normalised

        ################################################################################################################
        # 3. Write data
        ################################################################################################################

        self.save(output_combine_dir, res=True)

    def get_parameters_mc(self, mc_dirs):
        """
        Return array of parameter outputs for each MC trial

        :param mc_dirs: list: str
            List of directories of individual MC trial outputs

        :return:
            :a: numpy.ndarray
                Parameters from each MC trial
        """

        # Get number of parameters
        harm_out_path, harm_res_path = get_harm_paths(mc_dirs[0])
        n_a = self.open_harmonisation_output_file(harm_out_path)[0].shape[0]

        # Initialise array
        a = zeros((len(mc_dirs), n_a))

        for i, mc_dir in enumerate(mc_dirs):
            harm_out_path, harm_res_path = get_harm_paths(mc_dir)

            # Ignore any failed jobs
            if harm_out_path is not None:
                a[i, :] = self.open_harmonisation_output_file(harm_out_path)[0]

        return a[~nall(a == 0, axis=1)]

if __name__ == "__main__":
    pass

