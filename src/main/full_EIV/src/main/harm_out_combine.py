"""
Combine job MC outputs

Created on Thurs 18 May  2017 17:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Harmonisation Modules___'''
from config_functions import *
from harm_data_writer import HarmOutput

'''___Python Modules___'''
from sys import argv
from os.path import split
from os.path import join as pjoin
import os
import numpy as np
from numpy import savetxt, append, zeros, cov, mean, std
from glob import glob


class HarmOutputCombine(HarmOutput):
    """
    Class to statistically combine data from subdirectories of MC trials on harmonisation runs and write to
    harmonisation output files

    Sample Code:
    H = HarmOutputCombine("/path/to/output/mc/directories")
    H.run()

    """

    def __init__(self, output_dir):
        """
        Initilise class

        :param output_mc_dir: str
            path of directory containing individual MC trial output directories
        """

        self.output_dir = output_dir

    def run(self):
        """
        Run routine to combine and save MC trial output data
        """

        # Initialise
        output_dir = self.output_dir
        harm_out_path_best_estimate, harm_res_paths_best_estimate = get_harm_paths(output_dir)  # Get path
        output_mc_dir = pjoin(output_dir, 'mc')
        output_combine_dir = pjoin(output_mc_dir, 'combine')

        try:
            os.makedirs(output_combine_dir)
        except OSError:
            pass

        ################################################################################################################
        # 1. Compute MC statistics
        ################################################################################################################

        # a. Get individual mc trial output directories
        mc_dirs = [mc_dir for mc_dir in glob(pjoin(output_mc_dir, "*/")) if split(split(mc_dir)[0])[1] != "combine"]

        # Compute parameter statistics
        a = self.get_parameters_mc(mc_dirs)  # Get parameters from individual mc trials
        parameter_mean = mean(a, axis=0)
        parameter_std = std(a, axis=0)
        parameter_covariance_matrix = cov(a, rowvar=0)

        # b. Compute normalised residuals
        # i. compute residual variance

        # Check if monte carlo trials have residual data
        mc_res = False
        if get_harm_paths(mc_dir)[1] is not None:
            mc_res = True

        # Incrementally calculate variance
        H_res_00 = []
        if mc_res:
            for i, mc_dir in enumerate(mc_dirs):
                harm_res_paths_trial = get_harm_paths(mc_dir)[1]
                n_trial = i+1

                if n_trial == 1:
                    _, k_res_e, H_res_e, _, _, _, _, _ = self.open_harmonisation_residual_files(harm_res_paths_trial)
                    H_res_00.append(H_res_e[0,0])
                    k_res_variance = zeros(k_res_e.shape)
                    H_res_variance = zeros(H_res_e.shape)
                else:
                    _, k_res_i, H_res_i, _, _, _, _, _ = self.open_harmonisation_residual_files(harm_res_paths_trial)
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
            = self.open_harmonisation_output_file(harm_out_path_best_estimate)

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

        return a[~np.all(a == 0, axis=1)]

if __name__ == "__main__":

    def main():

        ################################################################################################################
        # Process configuration data
        ################################################################################################################

        # 1. Get configuration filename
        if len(argv) == 2:
            # else have usage:
            # argv[1] - path of job config file
            job_cfg_fname = os.path.abspath(argv[1])
        else:
            # EXAMPLE CONFIGURATION FILE
            job_cfg_fname = "Data/AVHRR_RSIM_3_3s/AVHRR_RSIM_3_sample3_sample.cfg"

        # 2. Read configuration data
        conf = {}   # dictionary to store data

        #  a. Read software config file
        software_cfg_fname = "software.cfg"
        conf['software'], conf['version'], conf['tag'] = read_software_cfg(software_cfg_fname)

        # b. Read job config file
        conf['job_id'], conf['matchup_dataset'], dataset_dir, parameter_path, output_dir, \
            sensor_functions_path, data_reader_path, _ = read_job_cfg(job_cfg_fname)

        ################################################################################################################
        # Get parameters from MC files
        ################################################################################################################

        HOuts = HarmOutputCombine(output_dir)
        HOuts.run()

        return 0

    main()
