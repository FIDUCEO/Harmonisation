"""
Created on Thurs April 13 2017 16:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Harmonisation Modules___'''
from harm_data_writer import HarmOutput
from plotting_functions import *
from radiance_functions import *
from config_functions import *

'''___Python Modules___'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import append, arange, ones, where, savetxt, loadtxt
from sys import argv
from os.path import join as pjoin
from os.path import abspath, split, join
from os import mkdir
from datetime import datetime
from numpy.random import rand

'''__Constants___'''

# Number of samples to include in all match-up scatter plots (-1 => plot all match-up data
N_SAMPLE = 10000

# Number of sample to include in simulated data plots
N_SIM = 20


class HarmDiag:

    def __init__(self, dataset_paths=None, parameter_path=None, harm_output_path=None, harm_res_paths=None,
                 output_dir=None, sensor_model=None, adjustment_model=None, software_cfg=None, data_reader=None):
        """
        Initialise harmonisation algorithm

        :param dataset_paths: list: str
            Paths of matchup series files in matchup dataset directory

        :param parameter_path: str
            Path of parameters to be used as starting point for solver

        :param harm_output_path: str
            Path of harmonisation output file

        :param harm_res_paths: str
            Path of harmonisation residual file

        :param output_dir: str
            Path of directory to store output data files in

        :param sensor_model: func
            Python function to calculate radiance and derivatives given input sensor state data

        :param adjustment_model: func
            Python function to calculate spectral adjustment factor between two sensors, k, and derivatives given
            the two sensor radiances

        :param software_cfg: dict:str
            dictionary of software configuration information

        :param data_reader: Class
            Python class to open harmonisation data. If none given default data reader used.
        """

        # Set required paths
        self.dataset_paths = dataset_paths
        self.parameter_path = parameter_path
        self.harm_output_path = harm_output_path
        self.harm_res_paths = harm_res_paths
        self.outDir = output_dir

        # Software info (currently unused)
        self.software_cfg = software_cfg

        # Functions
        self.sensor_model = sensor_model
        self.adjustment_model = adjustment_model

        # If none default data reader required add to class, else import default
        if data_reader is not None:
            self.HarmData = data_reader
        else:
            from harm_data_reader import HarmData
            self.HarmData = HarmData

    def run(self):
        """
        Generate the diagnostic plots
        """

        # Initialise

        # 1. Directories
        dataset_paths = self.dataset_paths
        parameter_path = self.parameter_path
        harm_output_path = self.harm_output_path
        harm_res_paths = self.harm_res_paths

        # 2. Functions
        sensor_model = self.sensor_model
        adjustment_model = self.adjustment_model

        ################################################################################################################
        # 1. Open Data
        ################################################################################################################

        print("Opening Data...")

        # Open output harmonisation data
        HOut = HarmOutput(harm_output_path, harm_res_paths)
        HOut.software_fullname = "_".join((HOut.software, HOut.software_version, HOut.software_tag, HOut.job_id))

        ################################################################################################################
        # 2. Calculate residuals
        ################################################################################################################

        # Open data match-up data match-up by match-up
        print("Calculating Residuals...")

        nodata = True

        r_sim = {}
        r_res_sim = {}
        r_res_sim_err = {}

        # Have to open all input parameters separately, as opening HData one file at a time

        # Make dummy a for individual file opening
        a_temp = ones(3)
        temp_parameter_path = join(split(parameter_path)[0], "a_temp.csv")
        savetxt(temp_parameter_path, a_temp, delimiter=",")

        # Open all a
        a = loadtxt(parameter_path, delimiter=',')
        n_s = a.shape[0]
        n_a = a.shape[1]
        a = a.flatten()
        Ia = zeros(a.shape[0])
        i_next = 0

        # Last times per sensor
        last_times_s = asarray([datetime.utcfromtimestamp(0) for i in range(n_s)], dtype=datetime)
        sensor_IDs = zeros(n_s)
        i_nextID = 0


        for dataset_path in dataset_paths:

            with self.HarmData([dataset_path], temp_parameter_path, sensor_model, adjustment_model) as HData:

                # Build Ia match-up by match-uo
                for Im in HData.idx['lm']:
                    for sensor in Im[:2]:
                        if (sensor not in Ia) and (sensor != -1):
                            Ia[i_next:i_next + n_a] = sensor
                            i_next += n_a
                # 1. All data

                # Calculate radiance from data with input coefficients
                temp_r = calc_R(HData, a=a, Ia=Ia, V=None)

                # Calculate residual between this and radiance from data with
                r_est = calc_R(HData, a=HOut.parameter, Ia=HOut.parameter_sensors, V=None)
                temp_r_res = calc_R_res(temp_r, r_est, uR=None)

                if nodata:
                    times = HData.times
                    last_times = array([HData.times[-1]])
                    lm = HData.idx['lm']
                    r = temp_r[:, :]
                    r_res = temp_r_res[:, :]
                    nodata = False

                else:
                    times = append(times, HData.times)
                    last_times = append(last_times, HData.times[-1])
                    lm = append(lm, HData.idx['lm'], axis=0)
                    r = append(r, temp_r, axis=0)
                    r_res = append(r_res, temp_r_res, axis=0)

                # Get last times by sensor
                for sensor_ID in HData.idx['lm'][0][:2]:
                    if sensor_ID != -1:
                        if sensor_ID not in sensor_IDs:
                            sensor_IDs[i_nextID] = sensor_ID
                            i_nextID += 1
                        i_sensor = where(sensor_IDs == sensor_ID)[0][0]
                        if last_times_s[i_sensor] < HData.times[-1]:
                            last_times_s[i_sensor] = HData.times[-1]

                # 2. Simulated data
                # Calculate radiance residual between input parameters and harmonised parameters for simulated
                # sampling of the data space
                values_sim = sim_sensor_values(HData, n_sim=N_SIM)

                for sensor in values_sim.keys():
                    i_s_in = asarray([True if s == sensor else False for s in Ia])
                    i_s_out = asarray([True if s == sensor else False for s in HOut.parameter_sensors])
                    if sensor not in r_sim.keys():
                        r_sim[sensor] = calc_R_sensor(values_sim[sensor], sensor_model=HData.sensor_model, a=a[i_s_in])
                        r_sim_est, uR_sim = calc_R_sensor(values_sim[sensor], sensor_model=HData.sensor_model,
                                                          a=HOut.parameter[i_s_out],
                                                          V=HOut.parameter_covariance_matrix[outer(i_s_out, i_s_out)].
                                                          reshape((sum(i_s_out), sum(i_s_out))))

                        r_res_sim[sensor], r_res_sim_err[sensor] = calc_R_res(r_sim[sensor], r_sim_est, uR_sim)

        # 3. Calculate per sensor and per match-up series statistics

        # a. radiance residual, r_res, statistics
        r_res_ave_s, r_res_sd_s = average_sensor(r_res, lm)   # per sensor
        r_res_ave_m, r_res_sd_m = average_mu(r_res, lm)       # per match-up series

        # b. match-up adjustment factor residual, k_res, statistics
        k_res_ave_m, k_res_sd_m = average_mu(HOut.k_res, lm)  # per match-up series

        # Select sample indices
        if (N_SAMPLE == -1) or (N_SAMPLE > r_res.shape[0]):
            i_sels = arange(sum(lm[:, 2]))
        else:
            i_sels = sort((sum(lm[:, 2]) * rand(N_SAMPLE)).astype(int))

        ################################################################################################################
        # 3. Make plots
        ################################################################################################################

        print("Making plots...")

        directories = [self.outDir,
                       pjoin(self.outDir, "cov"),
                       pjoin(self.outDir, "res"),
                       pjoin(self.outDir, "series"),
                       pjoin(self.outDir, "sensor")]

        for directory in directories:
            try:
                mkdir(directory)
            except OSError:
                pass

        # shared text box contents
        txt = ["Software Name: " + HOut.software_fullname, "Dataset: " + HOut.matchup_dataset]

        # shared axis labels
        xlbl_r = "Data Radiance, $L_{DATA}$"
        xlbl_t = "Match-up time"
        xlbl_m = "Match-up series"
        xlbl_s = "Sensor"
        ylbl_r_res = "Radiance Residual, $\Delta L = L_{EST} - L_{DATA}$"
        ylbl_k_res = "K Residual, $\Delta K = K_{EST} - K_{DATA}$"

        title_sample = " (Sample of "+str(N_SAMPLE)+" match-ups)"

        # shared axis tick labels
        albls = self.get_albls(HOut.parameter_sensors)
        slbls = self.get_slbls(HOut.parameter_sensors)
        mlbls = self.get_mlbls(lm)

        # time sort
        mlbls = [mlbl for last_time, mlbl in sorted(zip(last_times, mlbls), reverse=True)]
        r_res_ave_m = asarray([r_ave for last_time, r_ave in sorted(zip(last_times, r_res_ave_m), reverse=True)])
        r_res_sd_m = asarray([r_sd for last_time, r_sd in sorted(zip(last_times, r_res_sd_m), reverse=True)])
        k_res_ave_m = asarray([k_ave for last_time, k_ave in sorted(zip(last_times, k_res_ave_m), reverse=True)])
        k_res_sd_m = asarray([k_sd for last_time, k_sd in sorted(zip(last_times, k_res_sd_m), reverse=True)])

        # time sort
        slbls = [slbl for last_time, slbl in sorted(zip(last_times_s, slbls), reverse=True)]
        r_res_ave_s = asarray([r_ave for last_time, r_ave in sorted(zip(last_times_s, r_res_ave_s), reverse=True)])
        r_res_sd_s = asarray([r_sd for last_time, r_sd in sorted(zip(last_times_s, r_res_sd_s), reverse=True)])

        # 1. Covariance Matrix Plot
        plot_grid_heatmap(pjoin(self.outDir, "cov", "cov_heatmap.pdf"), HOut.parameter_covariance_matrix,
                          title="Parameter Covariance Matrix",
                          labels=albls, dividers=slbls)

        # 2. Plot variable residuals
        # a. radiance residual
        # i. space sampled data per sensor
        for sensor in r_res_sim.keys():
            plot_scatter(pjoin(self.outDir, "res", "r_res_vs_r_errbars_"+str(sensor)+".pdf"),
                         r_res_sim[sensor], r_sim[sensor], yerr=r_res_sim_err[sensor],
                         title="Sampled Radiance Residual - "+str(sensor), ylbl=ylbl_r_res, xlbl=xlbl_r, txt=txt, dash_ylines=[0])

        # ii. all data

        title_r_res = "Radiance Residual, $\Delta R$"
        plot_2dhist(pjoin(self.outDir, "res", "r_res_vs_r_2dhist.pdf"), r_res.flatten('F'), r.flatten('F'), bins=50,
                    title=title_r_res, xlbl=xlbl_r, ylbl=ylbl_r_res, txt=txt, solid_ylines=[0])
        plot_2dhist(pjoin(self.outDir, "res", "r_res_vs_time_2dhist.pdf"), r_res.flatten('F'), append(times, times), bins=50,
                    title=title_r_res, xlbl=xlbl_t, ylbl=ylbl_r_res, txt=txt, solid_ylines=[0])
        plot_scatter(pjoin(self.outDir, "res", "r_res_vs_r_scatter.png"), r_res[i_sels, :].flatten('F'), r[i_sels, :].flatten('F'),
                     title=title_r_res+title_sample, ylbl=ylbl_r_res, xlbl=xlbl_r, txt=txt, solid_ylines=[0])
        plot_scatter(pjoin(self.outDir, "res", "r_res_vs_time_scatter.png"), r_res[i_sels, :].flatten('F'), append(times[i_sels], times[i_sels]),
                     title=title_r_res+title_sample, ylbl=ylbl_r_res, xlbl=xlbl_t, txt=txt, solid_ylines=[0])

        # b. match-up adjustment factor residual

        title_k_res = "Match-up Adjustment Factor Residual, $\Delta K$"
        plot_2dhist(pjoin(self.outDir, "res", "k_res_vs_r_2dhist.pdf"), HOut.k_res, r[:, 1], bins=50,
                    title=title_k_res, xlbl=xlbl_r, ylbl=ylbl_k_res, txt=txt, solid_ylines=[0])
        plot_2dhist(pjoin(self.outDir, "res", "k_res_vs_time_2dhist.pdf"), HOut.k_res, times, bins=50,
                    title=title_k_res, xlbl=xlbl_r, ylbl=ylbl_k_res, txt=txt, solid_ylines=[0])
        plot_scatter(pjoin(self.outDir, "res", "k_res_vs_r_scatter.png"), HOut.k_res[i_sels], r[i_sels, 1],
                     title=title_k_res+title_sample, ylbl=ylbl_k_res, xlbl=xlbl_r, txt=txt, solid_ylines=[0])
        plot_scatter(pjoin(self.outDir, "res", "k_res_vs_time_scatter.png"), HOut.k_res[i_sels], times[i_sels],
                     title=title_k_res+title_sample, ylbl=ylbl_k_res, xlbl=xlbl_r, txt=txt, solid_ylines=[0])

        # c. variable residual (EIV only)
        if HOut.H_res is not None:
            m = HOut.H_res.shape[1]/2
            for i in range(m):
                title_xi_res = "$X_{"+str(i+1)+"}$ Residual, $\Delta X_{"+str(i+1)+"}$"
                ylbl_xi_res = "$X_{"+str(i+1)+"}$ Residual, $\Delta X_{"+str(i+1)+"} = X_{"+str(i+1)+", EST} - X_{"+str(i+1)+", DATA}$"
                plot_2dhist(pjoin(self.outDir, "res", "x"+str(i+1)+"_res_vs_r_2dhist.pdf"),
                            append(HOut.H_res[:, i], HOut.H_res[:, i+m]), r.flatten('F'), bins=50,
                            title=title_xi_res, xlbl=xlbl_r, ylbl=ylbl_xi_res, txt=txt, solid_ylines=[0])
                plot_2dhist(pjoin(self.outDir, "res", "x"+str(i+1)+"_res_vs_time_2dhist.pdf"),
                            append(HOut.H_res[:, i], HOut.H_res[:, i+m]), append(times, times), bins=50,
                            title=title_xi_res, xlbl=xlbl_t, ylbl=ylbl_xi_res, txt=txt, solid_ylines=[0])
                plot_scatter(pjoin(self.outDir, "res", "x"+str(i+1)+ "_res_vs_r_scatter.png"),
                             append(HOut.H_res[i_sels, i], HOut.H_res[i_sels, i+m]), r[i_sels, :].flatten('F'),
                             title=title_xi_res+title_sample, ylbl=ylbl_xi_res, xlbl=xlbl_r, txt=txt, solid_ylines=[0])
                plot_scatter(pjoin(self.outDir, "res", "x"+str(i+1)+"_res_vs_time_scatter.png"),
                             append(HOut.H_res[i_sels, i], HOut.H_res[i_sels, i+m]), append(times[i_sels], times[i_sels]),
                             title=title_xi_res+title_sample, ylbl=ylbl_xi_res, xlbl=xlbl_t, txt=txt, solid_ylines=[0])

        # 3. Plot radiance residual distribution
        # Todo - Properly implement radiance residual distribution graph
        #self.plot_hist(pjoin(self.outDir, "res_dist", "r_res_scaled_hist.pdf"))

        # 4. Plot match-up series
        title_r_res_m = r"Average Radiance Residual, $\overline{\Delta R}$, per Match-up Series"
        plot_scatter(pjoin(self.outDir, "series", "r_res_ave_vs_series_errbars.pdf"), r_res_ave_m, mlbls, yerr=r_res_sd_m,
                     title=title_r_res_m, ylbl=ylbl_r_res, xlbl=xlbl_m, txt=txt, solid_ylines=[0])

        title_k_res_m = r"Average Match-Up Adjustment Factor Residual, $\overline{\Delta K}$, per Match-up Series"
        plot_scatter(pjoin(self.outDir, "series", "k_res_ave_vs_series_errbars.pdf"), k_res_ave_m, mlbls, yerr=k_res_sd_m,
                      title=title_k_res_m, ylbl=ylbl_k_res, xlbl=xlbl_m, txt=txt, solid_ylines=[0])

        # 5. Plot sensors
        title_r_res_s = r"Average Radiance Residual, $\overline{\Delta R}$, per Sensor"
        plot_scatter(pjoin(self.outDir, "sensor", "r_res_ave_vs_sensor_errbars.pdf"), r_res_ave_s, slbls, yerr=r_res_sd_s,
                     title=title_r_res_s, ylbl=ylbl_r_res, xlbl=xlbl_s, txt=txt, solid_ylines=[0])

        return 0

    def get_albls(self, parameter_sensor_names):
        """
        Return labels for each harmonisation parameter from input Ia array from HarmOutput object

        :param parameter_sensor_names: numpy.ndarray
                harmonised parameter sensor names

        :return:
            :a_lbls: list:str
                harmonised parameter name of form "$a_{N, SS}$"
                where:
                - N is the parameter number
                - SS is the sensor number
        """

        a_nums = []
        sensor_prev = None
        for sensor in parameter_sensor_names:
            if sensor != sensor_prev:
                a_num = 0
            else:
                a_num += 1

            a_nums.append(a_num)
            sensor_prev = sensor

        a_lbls = ["$a_{" + str(n) + "}$" for n, s in zip(a_nums, parameter_sensor_names)]

        return a_lbls

    def get_slbls(self, parameter_sensor_names):
        """
        Return axis tick labels for each sensor

        :param parameter_sensor_names: numpy.ndarray
                harmonised parameter sensor names

        :return:
            :s_lbls: list:str
                match-up series name of form "SS"
        """

        # find unique sensors
        s_lbls = []
        sensors = ["a"]
        for i, sensor in enumerate(parameter_sensor_names):
            if sensor != sensors[-1]:

                sensors.append(sensor)

                if type(sensor) == float:
                    sensor = int(sensor)

                s_lbls.append(str(sensor))



        return s_lbls

    def get_mlbls(self, lm):
        """
        Return axis tick labels for each match-up series

        :param lm:
                match-up series description data

        :return:
            :m_lbls: list:str
                match-up series name of form "Si_Sj",
                where:
                > Si in the name of sensor i in the match-up series
                > Sj in the name of sensor j in the match-up series
        """

        # initialise list
        m_lbls = []

        # step through data array match-up series by match-up series generating each respective match-up name
        for row in lm[:, :2]:
            Si = str(row[0])
            Sj = str(row[1])
            m_lbls.append("_".join((Si, Sj)))

        return m_lbls


if __name__ == "__main__":

    def main():

        ################################################################################################################
        # Process configuration data
        ################################################################################################################

        # 1. Get configuration filename
        if len(argv) == 1:
            # EXAMPLE CONFIGURATION FILE
            job_cfg_fname = "Data/AVHRR_REAL_3_sample3/AVHRR_REAL_3_sample3.cfg"

        else:
            # else have usage:
            # argv[1] - path of job config file

            job_cfg_fname = os.path.abspath(argv[1])

        # 2. Read configuration data
        conf = {}   # dictionary to store data

        #  a. Read software config file
        software_cfg_fname = "software.cfg"
        conf['software'], conf['version'], conf['tag'] = read_software_cfg(software_cfg_fname)

        # b. Read job config file
        conf['job_id'], conf['matchup_dataset'], dataset_dir, parameter_path, output_dir, \
            sensor_functions_path, data_reader_path = read_job_cfg(job_cfg_fname)

        # 3. Get matchup data paths from directory
        dataset_paths = get_dataset_paths(dataset_dir)
        harm_output_path, harm_res_paths = get_harm_paths(output_dir)

        # 4. Import required specified functions
        sensor_functions = import_file(sensor_functions_path)
        harm_data_reader = import_file(data_reader_path)

        ################################################################################################################
        # Run diagnostics
        ################################################################################################################

        HDiag = HarmDiag(dataset_paths=dataset_paths,
                         parameter_path=parameter_path,
                         harm_output_path=harm_output_path,
                         harm_res_paths=harm_res_paths,
                         output_dir=pjoin(output_dir, "plots"),
                         sensor_model=sensor_functions.sensor_model,
                         adjustment_model=sensor_functions.adjustment_model,
                         software_cfg=conf,
                         data_reader=harm_data_reader.HarmData)
        HDiag.run()

        return 0

    main()
