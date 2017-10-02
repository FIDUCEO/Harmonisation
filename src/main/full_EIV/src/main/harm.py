"""
Main Operator of NPL Harmonsation Implementation
"""

'''___Python Modules___'''
import os.path
from os import makedirs
from sys import argv

'''___Third Party Modules___'''
from numpy import array_equal

'''___Harmonisation Modules___'''
from config_functions import *
from harm_data_writer import HarmOutput
from harm_data_errors import gen_errors
from harm_algo_EIV import HarmAlgo

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


class HarmOp:
    """
    This class runs the harmonisation process of Peter Harris for satellite sensor match-up data contained in the
    input directory

    Sample Code:

    .. code-block:: python

        H = HarmOp(directory)
        H.run()

    :Attributes:
        .. py:attribute:: dataset_paths

            *list:str*

            Paths of matchup series files in matchup dataset directory

        .. py:attribute:: parameter_path

            *str*

            Path of parameters to be used as starting point for solver

        .. py:attribute:: output_dir

            *str*

            Path of directory to store output data files in

        .. py:attribute:: sensor_model

            *func*

            Python function to calculate radiance and derivatives given input sensor state data

        .. py:attribute:: adjustment_model

            *func*

            Python function to calculate spectral adjustment factor between two sensors, k, and derivatives given the
            two sensor radiances

        .. py:attribute:: software

            *str*

            software implementation name

        .. py:attribute:: software_version

            *str*

            software implementation version

        .. py:attribute:: software_tag

            *str*

            software implementation vcs

        .. py:attribute:: job_id

            *str*

            job configuratoin ID

        .. py:attribute:: matchup_dataset

            *str*

            harmonisation dataset set name

        .. py:attribute:: HarmData

            *cls*

            Harmonisation data reader

        .. py:attribute:: hout_path

            *str*

            path to store harmonisation output file

        .. py:attribute:: hres_paths

            *str*

            path to store harmonisation residual files

    :Methods:
        .. py:method:: run(...):

            This function runs the harmonisation of satellite instrument calibration parameters for group of sensors
            with a reference sensor from the match-up data located in the input directory

    """

    def __init__(self, dataset_paths=None, parameter_path=None, output_dir=None, sensor_model=None,
                 adjustment_model=None, software_cfg=None, data_reader=None, hout_path=None, hres_paths=None):
        """
        Initialise harmonisation algorithm class

        :type dataset_paths: list:str
        :param dataset_paths: Paths of matchup series files in matchup dataset directory

        :type parameter_path: str
        :param parameter_path: Path of parameters to be used as starting point for solver

        :type output_dir: str
        :param output_dir: Path of directory to store output data files in

        :type sensor_model: func
        :param sensor_model: Python function to calculate radiance and derivatives given input sensor state data

        :type adjustment_model: func
        :param adjustment_model: Python function to calculate spectral adjustment factor between two sensors, k, and
        derivatives given the two sensor radiances

        :type software_cfg: dict:str
        :param software_cfg: dictionary of software configuration information

        :type data_reader: cls
        :param data_reader: Python class to open harmonisation data. If none given default data reader used.

        :type hout_path: str
        :param hout_path: path of harmonisation output file

        :type hres_paths: str
        :param hres_paths: path of harmonisation residual files
        """

        self.dataset_paths = None
        self.parameter_path = None
        self.output_dir = None
        self.sensor_model = None
        self.adjustment_model = None
        self.software = None
        self.software_version = None
        self.software_tag = None
        self.job_id = None
        self.matchup_dataset = None
        self.HarmData = None
        self.hout_path = None
        self.hres_paths = None

        if dataset_paths is not None:
            self.dataset_paths = dataset_paths
            self.parameter_path = parameter_path
            self.output_dir = output_dir

        if (sensor_model is not None) and (adjustment_model is not None):
            self.sensor_model = sensor_model
            self.adjustment_model = adjustment_model

        if software_cfg is not None:
            self.software = software_cfg['software']
            self.software_version = software_cfg['version']
            self.software_tag = software_cfg['tag']
            self.job_id = software_cfg["job_id"]
            self.matchup_dataset = software_cfg['matchup_dataset']

        else:
            self.software = "MM"
            self.software_version = "V.V"
            self.software_tag = "TTTTTTT"
            self.job_id = "CC"
            self.matchup_dataset = "TEST"

        if data_reader is not None:
            self.HarmData = data_reader

        else:
            import harm_data_reader
            self.HarmData = harm_data_reader.HarmData

        if hout_path is not None:
            self.hout_path = hout_path
            self.hres_paths = hres_paths

    def run(self, tolPC=TOLPC, tol=TOL, tolA=TOLA, tolB=TOLB, tolU=TOLU, show=False):
        """
        This function runs the harmonisation of satellite instrument calibration parameters for group of sensors with a
        reference sensor from the match-up data located in the input directory.

        It first reads the match-up data, computes a pre-conditioned solution with a sample of the data and then runs
        a Gauss-Newton iteration algorithm to perform the harmonisation.

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
        parameter_path = self.parameter_path
        output_dir = self.output_dir

        # 2. Functions
        sensor_model = self.sensor_model
        adjustment_model = self.adjustment_model

        # 3. Software Info
        software = self.software
        software_version = self.software_version
        software_tag = self.software_tag
        job_id = self.job_id
        matchup_dataset = self.matchup_dataset

        # 4. Residual Data
        hout_path = self.hout_path
        hres_paths = self.hres_paths

        # Default to save residual data
        res = True

        ################################################################################################################
        # 1.	Read Harmonisation Matchup Data
        ################################################################################################################

        print("Opening Data...")
        HData = self.HarmData(dataset_paths, parameter_path, sensor_model, adjustment_model)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # MC Trial Test Routine
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #
        # Generate errors for the data if MC trials
        if (hout_path is not None) and (hres_paths is not None):

            print("Adding Errors for MC Trial...")
            # # saving residuals not required for MC trial results
            # res = False

            # a. add residuals of previous run to find best estimates of data values
            with HarmOutput(hout_path, hres_paths) as HOut:
                HData.values[:, :] += HOut.H_res[:, :]
                HData.ks[:] += HOut.k_res[:]

                if array_equal(HData.idx['Ia'], HOut.parameter_sensors):
                    HData.a[:] = HOut.parameter[:]
                else:
                    raise Exception("Parameter mismatch: Ordering of match-up data and residual parameters different")

            # b. adjust best estimates with errors
            HData = gen_errors(HData)
        #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ################################################################################################################
        # 2.	Perform harmonisation
        ################################################################################################################

        Harmonisation = HarmAlgo(HData)

        HOut = HarmOutput()
        HOut.parameter, HOut.parameter_covariance_matrix, HOut.cost, \
        HOut.cost_dof, HOut.cost_p_value, HOut.H_res, HOut.k_res = Harmonisation.run()

        print "Final Solution:"
        print HOut.parameter
        print HOut.parameter_covariance_matrix

        ################################################################################################################
        # 3.	Write data to file
        ################################################################################################################

        print 'Writing data to file...'

        # Get metadata
        HOut.parameter_sensors = HData.idx["Ia"]
        HOut.lm = HData.idx['lm']
        HOut.software = software
        HOut.software_version = software_version
        HOut.software_tag = software_tag
        HOut.job_id = job_id
        startDate = str(HData.times[0].year) + '{:02d}'.format(HData.times[0].month) + str(HData.times[0].day)
        endDate = str(HData.times[-1].year) + '{:02d}'.format(HData.times[-1].month) + str(HData.times[-1].day)
        HOut.matchup_dataset = "_".join((matchup_dataset, startDate, endDate))

        HOut.save(output_dir, res=res)

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
            hout_dir = None
            n_trial = None

        elif len(argv) == 4:
            job_cfg_fname = os.path.abspath(argv[1])
            hout_dir = os.path.abspath(argv[2])
            n_trial = int(argv[3])

        # 2. Read configuration data
        conf = {}   # dictionary to store data

        #  a. Read software config file
        software_cfg_fname = "software.cfg"
        conf['software'], conf['version'], conf['tag'], conf['software_text'] = read_software_cfg(software_cfg_fname)

        # b. Read job config file
        conf['job_id'], conf['matchup_dataset'], dataset_dir, parameter_path, output_dir, \
            sensor_functions_path, data_reader_path, conf['job_text'] = read_job_cfg(job_cfg_fname)

        # 3. Get matchup data paths from directory
        dataset_paths = get_dataset_paths(dataset_dir)

        # 4. Import required specified functions
        sensor_functions = import_file(sensor_functions_path)
        harm_data_reader = import_file(data_reader_path)

        # 5. Get harmonisation output files paths if Monte Carlo run
        hout_path = None
        hres_paths = None
        if hout_dir is not None:

            # get paths of individual harmonisation files in directory
            hout_path, hres_paths = get_harm_paths(hout_dir)

            # set new output directory to one specific to MC trial
            output_dir = pjoin(output_dir, "mc", '{:03d}'.format(n_trial))

        # 6. Make output directory if it doesn't exist
        try:
            makedirs(output_dir)
        except OSError:
            pass

        ################################################################################################################
        # Run harmonisation
        ################################################################################################################

        # Initialise object
        H = HarmOp(dataset_paths=dataset_paths,
                   parameter_path=parameter_path,
                   output_dir=output_dir,
                   sensor_model=sensor_functions.sensor_model,
                   adjustment_model=sensor_functions.adjustment_model,
                   software_cfg=conf,
                   data_reader=harm_data_reader.HarmData,
                   hout_path=hout_path,
                   hres_paths=hres_paths)

        # Run algorithm
        H.run(tolPC=TOLPC, tol=TOL, tolA=TOLA, tolB=TOLB, tolU=TOLU, show=True)

        return 0

    main()
