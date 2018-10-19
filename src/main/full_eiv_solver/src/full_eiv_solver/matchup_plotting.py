"""
Harmonsation Output Plotting Implementation
"""

'''___Python Modules___'''
from os import makedirs
from os.path import basename, dirname
from os.path import join as pjoin
import sys
from sys import argv

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''
from common import *

main_directory = dirname(dirname(__file__))
sys.path.append(main_directory)
from nplcore import MatchUpVis

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "23/12/2017"
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

def try_makedirs(directory):
    try:
        makedirs(directory)
    except OSError:
        pass
    return 0

class MatchUpPlottingOp:
    """
    This class runs the matchup dataset plotting routiune
    """

    def __init__(self, dataset_paths=None, sensor_data_path=None, output_dir=None,
                 software_cfg=None, data_reader=None):
        """
        Initialise object

        :type dataset_paths: list:str
        :param dataset_paths: Paths of matchup series files in matchup dataset directory

        :type sensor_data_path: str
        :param sensor_data_path: Path of sensor data for of matchup sensors

        :type output_dir: str
        :param output_dir: Path of directory to store output data files in

        :type software_cfg: dict:str
        :param software_cfg: dictionary of software configuration information

        :type data_reader: cls
        :param data_reader: Python class to open harmonisation data. If none given default data reader used.
        """

        self.dataset_paths = None
        self.sensor_data_path = None
        self.output_dir = None
        self.plot_dir = None
        self.software = None
        self.software_version = None
        self.software_tag = None
        self.job_id = None
        self.matchup_dataset = None
        self.HarmData = None

        if (dataset_paths is not None) and (sensor_data_path is not None):
            self.dataset_paths = dataset_paths
            self.sensor_data_path = sensor_data_path
            self.output_dir = output_dir
            self.plots_dir = pjoin(dirname(dataset_paths[0]), "plots")
            try:
                makedirs(self.plots_dir)
            except OSError:
                pass

        if software_cfg is not None:
            #self.software = software_cfg['software']
            #self.software_version = software_cfg['version']
            #self.software_tag = software_cfg['tag']
            self.job_id = software_cfg["job_id"]
            self.matchup_dataset = software_cfg['matchup_dataset']

        else:
            self.software = "MM"
            self.software_version = "V.V"
            self.software_tag = "TTTTTTT"
            self.job_id = "CC"
            self.matchup_dataset = "TEST"

        if data_reader is not None:
            self.data_reader = data_reader

        else:
            from nplcore import MatchUp
            self.data_reader = MatchUp

    def run(self):
        """
        Generates plots for dataset specified by paths defined object attributes
        """

        # Initialise
        dataset_paths = self.dataset_paths
        sensor_data_path = self.sensor_data_path
        plots_dir = self.plots_dir

        ################################################################################################################
        # 1.	Read Match-up Data Data
        ################################################################################################################

        print"Match-up Dataset:"
        for path in dataset_paths:
            print ">", path

        print"\nSensor Data:"
        print ">", sensor_data_path

        print("\nOpening Files...")
        MatchUpData = self.data_reader(dataset_paths, sensor_data_path, open_uncertainty=False)

        ################################################################################################################
        # 2.	Make plots
        ################################################################################################################

        print("Generating Plots...")

        MUVisOp = MatchUpVis(MatchUpData)

        # Define plot directories
        overview_plots_dir = pjoin(plots_dir, "overview")

        # Make plot directories
        try_makedirs(overview_plots_dir)

        # 1. Make overview plots

        # a. Nominal kres vs. time
        MUVisOp.plot_kres_nom_nom_scatter(overview_plots_dir)
        MUVisOp.plot_kres_nom_nom_scatter_monthlymean(overview_plots_dir)

        # b. Nominal L vs. Nominal L1
        MUVisOp.plot_L1_nom_v_L2_nom_scatter(overview_plots_dir)

        # c. Satellite angles polar plot
        # todo - angles plots

        # d. Longitude-Latitude Polar Plot
        # todo - lon-lat plots

        print "\nPlots written to:", plots_dir


def main(job_cfg_fname):

    ################################################################################################################
    # Process configuration data
    ################################################################################################################

    print "Match-Up Output Plotting \n"

    print "Reading Job Config:", job_cfg_fname, "\n"

    # 1. Read configuration data
    conf = {}   # dictionary to store data

    # b. Read job config file
    conf['job_id'], conf['matchup_dataset'], dataset_dir, sensor_data_path,\
        output_dir, data_reader_path, conf['job_text'] = read_job_cfg(job_cfg_fname)

    # 2. Get matchup data paths from directory
    dataset_paths = get_dataset_paths(dataset_dir)

    # 3. Import required specified functions
    if basename(data_reader_path) == "DEFAULT":
        harm_data_reader = None
    else:
        harm_data_reader = import_file(data_reader_path).MatchUp

    ################################################################################################################
    # Run harmonisation
    ################################################################################################################

    # Initialise object
    M = MatchUpPlottingOp(dataset_paths=dataset_paths,
                          sensor_data_path=sensor_data_path,
                          output_dir=output_dir,
                          software_cfg=conf,
                          data_reader=harm_data_reader)

    # Run algorithm
    M.run()

    return 0

if __name__ == "__main__":
    main(os.path.abspath(argv[1]))
    # main("/home/seh2/src/DeployedProjects/full_EIV/src/test/test_configs/AVHRR_RSIM_3_test_newdata/AVHRR_RSIM_3_test_newdata.cfg")
