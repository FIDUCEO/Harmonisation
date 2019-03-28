"""
Tests for MatchUpReader Class
"""

'''___Built-In Modules___'''
import unittest
from os.path import join as pjoin
from os import getcwd
from os.path import dirname, abspath
from datetime import datetime

'''___Third-Party Modules___'''
from numpy.ma import MaskedArray
from numpy import ndarray, nan, isnan, array
from netCDF4 import Dataset

'''___harmonisation Modules___'''
from harmonisation.sensor_data.test_sim.test_sim import SENSOR_DATA


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "18/11/2017"
__credits__ = ["Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

# Constant
temp_data_directory = pjoin(getcwd(), "temp")
DATASET_PATHS = [abspath("../../../../../data/simulated_matchup/0_1.nc"),
                 abspath("../../../../../data/simulated_matchup/1_2.nc"),
                 abspath("../../../../../data/simulated_matchup/2_3.nc")]
IDX_SINGLE = {"Nm": [10000],
              "cNm": [0, 10000],
              "sensors": ["0", "1"],
              "Im": [[0, 1]],
              "sensor_ms": [1, 4],
              "n_sensor": [0, 1, 1, 1, 1],
              "n_mu": [1, 1, 1, 1, 1],
              "n_cov": [1, 1, 2, 3, 4],
              "N_var": [10000, 10000, 10000, 10000, 10000],
              "idx": [0, 10000, 20000, 30000, 40000, 50000]}
IDX_SINGLE_ADDITIONAL_VALUES = {"Nm": [10000],
                                "cNm": [0, 10000],
                                "sensors": ["0", "1"],
                                "Im": [[0, 1]],
                                "sensor_ms": [1, 4],
                                "n_sensor": [0, 1, 1, 1, 1],
                                "n_mu": [1, 1, 1, 1, 1],
                                "n_cov": [1, 1, 2, 3, 4],
                                "N_var": [10000, 10000, 10000, 10000, 10000],
                                "idx": [0, 10000, 20000, 30000, 40000, 50000],
                                "additional_values_name": ["nominal_measurand1", "nominal_measurand2",
                                                           "extra_variable"]}
IDX_MULTI = {"Nm": [10000, 10000, 10000],
             "cNm": [0, 10000, 20000, 30000],
             "sensors": ["0", "1", "2", "3"],
             "sensor_ms": [1, 4, 4, 4],
             "Im": [[0, 1], [1, 2], [2, 3]],
             "n_sensor": [0, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3],
             "n_mu": [1, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3],
             "n_cov": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
             "N_var": [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,
                       10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000],
             "idx": [0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000,
                      140000, 150000, 160000, 170000, 180000, 190000, 200000, 210000]}
IDX_MULTI_ADDITIONAL_VALUES = {"Nm": [10000, 10000, 10000],
                               "cNm": [0, 10000, 20000, 30000],
                               "sensors": ["0", "1", "2", "3"],
                               "sensor_ms": [1, 4, 4, 4],
                               "Im": [[0, 1], [1, 2], [2, 3]],
                               "n_sensor": [0, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3],
                               "n_mu": [1, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3],
                               "n_cov": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
                               "N_var": [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,
                                         10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000],
                               "idx": [0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000,
                                       120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000, 210000],
                               "additional_values_name": ["nominal_measurand1", "nominal_measurand2", "extra_variable"]}
TEST_UNC_TYPEID = [1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
TEST_UNC_UR = [0.0750426832, None, None, None, None, None, None, None, None, None, None,
           0.251087305, 0.251087305, 0.251087305, 0.251087305, 0.251087305,
           0.122101448, 0.122101448, 0.122101448, 0.122101448, 0.122101448]
TEST_UNC_W_I = [None, 0, 1, 1, 2, 2, 0, 1, 1, 2, 2, None, None, None, None, None, None, None, None, None]
TEST_UNC_U_I = [None, 0, 1, 1, 2, 2, 0, 1, 1, 2, 2, None, None, None, None, None, None, None, None, None]
TEST_UNC_US = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
           0.204141772, 0.204141772, 0.204141772, 0.204141772, 0.204141772]
TEST_UNC_US_I = [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
             1, 1, 2, 2, 3]

TEST_VALUES_BLOCK_FIRSTS = [40.052811, 951.55823, 1002.58868, 990.082462, 992.032639, 986.811812, 457.918812,
                            456.004804, 442.620179, 441.242585, 415.619208, 737.723959, 761.736216, 790.987766,
                            760.634506, 790.361844, 102.412389, 103.398928, 89.9778417, 112.579682, 102.499891]

TEST_KS_MU_FIRSTS = [5.18839393, -12.7125831, -12.2558266]

TEST_UNCK_TYPEID = [1, 1, 1]
TEST_UNCK_UR = [0.351025364, 0.351025364, 0.351025364]

TEST_ADDITIONAL_VALUES_BLOCK_FIRSTS = [[0, 45.5201530456543, 40.0528106689453],
                                       [nan, 34.1658058166504, 46.7008666992188],
                                       [0, 36.2559852600098, 48.6974601745605]]

TEST_ACROSS_TRACK_INDEX_FIRSTS = [1, nan, 1]
TEST_TIME_FIRSTS = [datetime(1970, 1, 1, 0, 0), datetime(1970, 1, 1, 2, 46, 40), datetime(1970, 1, 1, 5, 33, 20)]

SENSOR_DATA_PATH = pjoin(abspath("../../../.."), "sensor_data", "test_sim", "test_sim.py")
TEST_A = array([0., 0., 0., 0., 0., 0., 0., 0., 0.])
IDX_MULTI_SENSOR_DATA = {"Nm": [10000, 10000, 10000],
                         "cNm": [0, 10000, 20000, 30000],
                         "sensors": ["0", "1", "2", "3"],
                         "sensor_ms": [1, 4, 4, 4],
                         "Im": [[0, 1], [1, 2], [2, 3]],
                         "n_sensor": [0, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3, 1, 1, 2, 2, 3],
                         "n_mu": [1, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3, 1, 2, 2, 3, 3],
                         "n_cov": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
                         "N_var": [10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,
                                   10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000],
                         "idx": [0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000,
                                 120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000, 210000],
                         "parameter_sensor": ["1", "1", "1", "2", "2", "2", "3", "3", "3"],
                         "sensor_model_constant_sensor": []}

def setup():
    from harmonisation.core.matchup.matchupIO.MatchUpReader import MatchUpReader
    matchup = MatchUpReader()
    return matchup


def setup_open_nc():
    from harmonisation.core.matchup.matchupIO.MatchUpReader import MatchUpReader
    matchup = MatchUpReader()
    datasets = matchup.open_nc_datasets(DATASET_PATHS)
    return matchup, datasets


class TestMatchUpReader(unittest.TestCase):

    def test_open_nc_datasets___str_dir(self):
        dataset_directory = dirname(DATASET_PATHS[0])

        matchup = setup()

        test_datasets = matchup.open_nc_datasets(dataset_directory)

        for test_dataset, DATASET_PATH in zip(test_datasets, DATASET_PATHS):
            self.assertTrue(type(test_dataset) == Dataset)
            self.assertEqual(test_dataset.filepath(), DATASET_PATH)

    def test_open_nc_datasets___str_path(self):
        dataset_path = DATASET_PATHS[0]

        matchup = setup()

        test_dataset = matchup.open_nc_datasets(dataset_path)

        self.assertEqual(len(test_dataset), 1)
        self.assertTrue(type(test_dataset[0]) == Dataset)
        self.assertEqual(test_dataset[0].filepath(), dataset_path)

    def test_open_nc_datasets___list_path(self):
        matchup = setup()

        test_datasets = matchup.open_nc_datasets(DATASET_PATHS)

        for test_dataset, DATASET_PATH in zip(test_datasets, DATASET_PATHS):
            self.assertTrue(type(test_dataset) == Dataset)
            self.assertEqual(test_dataset.filepath(), DATASET_PATH)

    def test_get_idx_from_nc_datasets__single(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets([datasets[0]])
        self.assertItemsEqual(idxs.keys(), IDX_SINGLE.keys())

        for key in IDX_SINGLE.keys():
            self.assertItemsEqual(idxs[key], IDX_SINGLE[key])

    def test_get_idx_from_nc_datasets__multi(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        self.assertItemsEqual(idxs.keys(), IDX_MULTI.keys())

        for key in IDX_MULTI.keys():
            self.assertItemsEqual(idxs[key], IDX_MULTI[key])

    def test_get_idx_Nm_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_Nm = matchup.get_idx_Nm_from_nc_datasets([datasets[0]])

        self.assertItemsEqual(test_Nm, IDX_SINGLE["Nm"])

    def test_get_idx_Nm_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_Nm = matchup.get_idx_Nm_from_nc_datasets(datasets)

        self.assertItemsEqual(test_Nm, IDX_MULTI["Nm"])

    def test_get_idx_cNm_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_cNm = matchup.get_idx_cNm_from_nc_datasets([datasets[0]])

        self.assertItemsEqual(test_cNm, IDX_SINGLE["cNm"])

    def test_get_idx_cNm_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_cNm = matchup.get_idx_cNm_from_nc_datasets(datasets)

        self.assertItemsEqual(test_cNm, IDX_MULTI["cNm"])

    def test_get_idx_sensors_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_sensors, test_sensor_ms = matchup.get_idx_sensors_from_nc_datasets([datasets[0]])

        self.assertItemsEqual(test_sensors, IDX_SINGLE["sensors"])
        self.assertItemsEqual(test_sensor_ms, IDX_SINGLE["sensor_ms"])

    def test_get_idx_sensors_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_sensors, test_sensor_ms = matchup.get_idx_sensors_from_nc_datasets(datasets)

        self.assertItemsEqual(test_sensors, IDX_MULTI["sensors"])
        self.assertItemsEqual(test_sensor_ms, IDX_MULTI["sensor_ms"])

    def test_get_idx_Im_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_Im = matchup.get_idx_Im_from_nc_datasets([datasets[0]], IDX_SINGLE["sensors"])

        self.assertItemsEqual(test_Im, IDX_SINGLE["Im"])

    def test_get_idx_Im_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_Im = matchup.get_idx_Im_from_nc_datasets(datasets, IDX_MULTI["sensors"])

        self.assertItemsEqual(test_Im, IDX_MULTI["Im"])

    def test_get_idx_n_sensor_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_n_sensor = matchup.get_idx_n_sensor_from_nc_datasets(IDX_SINGLE["Im"], IDX_SINGLE["sensor_ms"])

        self.assertItemsEqual(test_n_sensor, IDX_SINGLE["n_sensor"])

    def test_get_idx_n_sensor_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_n_sensor = matchup.get_idx_n_sensor_from_nc_datasets(IDX_MULTI["Im"], IDX_MULTI["sensor_ms"])

        self.assertItemsEqual(test_n_sensor, IDX_MULTI["n_sensor"])

    def test_get_idx_n_mu_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_n_mu = matchup.get_idx_n_mu_from_nc_datasets(IDX_SINGLE["Im"], IDX_SINGLE["sensor_ms"])

        self.assertItemsEqual(test_n_mu, IDX_SINGLE["n_mu"])

    def test_get_idx_n_mu_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_n_mu = matchup.get_idx_n_mu_from_nc_datasets(IDX_MULTI["Im"], IDX_MULTI["sensor_ms"])

        self.assertItemsEqual(test_n_mu, IDX_MULTI["n_mu"])

    def test_get_idx_n_cov_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_n_cov = matchup.get_idx_n_cov_from_nc_datasets(IDX_SINGLE["Im"], IDX_SINGLE["sensor_ms"])

        self.assertItemsEqual(test_n_cov, IDX_SINGLE["n_cov"])

    def test_get_idx_n_cov_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_n_cov = matchup.get_idx_n_cov_from_nc_datasets(IDX_MULTI["Im"], IDX_MULTI["sensor_ms"])

        self.assertItemsEqual(test_n_cov, IDX_MULTI["n_cov"])

    def test_get_idx_N_var_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_N_var = matchup.get_idx_N_var_from_nc_datasets(IDX_SINGLE["Nm"], IDX_SINGLE["n_mu"])

        self.assertItemsEqual(test_N_var, IDX_SINGLE["N_var"])

    def test_get_idx_N_var_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_N_var = matchup.get_idx_N_var_from_nc_datasets(IDX_MULTI["Nm"], IDX_MULTI["n_mu"])

        self.assertItemsEqual(test_N_var, IDX_MULTI["N_var"])

    def test_get_idx_idxs_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_idxs = matchup.get_idx_idxs_from_nc_datasets(IDX_SINGLE["N_var"])

        self.assertItemsEqual(test_idxs, IDX_SINGLE["idx"])

    def test_get_idx_idxs_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_idxs = matchup.get_idx_idxs_from_nc_datasets(IDX_MULTI["N_var"])

        self.assertItemsEqual(test_idxs, IDX_MULTI["idx"])

    def test_get_get_idx_additional_values_name_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()

        test_additional_values_name = \
            matchup.get_idx_additional_values_name_from_nc_datasets([datasets[0]])

        self.assertItemsEqual(test_additional_values_name, IDX_SINGLE_ADDITIONAL_VALUES["additional_values_name"])

    def test_get_idx_additional_values_name_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()

        test_additional_values_name = \
            matchup.get_idx_additional_values_name_from_nc_datasets(datasets)

        self.assertItemsEqual(test_additional_values_name, IDX_MULTI_ADDITIONAL_VALUES["additional_values_name"])

    def test_get_idx_from_nc_datasets___single(self):
        matchup, datasets = setup_open_nc()
        test_idx = matchup.get_idx_from_nc_datasets([datasets[0]])

        self.assertItemsEqual(sorted(test_idx.keys()), sorted(IDX_SINGLE.keys()))

    def test_get_idx_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        test_idx = matchup.get_idx_from_nc_datasets(datasets)

        self.assertItemsEqual(sorted(test_idx.keys()), sorted(IDX_MULTI.keys()))

    def test_open_matrix_idxs_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        w_matrix_idxs, u_matrix_idxs = matchup.open_matrix_idxs_from_nc_datasets(datasets)

        self.assertItemsEqual(w_matrix_idxs, [[1, 1], [2, 1], [3, 1]])
        self.assertItemsEqual(u_matrix_idxs, [[1, 1], [2, 1], [3, 1]])

    def test_get_block_idxs(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        block_idxs = matchup.get_block_idxs(idxs)

        self.assertItemsEqual(block_idxs, [(0, 1, 1), (1, 1, 1), (1, 2, 1), (2, 2, 1), (2, 3, 1), (3, 3, 1), (1, 1, 2),
                                           (1, 2, 2), (2, 2, 2), (2, 3, 2), (3, 3, 2), (1, 1, 3), (1, 2, 3), (2, 2, 3),
                                           (2, 3, 3), (3, 3, 3), (1, 1, 4), (1, 2, 4), (2, 2, 4), (2, 3, 4), (3, 3, 4)])

    def test_open_unc_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        unc, w_matrices, u_matrices = matchup.open_unc_from_nc_datasets(datasets, idxs)

        # test unc

        self.assertItemsEqual([u.typeID for u in unc], TEST_UNC_TYPEID)

        for u, test_ur in zip(unc, TEST_UNC_UR):
            if test_ur is not None:
                self.assertEqual(type(u.uR), MaskedArray)
                self.assertEqual(u.uR.shape, (10000,))
                self.assertTrue(all(uR_i == test_ur for uR_i in u.uR))

        for u, test_w_i, test_u_i in zip(unc, TEST_UNC_W_I, TEST_UNC_U_I):
            if test_w_i is not None:
                self.assertEqual(u.w_i, test_w_i)
                self.assertEqual(u.u_i, test_u_i)

        for u, test_us, test_us_i in zip(unc, TEST_UNC_US, TEST_UNC_US_I):
            if test_us is not None:
                self.assertEqual(type(u.uS), MaskedArray)
                self.assertEqual(u.uS.shape, (10000,))
                self.assertTrue(all(uS == test_us for uS in u.uS))
                self.assertEqual(u.uS_i, test_us_i)

        for i, d in enumerate(datasets):
            self.assertTrue(all(a == b for a, b in zip(w_matrices[i].indices, matchup._open_w_matrix(d, 1).indices)))
            self.assertTrue(all(a == b for a, b in zip(w_matrices[i].indptr, matchup._open_w_matrix(d, 1).indptr)))
            self.assertTrue(all(a == b for a, b in zip(w_matrices[i].data, matchup._open_w_matrix(d, 1).data)))

        for i, d in enumerate(datasets):
            self.assertTrue(all(a == b for a, b in zip(u_matrices[i], matchup._open_u_matrix(d, 1))))

    def test_open_values_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        values = matchup.open_values_from_nc_datasets(datasets, idxs)

        self.assertEqual(type(values), ndarray)
        self.assertEqual(values.shape, (210000,))

        for i in range(21):
            self.assertAlmostEquals(values[i*10000], TEST_VALUES_BLOCK_FIRSTS[i], places=4)

    def test_open_ks_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        ks = matchup.open_ks_from_nc_datasets(datasets, idxs)

        self.assertEqual(type(ks), ndarray)
        self.assertEqual(ks.shape, (30000,))

        for i in range(3):
            self.assertAlmostEquals(ks[i * 10000], TEST_KS_MU_FIRSTS[i], places=4)

    def test_open_unck_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        unck = matchup.open_unck_from_nc_datasets(datasets, idxs)

        # test unc

        self.assertItemsEqual([u.typeID for u in unck], TEST_UNCK_TYPEID)

        for u, test_ur in zip(unck, TEST_UNCK_UR):
            if test_ur is not None:
                self.assertEqual(type(u.uR), MaskedArray)
                self.assertEqual(u.uR.shape, (10000,))
                self.assertTrue(all(uR_i == test_ur for uR_i in u.uR))

    def test_open_additional_values_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        additional_values, idxs = matchup.open_additional_values_from_nc_datasets(datasets, idxs)

        self.assertItemsEqual(idxs.keys(), IDX_MULTI_ADDITIONAL_VALUES.keys())

        for key in IDX_MULTI_ADDITIONAL_VALUES.keys():
            self.assertItemsEqual(idxs[key], IDX_MULTI_ADDITIONAL_VALUES[key])

        self.assertEqual(type(additional_values), ndarray)
        self.assertEqual(additional_values.shape, (30000,3))
        for i, additional_values_block_firsts in enumerate(TEST_ADDITIONAL_VALUES_BLOCK_FIRSTS):
            for j, additional_values_block_first_col in enumerate(additional_values_block_firsts):
                if isnan(additional_values_block_first_col):
                    self.assertTrue(isnan(additional_values[i*10000, j]))
                else:
                    self.assertAlmostEquals(additional_values[i*10000, j], additional_values_block_first_col, places=4)

    def test_open_time_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        time1, time2 = matchup.open_time_from_nc_datasets(datasets, idxs)

        self.assertEqual(type(time1), ndarray)
        self.assertEqual(time1.shape, (30000,))

        self.assertEqual(type(time2), ndarray)
        self.assertEqual(time2.shape, (30000,))

        for i, time_first in enumerate(TEST_TIME_FIRSTS):
            self.assertEquals(time1[i * 10000], time_first)
            self.assertEquals(time2[i * 10000], time_first)

    def test_open_optional_index_from_nc_datasets___multi(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        across_track_index1 = matchup.open_optional_index_from_nc_datasets("across_track_index1", datasets, idxs)

        self.assertEqual(type(across_track_index1), ndarray)
        self.assertEqual(across_track_index1.shape, (30000,))

        for i, across_track_index_first in enumerate(TEST_ACROSS_TRACK_INDEX_FIRSTS):
            if isnan(across_track_index_first):
                self.assertTrue(isnan(across_track_index1[i*10000]))
            else:
                self.assertEquals(across_track_index1[i*10000], across_track_index_first)

    def test_open_sensor_data(self):
        matchup, datasets = setup_open_nc()
        sensor_data = matchup.open_sensor_data(SENSOR_DATA_PATH)

        self.assertEqual(SENSOR_DATA.keys(), sensor_data.keys())

    def test_extract_sensor_data(self):
        matchup, datasets = setup_open_nc()
        idxs = matchup.get_idx_from_nc_datasets(datasets)
        a, sensor_model_constant, sensor_model, adjustment_model, idx = matchup.extract_sensor_data(SENSOR_DATA, idxs)

        self.assertItemsEqual(a, TEST_A)
        self.assertEqual(sensor_model_constant, [])

        self.assertEqual(len(sensor_model), 4)
        self.assertTrue((type(i) == function for i in sensor_model))

        self.assertEqual(len(adjustment_model), 4)
        self.assertTrue((type(i) == function for i in adjustment_model))

        for key in IDX_MULTI_SENSOR_DATA.keys():
            self.assertItemsEqual(idxs[key], IDX_MULTI_SENSOR_DATA[key])


if __name__ == '__main__':
    unittest.main()
