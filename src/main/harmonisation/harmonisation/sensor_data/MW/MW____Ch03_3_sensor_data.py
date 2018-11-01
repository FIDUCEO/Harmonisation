"""
AVHRR sensor data (for time-independent model)
"""

'''___Python Modules____'''
from collections import OrderedDict

'''___Third Party Modules____'''

'''___harmonisation Modules___'''
from MW_sensor_models import sensor_model_3 as MW_sensor_function
from MW_sensor_models import sensor_model_ref as reference_sensor_function
from MW_consts import W_AMSUB_T, W_MHS_T, const_AMSUB_3, const_MHS_3, antenna_corr_AMSUB_3, antenna_corr_MHS_3


'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "17/01/2018"
__credits__ = ["Imke Hans", "Martin Burgdorf"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


'''___Constants___'''
# sensor_model_parameter
# NB: 3 parameter model
sensor_model_parameter = [0.0, 0.0, 0.0]

# sensor_data_dict
# NB: same for all sensors (is different for NOAA-18 but is reference)
AMSUB_data_dict = {"sensor_model": MW_sensor_function,
                   "sensor_model_variables_num": 9,
                   "sensor_model_variable_labels": [r"$\bar{C}_{S}$", r"$\bar{C}_{IWCT}$",
                                                    r"\bar{T}_{IWCT, 1}", r"\bar{T}_{IWCT, 2}",
                                                    r"\bar{T}_{IWCT, 3}", r"\bar{T}_{IWCT, 4}",
                                                    r"\bar{T}_{IWCT, 5}", "$C_{E}$", "$T_{LO}$"],
                   "sensor_model_variable_names": ["C_S", "C_IWCT", "T_IWCT1", "T_IWCT2", "T_IWCT3",
                                                   "T_IWCT4", "T_IWCT5", "C_E", "T_LO"],
                   "nominal_measurand_label": "BT_{op}",
                   "sensor_model_parameter": sensor_model_parameter,
                   "sensor_model_constant": antenna_corr_AMSUB_3+W_AMSUB_T+const_AMSUB_3}

MHS_data_dict = {"sensor_model": MW_sensor_function,
                 "sensor_model_variables_num": 9,
                 "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{IWCT}",
                                                  r"\bar{T}_{IWCT, 1}", r"\bar{T}_{IWCT, 2}",
                                                  r"\bar{T}_{IWCT, 3}", r"\bar{T}_{IWCT, 4}",
                                                  r"\bar{T}_{IWCT, 5}", "C_E", "T_{LO}"],
                 "sensor_model_variable_names": ["C_S", "C_IWCT", "T_IWCT1", "T_IWCT2", "T_IWCT3",
                                                 "T_IWCT4", "T_IWCT5", "C_E", "T_LO"],
                 "nominal_measurand_label": "BT_{op}",
                 "sensor_model_parameter": sensor_model_parameter,
                 "sensor_model_constant": antenna_corr_MHS_3+W_MHS_T+const_MHS_3}

SENSOR_DATA = OrderedDict([("18",
                            {"sensor_model": reference_sensor_function,
                             "sensor_model_variable_labels": ["BT_{ref}"],
                             "sensor_model_variable_names": ["BTref"],
                             "nominal_measurand_label": "BT_{op}",
                             "sensor_model_variables_num": 1,
                             "sensor_model_parameter": [],
                             "sensor_model_constant": []}),
                           ("A", MHS_data_dict),
                           ("B", MHS_data_dict),
                           ("19", MHS_data_dict),
                           ("17", AMSUB_data_dict),
                           ("16", AMSUB_data_dict),
                           ("15", AMSUB_data_dict)])


class MWCh03_3:
    def __init__(self):
        self.name = "MW____Ch03_3"

    def get_sensor_data(self):
        return SENSOR_DATA

if __name__ == "__main__":
    pass
