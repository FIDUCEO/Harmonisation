"""
AVHRR sensor data for simulated match-up dataset (time-independent model)
"""

'''___Python Modules____'''
from collections import OrderedDict

'''___Third Party Modules____'''

'''___EIV Modules___'''
from AVHRR_sensor_models import sensor_model_4 as AVHRR_sensor_function
from AVHRR_sensor_models import sensor_model_ref as reference_sensor_function

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "17/01/2018"
__credits__ = ["Peter Harris", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"

'''___Constants___'''
SENSOR_DATA = OrderedDict([("m02",
                            {"sensor_model": reference_sensor_function,
                             "sensor_model_variables_num": 1,
                             "sensor_model_variable_labels": ["$L_{ref}$"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["Lref"],
                             "sensor_model_parameter": [],
                             "sensor_model_constant": []}),
                           ("n19",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-33.45, 0.0392, 0.0000078, 0.108],
                             "sensor_model_constant": []}),
                           ("n18",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [17.63, -0.001435, 0.0000218, -0.057],
                             "sensor_model_constant": []}),
                           ("n17",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-11.90, -0.0158, 0.000034, 0.05],
                             "sensor_model_constant": []}),
                           ("n16",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-14.24, -0.010093, 0.0000096, 0.0552],
                             "sensor_model_constant": []}),
                           ("n15",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-13.82, -7.31e-3, 1.90e-5, 0.0535],
                             "sensor_model_constant": []}),
                           ("n14",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-24.95, -0.01061, 0.000023, 0.0937],
                             "sensor_model_constant": []}),
                           ("n12",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-45.60, -0.00213, 0.000014, 0.1638],
                             "sensor_model_constant": []}),
                           ("n11",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-9.31, 0.04558, 0.0000099, 0.02041],
                             "sensor_model_constant": []}),
                           ("n10",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-10.44, 0.03732, 0.000023, 0.0293],
                             "sensor_model_constant": []}),
                           ("n09",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-34.41, -0.00582, 0.00001, 0.124],
                             "sensor_model_constant": []}),
                           ("n08",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-133.26, -0.02063, 0.0000165, 0.4626],
                             "sensor_model_constant": []}),
                           ("n07",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 5,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}", "T_o"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT", "T_o"],
                             "sensor_model_parameter": [-6.43, -0.0425, 0.000023, 0.0386],
                             "sensor_model_constant": []})])


class AVHRRSim4:
    def __init__(self):
        self.name = "AVHRR_sim_4"

    def get_sensor_data(self):
        return SENSOR_DATA


if __name__ == "__main__":
    pass
