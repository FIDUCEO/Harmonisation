"""
AVHRR sensor data for simulated, time-independent data
"""

'''___Python Modules____'''
from collections import OrderedDict

'''___Third Party Modules____'''

'''___NPL Modules___'''
from AVHRR_sensor_models_not_cython import sensor_model_3 as AVHRR_sensor_function
from AVHRR_sensor_models_not_cython import sensor_model_ref as reference_sensor_function

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "17/01/2018"
__credits__ = ["Peter Harris", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


sensor_data = OrderedDict([("m02",
                            {"sensor_model": reference_sensor_function,
                             "sensor_model_variables_num": 1,
                             "sensor_model_variable_labels": ["$L_{ref}$"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["Lref"],
                             "sensor_model_parameter": [],
                             "sensor_model_constant": []}),
                           ("n19",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [-2.45, 0.0392, 7.8e-6],
                             "sensor_model_constant": []}),
                           ("n18",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.227, -0.0014, 2.18e-5],
                             "sensor_model_constant": []}),
                           ("n17",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [2.479, -0.0157, 3.4e-5],
                             "sensor_model_constant": []}),
                           ("n16",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.91, -0.0095, 9.6e-6],
                             "sensor_model_constant": []}),
                           ("n15",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.802, -0.00631, 1.9e-5],
                             "sensor_model_constant": []}),
                           ("n14",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.802, -0.00631, 1.9e-5],
                             "sensor_model_constant": []}),
                           ("n12",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.802, -0.00631, 1.9e-5],
                             "sensor_model_constant": []}),
                           ("n11",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.802, -0.00631, 1.9e-5],
                             "sensor_model_constant": []}),
                           ("n10",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.802, -0.00631, 1.9e-5],
                             "sensor_model_constant": []}),
                           ("n09",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.802, -0.00631, 1.9e-5],
                             "sensor_model_constant": []}),
                           ("n08",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.802, -0.00631, 1.9e-5],
                             "sensor_model_constant": []}),
                           ("n07",
                            {"sensor_model": AVHRR_sensor_function,
                             "sensor_model_variables_num": 4,
                             "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                             "nominal_measurand_label": "L_{op}",
                             "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                             "sensor_model_parameter": [1.802, -0.00631, 1.9e-5],
                             "sensor_model_constant": []})])
