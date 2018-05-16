"""
AVHRR sensor data (for time-independent model)
"""

'''___Python Modules____'''
from collections import OrderedDict

'''___Third Party Modules____'''

'''___NPL Modules___'''
from AVHRR_sensor_models import sensor_model_3_scaled as AVHRR_sensor_function
from AVHRR_sensor_models import sensor_model_ref as reference_sensor_function

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "17/01/2018"
__credits__ = ["Peter Harris", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


# Build sensor_data dictionary

# 1. Initialise
sensor_data = OrderedDict()

# 2. Add Reference sensors
reference_sensors = ["aatsr", "atsr2", "atsr1"]
for reference_sensor in reference_sensors:
    sensor_data[reference_sensor] = {"sensor_model": reference_sensor_function,
                                     "sensor_model_variable_labels": ["L_{ref}"],
                                     "nominal_measurand_label": "L_{op}",
                                     "sensor_model_variable_names": ["Lref"],
                                     "sensor_model_variables_num": 1,
                                     "sensor_model_parameter": [],
                                     "sensor_model_constant": []}
# 3. Add Harmonisation Sensors
sensors = ["m02", "n19", "n18", "n17", "n16", "n15", "n14", "n12", "n11"]
for sensor in sensors:
    sensor_data[sensor] = {"sensor_model": AVHRR_sensor_function,
                           "sensor_model_variables_num": 4,
                           "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{ICT}", "C_E", "L_{ICT}"],
                           "nominal_measurand_label": "L_{op}",
                           "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "L_ICT"],
                           "sensor_model_parameter": [0.0, 0.0, 0.0],
                           "sensor_model_constant": []}