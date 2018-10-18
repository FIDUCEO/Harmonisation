"""
Function to read HIRS sensor data
"""

'''___Python Modules____'''
from collections import OrderedDict
from os.path import dirname
from os.path import join as pjoin

'''___Third Party Modules____'''

'''___NPL Modules___'''
from HIRS_sensor_models import sensor_model_3 as HIRS_sensor_function
from HIRS_sensor_models import sensor_model_ref as reference_sensor_function

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "21/01/2018"
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def read_coeff_data(fname):
    """
    :type fname: str
    :param fname: path of HIRS coeff file

    :return:
        :sensor_coeffs: *dict*

        HIRS coeff data per sensor
    """
    sensor_coeffs = {}
    with open(fname) as f:
        for i, line in enumerate(f):
            if i != 0:
                line_list = line.split()
                sensor_name = line_list[0]
                sensor_coeff = [float(v) for v in line_list[1:4]]
                sensor_coeffs[sensor_name] = sensor_coeff

    return sensor_coeffs

def return_sensor_data(channel_num):
    """
    :type channel_num: int
    :param channel_num: channel number

    :return:
        :sensor_data: *dict*

        Sensor data dictionary for given channel
    """
    # 1. Read coeff data
    sensor_coeffs = read_coeff_data(pjoin(dirname(__file__), "coeffs/coef_ch" + str(channel_num) + ".dat"))

    # 2. Define Sensor Data
    # a. Reference Sensor Data
    sensor_data = OrderedDict([("iasi",
                                {"sensor_model": reference_sensor_function,
                                 "sensor_model_variable_labels": ["L_{ref}"],
                                 "sensor_model_variable_names": ["Lref"],
                                 "sensor_model_variables_num": 1,
                                 "sensor_model_parameter": [],
                                 "sensor_model_constant": []})])

    sensors = ["ma", "mb", "n19", "n18", "n17", "n16", "n15", "n14", "n12", "n11", "n10",
               "n09", "n08", "n07", "n06", "tn"]

    # b. Harmonisation Sensor Data
    for sensor in sensors:
        sensor_data[sensor] = {"sensor_model": HIRS_sensor_function,
                               "sensor_model_variables_num": 5,
                               "sensor_model_variable_labels": [r"\bar{C}_{S}", r"\bar{C}_{IWCT}",
                                                                "C_E", "T_{ICT}", "R_{\mathrm{self}}"],
                               "sensor_model_variable_names": ["C_S", "C_ICT", "C_E", "T_ICT", "R_self"],
                               "nominal_measurand_label": "L_{op}",
                               "sensor_model_parameter": [0.0, 0.0, 0.0],
                               "sensor_model_constant": sensor_coeffs[sensor]}

    return sensor_data