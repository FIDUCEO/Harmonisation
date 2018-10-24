"""
Sensor Measurement functions for HIRS
"""

'''___Python Modules___'''

'''___Third Party Modules___'''
from numpy import vstack, ones, dot, column_stack, exp

'''___harmonisation Modules___'''


'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "16/01/2018"
__credits__ = ["Gerrit Holl", "Ralf Quast"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


'''___Natural Constants___'''
# Planck's Constant
H = 6.626070040e-34

# Speed of light
C = 299792458

# Boltzmann's Constant
KB = 1.38064852e-23

# Temperature of the Cosmic Microwave Background
TCMB = 2.72548

'''___Sensor Constants___'''
e = 0.98


def sensor_model_ref(sensor_model_variables, sensor_model_parameter, sensor_model_constant,
                     sensor_model_across_track_index, sensor_model_along_track_index, sensor_time):
    """
    :type sensor_model_variables: list:numpy.ndarray
    :param sensor_model_variables: list of sensor state variables

    :type sensor_model_parameter: list
    :param sensor_model_parameter: sensor calibration parameters

    :type sensor_model_constant: list
    :param sensor_model_constant: sensor model constants

    :return:
        :measurand: *numpy.ndarray*

        Determined sensor measurand

        :measurand_derivatives: *numpy.ndarray*

        Measurement model derivatives w.r.t. input quantities
    """

    R = sensor_model_variables[0]

    return R, None


def sensor_model_3(sensor_model_variables, sensor_model_parameter, sensor_model_constant,
                   sensor_model_across_track_index, sensor_model_along_track_index, sensor_time):
    """
    :type sensor_model_variables: *numpy.ndarray*
    :param sensor_model_variables: list of sensor state variables

    :type sensor_model_parameter: list
    :param sensor_model_parameter: sensor calibration parameters

    :type sensor_model_constant: list
    :param sensor_model_constant: sensor model constants

    :return:
        :measurand: *numpy.ndarray*

        Determined sensor measurand

        :measurand_derivatives: *numpy.ndarray*

        Measurement model derivatives w.r.t. input quantities
    """

    # Extract covariates
    Cs = sensor_model_variables[:, 0]
    Cict = sensor_model_variables[:, 1]
    Ce = sensor_model_variables[:, 2]
    Tict = sensor_model_variables[:, 3]
    Rself = sensor_model_variables[:, 4]
    m = len(Cs)

    v = sensor_model_constant[0]
    a = sensor_model_constant[1]
    b = sensor_model_constant[2]

    # Evaluate ICT radiance
    planck_exp = exp((v * H) / (KB * (a + b * Tict)))
    planck_c = (2.0 * H * v ** 3 / C ** 2)
    Rict = planck_c / (planck_exp - 1)

    # Evaluate radiance using sensor model
    Ja = vstack((ones(m), Rict * (Ce - Cs) / (Cict - Cs), (Ce - Cict) * (Ce - Cs))).T
    at = [sensor_model_parameter[0], (e + sensor_model_parameter[1]), sensor_model_parameter[2]]
    R = dot(Ja, at) - Rself

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dT_ICT
    # > dR/dRself
    # > dR/da1
    # > dR/da2
    # > dR/da3

    J = column_stack((at[1] * Rict * (Ce - Cict) / (Cict-Cs) ** 2 - at[2] * (Ce - Cict),
                      -at[1] * Rict * (Ce - Cs) / (Cict - Cs) ** 2 - at[2] * (Ce - Cs),
                      at[1] * Rict / (Cict - Cs) + at[2] * (2 * Ce - Cs - Cict),
                      (at[1] * (Ce - Cs) / (Cict - Cs)) *
                      Rict * planck_exp / (planck_exp-1.0) * (2.0*H*v*b) / (KB*(a+b*Tict)**2.0),
                      -1.0 * ones(m),
                      Ja))

    return R, J


def planck(v, T):
    return 2.0 * H * v**3 / C**2 / (exp(v*H/(KB*T))-1)
