"""
AVHRR sensor models
"""

'''___Python Modules____'''

'''___Third Party Modules____'''
from numpy import vstack, ones, dot, column_stack

'''___NPL Modules___'''

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris", "Jon Mittaz"]
__created__ = "17/01/2018"
__credits__ = []
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def sensor_model_ref(sensor_model_variables,
                     sensor_model_parameter,
                     sensor_model_constant,
                     sensor_model_across_track_index,
                     sensor_model_along_track_index, time):
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

    R = sensor_model_variables[:, 0]

    return R, None


def sensor_model_3(sensor_model_variables,
                   sensor_model_parameter,
                   sensor_model_constant,
                   sensor_model_across_track_index,
                   sensor_model_along_track_index, time):
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

    # Extract covariates
    Cs = sensor_model_variables[:, 0]
    Cict = sensor_model_variables[:, 1]
    Ce = sensor_model_variables[:, 2]
    Rict = sensor_model_variables[:, 3]
    m = len(Cs)

    # Evaluate radiance using sensor model
    Ja = vstack((ones(m), Rict * (Cs - Ce) / (Cs - Cict), (Cict - Ce) * (Cs - Ce))).T
    at = [sensor_model_parameter[0], 0.98514 + sensor_model_parameter[1], sensor_model_parameter[2]]
    R = dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/da0
    # > dR/da1
    # > dR/da2

    J = column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2 + at[2] * (Cict - Ce),
                      at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2 + at[2] * (Cs - Ce),
                      -at[1] * Rict / (Cs - Cict) - at[2] * (Cs + Cict - 2 * Ce),
                      at[1] * (Cs - Ce) / (Cs - Cict),
                      Ja))

    return R, J


def sensor_model_4((sensor_model_variables,
                   sensor_model_parameter,
                   sensor_model_constant,
                   sensor_model_across_track_index,
                   sensor_model_along_track_index, time)):
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

    # Extract covariates
    Cs = sensor_model_variables[:, 0]
    Cict = sensor_model_variables[:, 1]
    Ce = sensor_model_variables[:, 2]
    Rict = sensor_model_variables[:, 3]
    T = sensor_model_variables[:, 4]
    m = len(Cs)

    # Evaluate radiance using sensor model
    Ja = vstack((ones(m), Rict * (Cs - Ce) / (Cs - Cict), (Cict - Ce) * (Cs - Ce), T)).T
    at = [sensor_model_parameter[0], 0.98514 + sensor_model_parameter[1],
          sensor_model_parameter[2], sensor_model_parameter[3]]
    R = dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/dT_orb
    # > dR/da0
    # > dR/da1
    # > dR/da2
    # > dR/da3

    J = column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2 + at[2] * (Cict - Ce),
                      at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2 + at[2] * (Cs - Ce),
                      -at[1] * Rict / (Cs - Cict) - at[2] * (Cs + Cict - 2 * Ce),
                      at[1] * (Cs - Ce) / (Cs - Cict),
                      at[3] * ones(m),
                      Ja))

    return R, J

def sensor_model_3_scaled(sensor_model_variables,
                          sensor_model_parameter,
                          sensor_model_constant,
                          sensor_model_across_track_index,
                          sensor_model_along_track_index, time):
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

    # Extract covariates
    Cs = sensor_model_variables[:, 0]
    Cict = sensor_model_variables[:, 1]
    Ce = sensor_model_variables[:, 2]
    Rict = sensor_model_variables[:, 3]
    m = len(Cs)

    # Evaluate radiance using sensor model
    Ja = vstack((ones(m), Rict * (Cs - Ce) / (Cs - Cict), (Cict - Ce) * (Cs - Ce))).T
    at = [sensor_model_parameter[0], 0.98514 + sensor_model_parameter[1]/1000.0, sensor_model_parameter[2]/1000.0/1000.0]
    R = dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/da0
    # > dR/da1
    # > dR/da2

    J = column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2 + at[2] * (Cict - Ce),
                      at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2 + at[2] * (Cs - Ce),
                      -at[1] * Rict / (Cs - Cict) - at[2] * (Cs + Cict - 2 * Ce),
                      at[1] * (Cs - Ce) / (Cs - Cict),
                      Ja))

    return R, J


def sensor_model_4_scaled(sensor_model_variables, sensor_model_parameter, sensor_model_constant, sensor_model_index):
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

    # Extract covariates
    Cs = sensor_model_variables[:, 0]
    Cict = sensor_model_variables[:, 1]
    Ce = sensor_model_variables[:, 2]
    Rict = sensor_model_variables[:, 3]
    T = sensor_model_variables[:, 4]
    m = len(Cs)

    # Evaluate radiance using sensor model
    Ja = vstack((ones(m), Rict * (Cs - Ce) / (Cs - Cict), (Cict - Ce) * (Cs - Ce), T)).T
    at = [sensor_model_parameter[0], 0.98514 + sensor_model_parameter[1]/1000.0,
          sensor_model_parameter[2]/1000.0/1000.0, sensor_model_parameter[3]/1000.0]
    R = dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/dT_orb
    # > dR/da0
    # > dR/da1
    # > dR/da2
    # > dR/da3

    J = column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2 + at[2] * (Cict - Ce),
                      at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2 + at[2] * (Cs - Ce),
                      -at[1] * Rict / (Cs - Cict) - at[2] * (Cs + Cict - 2 * Ce),
                      at[1] * (Cs - Ce) / (Cs - Cict),
                      at[3] * ones(m),
                      Ja))

    return R, J


def sensor_model_37_scaled(sensor_model_variables, sensor_model_parameter, sensor_model_constant, sensor_model_index):
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

    # Extract covariates
    Cs = sensor_model_variables[:, 0]
    Cict = sensor_model_variables[:, 1]
    Ce = sensor_model_variables[:, 2]
    Rict = sensor_model_variables[:, 3]
    T = sensor_model_variables[:, 4]
    m = len(Cs)

    # Evaluate radiance using sensor model
    Ja = vstack((ones(m), Rict * (Cs - Ce) / (Cs - Cict), T-285.0)).T
    at = [sensor_model_parameter[0], 0.98514 + sensor_model_parameter[1]/1000.0, sensor_model_parameter[2]/1000.0]
    R = dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/dT_orb
    # > dR/da0
    # > dR/da1
    # > dR/da2
    # > dR/da3

    J = column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2,
                      at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2,
                      -at[1] * Rict / (Cs - Cict),
                      at[1] * (Cs - Ce) / (Cs - Cict),
                      at[2] * ones(m),
                      Ja))

    return R, J