"""
AVHRR sensor models
"""

'''___Python Modules____'''

'''___Third Party Modules____'''
from numpy import ones, vstack, column_stack, dot
import numpy as np
cimport numpy as np
DTYPEi = np.int32
DTYPEf = np.float32
DTYPEd = np.float64
ctypedef np.int32_t DTYPEi_t
ctypedef np.float32_t DTYPEf_t
ctypedef np.float64_t DTYPEd_t

'''___NPL Modules___'''

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris", "Jon Mittaz"]
__created__ = "17/01/2018"
__credits__ = []
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def sensor_model_ref(np.ndarray[DTYPEf_t, ndim=2] sensor_model_variables,
                     sensor_model_parameter,
                     sensor_model_constant,
                     sensor_model_across_track_index,
                     sensor_model_along_track_index,
                     time):
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

    cdef float[:] R = sensor_model_variables[:,0]

    return R, None


def sensor_model_3(np.ndarray[DTYPEf_t, ndim=2] sensor_model_variables,
                   np.ndarray[DTYPEf_t, ndim=1] sensor_model_parameter,
                   np.ndarray[DTYPEd_t, ndim=1] sensor_model_constant,
                   np.ndarray[DTYPEi_t, ndim=1] sensor_model_across_track_index,
                   np.ndarray[DTYPEi_t, ndim=1] sensor_model_along_track_index,
                   time):
    """
    :type sensor_model_variables: numpy.ndarray
    :param sensor_model_variables: sensor state variables

    :type sensor_model_parameter: numpy.ndarray
    :param sensor_model_parameter: sensor calibration parameters

    :type sensor_model_constant: numpy.ndarray
    :param sensor_model_constant: sensor model constants

    :type sensor_model_across_track_index: numpy.ndarray
    :param sensor_model_across_track_index: sensor model index in across track direction (unused)

    :type sensor_model_along_track_index: numpy.ndarray
    :param sensor_model_along_track_index: sensor model index in along track direction (unused)

    :type time: numpy.ndarray
    :param time: time of observation (unused)

    :return:
        :measurand: *numpy.ndarray*

        Determined sensor measurand

        :measurand_derivatives: *numpy.ndarray*

        Measurement model derivatives w.r.t. input quantities
    """

    # Extract covariates
    cdef:
        np.ndarray[DTYPEf_t, ndim=1] Cs = sensor_model_variables[:,0]
        np.ndarray[DTYPEf_t, ndim=1] Cict = sensor_model_variables[:,1]
        np.ndarray[DTYPEf_t, ndim=1] Ce = sensor_model_variables[:,2]
        np.ndarray[DTYPEf_t, ndim=1] Rict = sensor_model_variables[:,3]
        int m = len(Cs)

    # Evaluate radiance using sensor model
    cdef:
        np.ndarray[DTYPEf_t, ndim=2] Ja = np.vstack((np.ones(m, dtype=np.float32), Rict * (Cs - Ce) / (Cs - Cict),
                                                     (Cict - Ce) * (Cs - Ce))).T
        np.ndarray[DTYPEf_t, ndim=1] at = np.array([sensor_model_parameter[0], 0.98514 + sensor_model_parameter[1],
                                                    sensor_model_parameter[2]],
                                                   dtype=np.float32)
        np.ndarray[DTYPEf_t, ndim=1] R = np.dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/da0
    # > dR/da1
    # > dR/da2

    cdef np.ndarray[DTYPEf_t, ndim=2] J

    J = np.column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2 + at[2] * (Cict - Ce),
                         at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2 + at[2] * (Cs - Ce),
                         -at[1] * Rict / (Cs - Cict) - at[2] * (Cs + Cict - 2 * Ce),
                        at[1] * (Cs - Ce) / (Cs - Cict),
                        Ja))

    return R, J


def sensor_model_4(np.ndarray[DTYPEf_t, ndim=2] sensor_model_variables,
                   np.ndarray[DTYPEf_t, ndim=1] sensor_model_parameter,
                   np.ndarray[DTYPEd_t, ndim=1] sensor_model_constant,
                   np.ndarray[DTYPEi_t, ndim=1] sensor_model_across_track_index,
                   np.ndarray[DTYPEi_t, ndim=1] sensor_model_along_track_index,
                   time):
    """
    :type sensor_model_variables: numpy.ndarray
    :param sensor_model_variables: sensor state variables

    :type sensor_model_parameter: numpy.ndarray
    :param sensor_model_parameter: sensor calibration parameters

    :type sensor_model_constant: numpy.ndarray
    :param sensor_model_constant: sensor model constants

    :type sensor_model_across_track_index: numpy.ndarray
    :param sensor_model_across_track_index: sensor model index in across track direction (unused)

    :type sensor_model_along_track_index: numpy.ndarray
    :param sensor_model_along_track_index: sensor model index in along track direction (unused)

    :type time: numpy.ndarray
    :param time: time of observation (unused)

    :return:
        :measurand: *numpy.ndarray*

        Determined sensor measurand

        :measurand_derivatives: *numpy.ndarray*

        Measurement model derivatives w.r.t. input quantities
    """

    # Extract covariates
    cdef:
        np.ndarray[DTYPEf_t, ndim=1] Cs = sensor_model_variables[:,0]
        np.ndarray[DTYPEf_t, ndim=1] Cict = sensor_model_variables[:,1]
        np.ndarray[DTYPEf_t, ndim=1] Ce = sensor_model_variables[:,2]
        np.ndarray[DTYPEf_t, ndim=1] Rict = sensor_model_variables[:,3]
        np.ndarray[DTYPEf_t, ndim=1] To = sensor_model_variables[:,4]
        int m = len(Cs)

    # Evaluate radiance using sensor model
    cdef:
        np.ndarray[DTYPEf_t, ndim=2] Ja = np.vstack((np.ones(m, dtype=np.float32), Rict * (Cs - Ce) / (Cs - Cict),
                                                     (Cict - Ce) * (Cs - Ce), To)).T
        np.ndarray[DTYPEf_t, ndim=1] at = np.array([sensor_model_parameter[0], 0.98514 + sensor_model_parameter[1],
                                                    sensor_model_parameter[2], sensor_model_parameter[3]],
                                                   dtype=np.float32)
        np.ndarray[DTYPEf_t, ndim=1] R = np.dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/dTo
    # > dR/da0
    # > dR/da1
    # > dR/da2
    # > dR/da3

    cdef np.ndarray[DTYPEf_t, ndim=2] J

    J = np.column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2 + at[2] * (Cict - Ce),
                         at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2 + at[2] * (Cs - Ce),
                         -at[1] * Rict / (Cs - Cict) - at[2] * (Cs + Cict - 2 * Ce),
                        at[1] * (Cs - Ce) / (Cs - Cict),
                        at[3] * ones(m, dtype=np.float32),
                        Ja))

    return R, J


def sensor_model_3_scaled(np.ndarray[DTYPEf_t, ndim=2] sensor_model_variables,
                          np.ndarray[DTYPEf_t, ndim=1] sensor_model_parameter,
                          np.ndarray[DTYPEd_t, ndim=1] sensor_model_constant,
                          np.ndarray[DTYPEi_t, ndim=1] sensor_model_across_track_index,
                          np.ndarray[DTYPEi_t, ndim=1] sensor_model_along_track_index,
                          time):
    """
    :type sensor_model_variables: numpy.ndarray
    :param sensor_model_variables: sensor state variables

    :type sensor_model_parameter: numpy.ndarray
    :param sensor_model_parameter: sensor calibration parameters

    :type sensor_model_constant: numpy.ndarray
    :param sensor_model_constant: sensor model constants

    :type sensor_model_across_track_index: numpy.ndarray
    :param sensor_model_across_track_index: sensor model index in across track direction (unused)

    :type sensor_model_along_track_index: numpy.ndarray
    :param sensor_model_along_track_index: sensor model index in along track direction (unused)

    :type time: numpy.ndarray
    :param time: time of observation (unused)

    :return:
        :measurand: *numpy.ndarray*

        Determined sensor measurand

        :measurand_derivatives: *numpy.ndarray*

        Measurement model derivatives w.r.t. input quantities
    """

    # Extract covariates
    cdef:
        np.ndarray[DTYPEf_t, ndim=1] Cs = sensor_model_variables[:,0]
        np.ndarray[DTYPEf_t, ndim=1] Cict = sensor_model_variables[:,1]
        np.ndarray[DTYPEf_t, ndim=1] Ce = sensor_model_variables[:,2]
        np.ndarray[DTYPEf_t, ndim=1] Rict = sensor_model_variables[:,3]
        int m = len(Cs)

    # Evaluate radiance using sensor model
    cdef:
        np.ndarray[DTYPEf_t, ndim=2] Ja = np.vstack((np.ones(m, dtype=np.float32), Rict * (Cs - Ce) / (Cs - Cict),
                                                     (Cict - Ce) * (Cs - Ce))).T
        np.ndarray[DTYPEf_t, ndim=1] at = np.array([sensor_model_parameter[0],
                                                    0.98514 + sensor_model_parameter[1]/1000.0,
                                                    sensor_model_parameter[2]/1000.0/1000.0],
                                                   dtype=np.float32)
        np.ndarray[DTYPEf_t, ndim=1] R = np.dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/da0
    # > dR/da1
    # > dR/da2

    cdef np.ndarray[DTYPEf_t, ndim=2] J

    J = np.column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2 + at[2] * (Cict - Ce),
                         at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2 + at[2] * (Cs - Ce),
                         -at[1] * Rict / (Cs - Cict) - at[2] * (Cs + Cict - 2 * Ce),
                        at[1] * (Cs - Ce) / (Cs - Cict),
                        Ja))

    return R, J


def sensor_model_4_scaled(np.ndarray[DTYPEf_t, ndim=2] sensor_model_variables,
                          np.ndarray[DTYPEf_t, ndim=1] sensor_model_parameter,
                          np.ndarray[DTYPEd_t, ndim=1] sensor_model_constant,
                          np.ndarray[DTYPEi_t, ndim=1] sensor_model_across_track_index,
                          np.ndarray[DTYPEi_t, ndim=1] sensor_model_along_track_index,
                          time):
    """
    :type sensor_model_variables: numpy.ndarray
    :param sensor_model_variables: sensor state variables

    :type sensor_model_parameter: numpy.ndarray
    :param sensor_model_parameter: sensor calibration parameters

    :type sensor_model_constant: numpy.ndarray
    :param sensor_model_constant: sensor model constants

    :type sensor_model_across_track_index: numpy.ndarray
    :param sensor_model_across_track_index: sensor model index in across track direction (unused)

    :type sensor_model_along_track_index: numpy.ndarray
    :param sensor_model_along_track_index: sensor model index in along track direction (unused)

    :type time: numpy.ndarray
    :param time: time of observation (unused)

    :return:
        :measurand: *numpy.ndarray*

        Determined sensor measurand

        :measurand_derivatives: *numpy.ndarray*

        Measurement model derivatives w.r.t. input quantities
    """

    # Extract covariates
    cdef:
        np.ndarray[DTYPEf_t, ndim=1] Cs = sensor_model_variables[:,0]
        np.ndarray[DTYPEf_t, ndim=1] Cict = sensor_model_variables[:,1]
        np.ndarray[DTYPEf_t, ndim=1] Ce = sensor_model_variables[:,2]
        np.ndarray[DTYPEf_t, ndim=1] Rict = sensor_model_variables[:,3]
        np.ndarray[DTYPEf_t, ndim=1] To = sensor_model_variables[:,4]
        int m = len(Cs)

    # Evaluate radiance using sensor model
    cdef:
        np.ndarray[DTYPEf_t, ndim=2] Ja = np.vstack((np.ones(m, dtype=np.float32), Rict * (Cs - Ce) / (Cs - Cict),
                                                     (Cict - Ce) * (Cs - Ce), To)).T
        np.ndarray[DTYPEf_t, ndim=1] at = np.array([sensor_model_parameter[0],
                                                    0.98514 + sensor_model_parameter[1]/1000.0,
                                                    sensor_model_parameter[2]/1000.0/1000.0,
                                                    sensor_model_parameter[3]/1000.0],
                                                   dtype=np.float32)
        np.ndarray[DTYPEf_t, ndim=1] R = np.dot(Ja, at)

    # Evaluate derivatives
    # In the following order:
    # > dR/dC_S
    # > dR/dC_ICT
    # > dR/dC_E
    # > dR/dR_ICT
    # > dR/dTo
    # > dR/da0
    # > dR/da1
    # > dR/da2
    # > dR/da3

    cdef np.ndarray[DTYPEf_t, ndim=2] J

    J = np.column_stack((at[1] * Rict * (Ce - Cict) / (Cs - Cict) ** 2 + at[2] * (Cict - Ce),
                         at[1] * Rict * (Cs - Ce) / (Cs - Cict) ** 2 + at[2] * (Cs - Ce),
                         -at[1] * Rict / (Cs - Cict) - at[2] * (Cs + Cict - 2 * Ce),
                        at[1] * (Cs - Ce) / (Cs - Cict),
                        at[3] * ones(m, dtype=np.float32),
                        Ja))

    return R, J