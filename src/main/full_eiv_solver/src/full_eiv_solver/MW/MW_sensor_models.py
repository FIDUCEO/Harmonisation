"""
Sensor Measurement function for MW instruments aboard MHS and AMSU
"""

'''___Python Modules___'''

'''___Third Party Modules___'''
from numpy import array, zeros, cos, exp, float32
from math import pi

'''___Harmonisation Modules___'''

'''___Authorship___'''
__author__ = ["Sam Hunt"]
__created__ = "16/01/2018"
__credits__ = ["Martin Burgdorf", "Imke Hans", "Ralf Quast"]
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
dTw = 0
thetaS = pi*73.6/180.0


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

    R = sensor_model_variables[:, 0]

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

    # Sensor state variables
    CS = sensor_model_variables[:, 0]
    CICT = sensor_model_variables[:, 1]
    TICT1 = sensor_model_variables[:, 2]
    TICT2 = sensor_model_variables[:, 3]
    TICT3 = sensor_model_variables[:, 4]
    TICT4 = sensor_model_variables[:, 5]
    TICT5 = sensor_model_variables[:, 6]
    CE = sensor_model_variables[:, 7]

    # Sensor constants
    w1 = sensor_model_constant[-9]
    w2 = sensor_model_constant[-8]
    w3 = sensor_model_constant[-7]
    w4 = sensor_model_constant[-6]
    w5 = sensor_model_constant[-5]
    wtot = w1+w2+w3+w4+w5
    v = sensor_model_constant[-4]
    a = sensor_model_constant[-3]
    b = sensor_model_constant[-2]
    dTw = sensor_model_constant[-1]
    g = array([sensor_model_constant[i] for i in sensor_model_across_track_index])
    thetaE = array([pi*1.111*(i-45.5)/180.0 for i in sensor_model_across_track_index])

    # Harmonisation parameters
    a1 = sensor_model_parameter[0]
    a2 = sensor_model_parameter[1]
    a3 = sensor_model_parameter[2]

    # Compute Radiances
    TICT = (w1*TICT1 + w2*TICT2 + w3*TICT3 + w4*TICT4 + w5*TICT5)/wtot
    LICT = planck(v, a+b*(TICT+dTw))
    LS = planck(v, a+b*(TCMB+a3))
    LCMB = planck(v, a+b*TCMB)

    # Compute measurand
    LMEprime = LICT + (LICT-LS)/(CICT-CS)*(CE-CICT) + a1*(LICT-LCMB)**2/(CICT-CS)**2*(CE-CICT)*(CE-CS)
    cosangs = cos(2*thetaE)-cos(2*thetaS)
    dLpol = 0.5*a2*cosangs*(LICT-LMEprime)
    LME = LMEprime + dLpol

    L = 1.0/g*(LME-(1.0-g)*LCMB)

    # Compute measurand derivatives
    J = zeros((len(L), 12), dtype=float32)

    # dL/dCS
    J[:, 0] = 1.0/g * (CE-CICT)*(LICT-LS)/(CICT-CS)**2 * (1 + a1*(2*CE-CICT-CS)/(CICT-CS))

    # dL/dCICT
    J[:, 1] = 1.0/g * (CE-CS)*(LICT-LS)/(CICT-CS)**2 * (1 + a1*(LICT-LS)*(2*CE-CICT-CS)/(CICT-CS))

    # dL/dTICTn
    dLdTICT = 1.0/g * (1 + (CE-CICT)/(CICT-CS) + 2*a1*(CE-CS)*(CE-CICT)*(LICT-LS)/(CICT-CS)**2) * \
                 (LICT * (2*H*v*KB*b)/(KB*(a+b*(TICT+dTw)))**2 *
                      exp(H*v/(KB*(a+b*(TICT+dTw)))) / (exp(H*v/(KB*(a+b*(TICT+dTw))))-1))
    J[:, 2] = dLdTICT * w1 / wtot
    J[:, 3] = dLdTICT * w2 / wtot
    J[:, 4] = dLdTICT * w3 / wtot
    J[:, 5] = dLdTICT * w4 / wtot
    J[:, 6] = dLdTICT * w5 / wtot

    # dL/dCE
    J[:, 7] = 1.0/g * ((LICT-LS)/(CICT-CS)+a1*(LICT-LS)**2/(CICT-CS)**2)*(2*CE-CICT-CS)

    # dL/da1
    J[:, -3] = 1.0/g * (CE-CS)*(CE-CICT)*(LICT-LS)**2/(CICT-CS)**2

    # dL/da2
    J[:, -2] = 1.0/(2.0*g) * (LICT-LMEprime) * cosangs

    # dL/da3
    J[:, -1] = -1.0/g * (CE-CICT)/(CICT-CS) * (1.0 + 2.0*a1*(CE-CS)*(LICT-LS)/(CICT-CS)) * (1.0 - a2/(2.0*g)*cosangs)

    # Convert to BT (via RJ)
    RJ_const = C**2 / (2.0 * v**2 * KB)
    L = RJ_const * L
    J = RJ_const * J

    return L, J


def planck(v, T):
    return 2.0 * H * v**3 / C**2 / (exp(v*H/(KB*T))-1)
