"""
Created on Mon Feb 21  2017 10:00:00

Functions specific to describe the particular instrument

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
from numpy import vstack, ones, dot, column_stack


def sensor_model(a, X):
    """
    Function to return Radiances and Derivative for input calibration parameters and covariates

    :param a: numpy.ndarray, dtype=float
        array of sensor calibration parameters
    :param X: list: numpy.ndarray
        list of covariate data values in arrays

    :return:
        :R: numpy.ndarray
            Radiances
        :J: numpy.ndarray
            Derivatives
    """

    # Extract covariates
    Cs = X[0]
    Cict = X[1]
    Ce = X[2]
    Rict = X[3]
    T = X[4]
    m = len(Cs)

    # Evaluate radiance using sensor model
    Ja = vstack((ones(m), Rict * (Cs - Ce) / (Cs - Cict), (Cict - Ce) * (Cs - Ce), T)).T
    at = [a[0], 0.98514 + a[1], a[2], a[3]]
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


def adjustment_model(R):
    """
    Return adjustment factor array and derivatives from radiance values

    :param R: numpy.ndarray
        radiance values in an array

    :return:
        :B: numpy.ndarray
            adjustment factor array
        :J: numpy.ndarray
            derivatives
    """

    # evaluate adjustment model
    B = R

    # evaluate derivatives
    J = ones(len(R))

    return B, J