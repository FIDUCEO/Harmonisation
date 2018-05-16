"""
Implements the algorithm to calculate the pre-conditioner solution of the full algorithm to harmonise match-up data
"""

'''___Python Modules____'''
import os.path
from math import ceil
from copy import deepcopy

'''___Third Party Modules___'''
from numpy import loadtxt, append, eye, dot, zeros, ones, outer, diag, set_printoptions, nan, float32, asarray, count_nonzero, arange, int32
from numpy.linalg import norm, cholesky, solve, lstsq
from scipy.linalg import qr
from scipy.optimize import least_squares
set_printoptions(threshold=nan)

'''___Harmonisation Modules___'''
from GN_algo import GNAlgo

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class PCAlgo:
    """
    Class containing functionality to calculate the pre-conditioner solution for the full EIV algorithm to harmonise
    match-up data

    Sample Code:

    .. code-block:: python

        PC = PCAlgo(HData)
        PC.runGN()

    :Attributes:
        .. py:attribute:: calc_R

        *func*

        Function to calculate R

        .. py:attribute:: HData

        *harm_data_reader.HarmInputData*

        Harmonisation data (should be sampled by convert_data.sample4PC)

        .. py:attribute:: xyza

        *numpy.ndarray*

        Cost function

        .. py:attribute:: LK

        *numpy.ndarray*

        Cholesky factor of covariance matrix for measurements K

    :Methods:
        .. py:method:: runPC(...):

            Run algorithm to calculate the full harmonisation algorithm pre-conditioner solution

        .. py:method:: calc_LK(...):

            Evaluate Cholesky factor of covariance matrix for measurements K to include contributions from other
            measured variables (considered to be fixed)

        .. py:method:: get_f(...):

            Get weighted Jacobian matrix for measurements of K with respect to calibration parameters

        .. py:method:: get_Ja(...):

            Get weighted Jacobian matrix for measurements of K with respect to calibration parameters

        .. py:method:: calc_fJa(...):

            Evaluate weighted residual vector and weighted Jacobian matrix for measurements of K with respect to
            calibration parameters
    """

    def __init__(self, HData=None):
        """
        Initialise class

        :type HData: harm_data_reader.HarmInputData
        :param HData: Harmonisation data (should be sampled by convert_data.sample4PC)
        """

        # Initialise class
        self.HData = None
        self.xyza = None
        self.calc_R = None
        self.LK = None

        if HData is not None:
            self.HData = HData
            self.xyza = append(HData.values[:], HData.a[:])

            GN = GNAlgo()
            self.calc_R = GN.calc_R

    def runPC(self, tol=1e-6):
        """
        Run algorithm to calculate the full harmonisation algorithm pre-conditioner solution

        :type tol: float
        :param tol: tolerance of convergence

        :return:
            :a: *numpy.ndarray*

            Pre-conditioner estimated calibration parameters

            :S: *numpy.ndarray*

            Pre-conditioner solution
        """

        # The pre-condition solution to the harmonisation runs a simplified version of the full problem. All data are
        # considered variables in the full approach, in this case however we only vary the satellite sensor calibration
        # parameters to find the solution.

        HData = self.HData

        # initialise parameters
        N_var = HData.idx['idx'][-1]                                               # total number of variables
        N_a = len(HData.a)                                                         # total number of parameters for
                                                                                   # all sensors combined
        ncol = len(HData.values) + len(HData.a)                                    # number of columns
        niter = 0                                                                  # iteration counter
        mxiter = 20 + ceil(ncol**0.5)                                              # maximum allowed iterations
        conv = False                                                               # convergence boolean
        a0 = deepcopy(HData.a)                                                     # initial parameter estimates

        ################################################################################################################
        # 1. Determine parameter estimates
        ################################################################################################################

        while (conv is False) and (niter < mxiter):

            # Determine new estimates of calibration parameters a by solving nonlinear least-squares problem defined
            # by the data Kdatac with covariance matrix VK

            self.LK = self.calc_LK(self.HData)

            a = least_squares(self.get_f, a0, self.get_Ja, ftol=tol, xtol=tol, verbose=1)['x']

            # Check on convergence
            norma = norm(a - a0)
            if norma < tol:
                conv = True

            # Update current approximation and number of iterations
            a0 = a[:]

            niter += 1

        print 'Determined a:'
        print a

        ################################################################################################################
        # 2. Evaluate uncertainty
        ################################################################################################################

        J = self.get_Ja(a.flatten('F'))
        R = qr(J)[1][:N_a, :N_a]
        U = lstsq(R, eye(N_a, dtype=float32))[0]
        V = dot(U, U.T)

        ################################################################################################################
        # 3. Evaluate pre-conditioner
        ################################################################################################################

        S = cholesky(V)

        return a.astype(dtype=float32), S

    def calc_LK(self, HData):
        """
        Evaluate Cholesky factor of covariance matrix for measurements K to include contributions from other measured
        variables (considered to be fixed)

        :type HData: harm_data_reader.HarmInputData
        :param HData: Harmonisation data (should be sampled by convert_data.sample4PC)

        :return:
            :LK: *numpy.ndarray*

            Cholesky factor of covariance matrix for measurements K
        """

        # initialise parameters
        mc = HData.idx['cNm']                                                       # cumulative
        N_mu = HData.idx['cNm'][-1]                                                 # total match-ups (= number of ks)
        N_var = HData.idx['idx'][-1]                                                # total variables
        N_sensors = len(set(self.HData.idx['parameter_sensor']))                    # total number of sensors
        N_a = len(self.HData.a)                                                     # total number of parameters for
                                                                                    # all sensors combined
        N_cov = HData.idx['n_cov'][-1]                                              # total number of covariates

        # initialise arrays
        JT = zeros((N_mu, N_sensors), dtype=float32)
        VK = eye(N_mu, dtype=float32)

        # loop through match-up series
        for i, n_sensors in enumerate(HData.idx['Im']):

            n_mu = i + 1           # match-up series number
            uK = HData.unck[i].uR  # k uncertainty for match-up series

            # indices for data
            istart = mc[n_mu - 1]
            iend = mc[n_mu]

            # get radiances and adjustment model values for sensor 1 and sensor 2 of match-up series and build
            # Jacobian
            for j, n_sensor in enumerate(n_sensors):

                sensor_name = HData.idx['sensors'][n_sensor]

                R, JR = self.calc_R(self.xyza, HData.unc, HData.w_matrices, HData.u_matrices, HData.idx, HData.idx,
                                    HData.sensor_model,
                                    HData.sensor_model_constant,
                                    HData.getSensorAcrossTrackIndex(n_mu, sensor_name),
                                    HData.getSensorAlongTrackIndex(n_mu, sensor_name),
                                    HData.getSensorTime(n_mu, sensor_name),
                                    n_sensor, n_mu)

                # evaluate adjustment model for sensor radiances
                JB = HData.adjustment_model[n_sensor](R)[1]

                # first sensor or second sensor factor
                s = -1
                if j == 1:
                    s = 1

                # > if reference sensor
                if JR is None:

                    # find uncertainty data
                    indices = [(j, k) for j, k in zip(HData.idx['n_sensor'], HData.idx['n_mu'])]
                    im = indices.index((n_sensor, n_mu))

                    # undo conversion of data from ConvertData.convert4GN and add to radiances list
                    block_unc = HData.unc[im]

                    JK = s * diag(JB) * block_unc.uR / uK
                    VK[istart:iend, istart:iend] = VK[istart:iend, istart:iend] + dot(JK, JK.T)

                # > if sensor
                else:
                    # build covariate by covariate, in a way depending on correlation form
                    # NB: data should be sample in a way that there are no correlations from averaging
                    for n_cov in xrange(1, N_cov + 1):

                        # find block uncertainty data
                        indices = [(i1, i2, i3) for i1, i2, i3 in zip(HData.idx['n_sensor'],
                                                                      HData.idx['n_mu'],
                                                                      HData.idx['n_cov'])]
                        im = indices.index((n_sensor, n_mu, n_cov))
                        block_unc = HData.unc[im]

                        # a. random correlation
                        if block_unc.typeID == 1:
                            JK = s * diag(JB*JR[:, n_cov-1]) * block_unc.uR / uK

                        # b. random+systematic correlation
                        if block_unc.typeID == 2:
                            JK = s * diag(JB*JR[:, n_cov-1]) * block_unc.uR / uK
                            JT[istart:iend, n_sensor-1] = s * JB*JR[:, n_cov-1] * block_unc.uS / uK

                        VK[istart:iend, istart:iend] = VK[istart:iend, istart:iend] + dot(JK, JK.T)

        VK = VK + dot(JT, JT.T)

        # form Cholesky decomposition
        LK = cholesky(VK)

        return LK

    def get_f(self, a, *args):
        """
        Get weighted Jacobian matrix for measurements of K with respect to calibration parameters

        :type a: numpy.ndarray
        :param a: calibration parameter estimates

        :return:
            :f: *numpy.ndarray*

            weighted residual vector of estimates of k
        """

        f = self.calc_fJa(a)

        return f

    def get_Ja(self, a):
        """
        Get weighted Jacobian matrix for measurements of K with respect to calibration parameters

        :type a: numpy.ndarray
        :param a: Calibration parameter estimates

        :return:
            :Ja: *numpy.ndarray*

            Weighted Jacobian matrix

        """

        Ja = self.calc_fJa(a, ret="Ja")

        return Ja

    def calc_fJa(self, a, ret="f"):
        """
        Evaluate weighted residual vector and weighted Jacobian matrix for measurements of K with respect to
        calibration parameters

        :type a: numpy.ndarray
        :param a: Calibration parameter estimates

        :type ret: str
        :param ret: choose returned parameter. If:
            * ret = "f" (default) - return f
            * ret = "Ja" - return Ja

        :return:
            :f: *numpy.ndarray*

            Weighted residual vector of estimates of k

            *or, depending on ret,*

            :Ja: *numpy.ndarray*

            Weighted Jacobian matrix
        """

        HData = self.HData

        # initialise parameters
        mc = HData.idx['cNm']                                                       # cumulative
        N_mu = HData.idx['cNm'][-1]                                                 # total match-ups (= number of ks)
        N_var = HData.idx['idx'][-1]                                                # total variables
        N_a = len(self.HData.a)                                                     # total number of parameters for
                                                                                    # all sensors combined
        N_cov = HData.idx['n_cov'][-1]                                              # total number of covariates

        # update a estimate in xyza
        self.xyza[N_var:N_var + N_a] = a

        # initialise arrays
        Ja = zeros((N_mu, N_a), dtype=float32)
        k_est = zeros(N_mu, dtype=float32)  # initialise k_est for estimates of Ks

        # loop through match-up series
        for i, n_sensors in enumerate(HData.idx['Im']):

            n_mu = i + 1           # match-up series number
            uK = HData.unck[i].uR  # k uncertainty for match-up series

            # indices for data
            istart = mc[n_mu - 1]
            iend = mc[n_mu]

            # get radiances and adjustment model values for sensor 1 and sensor 2 of match-up series and build
            # Jacobian
            Rs = []
            Bs = []
            for j, n_sensor in enumerate(n_sensors):

                sensor_name = HData.idx['sensors'][n_sensor]

                R, JR = self.calc_R(self.xyza, HData.unc, HData.w_matrices, HData.u_matrices, HData.idx, HData.idx,
                                    HData.sensor_model,
                                    HData.sensor_model_constant,
                                    HData.getSensorAcrossTrackIndex(n_mu, sensor_name),
                                    HData.getSensorAlongTrackIndex(n_mu, sensor_name),
                                    HData.getSensorTime(n_mu, sensor_name),
                                    n_sensor, n_mu)
                Rs.append(R)

                # evaluate adjustment model for sensor radiances
                B, JB = HData.adjustment_model[n_sensor](R)
                Bs.append(B)  # Evaluate adjustment model

                # first sensor or second sensor factor
                s = -1
                if j == 1:
                    s = 1

                # build Jacobian matrix with respect to calibration parameters
                if JR is not None:
                    i_sensor_p = asarray([i for i, sn in enumerate(HData.idx["parameter_sensor"]) if sn == sensor_name], dtype=int32)
                    N_sensor_p = len(i_sensor_p)
                    Ja[istart:iend, i_sensor_p] = s * outer(JB, ones(N_sensor_p)) * JR[:, N_cov:N_cov+N_sensor_p] / uK[:, None]
            k_est[istart:iend] = (Bs[1] - Bs[0]) / uK  # Evaluate estimate of adjustment factor

        # Evaluate and return required parameter
        if ret == "f":
            # Evaluate weighted residual vector
            f = solve(self.LK, (k_est - HData.ks))
            return f

        if ret == "Ja":
            # Evaluate weighted Jacobian matrix
            Ja = solve(self.LK, Ja)
            return Ja


if __name__ == "__main__":
    pass
