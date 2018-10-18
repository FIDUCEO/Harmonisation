"""
Implements the Gauss-Newton iteration algorithm for match-up data with a pre-conditioned solution.
"""

'''___Python Modules____'''
import numpy
import pyximport
pyximport.install(setup_args={"include_dirs": numpy.get_include()},
                  reload_support=True)

import sys
from os import makedirs
from os.path import dirname
from os.path import join as pjoin
from copy import deepcopy
from multiprocessing import Pool
import cPickle
from time import time

'''___Third Party Modules____'''
from numpy import zeros, append, ones, dot, outer, hstack, array, eye, inf, asarray, int32, float32, float64, column_stack, count_nonzero, where
from numpy.linalg import norm, solve, cholesky
from scipy.sparse.linalg import LinearOperator
from math import ceil
from netCDF4 import Dataset
from pykrylov_lsmr import LSMRFramework
from pykrylov_minres import Minres

'''___NPL Modules___'''
from cython_lib_GN_algo import unconvert_Xs

matchupIO_directory = pjoin(dirname(dirname(dirname(dirname(dirname(__file__))))), "matchupIO")
sys.path.append(matchupIO_directory)
from MatchUp import MatchUp

matchupProcessing_directory = pjoin(dirname(dirname(dirname(dirname(dirname(__file__))))), "matchupProcessing")
sys.path.append(matchupProcessing_directory)
from transform2normind.Transform2NormInd import Transform2NormInd

matchupToolbox_directory = pjoin(dirname(dirname(dirname(dirname(dirname(__file__))))), "matchupToolbox")
sys.path.append(matchupToolbox_directory)
from utils.matchup_maths import evaluate_measurand, evaluate_adjusted_measurand, evaluate_K

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Arta Dillo", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


class GNAlgo(object):
    """
    This class implements the Gauss-Newton iteration algorithm of Peter Harris to harmonise satellite instrument
    calibration parameters for group of sensors with a reference sensor from match-up data, HData, and a
    pre-conditioned solution S.

    Sample Code:

    .. code-block::python

        GN = GNAlgo(HData, S)
        GN.runGN()

    :Attributes:
        .. py:attribute:: HData

        *harm_data_reader.HarmInputData*

        Input data to be harmonised

        .. py:attribute:: S

        *numpy.ndarray*

        Pre-conditioner solution

        .. py:attribute:: xyza

        *numpy.ndarray*

        Array containing variables and parameters

    :Methods:
        .. py:method:: runGN(...):

            Run Gauss-Newton Algorithm to perform harmonisation

        .. py:method:: calc_f(...):

            Return value for f, array containing the residual between the the current and original estimates of
            radiances, variables, and ks

        .. py:method:: get_JPx(...):

            Return array containing the product of JP and x for a given x

        .. py:method:: get_JPTx(...):

            Return array containing the product of JP transpose and x for a given x

        .. py:method:: calc_R(...):

            Return calculated radiance R and derivatives for given input data

        .. py:method:: unconvert_Xs(...):

            Return variable data for each covariate in the original form for a given sensor and match-up, undoing the
            reparameterisation performed in ConvertData.convert2ind()

        .. py:method:: calc_prod_JPx(...):

            Return the product of JP (or JP transpose) for a given x

        .. py:method:: calc_Px(...):

            Return value of x multiplied by preconditioner solution P (or transpose)

        .. py:method:: calc_unc(...):

            Return array containing the covariance matrices for the retrieved parameters for each sensor, derived
            analytically

        .. py:method:: calc_Hx(...):

            Return the product of H (P'*J'*J*P) with a given x
    """

    def __init__(self, HData=None, V_initial=None):
        """
        Initialise algorithm

        :type HData: HarmInputData
        :param HData: Input data to be harmonised

        :type V_initial: numpy.ndarray
        :param V_initial: Parameter covariance matrix estimation (from pre-conditioner)
        """

        # Initialise class
        self.HData = None
        self.S = None
        self.xyza = None
        self.Transform2NormIndOp = Transform2NormInd()
        if HData is not None:
            self.HData = HData

            # initialise current variable and parameter estimates
            self.xyza = append(self.HData.values, self.HData.a)
            pass

        if (HData is not None) and (V_initial is not None):
            self.S = cholesky(V_initial)
        elif (HData is not None) and (V_initial is None):
            self.S = eye(len(self.HData.a))

    def run(self, tol=1e-6, tolA=1e-8, tolB=1e8, tolU=1e-8, show=False, return_covariance=True):
        """
        Run Gauss-Newton Algorithm to perform harmonisation

        :type tol: float
        :param tol: Tolerance for convergance of GN algorithm

        :type tolA: float
        :param tolA: tolerance tolA for LSMR in GN algorithm

        :type tolB: float
        :param tolB: tolerance tolB for LSMR in GN algorithm

        :type tolU: float
        :param tolU: tolerance for uncertainty calculation convergence (rtol in Minres)

        :type show: bool
        :param show: boolean to decide if stdout output of algorithm

        :return:
            :a: *numpy.ndarray*

            Estimate of parameters

            :Va: *numpy.ndarray*

            Covariance matrix for the parameter estimates
        """

        if show:
            print "Initial Parameter Estimates:"
            print self.HData.a
            print "Determining Parameters..."

        # Useful parameters
        N_mu = self.HData.idx['cNm'][-1]              # total match-ups
        N_var = self.HData.idx['idx'][-1]             # total number of variables
        N_a = len(self.HData.a)                       # number calibration parameters

        ################################################################################################################
        # 1. Gauss Newton Solver
        ################################################################################################################

        # a. Preparation -----------------------------------------------------------------------------------------------
        # i. Iteration parameters
        niter = 0                                                                     # counter of iterations
        mxiter = ceil(N_var)                                                          # max number of iterations of GN
        mxiter_lsmr = ceil(N_var)                                                     # max number of iterations of LSMR
        conv = False                                                                  # convergence boolean
        GNlog = []

        # ii. Initialise Operators
        # - J LinearOperator
        J = LinearOperator((N_var+N_mu, N_var+N_a), matvec=self.get_JPx, rmatvec=self.get_JPTx)
        # - LSMR Operators
        LSMROp = LSMRFramework(J)

        # Calculate initial cost
        residuals = self.calc_f(self.xyza, self.HData)
        cost_previous = norm(residuals)**2
        # --------------------------------------------------------------------------------------------------------------

        # b. Gauss Newton Iterations -----------------------------------------------------------------------------------
        while (conv is False) and (niter < mxiter):
            # i. GN Step ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            niter += 1  # Update iteration counter

            # Determine Gauss-Newton step d as solution to linear least-squares
            # problem J*d = -f with the pre-conditioner applied.
            LSMROp.solve(-residuals, damp=0, atol=tolA, btol=tolA, conlim=tolB, itnlim=mxiter_lsmr, show=show)
            d = self.calc_Px(LSMROp.x)

            # Update parameter estimates, as well as residuals, cost and gradient
            self.xyza += d
            residuals = self.calc_f(self.xyza, self.HData)
            cost = norm(residuals)**2
            gradient = 2 * self.get_JPTx(residuals)
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # ii. Test convergence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            cost_reduction = cost_previous - cost
            cost_reduction_tol = tol*(1+cost)
            norm_d = norm(d, inf)
            norm_d_tol = (tol**0.5) * (1 + norm(self.xyza, inf))
            norm_g = norm(gradient, inf)
            norm_g_tol = (tol**(1./3.))*(1+cost)

            # Store cost at this iteration
            cost_previous = cost

            # Check for convergence
            if (cost_reduction >= 0) and (cost_reduction < cost_reduction_tol)\
                    and (norm_d < norm_d_tol) and (norm_g <= norm_g_tol):
                conv = True

            # Write log
            GNlog.append([niter, cost_reduction, cost_reduction_tol, norm_d, norm_d_tol, norm_g, norm_g_tol])
            if show:
                print "\n\t\t\t\tGNlog"
                print "niter\tU1\t\ttol1\t\tU2\t\ttol2\t\tU3\t\ttol3"
                for GN in GNlog:
                    print "{0:2d}\t{1:.2e}\t{2:.2e}\t{3:.2e}\t{4:.2e}\t{5:.2e}\t{6:.2e}"\
                          .format(GN[0], GN[1], GN[2], GN[3], GN[4], GN[5], GN[6])
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # --------------------------------------------------------------------------------------------------------------

        # Unpack solution
        a = self.xyza[N_var:]

        if show:
            print "Determined Parameter Estimates:"
            print a

        ################################################################################################################
        # 2. Uncertainty evaluation
        ################################################################################################################

        # Uncertainty evaluation
        parameter_covariance_matrix = zeros((len(a), len(a)))
        if return_covariance:
            if show:
                print 'Determining uncertainty...'

            parameter_covariance_matrix = self.calculate_parameter_covariance_matrix(tolU=tolU, show=show)

            if show:
                print "Determined Parameter Covariance Matrix:"
                print parameter_covariance_matrix

        ################################################################################################################
        # 3. Prepare Solution
        ################################################################################################################

        if show:
            print 'Preparing output...'

        MatchUpRes = MatchUp()
        MatchUpRes.idx = self.HData.idx
        MatchUpRes._original_idx = self.HData._original_idx
        MatchUpRes.unc = self.HData.unc
        MatchUpRes.unck = self.HData.unck
        MatchUpRes.w_matrices = self.HData.w_matrices
        MatchUpRes.u_matrices = self.HData.u_matrices
        MatchUpRes.values = residuals[:N_var]
        MatchUpRes.ks = residuals[N_var:]
        MatchUpRes = self.Transform2NormIndOp.reverse(MatchUpRes)

        cost_dof = N_var - N_mu - N_a
        cost_p_value = 0

        # Return fitted systematic errors
        n_uS = max([0] + [unc_i.uS_i for unc_i in self.HData.unc if (unc_i.typeID == 2) or (unc_i.typeID == 4)])
        systematic_errors = None
        systematic_error_sensors = None
        if n_uS != 0:
            systematic_errors = self.xyza[N_var - n_uS:N_var]

            n_uSs = [unc_i.uS_i if (unc_i.typeID == 2) or (unc_i.typeID == 4) else 0 for unc_i in self.HData.unc]
            systematic_error_sensors = [self.HData.idx['sensors'][self.HData.idx['n_sensor'][n_uSs.index(i)]]
                                               for i in range(1, n_uS+1)]

        return a, parameter_covariance_matrix, cost, cost_dof, cost_p_value, MatchUpRes.values, MatchUpRes.ks, \
               systematic_errors, systematic_error_sensors

    def calc_f(self, xyza, HData):
        """
        Return value for f, array containing the residual between the the current and original estimates of radiances,
        variables, and ks

        :type xyza: numpy.ndarray
        :param xyza: array containing the current estimates of variables and parameters

        :return:
            :f: *numpy.ndarray*

            array containing the differences between original values of *R*s, *X*s, and *K*s with the current estimates,
            structured as:
                    *f = [ dR | dX | dK ]*
        """

        # initialise parameters
        mc = HData.idx['cNm']                                                    # cumulative
        N_mu = HData.idx['cNm'][-1]                                              # total match-ups (= number of ks)
        N_var = HData.idx['idx'][-1]                                             # total variables

        # initialise f (length number of variables + number of ks)
        f = zeros(N_var + N_mu, dtype=float32)

        ################################################################################################################
        # 1. Calculate f for values
        ################################################################################################################

        f[0:N_var] = xyza[0:N_var] - HData.values[0:N_var]

        ################################################################################################################
        # 2. Calculate f for ks
        ################################################################################################################

        # Estimate of k, k_est, determined as,
        #
        # k_est =  B(R_2) - B(R_1),
        #
        # where:
        # - R_1/2 - Radiances from sensor 1 and sensor 2 respectively
        # - B - adjustment model

        MatchUpEstimate_NormInd = MatchUp()
        MatchUpEstimate_NormInd.idx = self.HData.idx
        MatchUpEstimate_NormInd._original_idx = self.HData._original_idx
        MatchUpEstimate_NormInd.a = deepcopy(xyza[N_var:])
        MatchUpEstimate_NormInd.unc = self.HData.unc
        MatchUpEstimate_NormInd.unck = self.HData.unck
        MatchUpEstimate_NormInd.w_matrices = self.HData.w_matrices
        MatchUpEstimate_NormInd.u_matrices = self.HData.u_matrices
        MatchUpEstimate_NormInd.sensor_model = self.HData.sensor_model
        MatchUpEstimate_NormInd.sensor_model_constant = self.HData.sensor_model_constant
        MatchUpEstimate_NormInd.adjustment_model = self.HData.adjustment_model
        MatchUpEstimate_NormInd.values = deepcopy(xyza[:N_var])
        MatchUpEstimate_NormInd.ks = zeros(N_mu, dtype=float32)

        MatchUpEstimate = self.Transform2NormIndOp.reverse(MatchUpEstimate_NormInd)
        MatchUpEstimate.ks = evaluate_K(MatchUpEstimate)

        # Normalise by uncertainty
        for i in xrange(len(MatchUpEstimate.idx['Im'])):
            istart = MatchUpEstimate.idx['cNm'][i]
            iend = MatchUpEstimate.idx['cNm'][i + 1]

            MatchUpEstimate.ks[istart:iend] /= MatchUpEstimate.unck[i].uR

        f[N_var:] = MatchUpEstimate.ks - HData.ks

        return f

    def get_JPx(self, x):
        """
        Return array containing the product of JP and x for a given x

        :type x: numpy.ndarray
        :param x: Input array to be multiplied by JP

        :return:
            :JPTx: *numpy.ndarray*

            Array containing the product of JP and x
        """

        # Call to function to calculate product
        JPx = self.calc_prod_JPx(x, transpose=False)

        return JPx

    def get_JPTx(self, x):
        """
        Return array containing the product of JP transpose and x for a given x

        :type x: numpy.ndarray
        :param x: Input array to be multiplied by the transpose of JP

        :return:
            :JPTx: *numpy.ndarray*

            Array containing the product of JP transpose and x
        """

        # Call to function to calculate product
        JPTx = self.calc_prod_JPx(x, transpose=True)

        return JPTx

    def calc_R(self, xyza, unc, w_matrices, u_matrices, idx, original_idx, sensor_model,
                     sensor_model_constant, sensor_across_track_index, sensor_along_track_index, sensor_time,
                     n_sensor, n_mu):
        """
        Return calculated radiance R and derivatives for given input data

        :type xyza: numpy.ndarray
        :param xyza: Array containing variables and parameters

        :type unc: numpy.ndarray
        :param unc: array containing uncertainty information for variables

        :type w_matrices: list
        :param w_matrices: list of dataset w matrices (in sparse format)

        :type u_matrices: list
        :param unc: list of dataset u matrices (diagonals)

        :type idx: dict
        :param idx: dictionary with entries describing the data structures

        :type original_idx: dict
        :param original_idx: dictionary with entries describing the data structures

        :type sensor_across_track_index: numpy.ndarray
        :param sensor_across_track_index: sensor across track indices per match-up

        :type sensor_along_track_index: numpy.ndarray
        :param sensor_along_track_index: sensor along track indices per match-up

        :type sensor_time: numpy.ndarray
        :param sensor_time: sensor time per match-up

        :type sensor_model: func
        :param sensor_model: function to calculate radiance from input variables and parameters

        :type n_sensor: int
        :param n_sensor: number of sensor to calculate radiances for

        :type n_mu: int
        :param n_mu: number of match-up series to calculate radiances for

        :return:
            :measurand: *numpy.ndarray*

            Sensor measurand

            :measaurand derivatives: *numpy.ndarray*

            Sensor measurand derivatives
        """

        # self.HData.sensor_model_constant

        sensor = idx['sensors'][n_sensor]

        # 1. Find sensor model constants
        sensor_constants = asarray([const for const, const_sensor in
                            zip(sensor_model_constant, idx['sensor_model_constant_sensor']) if const_sensor == sensor], dtype=float64)

        # 2. Find sensor model parameter estimates
        parameter_estimates = xyza[-len(idx['parameter_sensor']):]
        sensor_parameter_estimates = asarray([p for p, p_sensor in zip(parameter_estimates, idx['parameter_sensor'])
                                            if p_sensor == sensor], dtype=float32)

        # 3. Get sensor state data
        sensor_model_variables = self.unconvert_Xs(xyza, unc, idx, original_idx, n_sensor, n_mu)

        # sensor_model_variables = unconvert_Xs(xyza,
        #                                       asarray(idx["n_sensor"], dtype=int32),
        #                                       asarray(idx["idx"], dtype=int32),
        #                                       asarray(original_idx["idx"], dtype=int32),
        #                                       asarray(idx["n_cov"], dtype=int32),
        #                                       asarray(idx["n_mu"], dtype=int32),
        #                                       asarray(idx["N_var"], dtype=int32),
        #                                       unc,
        #                                       w_matrices,
        #                                       u_matrices,
        #                                       n_sensor,
        #                                       n_mu,
        #                                       idx["sensor_ms"][n_sensor])

        # 3. Compute sensor measurand
        measurand, measurand_derivatives = sensor_model[n_sensor](sensor_model_variables,
                                                                  sensor_parameter_estimates,
                                                                  sensor_constants,
                                                                  sensor_across_track_index,
                                                                  sensor_along_track_index,
                                                                  sensor_time)

        return measurand, measurand_derivatives

    def unconvert_Xs(self, xyza, unc, idx, original_idx, n_sensor, n_mu):
        """
        Return variable data for each covariate in the original form for a given sensor and match-up, undoing the
        reparameterisation performed in ConvertData.convert2ind()

        :type xyza: numpy.ndarray
        :param xyza: Converted variable and parameter data

        :type unc: list
        :param unc: Uncertainties associated with blocks of variable data

        :type idx: dict
        :param idx: Dictionary describing structure of variable data

        :type n_sensor: int
        :param n_sensor: number of sensor

        :type n_mu: int
        :param n_mu: Number of match-up

        :return:
            :X: *list:numpy.darray*

            Unconverted sensor state data
        """

        N_cov = idx["sensor_ms"][n_sensor]  # total number of covariates
        mcxyz = idx['idx']  # cumulative variables by block
        n_uS = max([0]+[unc_i.uS_i for unc_i in unc if (unc_i.typeID == 2) or (unc_i.typeID == 4)])

        Xs = zeros((idx['Nm'][n_mu-1], N_cov), dtype=float32)  # initialise list for covariate data

        # get covariate data covariate by covariate
        for n_cov in xrange(1, N_cov + 1):

            # find location of data in xyza
            indices = [(i1, i2, i3) for i1, i2, i3 in zip(idx['n_sensor'], idx['n_mu'], idx['n_cov'])]
            im = indices.index((n_sensor, n_mu, n_cov))
            ib = mcxyz[im]
            ie = ib + int(idx['N_var'][im])

            # undo conversion of data from ConvertData.convert4GN depending on correlation form
            # and add to radiances list
            block_unc = unc[im]

            # a. random correlation - unscale and add to covariate list
            if block_unc.typeID == 1:
                Xs[:, n_cov-1] = xyza[ib:ie] * block_unc.uR

            # b. random+systematic correlation - unscale components and recombine
            if block_unc.typeID == 2:
                # get index of required systematic value
                isys = idx['idx'][-1] - n_uS - 1 + block_unc.uS_i
                Xs[:, n_cov-1] = xyza[ib:ie]*block_unc.uR + xyza[isys]*block_unc.uS

            # c. structured correlation - transform independent counts to counts
            if block_unc.typeID == 3:

                # Retrieve required W and U matrices
                w = self.HData.w_matrices[block_unc.w_i]
                u = self.HData.u_matrices[block_unc.u_i]

                Xs[:, n_cov-1] = w.dot(u*xyza[ib:ie])

            # d. structured+systematic correlation - add systematic, unscale and transform independent counts to counts
            if block_unc.typeID == 4:
                # Retrieve required W and U matrices
                isys = self.HData.idx['idx'][-1] - n_uS - 1 + block_unc.uS_i
                w = self.HData.w_matrices[block_unc.w_i]
                u = self.HData.u_matrices[block_unc.u_i]

                Xs[:, n_cov-1] = w.dot(u * xyza[ib:ie] + xyza[isys] * block_unc.uS)

        return Xs

    def calc_prod_JPx(self, x, transpose=False):
        """
        Return the product of JP (or JP transpose) for a given x

        :globals:
            :self.HData: HarmInputData
                Harmonisation data object

        :type x: numpy.ndarray
        :param x: vector to multiply by JP (or JP transpose)

        :type transpose: bool
        :param transpose: Boolean to decide whether to multiply x by JP or (JP)T

        :return:
            :JPx: *numpy.ndarray*

            Array containing the product of JP (or JP transpose) with x
        """

        # Jacobian structured as,
        #
        #        d/d  Rs  X1.....XN   a
        #        Rs |               |    |
        #        X_1|       I       | 0  |
        #        ...|               |    |
        #  J =   X_N|               |    |
        #           |---------------+----|
        #        Ks |       U       | Ja |
        #           |               |    |
        #
        # with x vector as,
        #
        #     |x1|
        # x = |--|
        #     |x2|
        #
        # So that their product is structured as,
        #
        #      | I  | 0  | |x1|   |     x1      |
        # Jx = |----+----| |--| = |-------------|
        #      | U  | Ja | |x2|   |U*x1 + Ja*x2 |
        #
        # Can use this structure to in algorithm to speed up calculating product
        #
        # Algorithm structured into sections:
        # 1. First apply preconditioner to x
        # 2. Store x1 as first part of product array
        # 3. Calculate and store U*x1 + Ja*x2 in product array
        #
        # If multiplying by the transpose of J have,
        #
        #          | I  | UT | |x1|   | x1 + UT*x2  |
        # (J)T x = |----+----| |--| = |-------------|
        #          | 0  |JaT | |x2|   |U*x1 + Ja*x2 |
        #
        # In this case algorithm structured as:
        # 1. Store x1 as first part of product array
        # 2. Calculate other terms of product array
        # 3. Apply transpose of preconditioner
        #
        # This algorithm calculates JPx or (JP)T x depending on parameter transpose boolean

        # initialise parameters
        mc = self.HData.idx['cNm']                                                    # cumulative
        N_mu = self.HData.idx['cNm'][-1]                                              # total match-ups (= number of ks)
        N_a = len(self.HData.a)                                                       # total number of parameters for
                                                                                      # all sensors combined
        mcxyz = self.HData.idx['idx']                                                 # cumulative variables by block
        N_var = self.HData.idx['idx'][-1]                                             # total variables
        N_sensors = len(set(self.HData.idx['parameter_sensor']))                      # total number of sensors
        N_cov = self.HData.idx['n_cov'][-1]                                           # total number of covariates
        N_mu_s = len(self.HData.idx['Im'])                                            # total number of match-up series
        n_uS = max([0]+[unc_i.uS_i for unc_i in self.HData.unc if (unc_i.typeID == 2) or (unc_i.typeID == 4)])

        # initialise array
        if not transpose:
            JPx = zeros(N_var + N_mu, dtype=float32)         # initialise JPx (length number of variables + number of ks)
        elif transpose:
            JPx = zeros(N_var + N_a, dtype=float32)          # initialise JPx (length number of variables + number of as)

        ################################################################################################################
        # 1. Apply preconditioner if not transpose
        ################################################################################################################

        if not transpose:
            x = self.calc_Px(x)

        ################################################################################################################
        # 2. Evaluate rows of J corresponding to the reference radiances and covariates
        ################################################################################################################

        JPx[0:N_var] = x[0:N_var]

        ################################################################################################################
        # 3. Evaluate rows of J corresponding to the adjustment factors
        ################################################################################################################
        for i, n_sensors in enumerate(self.HData.idx['Im']):

            n_mu = i + 1

            # indices for match-up data
            istart = N_var + mc[n_mu - 1]
            iend = N_var + mc[n_mu]

            # match-up k uncertainty
            uK = self.HData.unck[i].uR

            # a. get radiances for sensor by sensor

            for j, n_sensor in enumerate(n_sensors):

                # First sensor or second sensor factor
                sensor_name = self.HData.idx['sensors'][n_sensor]
                s = -1
                if j == 1:
                    s = 1

                R, JR = self.calc_R(self.xyza, self.HData.unc, self.HData.w_matrices, self.HData.u_matrices,
                                    self.HData.idx, self.HData._original_idx, self.HData.sensor_model,
                                    self.HData.sensor_model_constant,
                                    self.HData.getSensorAcrossTrackIndex(n_mu, sensor_name),
                                    self.HData.getSensorAlongTrackIndex(n_mu, sensor_name),
                                    self.HData.getSensorTime(n_mu, sensor_name),
                                    n_sensor, n_mu)

                # evaluate radiances and derivatives
                JB = self.HData.adjustment_model[n_sensor](R)[1]  # evaluate derivatives of adjustment model for each sensor

                # b. build product of (unweighted) Jacobian with vector

                # > if reference sensor
                if JR is None:

                    # i. find variables location of data
                    indices = [(j, k) for j, k in zip(self.HData.idx['n_sensor'], self.HData.idx['n_mu'])]
                    im = indices.index((n_sensor, n_mu))
                    ib = mcxyz[im]
                    ie = ib + int(self.HData.idx['N_var'][im])
                    block_unc = self.HData.unc[im]

                    # ii. add products of Jacobian and vector to product vector
                    if not transpose:
                        JPx[istart:iend] = JPx[istart:iend] + s * JB * x[ib:ie] * block_unc.uR / uK
                    elif transpose:
                        JPx[ib:ie] = JPx[ib:ie] + s * JB * x[istart:iend] * block_unc.uR / uK

                # > if sensor
                else:
                    for n_cov in xrange(1, N_cov + 1):

                        # i. find variables location of data
                        indices = [(i1, i2, i3) for i1, i2, i3 in zip(self.HData.idx['n_sensor'],
                                                                      self.HData.idx['n_mu'],
                                                                      self.HData.idx['n_cov'])]
                        im = indices.index((n_sensor, n_mu, n_cov))
                        ib = mcxyz[im]
                        ie = ib + int(self.HData.idx['N_var'][im])
                        block_unc = self.HData.unc[im]

                        # ii. add products of Jacobian and vector to product vector, computed product depends on
                        #     correlation form:

                        # ~ random correlation
                        if block_unc.typeID == 1:
                            if not transpose:
                                JPx[istart:iend] = JPx[istart:iend] + s * JB * JR[:, n_cov-1] * x[ib:ie] * block_unc.uR / uK
                            elif transpose:
                                JPx[ib:ie] = JPx[ib:ie] + s * JB * JR[:, n_cov - 1] * x[istart:iend] * block_unc.uR / uK

                        # ~ random+systematic correlation
                        if block_unc.typeID == 2:

                            # get index of required systematic value
                            isys = self.HData.idx['idx'][-1] - n_uS - 1 + block_unc.uS_i

                            if not transpose:
                                # > random component
                                JPx[istart:iend] = JPx[istart:iend] + s*JB*JR[:, n_cov-1]*x[ib:ie]*block_unc.uR/uK
                                # > systematic component
                                JPx[istart:iend] = JPx[istart:iend] + s*(JB*JR[:, n_cov-1])*x[isys]*block_unc.uS/uK
                                pass
                            elif transpose:
                                # > random component
                                JPx[ib:ie] = JPx[ib:ie] + s*JB*JR[:, n_cov-1]*x[istart:iend]*block_unc.uR/uK
                                # > systematic component
                                JPx[isys] = JPx[isys] + s * dot(JB*JR[:, n_cov-1]*block_unc.uS/uK, x[istart:iend])
                                pass

                        # ~ structured correlation
                        if block_unc.typeID == 3:
                            w = self.HData.w_matrices[block_unc.w_i]
                            u = self.HData.u_matrices[block_unc.u_i]
                            if not transpose:
                                JPx[istart:iend] = JPx[istart:iend] + s*JB*JR[:, n_cov-1]*w.dot(u*x[ib:ie])/uK
                            elif transpose:
                                JPx[ib:ie] = JPx[ib:ie] + s*u*w.T.dot(JB*JR[:, n_cov-1]*x[istart:iend]/uK)

                        # ~ structured+systematic correlation
                        if block_unc.typeID == 4:
                            # get index of required systematic value
                            isys = self.HData.idx['idx'][-1] - n_uS - 1 + block_unc.uS_i

                            w = self.HData.w_matrices[block_unc.w_i]
                            ur = self.HData.u_matrices[block_unc.u_i]
                            us = block_unc.uS*ones(len(ur), dtype=float32)
                            if not transpose:
                                # > random component
                                JPx[istart:iend] = JPx[istart:iend] + s*JB*JR[:, n_cov-1]*w.dot(ur*x[ib:ie])/uK
                                # > systematic component
                                JPx[istart:iend] = JPx[istart:iend] + s*(JB*JR[:, n_cov-1])*w.dot(us*x[isys])/uK

                            elif transpose:
                                # > random component
                                JPx[ib:ie] = JPx[ib:ie] + s*ur*w.T.dot(JB*JR[:, n_cov-1]*x[istart:iend]/uK)
                                # > systematic component
                                JPx[isys] = JPx[isys] + s * dot(JB*JR[:, n_cov-1]*w.dot(us)/uK, x[istart:iend])

                    # iii. add terms for as
                    i_sensor_p = asarray([N_var+i for i, sn in enumerate(self.HData.idx["parameter_sensor"]) if sn == sensor_name])
                    N_sensor_p = len(i_sensor_p)

                    if not transpose:
                        JPx[istart:iend] = JPx[istart:iend] \
                                           + s * dot(outer(JB / uK, ones(N_sensor_p)) * JR[:, N_cov:N_cov + N_sensor_p], x[i_sensor_p])
                    elif transpose:
                        JPx[i_sensor_p] = JPx[i_sensor_p] \
                                     + s * dot((outer(JB, ones(N_sensor_p)) * JR[:, N_cov:N_cov + N_sensor_p]).T, x[istart:iend] / uK)

        ################################################################################################################
        # 4. Apply preconditioner transpose if transpose
        ################################################################################################################

        if transpose:
            JPx = self.calc_Px(JPx, transpose=True)

        return JPx

    def calc_Px(self, x, transpose=False):
        """
        Return value of x multiplied by preconditioner solution P (or transpose)

        :type x: numpy.ndarray
        :param x: Vector to be multiplied by x

        :type transpose: bool
        :param transpose: Parameter to decide if to calculate Px or PT x

        :return:
            :Px: *numpy.ndarray*

            Product of preconditioner solution with x
        """

        # Preconditioner solution matrix structured as,
        #
        #     | I | 0 |
        # P = |---+---|
        #     | 0 | S |
        #
        # with x vector as,
        #
        #     |x1|
        # x = |--|
        #     |x2|
        #
        # So that their product is structured as,
        #
        #      | I | 0 | |x1|   | x1 |
        # Px = |---+---| |--| = |----|
        #      | 0 | S | |x2|   |S x2|
        #
        # If multiplying by the transpose of P have
        #
        #      | I | 0  |  |x1|   |  x1 |
        # Px = |---+----|  |--| = |-----|
        #      | 0 | ST |  |x2|   |ST x2|

        # parameters
        N_var = self.HData.idx['idx'][-1]   # total number of variables
        N_tot = N_var + len(self.HData.a)   # total amount of data (number of variables + number of parameters)

        # calculate product
        if not transpose:
            Px = hstack((x[0:N_var], dot(self.S, x[N_var:N_tot])))
            return Px

        if transpose:
            PTx = hstack((x[0:N_var], dot(self.S.T, x[N_var:N_tot])))
            return PTx

    def calculate_parameter_covariance_matrix(self, tolU=1e-8, show=False):
        """
        Return array containing the covariance matrices for the retrieved parameters for each sensor, derived
        analytically

        :return:
            :Va: *numpy.ndarray*

            Array containing the covariance matrices for the retrieved parameters of each sensors
        """

        # initialise parameters
        N_a = len(self.HData.a)

        # initialise array
        parameter_covariance_matrix = zeros((N_a, N_a), dtype=float32)

        ijs = [[i, j] for i in range(N_a) for j in range(N_a)]

        pool = Pool()
        parameter_covariance_matrix_ijs = pool.map(self.calculate_parameter_covariance_matrix_ij, ijs)
        for parameter_covariance_matrix_ij, ij in zip(parameter_covariance_matrix_ijs, ijs):
            parameter_covariance_matrix[ij[0], ij[1]] = parameter_covariance_matrix_ij

        parameter_covariance_matrix = (parameter_covariance_matrix + parameter_covariance_matrix.T)/2

        return parameter_covariance_matrix

    def calculate_parameter_covariance_matrix_ij(self, ij, tolU=1e-8, show=False):
        # initialise parameters
        N_a = len(self.HData.a)  # total number of parameters for
        # all sensors combined
        N_var = self.HData.idx['idx'][-1]  # total variables

        i = ij[0]
        j = ij[1]

        # initialise H LinearOperator
        H = LinearOperator((N_var + N_a, N_var + N_a), matvec=self.calc_Hx)

        # initialise MINRES object
        MinresOp = Minres(H)

        # Select element to solve for
        c = zeros(N_var + N_a, dtype=float32)
        d = zeros(N_var + N_a, dtype=float32)
        c[N_var + i] = 1
        d[N_var + j] = 1

        MinresOp.solve(self.calc_Px(d, transpose=True), rtol=tolU, itnlim=500, show=show)
        dt = MinresOp.x
        Vij = dot(self.calc_Px(c, transpose=True), dt)

        return Vij

    def calc_Hx(self, x):
        """
        Return the product of H (P'*J'*J*P) with a given x

        :type x: numpy.ndarray
        :param x: Vector to multiply by JP (or JP transpose)

        :return:
            :Hx: *numpy.ndarray*

            Array containing the product of H with x
        """

        # Evaluate matrix-vector product (P'*J'*J*P)*x
        Hx = self.calc_prod_JPx(self.calc_prod_JPx(x), transpose=True)

        return Hx

    def save(self, directory):

        # 1. Initialise file for data
        try:
            makedirs(directory)
        except OSError:
            pass

        dataset = Dataset(pjoin(directory, "GNOp_data.nc"), "w")

        # 2. Create dimensions

        M_T_dim = dataset.createDimension('M_T', len(self.HData.values))
        V_T_dim = dataset.createDimension('V_T', len(self.xyza))
        n_dim = dataset.createDimension('n', self.S.shape[0])

        # 3. Create variables

        # > values_t_initial variable
        values_t_initial_var = dataset.createVariable('values_t_initial', 'f4', ('M_T',), zlib=True, complevel=9)
        values_t_initial_var.description = "Initial values of match-up data " \
                                           "(transformed to normalised, independent parameterisation)"
        values_t_initial_var[:] = self.HData.values

        # > variables_t_estimate variable
        variables_t_estimate_var = dataset.createVariable('variables_t_estimate', 'f4', ('V_T',), zlib=True,
                                                          complevel=9)
        variables_t_estimate_var.description = "Estimated values of solver variables " \
                                               "(transformed to normalised, independent parameterisation)"
        variables_t_estimate_var[:] = self.xyza

        # > preconditioner solution variable
        preconditioner_solution_var = dataset.createVariable('preconditioner_solution', 'f4', ('n', 'n',), zlib=True,
                                                             complevel=9)
        preconditioner_solution_var.description = "Solution to problem preconditioner"
        preconditioner_solution_var[:] = self.S

        # 5. Close
        dataset.close()

        # 6. Pickle data idx
        with open(pjoin(directory, "GNOp_data.HData.idx.pickle"), "wb") as output_file:
            cPickle.dump(self.HData.idx, output_file)

    def open(self, GNData_directory):

        with open(pjoin(GNData_directory, "GNOp_data.HData.idx.pickle"), "rb") as input_file:
            self.HData.idx = cPickle.load(input_file)

        dataset = Dataset(pjoin(GNData_directory, "GNOp_data.nc"))

        self.S = dataset.variables['preconditioner_solution'][:]
        self.HData.values = dataset.variables['values_t_initial'][:]
        self.xyza = dataset.variables['variables_t_estimate'][:]
        dataset.close()

        return 0

if __name__ == "__main__":
    pass

