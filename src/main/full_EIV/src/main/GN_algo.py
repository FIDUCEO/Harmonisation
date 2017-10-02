"""
Implements the Gauss-Newton iteration algorithm for match-up data with a pre-conditioned solution.
"""

'''___Python Modules____'''
from numpy import zeros, append, ones, dot, outer, hstack, array, eye, inf
from numpy.linalg import norm, solve
from scipy.sparse.linalg import LinearOperator
from math import ceil

'''___Third Party Modules____'''
from pykrylov_lsmr import LSMRFramework
from pykrylov_minres import Minres

'''___Harmonisation Modules___'''

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "09/01/2017"
__credits__ = ["Arta Dillo", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class GNAlgo:
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

        *harm_data_reader.HarmData*

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

        .. py:method:: unconvert_values(...):

            Return variable data for each covariate in the original form for all variable data, undoing the
            reparameterisation performed in ConvertData.convert2ind()

        .. py:method:: unconvert_ks(...):

            Return variable data for each covariate in the original for reparameterisation performed in
            ConvertData.convert2ind()
    """

    def __init__(self, HData=None, S=None):
        """
        Initialise algorithm

        :type HData: HarmData
        :param HData: Input data to be harmonised

        :type S: numpy.ndarray
        :param S: Pre-conditioner solution
        """

        # Initialise class
        self.HData = None
        self.S = None
        self.xyza = None

        if (HData is not None) and (S is not None):

            self.HData = HData
            self.S = S

            # initialise current variable and parameter estimates
            self.xyza = append(self.HData.values, self.HData.a)

            print "Initial Parameter Estimates:"
            print self.HData.a

    def runGN(self, tol=1e-6, tolA=1e-8, tolB=1e8, tolU=1e-8, show=False):
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

        # Initialise parameters
        N_mu = self.HData.idx['cNm'][-1]  # total match-ups (= number of ks)
        N_sensors = len(set([i for pair in self.HData.idx['Im'] for i in pair])) - 1  # total number of sensors
        N_a = len(self.HData.a)                                                       # total number of parameters for
                                                                                      # all sensors combined
        N_p = N_a / N_sensors                                                         # total number of parameters for
                                                                                      # each sensor
        N_var = self.HData.idx['idx'][-1]                                             # total number of variables
        niter = 0                                                                     # counter of iterations
        mxiter = ceil(N_var)                                                          # max number of iterations of GN
        mxiter_lsmr = ceil(N_var)                                                     # max number of iterations of LSMR
        conv = False                                                                  # convergence boolean

        # Initialise J LinearOperator
        J = LinearOperator((N_var+N_mu, N_var+N_a), matvec=self.get_JPx, rmatvec=self.get_JPTx)

        # Initialise LSMR object
        K = LSMRFramework(J)

        # Gauss-Newton iterations
        while (conv is False) and (niter < mxiter):

            # Evaluate vector f of weighted residual deviations, sum of squares of
            # residual deviations F and gradient g.

            # Initialise for first iteration
            if niter == 0:
                f = self.calc_f(self.xyza, self.HData)
                F0 = norm(f)**2
                GNlog = []

            niter += 1

            # Determine Gauss-Newton step d as solution to linear least-squares
            # problem J*d = -f with the pre-conditioner applied.
            K.solve(-f, damp=0, atol=tolA, btol=tolA, conlim=tolB, itnlim=mxiter_lsmr, show=show)
            d = self.calc_Px(K.x)

            # Update parameter estimates, as well as f, F and g
            self.xyza += d
            f = self.calc_f(self.xyza, self.HData)
            F = norm(f)**2
            g = 2 * self.get_JPTx(f)

            # Test convergence
            U1 = F0 - F
            tol1 = tol*(1+F)
            U2 = norm(d, inf)
            tol2 = (tol**0.5) * (1 + norm(self.xyza, inf))
            U3 = norm(g, inf)
            tol3 = (tol**(1./3.))*(1+F)

            # Update F0
            F0 = F

            # Check for convergence
            if (U1 > 0) and (U1 < tol1) and (U2 < tol2) and (U3 <= tol3):
                conv = True

            # Write log
            GNlog.append([niter, U1, tol1, U2, tol2, U3, tol3])
            if show:
                print "\n\t\t\t\tGNlog"
                print "niter\tU1\t\ttol1\t\tU2\t\ttol2\t\tU3\t\ttol3"
                for GN in GNlog:
                    print "{0:2d}\t{1:.2e}\t{2:.2e}\t{3:.2e}\t{4:.2e}\t{5:.2e}\t{6:.2e}"\
                          .format(GN[0], GN[1], GN[2], GN[3], GN[4], GN[5], GN[6])

        # Unpack solution
        a = self.xyza[N_var:N_var + N_sensors * N_p]

        print "Determined Parameter Estimates:"
        print a

        # Uncertainty evaluation
        print 'Determining uncertainty...'
        V = self.calc_unc(tolU, show=show)

        print 'Preparing output...'

        values_res = self.unconvert_values(f[:N_var], self.HData.unc, self.HData.idx, self.HData.idx_orig)
        k_res = self.unconvert_ks(f[N_var:], self.HData.unck, self.HData.idx)

        v = N_var - N_mu - N_a
        p = 0

        print "Determined Parameter Covariance Matrix"
        print V

        return a, V, F, v, p, values_res, k_res

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
        f = zeros(N_var + N_mu)

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

        # initialise k_est for estimates of Ks
        k_est = zeros(N_mu)

        # determine k_est per match-up series
        for i, n_sensors in enumerate(HData.idx['Im']):

            n_mu = i + 1
            # indices for data
            istart = mc[n_mu-1]
            iend = mc[n_mu]

            # a. get radiances for sensor 1 and sensor 2 of match-up series
            Rs = []
            for n_sensor in n_sensors:
                Rs.append(self.calc_R(xyza, HData.unc, HData.idx, HData.sensor_model, n_sensor, n_mu)[0])

            # b. calculate k_est for match-up series
            Bs = [HData.adjustment_model(R)[0] for R in Rs]          # Evaluate adjustment model

            k_est[istart:iend] = (Bs[1] - Bs[0]) / HData.unck[i].uR  # Evaluate estimate of adjustment factor

        # c. add difference between k_est and k (original) for all match-up series to f
        f[N_var:N_var+N_mu] = k_est - HData.ks

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

    def calc_R(self, xyza, unc, idx, sensor_model, n_sensor, n_mu):
        """
        Return calculated radiance R and derivatives for given input data

        :type xyza: numpy.ndarray
        :param xyza: Array containing variables and parameters

        :type unc: numpy.ndarray
        :param unc: array containing uncertainty information for variables

        :type idx: dict
        :param idx: dictionary with entries describing the data structures

        :type sensor_model: func
        :param sensor_model: function to calculate radiance from input variables and parameters

        :type n_sensor: int
        :param n_sensor: number of sensor to calculate radiances for

        :type n_mu: int
        :param n_mu: number of match-up series to calculate radiances for

        :return:
            :R: *numpy.ndarray*

            Calculated radiances

            :J: *numpy.ndarray*

            Calculated derivatives
        """

        # initialise parameters
        mcxyz = idx['idx']                                                 # cumulative variables by block
        N_var = idx['idx'][-1]                                             # total variables
        N_sensors = len(set([i for pair in idx['Im'] for i in pair])) - 1  # total number of sensors
        N_mu_s = len(idx['Im'])                                            # total number of match-up series
        N_p = len(xyza[idx['idx'][-1]:]) / N_sensors                       # total number of parameters in model
        N_cov = idx['n_cov'][-1]                                           # total number of covariates

        # > if reference sensor - get R data
        if n_sensor == 0:

            # find location of data in xyza
            indices = [(j, k) for j, k in zip(idx['n_sensor'], idx['n_mu'])]
            im = indices.index((n_sensor, n_mu))
            ib = mcxyz[im]
            ie = ib + int(idx['N_var'][im])

            # undo conversion of data from ConvertData.convert4GN and add to radiances list
            block_unc = unc[im]
            R = xyza[ib:ie] * block_unc.uR
            JR = None

        # > if sensor - determine R from sensor model
        else:
            # 1. get covariate data
            Xs = self.unconvert_Xs(xyza, unc, idx, n_sensor, n_mu)

            # 2. get parameters
            a = xyza[N_var + (n_sensor - 1) * N_p:N_var + n_sensor * N_p]

            # 3. compute radiances and derivatives
            R, JR = sensor_model(a, Xs)

        return R, JR

    def unconvert_Xs(self, xyza, unc, idx, n_sensor, n_mu):
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

        N_cov = idx['n_cov'][-1]  # total number of covariates
        mcxyz = idx['idx']  # cumulative variables by block
        N_sensors = len(set([i for pair in idx['Im'] for i in pair])) - 1  # total number of sensors
        N_mu_s = len(idx['Im'])  # total number of match-up series

        Xs = []  # initialise list for covariate data

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
            if block_unc.form == "r":
                Xs.append(xyza[ib:ie] * block_unc.uR)

            # b. random+systematic correlation - unscale components and recombine
            if block_unc.form == 'rs':
                # get index of required systematic value
                indices = [(i1, i2, i3) for i1, i2, i3 in zip(idx['n_sensor'], idx['n_mu'], idx['n_cov'])]
                im = indices.index((N_sensors, N_mu_s, n_cov))
                isys = mcxyz[im + 1] - N_sensors + n_sensor - 1

                Xs.append(xyza[ib:ie] * block_unc.uR + xyza[isys] * block_unc.uS)

            # c. averaging correlation - average raw counts to counts
            if block_unc.form == 'ave':
                Xs.append(block_unc.W.dot(xyza[ib:ie]))

        return Xs

    def calc_prod_JPx(self, x, transpose=False):
        """
        Return the product of JP (or JP transpose) for a given x

        :globals:
            :self.HData: HarmData
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
        N_sensors = len(set([i for pair in self.HData.idx['Im'] for i in pair]))-1    # total number of sensors
        N_p = N_a/N_sensors                                                           # total number of parameters in
                                                                                      # each sensor model
        N_cov = self.HData.idx['n_cov'][-1]                                           # total number of covariates
        N_mu_s = len(self.HData.idx['Im'])                                            # total number of match-up series

        # initialise array
        if not transpose:
            JPx = zeros(N_var + N_mu)         # initialise JPx (length number of variables + number of ks)
        elif transpose:
            JPx = zeros(N_var + N_a)          # initialise JPx (length number of variables + number of as)

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
            istart = N_var+mc[n_mu-1]
            iend = N_var+mc[n_mu]

            # match-up k uncertainty
            uK = self.HData.unck[i].uR

            # a. get radiances for sensor by sensor

            for j, n_sensor in enumerate(n_sensors):

                # First sensor or second sensor factor
                s = -1
                if j == 1:
                    s = 1

                R, JR = self.calc_R(self.xyza, self.HData.unc, self.HData.idx, self.HData.sensor_model, n_sensor, n_mu)
                                                           # evaluate radiances and derivatives
                JB = self.HData.adjustment_model(R)[1]     # evaluate derivatives of adjustment model for each sensor

                # b. build product of (unweighted) Jacobian with vector

                # > if reference sensor
                if n_sensor == 0:

                    # i. find variables location of data
                    indices = [(j, k) for j, k in zip(self.HData.idx['n_sensor'], self.HData.idx['n_mu'])]
                    im = indices.index((n_sensor, n_mu))
                    ib = mcxyz[im]
                    ie = ib + int(self.HData.idx['N_var'][im])
                    block_unc = self.HData.unc[im]

                    # ii. add products of Jacobian and vector to product vector
                    if not transpose:
                        JPx[istart:iend] = JPx[istart:iend] + s*JB*x[ib:ie]*block_unc.uR/uK
                    elif transpose:
                        JPx[ib:ie] = JPx[ib:ie] + s*JB*x[istart:iend]*block_unc.uR/uK

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

                        # ii. add products of Jacobian and vector to product vector, computed producted depends on
                        #     correlation form:

                        # ~ random correlation
                        if block_unc.form == "r":
                            if not transpose:
                                JPx[istart:iend] = JPx[istart:iend] + s*JB*JR[:, n_cov-1]*x[ib:ie]*block_unc.uR/uK
                            elif transpose:
                                JPx[ib:ie] = JPx[ib:ie] + s*JB*JR[:, n_cov-1]*x[istart:iend]*block_unc.uR/uK

                        # ~ random+systematic correlation
                        if block_unc.form == 'rs':

                            # get index of required systematic value
                            indices = [(i1, i2, i3) for i1, i2, i3 in zip(self.HData.idx['n_sensor'],
                                                                          self.HData.idx['n_mu'],
                                                                          self.HData.idx['n_cov'])]
                            im = indices.index((N_sensors, N_mu_s, n_cov))
                            isys = mcxyz[im + 1] - N_sensors + n_sensor - 1

                            if not transpose:
                                # > random component
                                JPx[istart:iend] = JPx[istart:iend] + s*JB*JR[:, n_cov-1]*x[ib:ie]*block_unc.uR/uK
                                # > systematic component
                                JPx[istart:iend] = JPx[istart:iend] + s*(JB*JR[:, n_cov-1])*x[isys]*block_unc.uS/uK
                            elif transpose:
                                # > random component
                                JPx[ib:ie] = JPx[ib:ie] + s*JB*JR[:, n_cov-1]*x[istart:iend]*block_unc.uR/uK
                                # > systematic component
                                JPx[isys] = JPx[isys] + s*dot(JB*JR[:, n_cov-1]*block_unc.uS/uK, x[istart:iend])

                        # ~ averaging correlation
                        if block_unc.form == 'ave':
                            if not transpose:
                                JPx[istart:iend] = JPx[istart:iend] \
                                                        + s*JB*JR[:, n_cov-1]*block_unc.W.dot(x[ib:ie])/uK
                            elif transpose:
                                JPx[ib:ie] = JPx[ib:ie] \
                                                        + s*block_unc.W.T.dot(JB*JR[:, n_cov-1]*x[istart:iend]/uK)

                    # iii. add terms for as
                    ib = N_var + (n_sensor - 1) * N_p
                    ie = N_var + n_sensor * N_p

                    if not transpose:
                        JPx[istart:iend] = JPx[istart:iend] \
                                           + s *dot(outer(JB/uK, ones(N_p))*JR[:, N_cov:N_cov+N_p], x[ib:ie])
                    elif transpose:
                        JPx[ib:ie] = JPx[ib:ie] \
                                     + s *dot((outer(JB, ones(N_p))*JR[:, N_cov:N_cov+N_p]).T, x[istart:iend]/uK)

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

    def calc_unc(self, tolU=1e-8, show=False):
        """
        Return array containing the covariance matrices for the retrieved parameters for each sensor, derived
        analytically

        :return:
            :Va: *numpy.ndarray*

            Array containing the covariance matrices for the retrieved parameters of each sensors
        """

        # initialise parameters
        N_a = len(self.HData.a)                                                       # total number of parameters for
                                                                                      # all sensors combined
        N_var = self.HData.idx['idx'][-1]                                             # total variables
        N_sensors = len(set([i for pair in self.HData.idx['Im'] for i in pair])) - 1  # total number of sensors
        N_p = N_a / N_sensors                                                         # total number of parameters in
                                                                                      # model
        N_cov = self.HData.idx['n_cov'][-1]                                           # total number of covariates

        # initialise array
        Va = zeros((N_p, N_p, N_sensors))

        # initialise H LinearOperator
        H = LinearOperator((N_var + N_a, N_var + N_a), matvec=self.calc_Hx)

        # initialise MINRES object
        K = Minres(H)

        for n_sensor in xrange(1, N_sensors+1):

            for n_p1 in xrange(N_p):
                d = zeros(N_var + N_a)
                d[N_var + (n_sensor-1)*N_p + n_p1] = 1
                K.solve(self.calc_Px(d, transpose=True), rtol=tolU, itnlim=500, show=show)
                dt = K.x

                for n_p2 in xrange(N_p):
                    c = zeros(N_var + N_a)
                    c[N_var + (n_sensor-1)*N_p + n_p2] = 1
                    Va[n_p2, n_p1, n_sensor-1] = dot(self.calc_Px(c, transpose=True), dt)

            Va[:, :, n_sensor-1] = (Va[:, :, n_sensor-1]+Va[:, :, n_sensor-1].T)/2

        V = zeros((N_p*N_sensors, N_p*N_sensors))

        for i in range(N_sensors):
            V[N_p*i:N_p*(i+1), N_p*i:N_p*(i+1)] = Va[:, :, i]

        return V

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

    def unconvert_values(self, values_con, unc, idx, idx_orig):
        """
        Return variable data for each covariate in the original form for all variable data, undoing the
        reparameterisation performed in ConvertData.convert2ind()

        :type values_con: numpy.ndarray
        :param values_con: Converted variable data

        :type unc: list
        :param unc: Uncertainties associated with blocks of variable data

        :type idx: dict
        :param idx: Dictionary describing structure of variable data

        :return:
            :H: *numpy.darrays*

            values in original H format
        """

        N_cov = idx['n_cov'][-1]  # total number of covariates
        mcxyz = idx['idx']  # cumulative variables by block
        N_mu = idx['cNm'][-1]
        N_sensors = len(set([i for pair in idx['Im'] for i in pair])) - 1  # total number of sensors
        N_mu_s = len(idx['Im'])  # total number of match-up series

        values = array([])  # initialise list for covariate data

        # get covariate data covariate by covariate
        for i in xrange(len(idx['n_mu'])):

            # find location of data in xyza
            ib = mcxyz[i]
            ie = ib + int(idx['N_var'][i])

            # undo conversion of data from ConvertData.convert4GN depending on correlation form
            # and add to radiances list
            block_unc = unc[i]

            # a. random correlation - unscale and add to covariate list
            if block_unc.form == "r":
                values = append(values, values_con[ib:ie] * block_unc.uR)

            # b. random+systematic correlation - unscale components and recombine
            if block_unc.form == 'rs':
                # get index of required systematic value
                indices = [(i1, i2, i3) for i1, i2, i3 in zip(idx['n_sensor'], idx['n_mu'], idx['n_cov'])]
                im = indices.index((N_sensors, N_mu_s, idx['n_cov'][i]))
                isys = mcxyz[im + 1] - N_sensors + idx['n_sensor'][i] - 1

                values = append(values, values_con[ib:ie] * block_unc.uR + values_con[isys] * block_unc.uS)

            # c. averaging correlation - average raw counts to counts
            if block_unc.form == 'ave':
                values = append(values, block_unc.W.dot(values_con[ib:ie]))

        # Reformat into data arrays
        H = zeros((N_mu, 2*N_cov))

        for k in range(len(idx['n_mu'])):

            # indices defining first and last positions in data matrix
            istartm = idx['cNm'][idx['n_mu'][k]-1]
            iendm = idx['cNm'][idx['n_mu'][k]]

            # indices defining first and last positions in data vector
            istartv = idx_orig['idx'][k]
            iendv = idx_orig['idx'][k+1]

            # if the sensor is the first sensor in the match-up series
            if idx['Im'][idx['n_mu'][k]-1][0] == idx['n_sensor'][k]:
                H[istartm:iendm, idx['n_cov'][k] - 1] = values[istartv:iendv]

            # if the sensor is the second sensor in the match-up series:
            if idx['Im'][idx['n_mu'][k]-1][1] == idx['n_sensor'][k]:
                H[istartm:iendm, idx['n_cov'][k]+N_cov-1] = values[istartv:iendv]

        return H

    def unconvert_ks(self, ks_con, unck, idx):
        """
        Return variable data for each covariate in the original form for all variable data, undoing the
        reparameterisation performed in ConvertData.convert2ind()

        :type ks_con: numpy.ndarray
        :param ks_con: Converted k data

        :type unck: list
        :param unck: Uncertainties associated with k data

        :type idx: dict
        :param idx: Dictionary describing structure of data

        :return:
            :ks: *numpy.ndarray*

            Unconverted form of input k
        """

        mc = idx['cNm']  # cumulative variables by block

        ks = array([])  # initialise list for covariate data

        # scale ks
        for i in xrange(len(idx['Im'])):
            ib = mc[i]
            ie = mc[i+1]

            ks = append(ks, ks_con[ib:ie]*unck[i].uR)

        return ks

if __name__ == "__main__":

    def main():
        return 0

    main()
