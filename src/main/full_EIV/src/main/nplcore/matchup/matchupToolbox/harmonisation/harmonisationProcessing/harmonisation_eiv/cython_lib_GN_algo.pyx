import numpy as np
cimport numpy as np
DTYPE = np.float32
ctypedef np.float32_t DTYPEf_t
ctypedef np.int32_t DTYPEi_t
ctypedef np.int8_t DTYPEis_t


cpdef np.ndarray[DTYPEf_t, ndim=2] unconvert_Xs(np.ndarray[DTYPEf_t, ndim=1] xyza,
                                                np.ndarray[DTYPEi_t, ndim=1] n_sensors,
                                                np.ndarray[DTYPEi_t, ndim=1] idxs,
                                                np.ndarray[DTYPEi_t, ndim=1] original_idxs,
                                                np.ndarray[DTYPEi_t, ndim=1] n_covs,
                                                np.ndarray[DTYPEi_t, ndim=1] n_mus,
                                                np.ndarray[DTYPEi_t, ndim=1] N_vars,
                                                list unc,
                                                list w_matrices,
                                                list u_matrices,
                                                int n_sensor,
                                                int n_mu,
                                                int N_cov):

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

    cdef:
        int n_cov, i1, i2, i3, im, ib, ie, ib_original, ie_original, isys
        int N_sensors = np.unique(n_sensors).shape[0]     # total number of sensors
        int N_mu_s = np.unique(n_mus).shape[0]     # total number of sensors
        np.ndarray[DTYPEf_t, ndim=1] u
        np.ndarray[DTYPEf_t, ndim=2] Xs
        list indices = [(i1, i2, i3) for i1, i2, i3 in zip(n_sensors, n_mus, n_covs)] # block data
        int n_uS = int(max([0]+[unc_i.uS_i for unc_i in unc if (unc_i.typeID == 2) or (unc_i.typeID == 4)]))

    for n_cov in xrange(1, N_cov + 1):

        # find location of data in xyza
        im = indices.index((n_sensor, n_mu, n_cov))
        ib = idxs[im]
        ie = idxs[im] + N_vars[im]
        ib_original = original_idxs[im]
        ie_original = original_idxs[im+1]

        if n_cov == 1:
            Xs = np.zeros((ie_original-ib_original, N_cov), dtype=np.float32)

        # undo conversion of data from ConvertData.convert4GN depending on correlation form
        # and add to radiances list
        block_unc = unc[im]

        # a. random correlation - unscale and add to covariate list
        if block_unc.typeID == 1:
            Xs[:, n_cov-1] = xyza[ib:ie] * block_unc.uR

        # b. random+systematic correlation - unscale components and recombine
        if block_unc.typeID == 2:
            # get index of required systematic value
            isys = idxs[-1] - n_uS - 1 + block_unc.uS_i

            Xs[:, n_cov-1] = xyza[ib:ie] * block_unc.uR + xyza[isys] * block_unc.uS

        # c. structured correlation - transform independent counts to counts
        if block_unc.typeID == 3:
            # Retrieve required W and U matrices
            w = w_matrices[block_unc.w_i]
            u = u_matrices[block_unc.u_i]

            Xs[:, n_cov-1] = w.dot(u * xyza[ib:ie])

        # d. structured+systematic correlation - add systematic, unscale and transform independent counts to counts
        if block_unc.typeID == 4:
            # Retrieve required W and U matrices and systematic uncertainty location
            isys = idxs[-1] - n_uS - 1 + block_unc.uS_i
            w = w_matrices[block_unc.w_i]
            u = u_matrices[block_unc.u_i]

            Xs[:, n_cov-1] = w.dot(u * xyza[ib:ie] + xyza[isys] * block_unc.uS)

    return Xs