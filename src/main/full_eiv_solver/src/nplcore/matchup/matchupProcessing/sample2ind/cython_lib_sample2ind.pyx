from random import sample
from numpy import asarray, zeros, bool_, array, append
import numpy as np
cimport numpy as np
DTYPE = np.float32
ctypedef np.float32_t DTYPEf_t
ctypedef np.int32_t DTYPEi_t
ctypedef np.int8_t DTYPEis_t


def get_sample_idxs(int[:] w_indices_mu, int mu_samples, list w_matrices):
    # Store truth value of if the transformed independent variables correlated with any samples per w matrix
    # Want to maintain that each transformed independent variables is only related to one sampled match up
    cdef:
        int w_i, i, j, col, row_idx, next_row_idx
        list sampled_cols, prev_sampled_cols
        np.ndarray[DTYPEi_t, ndim=1] idxs = array([], dtype=np.int32)
        np.ndarray[DTYPEi_t, ndim=1] row_cols, independence_test
        tuple row_idxs, next_row_idxs
        list row_cols_ws

    sampled_cols = [zeros(w_matrices[w_i].shape[1], dtype=bool_) for w_i in w_indices_mu]

    # initialise list to store match up series sample indices

    # Loop through match up series by match up
    # - i.e. by row for w matrices - using the new row index variable in the CSR sparse matrix format

    for i, (row_idxs, next_row_idxs) in enumerate(zip(zip(*[w_matrices[w_i].indptr[:-1] for w_i in w_indices_mu]),
                                                      zip(*[w_matrices[w_i].indptr[1:] for w_i in w_indices_mu]))):

        # Use new row index variable from the CSR sparse format to find the column indices of non-zero
        # in given row
        row_cols_ws = [asarray(w_matrices[w_i].indices[row_idx:next_row_idx])
                       for w_i, row_idx, next_row_idx in zip(w_indices_mu, row_idxs, next_row_idxs)]

        # Find if the row column indices per w matrix are all unused in previous samples - i.e. uncorrelated
        independence_test = array([], dtype=np.int32)     # test results, 1 - uncorrelated, 0 - correlated

        for j, row_cols in enumerate(row_cols_ws):

            # List of true/false if row cols previously used
            prev_sampled_cols = [bool(sampled_cols[j][row_col]) for row_col in row_cols]

            # All must be unused previously to be uncorrelated
            if all(col is False for col in prev_sampled_cols):
                independence_test = np.append(independence_test, array([1], dtype=np.int32))
            else:
                independence_test = np.append(independence_test, array([0], dtype=np.int32))

        # If row uncorrelated for all w matrices sample this match up
        if set(independence_test) == {1}:
            for j, row_cols in enumerate(row_cols_ws):
                sampled_cols[j][row_cols] = True
            idxs = append(idxs, array(i, dtype=np.int32))

    # If more match ups sampled than maximum randomly reduce to required amount
    if len(idxs) > mu_samples:
        return sample(idxs, mu_samples)
    return idxs