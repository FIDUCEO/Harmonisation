"""
Transform2NormInd Class - Processor to transform a MatchUp data object to set of variables with indpendent errors and
uncertainties of unity
"""

'''___Built-In Modules___'''
from copy import deepcopy
import sys
from os.path import dirname
from os.path import join as pjoin

'''___Third-Party Modules___'''
from numpy import zeros, bool_, float32, where
from scipy.sparse.linalg import lsqr

'''___NPL Modules___'''
matchupIO_directory = pjoin(dirname(dirname(dirname(__file__))), "matchupIO")
sys.path.append(matchupIO_directory)
from MatchUp import MatchUp

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "19/11/2017"
__credits__ = ["Jonathan Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class Transform2NormInd:
    """
    Transform2NormInd class operates on instances of ``eopy.matchup.matchupIO.MatchUp`` in-memory data products,
    returning a reparameterisation of the match-up such that the new set of variables has independent errors and
    uncertainties of unity

    :Methods:
        .. py:method:: run(...):

            Returns transformed instance of ``eopy.matchup.matchupIO.MatchUp``

        .. py:method:: reverse(...):

            Reverses transformation of ``eopy.matchup.matchupIO.MatchUp`` instance
    """

    def __init__(self):
        pass

    def run(self, MatchUpData):
        """
        Return a reparameterisation of the input data such that output data are independent quantities with
        uncertainties of unity

        :type MatchUpData: *eopy.matchup.matchupIO.MatchUp*
        :param MatchUpData: Input match-up data for transformation

        :return:
            :MatchUpData: *eopy.matchup.matchupIO.MatchUp*

            Transformed input data
        """

        # Convert data depending on it's correlation form to remove correlation:
        # 1. Random Form
        #    No correlation, so no action required. Scale by uncertainty.
        #
        # 2. Random+Systematic Form
        #    Separate random and systematic components:
        #    > random component - scale data by random uncertainty
        #    > systematic component - add 0 value for each block to the end of the final
        #                             block of the covariate
        #
        # 3. Average Form
        #    Simulate raw data used to compute averages (results in n_mu + n - 1 variables
        #    per block, where n_mu is the number of match-ups in the block and n is size of
        #    the averaging window)

        ################################################################################################################
        # 1. Determine idx for transformed data
        ################################################################################################################

        # Initialise copy of harmonisation idx to update for converted data
        # N.B. deepcopy() required to ensure copy of nested lists in dict
        independent_idx = deepcopy(MatchUpData.idx)

        # a. determine new number of variables per block ---------------------------------------------------------------
        for i, block_unc in enumerate(MatchUpData.unc):

            # i. number of variables remains the same for a block already with independent errors
            if block_unc.typeID == 1:
                pass

            # ii. number of variables for a block with error correlation defined by a w matrix changes to the number
            #     of independent variables it transforms to, or length of the corresponding u matrix
            elif (block_unc.typeID == 3) or (block_unc.typeID == 4):
                independent_idx['N_var'][i] = MatchUpData.u_matrices[block_unc.u_i].shape[0]

        # iii. Number of systematic error variables to introduce
        n_uS = max([0]+[unc_i.uS_i for unc_i in MatchUpData.unc if (unc_i.typeID == 2) or (unc_i.typeID == 4)])

        # --------------------------------------------------------------------------------------------------------------

        # b. determine new data block indices --------------------------------------------------------------------------
        idxs = [0]
        total = 0
        for N in independent_idx['N_var']:
            total += N
            idxs.append(int(total))

        # Add systematic error variables
        idxs[-1] += n_uS

        independent_idx['idx'] = idxs
        # --------------------------------------------------------------------------------------------------------------

        ################################################################################################################
        # 2. Determine transformed independent normalised data
        ################################################################################################################

        # Initialise transformed data product
        MatchUpNormInd = MatchUp()
        MatchUpNormInd.values = zeros(independent_idx['idx'][-1], dtype=float32)
        MatchUpNormInd.ks = zeros(independent_idx['cNm'][-1])
        MatchUpNormInd.a = MatchUpData.a[:]
        MatchUpNormInd.sensor_model = MatchUpData.sensor_model
        MatchUpNormInd.sensor_model_constant = MatchUpData.sensor_model_constant
        MatchUpNormInd.adjustment_model = MatchUpData.adjustment_model
        MatchUpNormInd.idx = independent_idx
        MatchUpNormInd._original_idx = MatchUpData._original_idx
        MatchUpNormInd.unc = MatchUpData.unc
        MatchUpNormInd.unck = MatchUpData.unck
        MatchUpNormInd.w_matrices = MatchUpData.w_matrices
        MatchUpNormInd.u_matrices = MatchUpData.u_matrices
        MatchUpNormInd.across_track_index1 = MatchUpData.across_track_index1
        MatchUpNormInd.across_track_index2 = MatchUpData.across_track_index2
        MatchUpNormInd.along_track_index1 = MatchUpData.along_track_index1
        MatchUpNormInd.along_track_index2 = MatchUpData.along_track_index2

        # Convert data block by depending on correlation form
        for i, block_unc in enumerate(MatchUpData.unc):
            istart = MatchUpData.idx['idx'][i]                        # start of original dataset values data block
            iend = istart + int(MatchUpData.idx['N_var'][i])          # number of variables in tranformed data block
            istart_i = independent_idx['idx'][i]                      # start of transformed dataset values data block
            iend_i = istart_i + int(MatchUpNormInd.idx['N_var'][i])   # number of variables in tranformed data block

            # a. independent type correlation --------------------------------------------------------------------------
            if (block_unc.typeID == 1) or (block_unc.typeID == 2):
                # scale data by uncertainty
                MatchUpNormInd.values[istart_i:iend_i] = MatchUpData.values[istart:iend]/block_unc.uR

            # c. structured type correlation ---------------------------------------------------------------------------
            elif (block_unc.typeID == 3) or (block_unc.typeID == 4):
                # Simulate independent data, X_ind, by determining solution to,
                #       X = W X_ind,
                # where W is the W matrix and X is the original data

                # Retrieve required W matrix and u matrix
                w_matrix = MatchUpData.w_matrices[block_unc.w_i]
                u_matrix = MatchUpData.u_matrices[block_unc.u_i]

                encountered_cols = zeros(w_matrix.shape[1], dtype=bool_)

                for i_row, i_values in enumerate(xrange(istart, iend)):

                    row_cols = w_matrix.indices[w_matrix.indptr[i_row]:w_matrix.indptr[i_row+1]]
                    row_encountered_cols = [bool(encountered_cols[row_col]) for row_col in row_cols]
                    row_unencountered_cols_num = row_encountered_cols.count(False)

                    w_row = w_matrix.data[w_matrix.indptr[i_row]:w_matrix.indptr[i_row+1]]
                    n_w = len(row_cols)

                    if row_unencountered_cols_num == n_w:
                        # averaged values
                        for col in row_cols:
                            MatchUpNormInd.values[istart_i+col] = MatchUpData.values[i_values] / u_matrix[col]
                            encountered_cols[col] = True
                        pass

                    elif 0 < row_unencountered_cols_num < n_w:

                        row_idx_e = [(i, col) for i, col in enumerate(row_cols) if encountered_cols[col] == True]
                        row_idx_une = [(i, col) for i, col in enumerate(row_cols) if encountered_cols[col] == False]

                        average_sofar = sum([w_row[col[0]] * MatchUpNormInd.values[istart_i+col[1]] *
                                             u_matrix[col[1]] for col in row_idx_e])
                        weight_remaining = sum([w_row[col[0]] for col in row_idx_une])

                        for col in row_idx_une:
                            MatchUpNormInd.values[istart_i+col[1]] = (MatchUpData.values[i_values] - average_sofar) / weight_remaining / u_matrix[col[1]]
                            encountered_cols[col[1]] = True
                        pass

                    elif row_unencountered_cols_num == 0:
                        pass

            # ----------------------------------------------------------------------------------------------------------

        # d. scale ks --------------------------------------------------------------------------------------------------
        for i in xrange(len(MatchUpData.idx['Im'])):
            istart = MatchUpData.idx['cNm'][i]
            iend = MatchUpData.idx['cNm'][i + 1]

            MatchUpNormInd.ks[istart:iend] = MatchUpData.ks[istart:iend] / MatchUpData.unck[i].uR
        # --------------------------------------------------------------------------------------------------------------

        return MatchUpNormInd

    def reverse(self, MatchUpNormInd):
        """
        Return reparameterised input match-up data in its original parameterisation

        :type MatchUpNormInd: *eopy.matchup.matchupIO.MatchUp*
        :param MatchUpNormInd: Transformed match-up data

        :return:
            :MatchUpData: *eopy.matchup.matchupIO.MatchUp*

            Input data with transformation reversed
        """

        # Initialise untransformed data product
        MatchUpData = MatchUp()
        MatchUpData.values = zeros(MatchUpNormInd._original_idx['idx'][-1])
        MatchUpData.ks = zeros(MatchUpNormInd._original_idx['cNm'][-1])
        MatchUpData.a = MatchUpNormInd.a[:]
        MatchUpData.sensor_model = MatchUpNormInd.sensor_model
        MatchUpData.sensor_model_constant = MatchUpNormInd.sensor_model_constant
        MatchUpData.adjustment_model = MatchUpNormInd.adjustment_model
        MatchUpData.idx = MatchUpNormInd._original_idx
        MatchUpData._original_idx = MatchUpNormInd._original_idx

        # todo - review how to better use memory here
        MatchUpData.unc = MatchUpNormInd.unc
        MatchUpData.unck = MatchUpNormInd.unck
        MatchUpData.w_matrices = MatchUpNormInd.w_matrices
        MatchUpData.u_matrices = MatchUpNormInd.u_matrices
        MatchUpData.across_track_index1 = MatchUpNormInd.across_track_index1
        MatchUpData.across_track_index2 = MatchUpNormInd.across_track_index2
        MatchUpData.along_track_index1 = MatchUpNormInd.along_track_index1
        MatchUpData.along_track_index2 = MatchUpNormInd.along_track_index2

        # Required to find systematic errors
        n_uS = max([0]+[unc_i.uS_i for unc_i in MatchUpData.unc if (unc_i.typeID == 2) or (unc_i.typeID == 4)])

        for i, block_unc in enumerate(MatchUpData.unc):

            istart = MatchUpData.idx['idx'][i]  # start of original dataset values data block
            iend = istart + int(MatchUpData.idx['N_var'][i])  # number of variables in tranformed data block
            istart_i = MatchUpNormInd.idx['idx'][i]  # start of transformed dataset values data block
            iend_i = istart_i + int(MatchUpNormInd.idx['N_var'][i])  # number of variables in tranformed data block

            # a. random correlation - rescale and add to covariate list
            if block_unc.typeID == 1:
                MatchUpData.values[istart:iend] = MatchUpNormInd.values[istart_i:iend_i] * block_unc.uR

            # b. random+systematic correlation - rescale components and recombine
            if block_unc.typeID == 2:

                # get index of required systematic value
                isys = MatchUpNormInd.idx['idx'][-1] - n_uS - 1 + block_unc.uS_i

                MatchUpData.values[istart:iend] = MatchUpNormInd.values[istart_i:iend_i]*block_unc.uR
                MatchUpData.values[istart:iend] += MatchUpNormInd.values[isys]*block_unc.uS

            # c. structured correlation - transform from independent to original variables
            if block_unc.typeID == 3:

                # Retrieve required W matrix and u matrix
                w = MatchUpData.w_matrices[block_unc.w_i]
                u = MatchUpData.u_matrices[block_unc.u_i]

                MatchUpData.values[istart:iend] = w.dot(u*MatchUpNormInd.values[istart_i:iend_i])

            # d. structured+systematic correlation - add sys error, then transform from independent to original variable
            if block_unc.typeID == 4:
                # Retrieve required W matrix and u matrix
                w = MatchUpData.w_matrices[block_unc.w_i]
                u = MatchUpData.u_matrices[block_unc.u_i]
                isys = MatchUpNormInd.idx['idx'][-1] - n_uS - 1 + block_unc.uS_i

                MatchUpData.values[istart:iend] = w.dot(MatchUpNormInd.values[istart_i:iend_i] * u
                                                        + MatchUpNormInd.values[isys] * block_unc.uS)

        # d. rescale ks ------------------------------------------------------------------------------------------------
        for i in xrange(len(MatchUpData.idx['Im'])):
            istart = MatchUpData.idx['cNm'][i]
            iend = MatchUpData.idx['cNm'][i + 1]

            MatchUpData.ks[istart:iend] = MatchUpNormInd.ks[istart:iend] * MatchUpData.unck[i].uR

        return MatchUpData

if __name__ == "__main__":
    pass
