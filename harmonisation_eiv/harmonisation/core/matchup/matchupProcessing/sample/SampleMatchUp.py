"""
Sample2Ind Class - Processor to sample a MatchUp data object such that the sample has independent errors
"""

'''___Built-In Modules___'''
from random import sample
from copy import deepcopy

'''___Third-Party Modules___'''
from numpy import zeros, arange
from numpy import sum as npsum

'''___harmonisation Modules___'''
from harmonisation import MatchUp, Uncertainty


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "23/11/2017"
__credits__ = ["Jonathan Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class SampleMatchUp:
    """
    SampleMatchUp class operates on instances of ``eopy.matchup.matchupIO.MatchUp`` in-memory data products, returning a
    random sample of the match-up

    :Methods:
        .. py:method:: run(...):

            Returns sampled instance of ``eopy.matchup.matchupIO.MatchUp``
    """

    def __init__(self):
        pass

    def run(self, MatchUpData, sf=1, samples=None, show=False):
        """
        Return a randomly sampled instance of ``eopy.matchup.matchupIO.MatchUp``

        :type MatchUpData: *eopy.matchup.matchupIO.MatchUp*
        :param MatchUpData: Input match-up data for sampling

        :type sf: int
        :param sf: Sampling factor

        :type samples: int
        :param samples: Number of samples to include per match-up series

        :return:
            :MatchUpSample: *eopy.matchup.matchupIO.MatchUp*

            Sampled harmonisation data
        """

        # initialise parameters
        mcxyz = MatchUpData.idx['idx']       # cumulative total of variables data block
        mc = MatchUpData.idx['cNm']          # cumulative total of match-ups by series

        ################################################################################################################
        # 1. Determine Match Ups to Include in Sample
        ################################################################################################################

        # Choose a sample of match ups such that data is left with independent errors.
        #
        # N.B. - This is not possible for a fully systematic effect, but data with error correlation structures defined
        #        by w matrices can be sampled to achieve this

        sampling_idxs = {}  # initialise dictionary of sampling indices per match up series

        n_mus = set(MatchUpData.idx['n_mu'])

        # Find sampling indices by match-up series
        for n_mu in n_mus:

            # total number of match up in series (should be the same for each variable so take the first one)
            mu_total = [MatchUpData.idx['N_var'][i] for i, n_mu_i in enumerate(MatchUpData.idx['n_mu']) if n_mu_i == n_mu][0]
            mu_samples = mu_total/sf                  # required number sample match ups (determined by scaling factor)

            idxs = sorted(sample(arange(mu_total), mu_samples))
            sampling_idxs[n_mu] = idxs       # add match up sample indices to dictionary

        ################################################################################################################
        # 2. Sample Data
        ################################################################################################################

        # a. Initialise sampled harmonisation data product -------------------------------------------------------------
        MatchUpSample = MatchUp()
        MatchUpSample.a = MatchUpData.a[:]
        MatchUpSample.sensor_model = MatchUpData.sensor_model
        MatchUpSample.sensor_model_constant = MatchUpData.sensor_model_constant
        MatchUpSample.adjustment_model = MatchUpData.adjustment_model
        # --------------------------------------------------------------------------------------------------------------

        # b. Update idx attribute of MatchUpSample to describe structure of sampled data -------------------------------

        # Start with copy of full dataset idx dictionary attribute
        # N.B. - deepcopy because of python memory issues with lists nested in dicts
        MatchUpSample.idx = deepcopy(MatchUpData.idx)

        # Formulate required replacement idx entries (several can remain the same as the full dataset)
        idxs = [0]
        total = 0
        for i, n_mu in enumerate(MatchUpData.idx['n_mu']):
            block_samples = len(sampling_idxs[n_mu])
            MatchUpSample.idx['N_var'][i] = block_samples
            total += block_samples
            idxs.append(int(total))
        MatchUpSample.idx['idx'] = idxs

        cNm = [0]
        total = 0
        for i, n_mu in enumerate(n_mus):
            n_mu_sample = len(sampling_idxs[n_mu])
            MatchUpSample.idx['Nm'][i] = n_mu_sample
            total += n_mu_sample
            cNm.append(total)
        MatchUpSample.idx['cNm'] = cNm

        if show:
            print "Sample Size: ", MatchUpSample.idx['Nm']
        # --------------------------------------------------------------------------------------------------------------

        # c. Sample variables and respective uncertainty ---------------------------------------------------------------

        # Initialise data arrays
        MatchUpSample.values = zeros(MatchUpSample.idx['idx'][-1])
        MatchUpSample.unc = [0] * len(MatchUpSample.idx['n_cov'])

        #  Sample data by data block
        for i, block_unc in enumerate(MatchUpData.unc):
            istart = mcxyz[i]                             # start of full dataset values data block
            iend = mcxyz[i+1]                             # end of full dataset values data block
            istart_s = MatchUpSample.idx['idx'][i]         # start of sampled values data block
            iend_s = MatchUpSample.idx['idx'][i+1]         # end of sampled values data block
            s_idx = sampling_idxs[MatchUpData.idx['n_mu'][i]]   # indices of match up to sample within values data block

            # i. Sample values
            MatchUpSample.values[istart_s:iend_s] = MatchUpData.values[istart:iend][s_idx]

            # ii. Sample values uncertainty data

            if block_unc.form == 'ave':
                # If block error correlation form was defined by a w matrix
                # - now simplified to random error correlation by sampling choice

                # Initialise uncertainty data array
                MatchUpSample.unc[i] = Uncertainty("r", zeros(len(s_idx)))

                # Retrieve required W matrix and uncertainty vector and hence determine new random uncertainties
                w = MatchUpData.w_matrices[block_unc.w_i]
                u = MatchUpData.uncertainty_vectors[block_unc.u_i]
                for j, s_i in enumerate(s_idx):
                    col_start = w.indices[w.indptr[j]]
                    col_end = w.indices[w.indptr[j+1]-1]+1
                    MatchUpSample.unc[i].uR[j] = npsum(u[col_start:col_end]**2) ** 0.5

            else:
                # If block error correlation form random or random and systematic simplify to random and sample
                MatchUpSample.unc[i] = Uncertainty("r", deepcopy(block_unc.uR[s_idx]))
        # --------------------------------------------------------------------------------------------------------------

        # d. sample k --------------------------------------------------------------------------------------------------

        # Initialise k data arrays
        MatchUpSample.ks = zeros(MatchUpSample.idx['cNm'][-1])
        MatchUpSample.unck = [0] * len(MatchUpSample.idx['Nm'])

        # Sample k and respective uncertainty data by match-up series
        for i, mu_unck in enumerate(MatchUpData.unck):
            istart = mc[i]                           # start of full dataset k data block
            iend = mc[i+1]                           # end of full dataset k data block
            istart_s = MatchUpSample.idx['cNm'][i]    # start of sampled dataset k data block
            iend_s = MatchUpSample.idx['cNm'][i+1]    # end of sampled dataset k data block
            s_idx = sampling_idxs[i+1]               # indices of match up to sample within k data block

            # i. Sample data
            MatchUpSample.ks[istart_s:iend_s] = MatchUpData.ks[istart:iend][s_idx]

            # ii. Sample uncertainties
            MatchUpSample.unck[i] = Uncertainty("r", deepcopy(mu_unck.uR[s_idx]))
        # --------------------------------------------------------------------------------------------------------------

        # d. sample times ----------------------------------------------------------------------------------------------

        # todo - write sampling of times

        # --------------------------------------------------------------------------------------------------------------

        # e. sample additional variables -------------------------------------------------------------------------------

        # todo - write sampling of additional variables

        # --------------------------------------------------------------------------------------------------------------

        return MatchUpSample

if __name__ == "__main__":
    pass
