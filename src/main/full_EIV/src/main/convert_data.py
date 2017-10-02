"""
Created on Fri Jan 13  2017 09:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules____'''
from copy import deepcopy

from numpy import array, append, zeros
from numpy import sum as npsum

'''___Harmonisation Modules___'''
from harm_data_reader import HarmData
from correl_forms import CorrelForm


class ConvertData:
    """
    Class contains methods to convert input data into forms suitable for the pre-conditioner algorithm and the
    Gauss-Netwon harmonisation algorithm

    :Methods:
        :convert2ind:
            reparameterises input data such that output data are independent quantities, required for the Gauss-Newton
            and pre-conditioner algorithm
        :sample4PC:
            sample data so that the only remaining correlations arise from systematic effects, required for the
            pre-conditioner algorithm

    """

    def convert2ind(self, HData):
        """
        Return a reparameterisation of the input data such that output data are independent quantities (suitable for
        use in the pre-conditioner and Gauss-Newton algorithms)

        :param HData: HarmData
            Input harmonisation data for conversion

        :return:
            :HData: HarmData
                Reparameterised harmonisation data suitable for pre-conditioner and GN algorithm
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

        # initialise empty array to store converted harmonisation data
        new_values = array([])

        # initialise copy of harmonisation idx to update for converted data
        new_idx = deepcopy(HData.idx)      # deep copy required to ensure copy of nested lists in dict

        # counter to track block number through covariate, e.g. 2*N_mu-1 blocks per covariate
        # (used to add systematic uncertainty data to the end of the all blocks of a given
        #  covariate)
        n_b_cov = 0
        N_bpcov = 2*HData.idx['n_mu'][-1] - 1    # total number of blocks per covariate

        # convert data block by depending on correlation form
        for i, block_unc in enumerate(HData.unc):

            # block parameters
            ib = HData.idx['idx'][i]             # start idx of data block
            ie = HData.idx['idx'][i+1]           # end idx of data block
            nm = int(HData.idx['N_var'][i])      # number of variables in block

            # 1. random type correlation - scale by uncertainty
            if block_unc.form == "r":
                ur = block_unc.uR                # uncertainty
                # store old values in new array (scaling by uncertainty)
                new_values = append(new_values, HData.values[ib:ie]/ur)

            # 2. random+systematic type correlation - separate components
            elif block_unc.form == "rs":
                # a. random component - scale by uncertainty
                uR = block_unc.uR                # random uncertainty

                # store old values in new array (scaling by random uncertainty)
                new_values = append(new_values, HData.values[ib:ie]/uR)

                # b. systematic component - at final block of given covariate add
                #                           systemic data for all blocks
                n_b_cov += 1           # add to count of blocks through covariate

                # when block number is total number of blocks of covariate
                if n_b_cov == N_bpcov:
                    # store systematic uncertainty for each block
                    # todo - changes needed for Jon's data?
                    N_sensors = len(set(HData.idx['n_sensor']))-1

                    new_values = append(new_values, zeros(N_sensors))

                    # reset counter
                    new_idx['idx'][i+1:] = [j+N_sensors for j in new_idx['idx'][i+1:]]
                    n_b_cov = 0

            # 3. averaging type correlation - simulate data without averaging
            elif block_unc.form == "ave":

                # initialise array
                Htemp = zeros(block_unc.W.shape[1])
                uRtemp = zeros(block_unc.W.shape[1])

                next_idx = 0
                col = 0
                first_idx_prev = 0

                for j in xrange(HData.idx['N_var'][i]):

                    i_value = j + ib

                    first_idx = block_unc.W.indices[col]
                    step = first_idx - first_idx_prev

                    uR = block_unc.uR[j][block_unc.uR[j] != 0]                # valid scanline uncertainties
                    n_w = len(uR)                                             # size of averaging window
                    w = block_unc.W[j, first_idx:first_idx+n_w].toarray()[0]  # W row

                    if (step == n_w) or (j == 0):
                        # averaged values
                        Htemp[next_idx:next_idx+n_w] = HData.values[i_value]/uR
                        uRtemp[next_idx:next_idx+n_w] = uR
                        next_idx = next_idx + n_w

                    elif 0 < step < n_w:

                        istartk = next_idx + step - n_w
                        iendk = next_idx + step - 1

                        try:
                            uRtemp[next_idx:iendk+1] = uR[-(iendk+1-next_idx):]

                            # fill all but last missing data of average with averaged value
                            Htemp[next_idx:iendk] = HData.values[i_value] / uR[-(iendk+1-next_idx):-1]

                            # compute missing final value in average
                            Htemp[iendk] = (HData.values[i_value]-sum(Htemp[istartk:iendk] * w[:-1])) / w[-1]

                        except ValueError:
                            print 'i', i
                            print 'j', j
                            print 'step', step
                            print 'n_w', n_w
                            print 'next_idx', next_idx
                            print 'istartk', istartk
                            print 'iendk', iendk
                            print 'uR[-(iendk+1-next_idx):]', uR[-(iendk+1-next_idx):]
                            pass

                        next_idx += step

                    elif step == 0:
                        pass

                    first_idx_prev = first_idx
                    col += n_w

                #For testing
                #Before = HData.values[ib:ie]
                #After = block_unc.W.dot(Htemp)

                # store new data values
                new_values = append(new_values, Htemp)
                block_unc.uR = uRtemp

                # update N_var of HData.idx to count new variables
                new_len = len(Htemp)
                new_idx['N_var'][i] = new_len
                new_idx['idx'][i+1:] = [old_idx+new_len-(ie-ib) for old_idx in new_idx['idx'][i+1:]]

        # replace old data values array with converted form
        HData.values = new_values
        HData.idx = new_idx

        # scale ks
        for i in xrange(len(HData.idx['Im'])):
            istart = HData.idx['cNm'][i]
            iend = HData.idx['cNm'][i + 1]

            HData.ks[istart:iend] /= HData.unck[i].uR

        return HData

    def sample4PC(self, HData, sf):
        """
        Return sample of data for which the only data correlations arise from systematic effects

        :param HData: HarmData
            Input harmonisation data for conversion
        :param sf: int
            Sampling factor

        :return:
            :HData_sample: HarmData
                Sampled harmonisation data
        """

        # initialise parameters
        mcxyz = HData.idx['idx']       # cumulative total of variables data block
        mc = HData.idx['cNm']          # cumulative total of match-ups by series

        # initialise sampled harmonisation data product
        HData_sample = HarmData()
        HData_sample.idx = deepcopy(HData.idx)
        HData_sample.unc = deepcopy(HData.unc[:])
        HData_sample.unck = deepcopy(HData.unck[:])
        HData_sample.a = HData.a[:]
        HData_sample.sensor_model = HData.sensor_model
        HData_sample.adjustment_model = HData.adjustment_model

        ################################################################################################################
        # 1. Sample Data
        ################################################################################################################

        # a. find sampling indices

        n_mus = set(HData.idx['n_mu'])
        sampling_idxs = {}

        # find sampling indices per match-up series
        for n_mu in n_mus:

            # find W for covariate with largest moving average window (i.e. responsible for the most correlation)
            n_w = 0
            W = 0
            for i, block_unc in enumerate(HData.unc):
                if HData.idx['n_mu'][i] == n_mu:
                    if block_unc.form == 'ave':
                        if block_unc.uR.shape[1] > n_w:
                            n_w = block_unc.uR.shape[1]
                            W = block_unc.W

            # produce sampling indices
            stop = False
            istartW = 0
            last_idx = 0
            idx = 0
            idxs = [idx]
            while stop is False:

                for j, first_idx in enumerate(W.indices[istartW::n_w]):

                    step = first_idx - last_idx

                    current_idx = idx + j
                    final_idx = len(W.indices[::n_w]) - 1

                    if current_idx == final_idx:
                        sampling_idxs[n_mu] = idxs
                        stop = True
                        break

                    elif step >= n_w:
                        # averaged values
                        idx += j
                        idxs.append(idx)
                        last_idx = first_idx
                        istartW += j*n_w
                        break

        # b. sample variables

        # update idx attribute of HData_sample to describe structure of sampled data
        idxs = [0]
        total = 0
        for i, n_mu in enumerate(HData.idx['n_mu']):
            block_samples = len(sampling_idxs[n_mu])
            HData_sample.idx['N_var'][i] = block_samples
            total += block_samples
            idxs.append(int(total))
        HData_sample.idx['idx'] = idxs

        # sample variables and respective uncertainty data by data block
        HData_sample.values = zeros(HData_sample.idx['idx'][-1])

        for i, block_unc in enumerate(HData.unc):

            # block indices
            istart = mcxyz[i]
            iend = mcxyz[i+1]

            istart_s = HData_sample.idx['idx'][i]
            iend_s = HData_sample.idx['idx'][i+1]

            s_idx = sampling_idxs[HData.idx['n_mu'][i]]

            HData_sample.values[istart_s:iend_s] = HData.values[istart:iend][s_idx]

            if block_unc.form == "ave":
                HData_sample.unc[i] = CorrelForm("r", zeros(len(s_idx)))
                for j, s_i in enumerate(s_idx):
                    HData_sample.unc[i].uR[j] = npsum(block_unc.W[s_i, :].toarray()[0]**2) ** 0.5
            else:
                HData_sample.unc[i].uR = deepcopy(block_unc.uR[s_idx])

        # c. sample ks
        cNm = [0]
        total = 0
        for i, n_mu in enumerate(n_mus):
            n_mu_sample = len(sampling_idxs[n_mu])
            HData_sample.idx['Nm'][i] = n_mu_sample
            total += n_mu_sample
            cNm.append(total)
        HData_sample.idx['cNm'] = cNm

        print "Sample Size: ", HData_sample.idx['Nm']

        # sample k and respective uncertainty data by match-up series
        HData_sample.ks = zeros(HData_sample.idx['cNm'][-1])

        for i, mu_unck in enumerate(HData.unck):

            n_mu = i+1

            # match-up series indices
            istart = mc[i]
            iend = mc[i+1]

            istart_s = HData_sample.idx['cNm'][i]
            iend_s = HData_sample.idx['cNm'][i+1]

            s_idx = sampling_idxs[n_mu]

            # sample data
            HData_sample.ks[istart_s:iend_s] = HData.ks[istart:iend][s_idx]
            HData_sample.unck[i].uR = deepcopy(mu_unck.uR[s_idx])

        ################################################################################################################
        # 2. Convert to Independent Data
        ################################################################################################################

        HData_sample = self.convert2ind(HData_sample)

        return HData_sample

if __name__ == "__main__":

    def main():
        return 0

    main()
