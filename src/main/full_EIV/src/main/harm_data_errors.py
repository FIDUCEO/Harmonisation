"""
Created on Tues April 18  2017 11:00:00
@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
from numpy import zeros
from numpy.random import normal


def gen_errors(HData):
    """
    Return input HData object with values and ks adjusted by errors respecting its uncertainty structure

    :param HData: harm_data_reader.HarmData
        input harmonisation match-up data object

    :return:
        :HData: harm_data_reader.HarmData
            Input HData object with values and ks adjusted with errors respecting its uncertainty structure
    """

    for i, (block_unc, cov, mu) in enumerate(zip(HData.unc, HData.idx['n_cov'], HData.idx['n_mu'])):

        # indices defining first and last positions in data matrix
        istart = HData.idx['cNm'][mu - 1]
        iend = HData.idx['cNm'][mu]

        # index defining column in data matrix
        if HData.idx['Im'][mu - 1][0] == HData.idx['n_sensor'][i]:  # if the sensor is the first sensor in the match-up series
            col = cov - 1
        if HData.idx['Im'][mu - 1][1] == HData.idx['n_sensor'][i]:  # if the sensor is the second sensor in the match-up series
            col = cov + HData.values.shape[1]/2 - 1

        if block_unc.form == 'r':
            HData.values[istart:iend, col] = normal(loc=HData.values[istart:iend, col], scale=block_unc.uR)

        elif block_unc.form == 'rs':
            # todo - write harmonisation
            print 'shouldn''t be here'

        elif block_unc.form == 'ave':
            HData.values[istart:iend, col] += block_unc.W.dot(normal(loc=zeros(block_unc.W.indices[-1] + 1)))

    for i in xrange(len(HData.idx['Im'])):
        istart = HData.idx['cNm'][i]
        iend = HData.idx['cNm'][i + 1]

        HData.ks[istart:iend] = normal(loc=HData.ks[istart:iend], scale=HData.unck[i].uR)

    return HData

if __name__ == "__main__":

    def main():
        return 0

    main()
