"""
Contains the functionality to run the NPL Full EIV harmonisation method.
"""

'''___Python Modules___'''

'''___Third Party Modules___'''

'''___Harmonisation Modules___'''
import convert_data
from pc_algo import PCAlgo
from GN_algo import GNAlgo

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "16/04/2017"
__credits__ = ["Arta Dillo", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class HarmAlgo:
    """
    Class to run the NPL Full EIV harmonisation method.

    Sample code:

    .. code-block:: python

        H = HarmAlgo(HData)
        H.run()

    where ``HData`` is a ``harm_data_reader.HarmData`` object, containing the match-up data to be harmonised

    :Attributes:
        ..py:attribute:: HData:

        *harm_data_reader.HarmData*

        Input harmonisation data object containing match-up data to be harmonised

        ..py:attribute:: convert_data:

        *obj*

        Object containing functionality to convert *harm_data_reader.HarmData* objects

    :Methods:
        .. py:method:: run(...):

            Return harmonised parameters and diagnostic data for input harmonisaton match-up data
    """

    def __init__(self, HData):
        """
        Initialise HarmAlgo class

        :type HData: harm_data_reader.HarmData
        :param HData: Input harmonisation data object containing match-up data to be harmonised
        """

        # Initialise class
        self.convert_data = convert_data.ConvertData()
        self.HData = HData

    def run(self):
        """
        Return harmonised parameters and diagnostic data for input harmonisaton match-up data

        :type HData: harm_data_writer.HarmData
        :param HData: Input harmonisation match-up data object

        :return:
            :a: *numpy.ndarray*

            Harmonised parameters

            :Ia: *numpy.ndarray*

            Harmonised parameter sensor names

            :V: *numpy.ndarray*

            Harmonised parameters convariance matrix

            :F: *float*

            Objective function final value

            :v: *float*

            Objective function degrees of freedom

            :p: *float*

            Chi-squared probability
        """

        ################################################################################################################
        # 1.	Prepare Data
        ################################################################################################################

        HData = self.HData

        # Flatten values into required 1d form
        HData.values = HData.flatten_values(HData.values, HData.idx)

        ################################################################################################################
        # 2.	Compute Approximate Solution to find Pre-conditioner to Full Problem
        ################################################################################################################

        print("Determine approximate solution to find pre-conditioner to full problem...")

        # a. sample data for preconditioning
        HData_sample = self.convert_data.sample4PC(HData, sf=1)

        # b. determine preconditioner solution
        PC = PCAlgo(HData_sample)
        a_PC, S = PC.runPC(tol=1e-6)
        HData.a = a_PC  # set PC output parameters as current parameter estimates

        ################################################################################################################
        # 3.	Compute Full Solution using Gauss-Newton Algorithm
        ################################################################################################################

        print("Computing full solution...")

        # a. reparameterise input data such that output data are independent quantities
        HData = self.convert_data.convert2ind(HData)

        # b. run GN algorithm on modified data
        GN = GNAlgo(HData, S)
        a, V, F, v, p, H_res, K_res = GN.runGN(show=True)

        return a, V, F, v, p, H_res, K_res

if __name__ == "__main__":
    pass

