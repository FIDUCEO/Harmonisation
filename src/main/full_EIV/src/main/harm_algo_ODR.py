"""
Created on Sun April 16  2017 11:00:00

@author: Arta Dillo, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''

'''___Harmonisation Modules___'''


class HarmAlgo:
    """
    Class to run the NPL ODR harmonisation method.

    Sample code:

    H = HarmAlgo(HData)
    H.run()

    where HData is a harm_data_reader.HarmData object, containing the match-up data to be harmonised
    """

    def __init__(self, HData):
        """
        Initialise HarmAlgo class

        :param HData: harm_data_reader.HarmData
            Input harmonisation data object containing match-up data to be harmonised
        """

        self.HData = HData

    def run(self):
        """
        Return harmonised parameters and diagnostic data for input harmonisaton match-up data

        :param HData: harm_data_writer.HarmData
                Input harmonisation match-up data object

        :return:

            :a: numpy.ndarray
                harmonised parameters

            :Ia: numpy.ndarray
                harmonised parameter sensor names

            :V: numpy.ndarray
                harmonised parameters convariance matrix

            :F: float
                objective function final value

            :v: float
                objective function degrees of freedom

            :p: float
                chi-squared probability
        """

        HData = self.HData

        a = None
        V = None
        F = None
        v = None
        p = None
        H_res = None
        K_res = None

        return a, V, F, v, p, H_res, K_res

if __name__ == "__main__":

    def main():
        return 0

    main()