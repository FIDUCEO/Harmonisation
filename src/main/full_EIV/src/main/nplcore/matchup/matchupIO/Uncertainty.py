"""
Uncertainty Class - Class to hold data uncertainty information
"""

'''___Built-In Modules___'''

'''___Third-Party Modules___'''

'''___EOPy Modules___'''

'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "18/11/2017"
__credits__ = ["Jonathan Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class Uncertainty:
    """
    Object containing information of a variables uncertainty values and form

    :Attributes:
        .. py:attribute:: form

            *str*

            Correlation form of uncertainty, can have values:

            * "r" - random form
            * "rs" - random + systematic form
            * "ave" - w matrix form

        + additional attribute/s containing uncertainty data depending on correlation form:

        * "r" form

            .. py:attribute:: uR

                *numpy.ndarray*

                Standard uncertainty per match-up caused by independent error effects

        * "rs" form

            .. py:attribute:: uR

                *numpy.ndarray*

                Standard uncertainty per match-up caused by independent error effects


            .. py:attribute:: uS

                *numpy.ndarray*

                Standard uncertainty per match-up caused by systematic error effects

        * "ave" form

            .. py:attribute:: w_i

                *int*

                w matrix index in ``eopy.matchup.matchupIO.MatchUp`` *w_matrices* attribute

            .. py:attribute:: u_i

                *int*

                u matrix index in ``eopy.matchup.matchupIO.MatchUp`` *u_matrices* attribute
    """

    def __init__(self, typeID, data_tuple):
        """
        Take user input covariate uncertainty correlation form and data and apply as attributes of the class

        :type form: str
        :param form: Correlation form, either "r", "rs" or "ave"

        :param data_tuple: tuple
        :param data_tuple: Uncertainty values - structure depends on form, as follows:

        * "r" - uR

            + uR(*float*) - standard uncertainty per match-up caused by independent error effects

        * "rs" - (uR, uS)

            + uR(*float*) - standard uncertainty per match-up caused by independent error effects
            + uS(*float*) - standard uncertainty per match-up caused by systematic error effects

        * "ave" - (uR, times, corrData)

            + w_i(*int*) - w matrix index in ``eopy.matchup.matchupIO.MatchUp`` *w_matrices* attribute
            + u_i(*int*) - u matrix index in ``eopy.matchup.matchupIO.MatchUp`` *u_matrices* attribute
        """

        # set form attribute
        self.typeID = typeID

        # set data value attributes from input data_tuple
        if typeID == 1:
            self.type = "independent"
            self.uR = data_tuple

        elif typeID == 2:
            self.type = "independent + systematic"
            self.uR = data_tuple[0]
            self.uS = data_tuple[1]

        elif typeID == 3:
            self.type = "structured"
            self.w_i = data_tuple[0]
            self.u_i = data_tuple[1]

        elif typeID == 4:
            self.type = "structured + systematic"
            self.w_i = data_tuple[0]
            self.u_i = data_tuple[1]
            self.uS = data_tuple[2]


if __name__ == "__main__":

    def main():
        return 0

    main()
