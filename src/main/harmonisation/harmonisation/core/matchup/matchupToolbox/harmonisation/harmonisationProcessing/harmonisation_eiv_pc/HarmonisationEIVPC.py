"""
Contains the functionality to run the NPL Full EIV harmonisation method preconditioner
"""

'''___Built-In Modules___'''
from copy import deepcopy

'''___Third Party Modules___'''


'''___harmonisation Modules___'''
from harmonisation import HarmonisationResult
from pc_algo import PCAlgo


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "22/05/2018"
__credits__ = []
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class HarmonisationEIVPC:
    """
    Class to run the NPL Full EIV harmonisation method preconditioner

    Sample code:

    .. code-block:: python

        H = HarmAlgo(HData)
        H.run()

    where ``HData`` is a ``harm_data_reader.HarmInputData`` object, containing the match-up data to be harmonised

    :Attributes:
        ..py:attribute:: HData:

            *harm_data_reader.HarmInputData*

            Input harmonisation data object containing match-up data to be harmonised

    :Methods:
        .. py:method:: run(...):

            Return harmonised parameters and diagnostic data for input harmonisaton match-up data
    """

    def __init__(self):
        pass

    def run(self, HData, show=True, return_covariance=True):

        PCAlgoOp = PCAlgo(HData)

        HarmonisationOutput = HarmonisationResult()
        HarmonisationOutput.parameter_sensors = HData.idx["parameter_sensor"]
        HarmonisationOutput.idx = deepcopy(HData._original_idx)

        HarmonisationOutput.parameter, HarmonisationOutput.parameter_covariance_matrix, \
            HarmonisationOutput.cost, HarmonisationOutput.cost_dof, \
                HarmonisationOutput.cost_p_value, HarmonisationOutput.values_res, \
                    HarmonisationOutput.ks_res = PCAlgoOp.runPC(tol=1e-6)

        return HarmonisationOutput

if __name__ == "__main__":
    pass
