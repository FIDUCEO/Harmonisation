"""
Contains the functionality to run the NPL Full EIV harmonisation method
"""

'''___Python Modules___'''
from copy import deepcopy
from time import time
from numpy.linalg import cholesky

'''___Third Party Modules___'''

'''___harmonisation Modules___'''
from harmonisation import Sample2Ind, Transform2NormInd, HarmonisationResult, HarmonisationEIVPC
from GN_algo import GNAlgo


'''___Authorship___'''
__author__ = ["Sam Hunt", "Peter Harris"]
__created__ = "16/04/2017"
__credits__ = ["Arta Dillo", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class HarmonisationEIV:
    """
    Class to run the NPL Full EIV harmonisation method.

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

    def run(self, HData,
                  pc_input=None, save_pc=None, gn_input=None, save_gn=None, save_directory_GNOp=None,
                  show=1):
        """
        Return harmonised parameters and diagnostic data for input harmonisaton match-up data

        :type HData: harm_data_writer.HarmInputData
        :param HData: Input match-up data object



        :return:
            :a: *numpy.ndarray*

            Harmonised parameters
        """

        Transform2NormIndOp = Transform2NormInd()

        ################################################################################################################
        # 1.	Compute Approximate Solution to find Pre-conditioner to Full Problem
        ################################################################################################################

        return_covariance = True
        if save_directory_GNOp is not None:
            return_covariance = False

        if pc_input is None:
            if show != 0:
                print "- Determine approximate solution to find pre-conditioner to full problem..."

            t1 = time()
            # a. sample data for preconditioning
            Sample2IndOp = Sample2Ind()
            HData_sample = Sample2IndOp.run(HData, sf=0.1, show=(show==2))

            HData_sample = Transform2NormIndOp.run(HData_sample)

            # b. determine preconditioner solution
            print "Beginning Solver..."
            PCOp = HarmonisationEIVPC()
            preconditioner = PCOp.run(HData_sample)

            del HData_sample
            t2 = time()
            print "t_PC:", str(t2-t1)

        else:
            preconditioner = HarmonisationResult(pc_input)

        HData.a = preconditioner.parameter  # set PC output parameters as current parameter estimates

        if save_pc:
            preconditioner.save(save_pc, save_residuals=False)

        ################################################################################################################
        # 2.	Compute Full Solution using EIV Gauss-Newton Algorithm
        ################################################################################################################

        if show != 0:
            print "Computing full solution..."

        if gn_input is None:
            if show != 0:
                print " - Transforming to Independent Variables..."
            # a. reparameterise input data such that output data are independent quantities
            HData = Transform2NormIndOp.run(HData)

            # b. run GN algorithm on modified data
            GNOp = GNAlgo(HData, preconditioner.parameter_covariance_matrix)
        else:
            if show != 0:
                print " - Opening Transformed Independent Variables..."
            GNOp = GNAlgo(HData)
            GNOp.open(gn_input)
            if show != 0:
                print " - Applying approximate solution to pre-conditioner to full problem..."
            GNOp.S = cholesky(preconditioner.parameter_covariance_matrix)

        HarmonisationOutput = HarmonisationResult()
        HarmonisationOutput.parameter_sensors = HData.idx["parameter_sensor"]
        HarmonisationOutput.idx = deepcopy(HData._original_idx)

        HarmonisationOutput.parameter, HarmonisationOutput.parameter_covariance_matrix, \
            HarmonisationOutput.cost, HarmonisationOutput.cost_dof, \
                HarmonisationOutput.cost_p_value, HarmonisationOutput.values_res,\
                    HarmonisationOutput.ks_res, systematic_errors, systematic_error_sensors \
                        = GNOp.run(show=(show == 2), return_covariance=return_covariance)

        if systematic_errors is not None:
            HarmonisationOutput.additional_variables['systematic_errors'] = {'data': systematic_errors,
                                                                             'dim': 's',
                                                                             'dtype': 'f4',
                                                                             'Description': 'Fitted systematic errors'}

            HarmonisationOutput.additional_variables['systematic_error_sensors'] = {'data': systematic_error_sensors,
                                                                                    'dim': 's',
                                                                                    'dtype': 'S1',
                                                                                    'Description': 'Sensors for fitted'
                                                                                                   'systematic errors'}

        if save_directory_GNOp is not None:
            GNOp.save(save_directory_GNOp)

        return HarmonisationOutput

if __name__ == "__main__":
    pass
