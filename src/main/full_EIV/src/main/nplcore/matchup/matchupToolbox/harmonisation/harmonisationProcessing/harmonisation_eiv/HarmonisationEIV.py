"""
Contains the functionality to run the NPL Full EIV harmonisation method
"""

'''___Python Modules___'''
import sys
from os.path import dirname
from os import makedirs
from os.path import join as pjoin
from copy import deepcopy
import cPickle
from time import time

'''___Third Party Modules___'''
from netCDF4 import Dataset

'''___Harmonisation Modules___'''
matchupProcessing_directory = pjoin(dirname(dirname(dirname(dirname(dirname(__file__))))), "matchupProcessing")
sys.path.append(matchupProcessing_directory)
from sample2ind.Sample2Ind import Sample2Ind
from transform2normind.Transform2NormInd import Transform2NormInd

harmonisationIO_directory = pjoin(dirname(dirname(dirname(__file__))), "harmonisationIO")
sys.path.append(harmonisationIO_directory)
from HarmonisationResult import HarmonisationResult

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

    def run(self, HData, show=True, S_PC=None, a_PC=None, open_directory_GNOp=None, save_directory_GNOp=None):
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
        sensor_data = HData.sensor_data
        return_covariance = True
        if save_directory_GNOp is not None:
            return_covariance = False

        if (S_PC is None) and (a_PC is None):
            if show:
                print "- Determine approximate solution to find pre-conditioner to full problem..."

            t1 = time()
            # a. sample data for preconditioning
            Sample2IndOp = Sample2Ind()
            HData_sample = Sample2IndOp.run(HData, sf=0.1, show=show)

            HData_sample = Transform2NormIndOp.run(HData_sample)

            # b. determine preconditioner solution
            print "Beginning Solver..."
            PC = PCAlgo(HData_sample)

            a_PC, S_PC = PC.runPC(tol=1e-6)
            del HData_sample
            t2 = time()
            print "t_PC:", str(t2-t1)

        HData.a = a_PC  # set PC output parameters as current parameter estimates

        ################################################################################################################
        # 2.	Compute Full Solution using EIV Gauss-Newton Algorithm
        ################################################################################################################

        if show:
            print "Computing full solution..."

        if open_directory_GNOp is None:
            if show:
                print " - Transforming to Independent Variables..."
            # a. reparameterise input data such that output data are independent quantities
            HData = Transform2NormIndOp.run(HData)

            # b. run GN algorithm on modified data
            GNOp = GNAlgo(HData, S_PC)
        else:
            if show:
                print " - Opening Transformed Independent Variables..."
            GNOp = GNAlgo(HData)
            GNOp.open(open_directory_GNOp)
            if show:
                print " - Applying approximate solution to pre-conditioner to full problem..."
            GNOp.S = S_PC

        HarmonisationOutput = HarmonisationResult()
        HarmonisationOutput.parameter_sensors = HData.idx["parameter_sensor"]
        HarmonisationOutput.idx = deepcopy(HData._original_idx)

        HarmonisationOutput.parameter, HarmonisationOutput.parameter_covariance_matrix, \
            HarmonisationOutput.cost, HarmonisationOutput.cost_dof, \
                HarmonisationOutput.cost_p_value, HarmonisationOutput.values_res,\
                    HarmonisationOutput.ks_res, systematic_errors, systematic_error_sensors \
                        = GNOp.run(show=show, return_covariance=return_covariance)

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
