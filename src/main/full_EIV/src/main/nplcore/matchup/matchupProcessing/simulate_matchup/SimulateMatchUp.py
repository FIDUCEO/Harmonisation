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
from numpy import zeros, linspace, float32

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


class SimulateMatchUp:
    """
    SimulateMatchUp class operates on instances of ``eopy.matchup.matchupIO.MatchUp`` in-memory data products,
    returning a space filling simulation of the match-ups based on the data ranges an input dataset

    :Methods:
        .. py:method:: run(...):

            Returns simulated instance of ``eopy.matchup.matchupIO.MatchUp``
    """

    def __init__(self):
        pass

    def run(self, MatchUpData, n_samples_mu=20, verbose=False):
        """
        Returning a space filling simulation of the match-ups based on the data ranges an input dataset

        :type MatchUpData: *eopy.matchup.matchupIO.MatchUp*
        :param MatchUpData: Input match-up data

        :type n_samples_mu: int
        :param n_samples_mu: Number of simulation samples per matchup series

        :type

        :return:
            :MatchUpData: *eopy.matchup.matchupIO.MatchUp*

            Simulated matchup data
        """

        # Initialise transformed data product
        MatchUpSimulation = MatchUp()
        MatchUpSimulation.idx = self.return_simulation_idx(MatchUpData.idx, n_samples_mu)
        MatchUpSimulation.values = self.return_simulation_x(MatchUpData.values, MatchUpData.idx['idx'],
                                                            MatchUpSimulation.idx['idx'])
        MatchUpSimulation.ks = self.return_simulation_x(MatchUpData.ks, MatchUpData.idx['cNm'],
                                                        MatchUpSimulation.idx['cNm'])
        MatchUpSimulation.a = MatchUpData.a[:]
        MatchUpSimulation.sensor_model = MatchUpData.sensor_model
        MatchUpSimulation.sensor_model_constant = MatchUpData.sensor_model_constant
        MatchUpSimulation.adjustment_model = MatchUpData.adjustment_model
        MatchUpSimulation._original_idx = MatchUpData._original_idx

        if verbose:
            sensors = MatchUpSimulation.idx['sensors']
            plotted_sensors = []

            i_mu = 0
            while set(plotted_sensors) != set(sensors):

                pair = MatchUpSimulation.idx['Im'][i_mu]
                for i_pair, n_sensor in enumerate(pair):
                    sensor_name = sensors[n_sensor]
                    print sensor_name
                    print "Variable\tMin\t\t\tMax"
                    for cov in range(1, MatchUpSimulation.idx['sensor_ms'][n_sensor]+1):
                        maximum = max(MatchUpSimulation.getVariableData(cov, sensor_name, i_mu+1))
                        minimum = min(MatchUpSimulation.getVariableData(cov, sensor_name, i_mu+1))
                        print str(cov) + "\t\t\t" + str(minimum) + "\t\t" + str(maximum)
                    print "\n"
                    plotted_sensors.append(sensor_name)
                i_mu += 1

        # todo - improve simplifications of simulation, does not handle unc/k, w/u_matrices, across/along_track_index1/2
        return MatchUpSimulation

    def return_simulation_idx(self, idx, n_samples_mu):
        # Initialise copy of harmonisation idx to update for converted data
        # N.B. deepcopy() required to ensure copy of nested lists in dict
        simulation_idx = deepcopy(idx)

        simulation_idx['N_var'] = [int(n_samples_mu) for x in simulation_idx['N_var']]
        simulation_idx['Nm'] = [int(n_samples_mu) for x in simulation_idx['Nm']]

        idxs = [0]
        total = 0
        for N in simulation_idx['N_var']:
            total += N
            idxs.append(int(total))

        simulation_idx['idx'] = idxs

        cNm = [0]
        total = 0
        for N in simulation_idx['Nm']:
            total += N
            cNm.append(int(total))

        simulation_idx['cNm'] = cNm
        return simulation_idx

    def return_simulation_x(self, x, idx, simulation_idx):

        sim = self.space_fill_sim

        simulation_x = zeros(simulation_idx[-1], dtype=float32)
        for i, (istart, iend, istart_sim, iend_sim) in enumerate(zip(idx[:-1], idx[1:], simulation_idx[:-1], simulation_idx[1:])):
            simulation_x[istart_sim:iend_sim] = sim(x[istart:iend], iend_sim-istart_sim)

        return simulation_x

    def space_fill_sim(self, values, n_sim):
        values_sim = linspace(min(values), max(values), n_sim)
        return values_sim
