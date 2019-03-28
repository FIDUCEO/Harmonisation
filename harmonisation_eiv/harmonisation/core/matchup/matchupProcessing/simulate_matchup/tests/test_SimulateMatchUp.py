"""
Tests for SimulateMatchUp
"""

'''___Built-In Modules___'''
import unittest

'''___Third-Party Modules___'''
from numpy import array, append
from scipy.sparse import csr_matrix

'''___harmonisation Modules___'''
from harmonisation import MatchUp, SimulateMatchUp


'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "23/11/2017"
__credits__ = ["Peter Harris"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


class TestSimulateMatchUp(unittest.TestCase):
    def test_space_fill_sim(self):

        # input
        x = array([2.0, 4.0, 1.0, 5.0, 3.0])
        n_sim = 5

        # expected output
        x_sim = array([1.0, 2.0, 3.0, 4.0, 5.0])

        # run test
        SimulateMatchUpOp = SimulateMatchUp()
        x_sim_test = SimulateMatchUpOp.space_fill_sim(x, n_sim)

        # check results
        for x_sim_test_i, x_sim_i in zip(x_sim_test, x_sim):
            self.assertEquals(x_sim_test_i, x_sim_i)

    def test_return_simulation_x(self):

        # input
        x = array([1.0, 5.0, 10.0, 0.0])
        idx = [0, 2, 4]
        idx_sim = [0, 5, 10]

        # expected output
        SimulateMatchUpOp = SimulateMatchUp()
        x_sim = append(SimulateMatchUpOp.space_fill_sim(x[idx[0]:idx[1]], idx_sim[1]-idx_sim[0]),
                       SimulateMatchUpOp.space_fill_sim(x[idx[1]:idx[2]], idx_sim[2]-idx_sim[1]))

        # run test
        x_sim_test = SimulateMatchUpOp.return_simulation_x(x, idx, idx_sim)

        # check results
        for x_sim_test_i, x_sim_i in zip(x_sim_test, x_sim):
            self.assertEquals(x_sim_test_i, x_sim_i)

    def test_return_simulation_idx(self):

        # input
        idx = {'N_var': [7, 7, 5, 5],
               'Nm': [7, 5],
               'cNm': [0, 7, 12],
               'idx': [0, 7, 14, 19, 24],
               'another': [1, 2, 3, 4]}
        n_sim = 2

        # expected output
        idx_sim = {'N_var': [2, 2, 2, 2],
                   'Nm': [2, 2],
                   'cNm': [0, 2, 4],
                   'idx': [0, 2, 4, 6, 8],
                   'another': [1, 2, 3, 4]}

        # run test
        SimulateMatchUpOp = SimulateMatchUp()
        idx_sim_test = SimulateMatchUpOp.return_simulation_idx(idx, n_sim)

        # check results
        for key in idx_sim.keys():
            self.assertSequenceEqual(idx_sim[key], idx_sim_test[key], key)

    def test_run(self):

        # input
        MatchUpData = MatchUp()
        MatchUpData.values = array([1.0, 4.0, 5.0, 6.0, 4.0, 7.0, 3.0, 3.0, 2.0, 1.0])
        MatchUpData.ks = array([1.0, 4.0, 3.0, 3.0, 3.0])
        MatchUpData.idx = {'N_var': [2, 2, 3, 3],
                           'Nm': [2, 3],
                           'cNm': [0, 2, 5],
                           'idx': [0, 2, 4, 7, 10],
                           'another': [1, 2, 3, 4]}
        MatchUpData.a = 'test'
        MatchUpData.sensor_model = 'test'
        MatchUpData.sensor_model_constant = 'test'
        MatchUpData.adjustment_model = 'test'
        MatchUpData._original_idx = 'test'

        n_samples_mu = 5

        # expected output
        SimulateMatchUpOp = SimulateMatchUp()
        idx_sim = SimulateMatchUpOp.return_simulation_idx(MatchUpData.idx, n_samples_mu)
        values_sim = SimulateMatchUpOp.return_simulation_x(MatchUpData.values,
                                                           MatchUpData.idx['idx'], idx_sim['idx'])
        ks_sim = SimulateMatchUpOp.return_simulation_x(MatchUpData.ks,
                                                           MatchUpData.idx['cNm'], idx_sim['cNm'])
        a_sim = MatchUpData.a
        sensor_model_sim = MatchUpData.sensor_model
        sensor_model_constant_sim = MatchUpData.sensor_model_constant
        adjustment_model_sim = MatchUpData.adjustment_model
        _original_idx_sim = MatchUpData._original_idx

        # run test
        MatchUpSimulation = SimulateMatchUpOp.run(MatchUpData, n_samples_mu)

        # check output

        for x_sim_test_i, x_sim_i in zip(values_sim, MatchUpSimulation.values):
            self.assertEquals(x_sim_test_i, x_sim_i)

        for x_sim_test_i, x_sim_i in zip(ks_sim, MatchUpSimulation.ks):
            self.assertEquals(x_sim_test_i, x_sim_i)

        for key in idx_sim.keys():
            self.assertSequenceEqual(idx_sim[key], MatchUpSimulation.idx[key], key)

        self.assertEqual(MatchUpSimulation.a, a_sim)
        self.assertEqual(MatchUpSimulation.sensor_model, sensor_model_sim)
        self.assertEqual(MatchUpSimulation.sensor_model_constant, sensor_model_constant_sim)
        self.assertEqual(MatchUpSimulation.adjustment_model, adjustment_model_sim)
        self.assertEqual(MatchUpSimulation._original_idx, _original_idx_sim)

if __name__ == "__main__":
    unittest.main()