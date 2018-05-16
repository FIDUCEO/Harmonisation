"""
Module containing functions to evaluate MatchUp data
"""

'''___Python Modules___'''

'''___Third Party Modules___'''
from numpy import zeros, asarray, float32, float64, array, outer, diag

'''___Harmonisation Modules___'''

'''___Authorship___'''
__author__ = "Sam Hunt"
__created__ = "21/11/2017"
__credits__ = ["Peter Harris", "Jon Mittaz"]
__version__ = "0.0"
__maintainer__ = "Sam Hunt"
__email__ = "sam.hunt@npl.co.uk"
__status__ = "Development"


def evaluate_measurand(MatchUpData, matchup_series=0, sensor_name=0, return_derivatives=False,
                       parameter_covariance_matrix=None):
    """
    Return evaluated sensor function measurand for match-up dataset

    :type MatchUpData: eopy.matchup.matchupIO.MatchUp
    :param MatchUpData: Match-up data to evaluate

    :type matchup_series: float/list
    :param matchup_series: Selected match-up series to compute the measurand for
    NB: Default 0 - all series

    :type return_derivatives: bool
    :param return_derivatives: switch to determine whether to return derivatives of measurands

    :return:
        :measurand: *numpy.ndarray*

        Evaluated sensor model measurand

        :measurand_derivatives: *numpy.ndarray*

        Derivatives of sensor model measurand
        NB: Only returned if return_derivatives is True
    """

    ####################################################################################################################
    # 1. Determine match-up series to retrieve measurand for
    ####################################################################################################################

    # a. Build list of measurement series to retreive measurand for
    if matchup_series == 0:
        matchup_series = list(set(MatchUpData.idx['n_mu']))
    elif type(matchup_series) == int:
        matchup_series = [matchup_series]

    # b. Determine idxs of selected measurands
    selection_Nm = [m for i, m in enumerate(MatchUpData.idx['Nm']) if i+1 in matchup_series]
    selection_cNm = [0]
    total = 0
    for m in selection_Nm:
        total += m
        selection_cNm.append(total)

    ####################################################################################################################
    # 2. Determine required measurand
    ####################################################################################################################

    # Initialise output arrays
    measurand = zeros((selection_cNm[-1], 2), dtype=float32)
    measurand_derivatives = None
    if return_derivatives:
        measurand_derivatives = [[zeros((selection_cNm[-1],
                                         MatchUpData.idx['sensor_ms'][Im[0]]+MatchUpData.idx['sensor_a'][Im[0]]),
                                        dtype=float32),
                                  zeros((selection_cNm[-1],
                                         MatchUpData.idx['sensor_ms'][Im[1]]+MatchUpData.idx['sensor_a'][Im[1]]),
                                        dtype=float32)]
                                  for Im in MatchUpData.idx['Im']]

    u_measurand = None
    if parameter_covariance_matrix is not None:
        u_measurand = zeros((selection_cNm[-1], 2), dtype=float32)

    # Compute measurand per selected match-up series
    block_idxs = [(j, k) for j, k in zip(MatchUpData.idx['n_sensor'], MatchUpData.idx['n_mu'])]
    for i_sel_matchup, n_mu in enumerate(matchup_series):

        # Find sensor in selected match-up
        sensor_pair = MatchUpData.idx['Im'][n_mu - 1]

        # Row location of measurand values in measurand array
        i_start_m = selection_cNm[i_sel_matchup]
        i_end_m = selection_cNm[i_sel_matchup+1]

        # Compute measurand per sensor of selected match-up series
        for i_sensor_pair, n_sensor in enumerate(sensor_pair):

            sensor_name = MatchUpData.idx["sensors"][n_sensor]

            # Find location of covariate data in 1D data structure
            cov_block_idxs = [i for i, block_idx in enumerate(block_idxs) if block_idx == (n_sensor, n_mu)]
            cov_idxs = [(MatchUpData.idx['idx'][cov_block_idx], MatchUpData.idx['idx'][cov_block_idx+1])
                       for cov_block_idx in cov_block_idxs]

            # Compute measurand and store in measurand array
            sensor_a = asarray([a for s, a in zip(MatchUpData.idx['parameter_sensor'], MatchUpData.a) if s == sensor_name], dtype=float32)
            sensor_c = asarray([c for s, c in zip(MatchUpData.idx['sensor_model_constant_sensor'], MatchUpData.sensor_model_constant) if s == sensor_name], dtype=float64)

            # Select Sensor State Variables
            sensor_X = zeros((cov_idxs[0][1] - cov_idxs[0][0], len(cov_idxs)), dtype=float32)
            for i_cov_idx, cov_idx in enumerate(cov_idxs):
                sensor_X[:, i_cov_idx] = MatchUpData.values[cov_idx[0]:cov_idx[1]]

            # Select Sensor Indices
            sensor_xt_i = MatchUpData.getSensorAcrossTrackIndex(n_mu, sensor_name)
            sensor_at_i = MatchUpData.getSensorAlongTrackIndex(n_mu, sensor_name)

            # Select Sensor Match-up Time
            sensor_t = MatchUpData.getSensorTime(n_mu, sensor_name)

            measurand[i_start_m:i_end_m, i_sensor_pair], J \
                = MatchUpData.sensor_model[n_sensor](sensor_X, sensor_a,
                                                     sensor_c, sensor_xt_i, sensor_at_i, sensor_t)

            if return_derivatives:
                measurand_derivatives[n_mu - 1][i_sensor_pair] = J

            if u_measurand is not None:
                if J is not None:
                    m = MatchUpData.idx['sensor_ms'][n_sensor]

                    loc_a = array([True if s == sensor_name else False for s in MatchUpData.idx['parameter_sensor']])
                    loc_a_arr = outer(loc_a.T, loc_a)
                    V = parameter_covariance_matrix[loc_a_arr]
                    V = V.reshape((len(V)**0.5, len(V)**0.5))

                    u_measurand[i_start_m:i_end_m, i_sensor_pair] = (diag(J[:, m:].dot(V).dot(J[:, m:].T))) ** 0.5

                    pass

    if return_derivatives:
        return measurand, measurand_derivatives

    if u_measurand is not None:
        if return_derivatives:
            return measurand, measurand_derivatives, u_measurand
        return measurand, u_measurand

    return measurand


def evaluate_adjusted_measurand(MatchUpData, matchup_series=0, return_derivatives=False, return_measurand_derivatives=False):
    """
    Return evaluated match-up adjusted measurand following adjustment model for match-up dataset

    :type MatchUpData: eopy.matchup.matchupIO.MatchUp
    :param MatchUpData: Match-up dataset to evaluate

    :type matchup_series: float/list
    :param matchup_series: Selected match-up series to compute the measurand for
    NB: Default 0 - all series

    :type return_derivatives: bool
    :param return_derivatives: switch to determine whether to return derivatives

    :return:
        :adjusted_measurand_derivatives: *numpy.ndarray*

        Evaluated match-up adjustment model

        :adjusted_measurand_derivatives:

        Derivatives of match-up adjustment model
        NB: Only returned if return_derivative is True
    """

    ####################################################################################################################
    # 1. Determine match-up series to retrieve measurand for
    ####################################################################################################################

    # a. Build list of measurement series to retreive measurand for
    if matchup_series == 0:
        matchup_series = list(set(MatchUpData.idx['n_mu']))
    elif type(matchup_series) == int:
        matchup_series = [matchup_series]

    # b. Determine idxs of selected measurands
    selection_Nm = [m for i, m in enumerate(MatchUpData.idx['Nm']) if i+1 in matchup_series]
    selection_cNm = [0]
    total = 0
    for m in selection_Nm:
        total += m
        selection_cNm.append(total)

    ####################################################################################################################
    # 2. Determine required adjustments
    ####################################################################################################################

    # todo - would need updating for a sensor specific adjustment model
    # NB: This bit is deliberately over complicated to allow for sensor specific adjustment model in future

    # a. Initialise output arrays

    # Start with untransformed measurand
    if return_measurand_derivatives:
        adjusted_measurand, measurand_derivatives = evaluate_measurand(MatchUpData, matchup_series=matchup_series,
                                                                       return_derivatives=True)
    else:
        adjusted_measurand = evaluate_measurand(MatchUpData, matchup_series=matchup_series, return_derivatives=False)

    # Initialise derivatives array if required
    adjusted_measurand_derivatives = None
    if return_derivatives:
        adjusted_measurand_derivatives = zeros((selection_cNm[-1], 2))

    # Compute measurand per selected match-up series
    for i_sel_matchup, n_mu in enumerate(matchup_series):

        # Find sensor in selected match-up
        sensor_pair = MatchUpData.idx['Im'][n_mu - 1]

        # Row location of measurand values in measurand array
        i_start_m = selection_cNm[i_sel_matchup]
        i_end_m = selection_cNm[i_sel_matchup+1]

        # Compute measurand per sensor of selected match-up series
        for i_sensor_pair, n_sensor in enumerate(sensor_pair):

            # Compute adjusted measurand and store in adjusted measurand array
            if return_derivatives:
                adjusted_measurand[i_start_m:i_end_m, i_sensor_pair],\
                    adjusted_measurand_derivatives[i_start_m:i_end_m, i_sensor_pair] \
                        = MatchUpData.adjustment_model[n_sensor](adjusted_measurand[i_start_m:i_end_m, i_sensor_pair])
            else:
                adjusted_measurand[i_start_m:i_end_m, i_sensor_pair] \
                    = MatchUpData.adjustment_model[n_sensor](adjusted_measurand[i_start_m:i_end_m, i_sensor_pair])[0]

    if return_derivatives and return_measurand_derivatives:
        return adjusted_measurand, adjusted_measurand_derivatives, measurand_derivatives
    if return_derivatives:
        return adjusted_measurand, adjusted_measurand_derivatives
    return adjusted_measurand


def evaluate_K(MatchUpData, matchup_series=0):
    """
    Return match-up adjustment factor per match-up for match-up dataset

    :type MatchUpData: eopy.matchup.matchupIO.MatchUp
    :param MatchUpData: Match-up data to evaluate

    :type matchup_series: float/list
    :param matchup_series: Selected match-up series to compute the measurand for
    NB: Default 0 - all series

    :return:
        :K: *numpy.ndarray*

        Match-up adjustment factor
    """

    adjusted_measurand = evaluate_adjusted_measurand(MatchUpData,
                                                     matchup_series=matchup_series, return_derivatives=False)

    return adjusted_measurand[:, 1] - adjusted_measurand[:, 0]


if __name__ == "__main__":
    pass
