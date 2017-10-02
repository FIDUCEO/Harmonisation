"""
Functions required to calculate radiance (and sample radiance) residual

Created on Wed May 03 2017 09:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Harmonisation Modules___'''
from harm_data_reader import HarmData

'''___Python Modules___'''
from numpy import zeros, linspace, diag, sort, cumsum, append, asarray, outer, sum, amax, amin, ones, mean, std, array, arange, nan, isnan
from numpy.random import rand
from copy import deepcopy


def sample_H(HData, n_sample):
    """
    Sample harmonisation data

    :param H: harm_data_reader.HarmData
        harmonisation data object
    :param n_sample: int
        number of samples per match-up series

    :return:
        :HData_sample: harm_data_reader
            sampled harmonisation data
        :i_sels: list: int
            indices of randomly selected data
    """

    # If sample size is larger than the data return the data and a i_sels of all the data indices
    if (n_sample >= HData.values.shape[0]) or (n_sample == -1):
        return HData, arange(0, HData.values.shape[0])

    else:
        HData_sample = HarmData()
        HData_sample.a = HData.a[:]
        HData_sample.sensor_model = HData.sensor_model
        HData_sample.adjustment_model = HData.adjustment_model

        ################################################################################################################
        # 1. Sample Data
        ################################################################################################################

        # Select random indices to sample from data
        i_sels = sort((HData.values.shape[0] * rand(n_sample)).astype(int))

        # Perform sampling
        HData_sample.values = HData.values[i_sels, :]
        HData_sample.times = HData.times[i_sels]

        ################################################################################################################
        # 2. Update indices
        ################################################################################################################

        # start with copy of old idx
        HData_sample.idx = deepcopy(HData.idx)

        # count how many match-ups are removed from each series
        n_mu = len(HData.idx['Im'])
        HData_sample.idx['Nm'] = zeros(n_mu, dtype=int)
        for i_sel in i_sels:
            for i in range(n_mu):
                i_start = HData.idx['cNm'][i]
                i_end = HData.idx['cNm'][i+1]
                if i_start <= i_sel < i_end:
                    HData_sample.idx['Nm'][i] += 1

        HData_sample.idx['cNm'] = append(zeros(1), cumsum(HData_sample.idx['Nm'])).astype(int)
        HData_sample.idx['N_var'] = [int(HData_sample.idx['Nm'][n - 1]) for n in HData_sample.idx['n_mu']]
        HData_sample.idx['idx'] = append(zeros(1), cumsum(HData_sample.idx['N_var'])).astype(int)

        return HData_sample, i_sels


def sim_sensor_values(HData, n_sim=20):
    """
    Simulate harmonisation data

    :param HData: harm_data_reader.HarmData
        input harmonisation data
    :param n_sim: int
        number of simulated match-ups per match-up series (default = 20)

    :return:
        :values_sim: dict:numpy.ndarray
            simulated harmonisation data per sensor
    """

    # get list of sensors removing reference if present
    sensors = filter(lambda a: a != -1, HData.idx['sensors'])

    # initialise data dict
    values_sim = {}

    for sensor in sensors:
        # get max/min values for given sensor across all match-up

        m = HData.values.shape[1]/2

        max_sensor = -1*ones(m)
        min_sensor = -1*ones(m)

        # step through sensor match-up by match-up checking for largest/smallest values
        for i, Im in enumerate(HData.idx['Im']):
            for j, sensor_m_id in enumerate(Im):

                sensor_m = HData.idx['sensors'][sensor_m_id]
                # if sensor of given match-up is the sensor we are looking for
                if sensor == sensor_m:

                    # pull out the data for this sensor in this match-up
                    istart = HData.idx['cNm'][i]
                    iend = HData.idx['cNm'][i+1]

                    jstart = j*m
                    jend =j*m + m

                    max_sensor_m = amax(HData.values[istart:iend, jstart:jend], axis=0)
                    min_sensor_m = amin(HData.values[istart:iend, jstart:jend], axis=0)

                    # and replace old values for max and
                    for k, (max_sensor_m_v, min_sensor_m_v) in enumerate(zip(max_sensor_m, min_sensor_m)):

                        # if match-up max values greater than previous match-up values replace
                        if max_sensor_m_v > max_sensor[k]:
                            max_sensor[k] = max_sensor_m_v

                        # if match-up max values greater than previous match-up values replace
                        if min_sensor_m_v < min_sensor[k]:
                            min_sensor[k] = min_sensor_m_v
                        # initialised value -1 so always replace on first run
                        elif min_sensor[k] == -1:
                            min_sensor[k] = min_sensor_m_v

        values_sim_s = zeros((n_sim, m))
        for l, (max_sensor_i, min_sensor_i) in enumerate(zip(max_sensor, min_sensor)):
            values_sim_s[:, l] = linspace(min_sensor_i, max_sensor_i, n_sim)

        values_sim[sensor] = deepcopy(values_sim_s)

    return values_sim


def calc_R(HData, a=None, Ia=None, V=None):
    """
    Calculate radiance from input harmonisation data

    :param HData: harm_data_reader.HarmData
        harmonisation data
    :param a: numpy.ndarray
        harmonisation parameters
    :param V: numpy.ndarray
        covariance matrix of parameter estimates, if a are parameter estimates not input values

    :return:
        :R: numpy.ndarray
            radiance
    """

    if a is None:
        a = HData.a[:]
    if Ia is None:
        Ia = HData.idx['Ia']

    m = HData.values.shape[1] / 2

    R = zeros((HData.values.shape[0], 2))
    if V is not None:
        uR = zeros((HData.values.shape[0], 2))

    for i, Im in enumerate(HData.idx['Im']):

        istart = HData.idx['cNm'][i]
        iend = HData.idx['cNm'][i + 1]

        for j, n_sensor in enumerate(Im):
            colHs = j * m
            colHe = j * m + m

            colR = j

            if isnan(HData.values[iend-1, colHs:colHe]).any():
                R[istart:iend, colR] = nan
                if V is not None:
                    uR[istart:iend, colR] = nan

            else:
                i_s = asarray([True if s == HData.idx['sensors'][n_sensor] else False for s in Ia], dtype=bool)

                R[istart:iend, colR], J = HData.sensor_model(a[i_s], HData.values[istart: iend, colHs:colHe].T)
                if V is not None:
                    uR[istart:iend, colR] = (diag(J[:, m:].dot(V[outer(i_s, i_s)].reshape((sum(i_s), sum(i_s)))).dot(J[:, m:].T))) ** 0.5

    if V is None:
        return R
    else:
        return R, uR


def calc_R_sensor(values, sensor_model, a, V=None):
    """
    Calculate radiance from input sensor data

    :param values: numpy.ndarray
        sensor data with variable by column and observation by row
    :param sensor_model: function
        sensor model to combine sensor data a parameters to evalulate radiance
    :param a: numpy.ndarray
        harmonisation parameters
    :param V: numpy.ndarray
        covariance matrix of parameter estimates, if a are parameter estimates not input values

    :return:
        :R: numpy.ndarray
            radiance
    """

    R, J = sensor_model(a, values.T)
    if V is not None:
        uR = (diag(J[:, values.shape[1]:].dot(V).dot(J[:, values.shape[1]:].T))) ** 0.5

    if V is None:
        return R
    else:
        return R, uR


def calc_R_res(R, R_est, uR=None):
    """
    Return R_res and errorbars

    :param R: numpy.ndarray
        Input radiance values
    :param R_est: numpy.narray
        Estimated radiance values

    :return:
        :R_res: numpy.ndarray
            Radiance residual between estimated and input radiance values
        :errorbars: numpy.ndarray
            Errorbars of radiance residual calculated in the GUM way
    """

    R_res = R_est - R
    if uR is not None:
        errorbars = uR * 2

        return R_res, errorbars

    else:
        return R_res

def average_sensor(values, lm):
    """
    Return averages of input values per sensor

    :param values: numpy.ndarray
        array of data to organised as H

    :param lm: numpy.ndarray
        array to describe match-up pairs. One row per match-up pair, column as [sensor1, sensor2, n_mu]

    :return:
        :mean: list:float
            mean per sensor of input values
        :sd: list:float
            standard deviation per sensor of input values
    """

    sensors_full = lm[:, :2].flatten()
    sensors = []
    for sensor in sensors_full:
        if (sensor != -1) and (sensor not in sensors):
            sensors.append(sensor)

    means = zeros(len(sensors))
    sd = zeros(len(sensors))
    for i, sensor in enumerate(sensors):
        values_sensor = get_values_sensor(values, lm, sensor)
        means[i] = mean(values_sensor)
        sd[i] = std(values_sensor)

    return means, sd


def get_values_sensor(values, lm, sensor):
    """
    Return values per sensor

    :param values: numpy.ndarray
        array of data to select per sensor data from

    :param lm: numpy.ndarray
        array to describe match-up pairs. One row per match-up pair, column as [sensor1, sensor2, n_mu]

    :param sensor: int
        sensor to select data from

    :return:
        :values_sensor: numpy.ndarray
            values from values array for just selected sensor

    """

    m = values.shape[1]/2
    cNm = append(zeros(1), cumsum(lm[:, 2])).astype(int)

    # initialise array
    n_mu = array([], dtype=int)
    for i, Im in enumerate(lm):
        if (sensor == Im[0]) or (sensor == Im[1]):
            n_mu = append(n_mu, array(Im[2])).astype(int)

    N_mu = int(sum(n_mu))
    cn_mu = append(zeros(1), cumsum(n_mu)).astype(int)
    values_sensor = zeros((N_mu, m))

    i_mu_sensor = 0
    # populate array
    for i, Im in enumerate(lm[:, :2]):
        for j, sensor_m in enumerate(Im):

            # if sensor of given match-up is the sensor we are looking for
            if sensor == sensor_m:
                # pull out the data for this sensor in this match-up
                istart = cNm[i]
                iend = cNm[i+1]
                jstart = j * m
                jend = j * m + m

                # and put it into the right place in values sensor
                istart_sensor = cn_mu[i_mu_sensor]
                iend_sensor = cn_mu[i_mu_sensor + 1]

                i_mu_sensor += 1

                values_sensor[istart_sensor:iend_sensor] = values[istart:iend, jstart:jend]

    return values_sensor


def average_mu(values, lm):
    """
    Return averages of input values per match-up

    :param values: harm_data_reader.HarmData
        array of data to

    :param lm: numpy.ndarray
        array to describe match-up pairs. One row per match-up pair, column as [sensor1, sensor2, n_mu]

    :return:
        :mean: list:float
            mean per sensor of input values
        :sd: list:float
            standard deviation per sensor of input values
    """

    mus = list(range(lm.shape[0]))

    means = zeros(len(mus))
    sd = zeros(len(mus))
    for i, mu in enumerate(mus):
        values_mu = get_values_mu(values, lm, mu)

        #remove nans
        values_mu = values_mu[~isnan(values_mu)]

        means[i] = mean(values_mu)
        sd[i] = std(values_mu)

    return means, sd


def get_values_mu(values, lm, mu):
    """
    Return values per sensor

    :param values: numpy.ndarray
        array of data to select per sensor data from

    :param idx: dict
        harm_data_reader.HarmData.idx attribute

    :param mu: int
        match-up number to select data for

    :return:
        :values_sensor: numpy.ndarray
            values from values array for just selected sensor

    """

    cNm = append(zeros(1), cumsum(lm[:, 2])).astype(int)
    istart = cNm[mu]
    iend = cNm[mu+1]
    values_mu = values[istart:iend]

    return values_mu

if __name__ == "__main__":

    def main():
        return 0

    main()
