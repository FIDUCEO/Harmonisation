import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta


def time_average(data, times, time_period="month"):

    if time_period == "hour":
        date_figs = 4
        time_delta = relativedelta(hours=1)
    if time_period == "day":
        date_figs = 3
        time_delta = relativedelta(days=1)
    if time_period == "month":
        date_figs = 2
        time_delta = relativedelta(months=1)
    if time_period == "year":
        date_figs = 1
        time_delta = relativedelta(years=1)

    averages = []
    average_times = []

    time_min = [i for i in times.min().timetuple()]
    time_min_period = [1970, 1, 1, 0, 0]
    for i in range(date_figs):
        time_min_period[i] = time_min[i]
    current_time = dt.datetime(*time_min_period)

    time_max = [i for i in times.max().timetuple()]
    time_max_period = [1970, 1, 1, 0, 0]
    for i in range(date_figs):
        time_max_period[i] = time_max[i]
    last_time = dt.datetime(*time_max_period)

    while current_time <= last_time:
        next_time = current_time + time_delta
        ind = np.where((times > current_time) & (times < next_time))
        vals = data[ind]

        try:
            average = np.average(vals)
        except:
            average = np.nan

        averages.append(average)
        average_times.append(current_time)
        current_time = next_time

    return np.asarray(averages), np.asarray(average_times).astype(dt.datetime)


def bin(Y, X, bins=20, X_range=None):

    if X_range is None:
        X_range = [min(X), max(X)]

    averageY = np.zeros(bins)
    bin_ranges = [X_range[0] + i*(X_range[1]-X_range[0])/bins for i in range(bins+1)]
    bins = np.asarray([(bin_max+bin_min)/2 for bin_min, bin_max in zip(bin_ranges[:-1], bin_ranges[1:])])

    for i, (bin_min, bin_max) in enumerate(zip(bin_ranges[:-1], bin_ranges[1:])):

        ind = np.where((X >= bin_min) & (X < bin_max))
        vals = Y[ind]

        try:
            averageY_i = np.average(vals)
        except:
            averageY_i = np.nan

        averageY[i] = averageY_i

    return averageY, bins

if __name__ == "__main__":
    pass

