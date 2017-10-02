"""
General plotting functions

Created on Wed May 03 2017 09:00:00

@author: Peter Harris, NPL\MM
@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import arange, where, zeros, linspace, nan, isnan, ndenumerate
from datetime import datetime

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)


def plot_grid_heatmap(savePath, data, title, labels=None, dividers=None, range=None):
    """
    Generate heatmap of input array and save to file

    :param savePath: str
        Path (including extension) of file to save plot to

    :param data: numpy.ndarray
        2d numpy array to plot

    :param labels: list:str
        List of axis labels for each element of the grib

    :param range: list: float
        Range of colour scale for heatmap

    """

    fig1, (ax1) = plt.subplots(1)
    ax1.imshow(data, cmap="bwr", interpolation='nearest', vmin=-1, vmax=1)
    ax1.set_title(title)

    # Add labels to each pixel
    fontsize_pixel = 14. - 1./3.*len(labels)
    for (j, i), label in ndenumerate(data):
        label = round(label, 3)
        x_pos = (i+1.0)*(1.0/data.shape[0])-(1.0/(2*data.shape[0]))
        y_pos = 1-(j+1.0)*(1.0/data.shape[1])+(1.0/(2*data.shape[1]))
        ax1.text(x_pos, y_pos, label, ha='center', va='center', transform=ax1.transAxes, fontsize=fontsize_pixel)

    # Add axis labels
    fontsize_label = 14. - 1./3.*len(labels)
    for i, label in enumerate(labels):
        x_pos = (i + 1.0) * (1.0 / data.shape[0]) - (1.0 / (2 * data.shape[0]))
        y_pos = 1 - (i + 1.0) * (1.0 / data.shape[1]) + (1.0 / (2 * data.shape[1]))
        ax1.text(x_pos, -0.05, label, ha='center', va='center', transform=ax1.transAxes, fontsize=fontsize_label)
        ax1.text(-0.05, y_pos, label, ha='center', va='center', transform=ax1.transAxes, fontsize=fontsize_label)

    # Add dividers
    if dividers is not None:
        n_d = len(dividers)
        for i, divider in enumerate(dividers):

            x_pos = (i + 1.0) * (1.0 / n_d) - (1.0 / (2 * n_d))
            y_pos = 1 - (i + 1.0) * (1.0 / n_d) + (1.0 / (2 * n_d))
            ax1.text(x_pos, -0.1, divider, ha='center', va='center', transform=ax1.transAxes)
            ax1.text(-0.1, y_pos, divider, ha='right', va='center', transform=ax1.transAxes)


        xlims = ax1.get_xlim()
        xrange = xlims[1] - xlims[0]

        ylims = ax1.get_ylim()
        yrange = ylims[0] - ylims[1]
        for i in arange(1, n_d):
            ax1.axvline(xrange/n_d*i+xlims[0], color='k')
            ax1.axhline(yrange/n_d*i+ylims[1], color='k')

    ax1.tick_params(bottom='off', top='off', labelbottom='off', right='off', left='off',
                    labelleft='off')

    ax1.set_xticks([])
    ax1.set_yticks([])

    #plt.tight_layout()
    plt.savefig(savePath)
    plt.close(fig1)

    return 0


def plot_scatter(savePath, Y, X, xlbl, ylbl, title, yerr=None, xerr=None, txt=None, *args, **kwargs):
    """
    Generate x-y plot of input data (optionally including error bars) and save to file

    :param savePath: str
        Path (including extension) of file to save plot to

    :param X: numpy.ndarray
        X parameter

    :param Y: numpy.ndarray
        Y parameter

    :param xlbl: str
        x axis label

    :param ylbl: str
        y axis label

    :param ylbl: str
        y axis label

    :param yerr: numpy.ndarray
        yerr parameter

    :param xerr: numpy.ndarray
        xerr parameter

    :param txt: list:str
        list of strings to stack in textbox

    """

    switch2labels = False
    if all(isinstance(x, str) for x in X):
        switch2labels = True
        X_labels = X
        X = arange(len(X), dtype=float)

    # initialise plot
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)

    # Ignore nans
    if (str(type(X)) == "<type 'numpy.ndarray'>"):
        if(X.dtype == float):
            Y = Y[~isnan(X)]
            X = X[~isnan(X)]
    if (str(type(Y)) == "<type 'numpy.ndarray'>"):
        if (Y.dtype == float):
            X = X[~isnan(Y)]
            Y = Y[~isnan(Y)]

    # set scales

    # ys
    if yerr is not None:
        Y_max = max(Y+yerr)
        Y_min = min(Y-yerr)
    else:
        Y_max = max(Y)
        Y_min = min(Y)

    # xs
    Y_width = Y_max - Y_min
    Y_inc = Y_width / Y.shape[0]
    Y_min -= Y_inc
    Y_max += Y_inc

    if xerr is not None:
        X_max = max(X+xerr)
        X_min = min(X+xerr)
    else:
        X_max = max(X)
        X_min = min(X)

    X_width = X_max - X_min
    X_inc = X_width / X.shape[0]
    X_min -= X_inc
    X_max += X_inc

    ax.set_xlim([X_min, X_max])
    ax.set_ylim([Y_min, Y_max])

    # add axis labels and title
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    ax.set_title(title)

    # add any required horizontal or vertical lines
    if "dash_ylines" in kwargs:
        for y in kwargs["dash_ylines"]:
            ax.axhline(y=y, linewidth=0.7, color='k', linestyle=':')

    if "dash_xlines" in kwargs:
        for x in kwargs["dash_xlines"]:
            ax.axvline(x=x, linewidth=0.7, color='k', linestyle=':')

    if "solid_ylines" in kwargs:
        for y in kwargs["solid_ylines"]:
            ax.axhline(y=y, linewidth=0.7, color='k')

    if "solid_xlines" in kwargs:
        for x in kwargs["solid_xlines"]:
            ax.axvline(x=x, linewidth=0.7, color='k')

    # add error bars to data points
    eb = ax.errorbar(X, Y, yerr=yerr, xerr=xerr,
                     marker='o',
                     color='#333333',
                     ecolor='#333333',
                     elinewidth=2,
                     markerfacecolor='red',
                     markeredgecolor='#333333',
                     markeredgewidth=1,
                     markersize=4,
                     capsize=0,
                     linestyle='None')

    if switch2labels is True:
        ax.set_xticks(X)
        ax.set_xticklabels(X_labels, rotation=90)

    if txt is not None:
        s = ""
        for string in txt:
            s += string + "\n"
        s = s[:-2]
        ax.text(0.05, 0.9, s, fontsize=8, bbox=dict(edgecolor='k', facecolor='w'), transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(savePath)
    plt.close(fig)
    return 0


def plot_2dhist(savePath, Y, X, bins, xlbl, ylbl, title, txt=None, *args, **kwargs):
    """
    Generate 2d histogram (heatmap) of input x-y data and save to file

    :param savePath: str
        Path (including extension) of file to save plot to

    :param X: numpy.ndarray
        X parameter

    :param Y: numpy.ndarray
        Y parameter

    :param bins: int
        total number of bins

    :param xlbl: str
        x axis label

    :param ylbl: str
        y axis label

    :param ylbl: str
        y axis label

    :param txt: list:str
        list of strings to stack in textbox

    """

    # initialise figure
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)

    # add axis labels and title
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    ax.set_title(title)

    # Change plotting variables to seconds since epoch if datetime as hist2d can't handle it
    X_dates = None
    if type(X[0]) is datetime:
        X_dates = X[:]
        epoch = datetime.utcfromtimestamp(0)
        X = zeros(X_dates.shape[0])
        for i, date in enumerate(X_dates):
            X[i] = (date - epoch).total_seconds()

    # Ignore nans
    X = X[~isnan(X)]
    Y = Y[~isnan(X)]
    X = X[~isnan(Y)]
    Y = Y[~isnan(Y)]

    h, x, y, p = ax.hist2d(X, Y, bins=bins, cmap="Reds")

    if txt is not None:
        s = ""
        for string in txt:
            s += string + "\n"
        s = s[:-2]
        ax.text(0.05, 0.9, s, fontsize=8, bbox=dict(edgecolor='k', facecolor='w'), transform=ax.transAxes)

    # add any required horizontal or vertical lines
    if "dash_ylines" in kwargs:
        for y in kwargs["dash_ylines"]:
            ax.axhline(y=y, linewidth=1, color='k', linestyle=':')

    if "dash_xlines" in kwargs:
        for x in kwargs["dash_xlines"]:
            ax.axvline(x=x, linewidth=1, color='k', linestyle=':')

    if "solid_ylines" in kwargs:
        for y in kwargs["solid_ylines"]:
            ax.axhline(y=y, linewidth=1, color='k')

    if "solid_xlines" in kwargs:
        for x in kwargs["solid_xlines"]:
            ax.axvline(x=x, linewidth=1, color='k')

    # Set the axis labels as dates if required
    if X_dates is not None:
        years = list(set([time.year for time in X_dates]))  # get unique years
        year_max = max(years)
        year_min = min(years)
        years = arange(year_min, year_max+1).astype(int)
        if years.shape[0] > 15:
            years = linspace(year_min, year_max, num=8).astype(int)
        epoch = datetime.utcfromtimestamp(0)
        seconds = [(datetime(day=1, month=1, year=year)-epoch).total_seconds() for year in years]
        ax.set_xticks(seconds)
        ax.set_xticklabels(years, rotation=90)

    plt.colorbar(p, ax=ax)  # , ticks=[-1, 0, 1])
    plt.tight_layout()
    plt.savefig(savePath)
    plt.close(fig)

    return 0

if __name__ == "__main__":

    def main():
        return 0

    main()
