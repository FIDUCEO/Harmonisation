'''___Built-In Modules___'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import arange, where, zeros, linspace, nan, isnan, ndenumerate, ndarray
from scipy.stats import kde
from datetime import datetime
from numpy import array
import matplotlib
#from mpl_toolkits.basemap import Basemap
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'


def plot_scatter(savePath, Y, X, xlbl=None, ylbl=None, title=None, yerr=None, xerr=None, txt=None, line='None',
                 legend=None, legend_loc="top-right", fig_size=(10,8), ylim=None, xlim=None, *args, **kwargs):
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
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot(111)
    # Ignore nans
    if (type(X) == ndarray) or (type(Y) == ndarray):
        if (type(X) == ndarray) and (X.dtype == float):
            Y = Y[~isnan(X)]
            yerr = yerr[~isnan(X)] if yerr is not None else None
            xerr = xerr[~isnan(X)] if xerr is not None else None
            X = X[~isnan(X)]
        if (type(Y) == ndarray) and Y.dtype == float:
            X = X[~isnan(Y)]
            yerr = yerr[~isnan(Y)] if yerr is not None else None
            xerr = xerr[~isnan(Y)] if xerr is not None else None
            Y = Y[~isnan(Y)]

    else:
        if X[0].dtype == float:
            Y = [Y_i[~isnan(X_i)] for X_i, Y_i in zip(X, Y)]
            yerr = [yerr_i[~isnan(X_i)] for X_i, yerr_i in zip(X, yerr)] if yerr is not None else None
            xerr = [xerr_i[~isnan(X_i)] for X_i, xerr_i in zip(X, xerr)] if xerr is not None else None
            X = [X_i[~isnan(X_i)] for X_i in X]
        if Y[0].dtype == float:
            X = [X_i[~isnan(Y_i)] for X_i, Y_i in zip(X, Y)]
            yerr = [yerr_i[~isnan(Y_i)] for Y_i, yerr_i in zip(Y, yerr)] if yerr is not None else None
            xerr = [xerr_i[~isnan(Y_i)] for Y_i, xerr_i in zip(Y, xerr)] if xerr is not None else None
            Y = [Y_i[~isnan(Y_i)] for Y_i in Y]

        for i in range(len(Y)):
            if yerr is not None:
                print Y[i].shape, X[i].shape, yerr[i].shape

    if legend == None:
        if (type(X) == list) and ((type(X[0]) == list) & (type(X[0]) == ndarray)):
            legend = ['__nolegend__' for i in range(len(X))]
        else:
            legend = ['__nolegend__']

    # set scales

    # # ys
    # if yerr is not None:
    #     Y_max = max(Y+yerr)
    #     Y_min = min(Y-yerr)
    # else:
    #     Y_max = max(Y)
    #     Y_min = min(Y)
    #
    # # xs
    # Y_width = Y_max - Y_min
    # Y_inc = Y_width / Y.shape[0]
    # Y_min -= Y_inc
    # Y_max += Y_inc
    #
    # if xerr is not None:
    #     X_max = max(X+xerr)
    #     X_min = min(X+xerr)
    # else:
    #     X_max = max(X)
    #     X_min = min(X)
    #
    # X_width = X_max - X_min
    # X_inc = X_width / X.shape[0]
    # X_min -= X_inc
    # X_max += X_inc

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

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
    if (type(X) != list) and ((type(X[0]) != list) & (type(X[0]) != ndarray)):
        X = [X]
        Y = [Y]

    colours = ["#1f78b4",
               "#a6cee3",
               "#b2df8a",
               "#33a02c",
               "#fb9a99",
               "#e31a1c",
               "#fdbf6f",
               "#ff7f00",
               "#cab2d6"]

    if len(X) > len(colours):
        cmap = plt.get_cmap('rainbow')
        colours = [cmap(i) for i in np.linspace(0, 1, len(X))]

    for i, (X_i, Y_i) in enumerate(zip(X, Y)):
        ax.errorbar(X_i, Y_i,
                    yerr=yerr[i] if yerr is not None else None,
                    xerr=xerr[i] if xerr is not None else None,
                    marker='o',
                    ecolor='#333333',
                    elinewidth=2,
                    markerfacecolor=colours[i % len(colours)],
                    markeredgecolor='#333333',
                    markeredgewidth=1,
                    markersize=4,
                    capsize=0,
                    linestyle=line,
                    color=colours[i % len(colours)],
                    label=legend[i])

    if set(legend) != {'__nolegend__'}:
        if legend_loc == "top_right":
            ax.legend()
        elif legend_loc == "outside":
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints=1)
        else:
            ax.legend()

    # if switch2labels is True:
    #     ax.set_xticks(X)
    #     ax.set_xticklabels(X_labels, rotation=90)

    if txt is not None:
        s = ""
        for string in txt:
            s += string + "\n"
        s = s[:-2]
        ax.text(0.05, 0.9, s, fontsize=8, bbox=dict(edgecolor='k', facecolor='w'), transform=ax.transAxes)

    #plt.tight_layout()
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
    Y = Y[~isnan(X)]
    X = X[~isnan(X)]
    X = X[~isnan(Y)]
    Y = Y[~isnan(Y)]

    h, x, y, p = ax.hist2d(X, Y, bins=bins, cmap="Reds")

    if txt is not None:
        s = ""
        for string in txt:
            s += string + "\n"
        s = s[:-2]
        ax.text(0.05, 0.9, s, fontsize=10, bbox=dict(edgecolor='k', facecolor='w'), transform=ax.transAxes)

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


def plot_density(path, y, x, bins, xlbl, ylbl, title, txt=None, *args, **kwargs):
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)

    # Remove nans from input data
    if x.dtype == float:
        y = y[~isnan(x)]
        x = x[~isnan(x)]
    if y.dtype == float:
        x = x[~isnan(y)]
        y = y[~isnan(y)]

    # add axis labels and title
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    ax.set_title(title)

    if txt is not None:
        s = ""
        for string in txt:
            s += string + "\n"
        s = s[:-2]
        ax.text(0.05, 0.9, s, fontsize=10, bbox=dict(edgecolor='k', facecolor='w'), transform=ax.transAxes)

    x_dates = None
    if x.dtype == datetime:
        x_dates = x
        x = array([(t - datetime(1970, 1, 1)).total_seconds() for t in x])

    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():bins * 1j, y.min():y.max():bins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    ax.set_xlim([x.min(), x.max()])
    ax.set_ylim([y.min(), y.max()])

    # Make the plot
    if x_dates is not None:
        years = list(set([time.year for time in x_dates]))  # get unique years
        year_max = max(years)
        year_min = min(years)
        years = arange(year_min, year_max+1).astype(int)
        if years.shape[0] > 15:
            years = linspace(year_min, year_max, num=8).astype(int)
        epoch = datetime.utcfromtimestamp(0)
        seconds = [(datetime(day=1, month=1, year=year)-epoch).total_seconds() for year in years]
        ax.set_xticks(seconds)
        ax.set_xticklabels(years, rotation=90)

    # add any required horizontal or vertical lines
    if "dash_ylines" in kwargs:
        for yline in kwargs["dash_ylines"]:
            ax.axhline(y=yline, linewidth=1, color='k', linestyle=':')

    if "dash_xlines" in kwargs:
        for xline in kwargs["dash_xlines"]:
            ax.axvline(x=xline, linewidth=1, color='k', linestyle=':')

    if "solid_ylines" in kwargs:
        for yline in kwargs["solid_ylines"]:
            ax.axhline(y=yline, linewidth=1, color='k')

    if "solid_xlines" in kwargs:
        for xline in kwargs["solid_xlines"]:
            ax.axvline(x=xline, linewidth=1, color='k')

    # Change color palette
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap="BuGn")
    plt.colorbar()
    plt.savefig(path)
    plt.show()


# def plot_map_density(path, lat, lon, bins, title):
#     fig = plt.figure()
#     ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
#
#     m = Basemap(projection='kav7', lon_0=0, lat_0=0)
#     m.drawcoastlines()
#
#     # draw parallels.
#     m.drawparallels(np.arange(-90., 99., 30.))
#     m.drawmeridians(np.arange(-180., 180., 60.))
#
#     k = kde.gaussian_kde([lon, lat])
#     loni, lati = np.mgrid[lon.min():lon.max():bins * 1j, lat.min():lat.max():bins * 1j]
#     zi = k(np.vstack([loni.flatten(), lati.flatten()]))
#     im = m.pcolormesh(loni, lati, zi.reshape(loni.shape), shading="Flat", cmap="BuGn", latlon=True)
#     cb = m.colorbar(im, "bottom", size="5%", pad="2%")
#     ax.set_title(title)
#
#     plt.savefig(path)
#     return 0