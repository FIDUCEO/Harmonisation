"""
Created 2017/02/15

@author: Sam Hunt, NPL\ENV
"""

'''___Python Modules___'''
import matplotlib.pyplot as plt
from numpy import ndenumerate
from numpy.random import random


def plot_corr_heatmap(A, title=None, labels=None, save_path=None):
    """
    Plot heatmap of input array

    :param A: numpy.ndarray
        input array to plot heatmap
    :param save_path: str
        path to save chart image to
    """

    fig1, (ax1) = plt.subplots(1)
    ax1.imshow(A, cmap="bwr", interpolation='nearest', vmin=-1, vmax=1)
    ax1.set_title(title)
    for (j, i), label in ndenumerate(A):
        label = round(label, 2)
        ax1.text(i, j, label, ha='center', va='center')

    ax1.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='on', right='off', left='off',
                    labelleft='on')

    tick_labels = []
    for l in labels:
        tick_labels.append("")
        tick_labels.append(l)
    ax1.set_xticklabels(tick_labels)
    ax1.set_yticklabels(tick_labels)

    if save_path is None:
        plt.show()
    else:
        plt.savefig(save_path)

    return 0

if __name__ == "__main__":
    def sample_code():
        """
        Usage example for plot_corr_heatmap
        """

        A = random((4,4))

        plot_corr_heatmap(A, title="Parameter Correlation Heat Map", labels=['a0', 'a1', 'a2', 'a3'], save_path="test.png")

        return 0

    sample_code()
