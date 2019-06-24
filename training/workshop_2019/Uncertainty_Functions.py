import warnings
warnings.simplefilter('ignore') # filter some warning messages
import numpy as np
import ipywidgets as widgets
import matplotlib.pyplot as plt
import xarray as xr



# ~~~~~~~  Key Functions  ~~~~~~~~~~#
def Scan(data, size, centre):
    Window = np.zeros(np.shape(data))
    mod = np.zeros(np.shape(data))
    im_size = np.shape(data)[0]
    x1 = centre[1] - int(size[1] / 2)
    y1 = centre[0] - int(size[0] / 2)
    top_left = np.array([x1, y1])
    for i in range(size[1]):
        for j in range(size[0]):
            Window[x1 + i - 1, y1 + j - 1] = data[x1 + i - 1, y1 + j - 1]
    mod = np.where(data == Window, 0, data)
    return np.array([Window, mod])


def CombinedDa(data, param, Size, Averaging=False, ShowStructured=False):
    N = Size ** 2
    sss = 0
    Window, Background = data
    if ShowStructured == True:
        sss = 1
    Bindependent, Bcommon, Bstructured = Background
    Windependent, Wcommon, Wstructured = Window
    combined = np.zeros(np.shape(Bindependent))
    if Averaging == False:
        return np.sqrt((1 - param) * (Bindependent ** 2) + (param) * Bcommon ** 2 + sss * Bstructured ** 2)
    if Averaging == True:
        idv_avg = np.sqrt(np.sum(Windependent[np.nonzero(Windependent)] ** 2) / (N ** 2))
        # print("Individual = ",idv_avg)
        cm_avg = np.mean(Wcommon[np.nonzero(Wcommon)])
        # print("Common = ",cm_avg)
        str_avg = np.sqrt((1 / N * 2) * (np.sum(np.nonzero(Windependent ** 2))))
        u_comb = np.sqrt((1 - param) * (idv_avg ** 2) + (param) * (cm_avg ** 2) + sss * str_avg ** 2)
        win_comb = np.where(Bindependent == combined, u_comb, 0)
        # print("Combined = ",u_comb)
        return win_comb


# ~~~~~~~~~ Display Functions ~~~~~~~~~#
def IndependentUFig(N_values, x_val, uncertainties, cbar_range, cmap, Names, Titlestart):
    plt.figure(figsize=(8, 5))
    plt.subplot(121)
    plt.imshow(uncertainties[0], vmin=cbar_range[0], vmax=cbar_range[1], cmap=cmap, aspect='auto')
    plt.title(Titlestart + Names[0] + " Error", size=16)
    plt.colorbar()
    plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
                    labelright='off', labelbottom='off')
    plt.subplot(122)
    plt.plot(1 / np.sqrt(N_values), 'x', color="red")
    plt.plot(x_val, 1 / np.sqrt(x_val), "-", color="black")
    plt.ylabel("Independent Uncertaity (%)", size=14)
    plt.xlabel("Value of N")
    plt.xlim(0, 15)
    plt.xticks(np.arange(0, 16, 1))
    plt.axhline(y=0.267, color='grey', linestyle='dashed')
    plt.text(16, 0.265, "~0.267")
    plt.axvline(x=14, color='grey', linestyle='dashed')
    plt.title("$ U_c \propto 1/ \sqrt{N}$")
    plt.tight_layout()
    plt.show()


def CommonUFig(Nvalues, x_val, uncertainties, cbar_range, cmap, Names, Titlestart):
    plt.figure(figsize=(8, 6))
    plt.subplot(121)
    plt.imshow(uncertainties[1], vmin=cbar_range[0], vmax=cbar_range[1], cmap=cmap, aspect='auto')
    plt.title(Titlestart + Names[1] + " Error", size=16)
    plt.colorbar()
    plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
                    labelright='off', labelbottom='off')
    plt.subplot(122)
    plt.plot((x_val / x_val) * 0.076, 'x', color="red")
    plt.plot(x_val, (x_val / x_val) * 0.076, "-", color="black")
    plt.ylabel("Common Uncertaity (%)", size=14)
    plt.xlabel("Value of N")
    plt.xlim(0, 15)
    plt.xticks(np.arange(0, 16, 1))
    # plt.axhline(y=0.076, color='grey',linestyle='dashed')
    plt.text(16, 0.076, "~0.076")
    plt.axvline(x=14, color='grey', linestyle='dashed')
    # plt.title("$ U_c \propto 1/ \sqrt{N}$")
    plt.tight_layout()
    plt.show()


def ThreeErrors(uncertainties, cbar_range, cmap, Titlestart, Names, v1=True):
    if v1==True:
        plt.figure(figsize=(6, 3))
        i = 0
        for elem in uncertainties:
            im = np.array([1, len(uncertainties), i + 1])
            plt.subplot(*im)
            plt.imshow(uncertainties[i], vmin=cbar_range[0], vmax=cbar_range[1], cmap=cmap, aspect='auto')
            plt.title(Titlestart + Names[i] + " Error", size=12)
            i += 1
            plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
                            labelright='off', labelbottom='off')
        plt.tight_layout()
        plt.colorbar()
        plt.show()
    if v1==False:
        plt.figure(figsize=(6, 3))
        for i in range(len(uncertainties)-1):
            im = np.array([1, len(uncertainties), i + 1])
            plt.subplot(*im)
            plt.imshow(uncertainties[i], vmin=cbar_range[0], vmax=cbar_range[1], cmap=cmap, aspect='auto')
            plt.title(Titlestart + Names[i] + " Error", size=12)
            plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
                            labelright='off', labelbottom='off')
        plt.tight_layout()
        plt.colorbar()
        plt.show()