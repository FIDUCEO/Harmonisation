# Lisbon Workshop, 25-27 June 2019

## How to fetch the interactive lecture files

### Getting the interactive lecture files from Dropbox

The interactive lecture files are available on Dropbox for download. Do this in the following way:

1.	Follow this [link](https://www.dropbox.com/sh/zowwp25rsj713on/AABYbsVKhu2yJpTkMjAaw5Qba?dl=0) to the interactive lecture folder on Dropbox.
2.	Download the contents of this folder to your computer, by clicking download then direct download.
3.	Unzip the files to a folder wherever you like on your computer.


### Getting the interactive lecture files from GitHub

Users familiar with git may prefer to access the files by cloning the FIDUCEO [Harmonisation repository](https://github.com/FIDUCEO/Harmonisation). The workshop files are located in `Harmonisation/training/workshop_2019`.

## How to get started on Linux

### Installing prerequisites

To use the interactive lecture you need to have Python 3 and [Jupyter](https://jupyter.org) installed. If you do not have Python 3
or Jupyter installed, install it. The most simple way may be to install [Anaconda](https://www.anaconda.com).

Either download the [Anaconda installer](https://www.anaconda.com/download/#macos) and follow the instructions on the Anaconda web
site or use your system's package manager to install Jupyter and Python 3.

In addition to python modules that come with Anaconda, cartopy and xarray must installed. To do this run the following in the command line:

    conda install -c conda-forge cartopy
    conda install xarray

### Opening the interactive lecture

Now you are ready to start the interactive lecture. Open a Terminal window, `cd` into the directory where the lecture files reside and type

    jupyter notebook

A browser window opens and shows the directory contents.  Click on an interactive lecture iPython notebook file and play with it!

### Uninstalling prerequisites

If you installed Anaconda from the Anaconda web site, the uninstallation procedure described for macOS below may work for Linux
as well.

## How to get started on macOS

### Installing prerequisites

To use the interactive lecture you need to have Python 3 and [Jupyter](https://jupyter.org) installed. If you do not have Python 3 or Jupyter installed, install it. The most simple way is to install [Anaconda](https://www.anaconda.com).

#### Anaconda

Download the [Anaconda installer](https://www.anaconda.com/download/#macos) and follow the instructions on the Anaconda web site.

In addition to python modules that come with Anaconda, cartopy and xarray must installed. To do this run the following in the command line:

    conda install -c conda-forge cartopy
    conda install xarray

### Opening the interactive lecture

Now you are ready to start the interactive lecture. Open a Terminal window, `cd` into the directory where the lecture files reside and type

    jupyter notebook

A browser window opens and shows the directory contents.  Click on an interactive lecture iPython notebook file and play with it!

### Uninstalling prerequisites

#### Anaconda

To uninstall Anaconda open a Terminal window and type

    conda install anaconda-clean

to install the Anaconda uninstaller. Then type

    anaconda-clean --yes

to execute the uninstaller. Then remove all Anaconda directories from your home directory

    rm -rf ~/anaconda[23]/
    rm -rf ~/.anaconda_backup/

and edit your `~/.profile` to remove the Anaconda directory from your `PATH` environment variable.

## How to get started on Windows

### Installing prerequisites

To use the interactive lecture you need to have Python 3 and [Jupyter](https://jupyter.org) installed. If you do not have Python 3 or Jupyter installed, install it. The most simple way is to install [Anaconda](https://www.anaconda.com). Download the
Anaconda installer and follow the instructions on the Anaconda web site.

In addition to python modules that come with Anaconda, cartopy and xarray must installed. To do this run the following in the command line:

    conda install -c conda-forge cartopy
    conda install xarray

### Opening the interactive lecture

On the folder that contains the interactive lecture files hold down the shift key and right-click. On the menu that appears
click “Open command window here”. Then get started with the following command:

    jupyter notebook

A browser window opens and shows the directory contents.  Click on an interactive lecture iPython notebook file and play with it!

### Uninstalling prerequisites

Navigate to "Uninstall a Program" in your Control Panel. Uninstall your Anaconda Python distribution from
there in the usual way. 
