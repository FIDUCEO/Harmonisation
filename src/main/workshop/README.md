# Lisbon Workshop, 17-19 April 2018

## How to fetch the interactive lecture from GitHub

Open a web browser and visit the [Lisbon Workshop on GitHub](https://github.com/FIDUCEO/Harmonisation/tree/master/src/main/workshop).
Right click on the [harmonisation_interactive_lecture.ipynb](https://github.com/FIDUCEO/Harmonisation/blob/master/src/main/workshop/harmonisation_interactive_lecture.ipynb)
file listed on top of the page to download it to your computer. You may download the interactive lecture to any location. 

## How to get started on macOS

### Installing prerequisites

#### Anaconda way

To use the interactive lecture you need to have Python 3 and [Jupyter](https://jupyter.org) installed. If you do not have Python 3
or Jupyter installed, install it. The most simple way is to install [Anaconda](https://www.anaconda.com). Download the
Anaconda installer and follow the instructions on the Anaconda web site.

#### Homebrew way

For more experienced users the more convenient way to install Python 3 and Jupyter may be the [Homebrew](https://brew.sh)
package manager, which allows you to install everything you need, and get rid of it when you do not need it anymore,
including itself.

To install Homebrew open a Terminal window and type

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

To install Python 3 type

    brew install python

To bring Python's Pip package manager up to date type

    pip3 install --upgrade pip setuptools wheel

To install Jupyter (and Python modules required by the interactive lecture) type

    pip3 install jupyter matplotlib numpy scipy

You have completed the installation of prerequisites.

### Opening the interactive lecture

Now you are ready to start the interactive lecture. Open a Terminal window, `cd` into the directory where the lecture resides
and type

    jupyter notebook

A browser window opens and shows the directory contents. Click on the interactive lecture and play with it!

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
    
#### Homebrew

To uninstall Jupyter (and Python modules required by the interactive lecture) open a Terminal window and type

    pip3 uninstall jupyter matplotlib numpy scipy

To uninstall Python 3 type

    brew uninstall python

To uninstall homebrew itself type

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"

## How to get started on Windows

### Installing prerequisites

To use the interactive lecture you need to have Python 3 and [Jupyter](https://jupyter.org) installed. If you do not have Python 3
or Jupyter installed, install it. The most simple way is to install [Anaconda](https://www.anaconda.com). Download the
Anaconda installer and follow the instructions on the Anaconda web site.

### Opening the interactive lecture

On the folder that contains the interactive lecture hold down the shift key and right-click. On the menu that appears
click “Open command window here”. Then get started with the following command:

    jupyter notebook

A browser window opens and shows the directory contents. Click on the interactive lecture and play with it!

### Uninstalling prerequisites

Navigate to "Uninstall a Program" in your Control Panel. Uninstall your Anaconda Python distribution from
there in the usual way. 
