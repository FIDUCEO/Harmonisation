# Lisbon Workshop, 17-19 April 2018


## How to get started (mac OS)

### Installing prerequisites

To use the notebook you need to have Python 3 and [Jupyter](https://jupyter.org) installed. If you do not have Python 3
installed, install it.

The most convenient way to install Python 3 is the [Homebrew](https://brew.sh) package manager, which allows you to install everything you need, and get rid of it when you do not need it anymore, including itself. If you do not have installed Homebrew, do so by opening a Terminal window and typing

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Then install Python 3 by typing

    brew install python

and bring the Python's Pip package manager up to date

    pip3 install --upgrade pip setuptools wheel

Then install Jupyter by typing

    pip3 install jupyter

You have completed the installation of prerequisites.

### Opening the notebook

Now you are ready to use the notebook. `cd` into the notebook directory and type

    jupyter notebook

A browser window opens and shows the directory contents. Click on the notebook and play with it!

### Uninstalling prerequisites
 
To uninstall Jupyter open a Terminal window and type

    pip3 uninstall jupyter

To uninstall Python 3 open a Terminal window and type

    brew uninstall python

To uninstall homebrew open a Terminal window and type

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"
