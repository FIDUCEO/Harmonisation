# Lisbon Workshop, 17-19 April 2018


## How to get started (mac OS)

### Installing prerequisites

To use the notebook you need to have Python 3 and [Jupyter](https://jupyter.org) installed. If you do not have Python 3
or Jupyter installed, install it. The most simple way is to install [Anaconda](https://www.anaconda.com).

For the more experienced users the more convenient way to install Python 3 and Jupyter may be the [Homebrew](https://brew.sh)
package manager, which allows you to install everything you need, and get rid of it when you do not need it anymore,
including itself.

To install Homebrew open a Terminal window and type

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

To install Python 3 type

    brew install python

To bring Python's Pip package manager up to date type

    pip3 install --upgrade pip setuptools wheel

To install Jupyter type

    pip3 install jupyter

You have completed the installation of prerequisites.

### Opening the notebook

Now you are ready to use the notebook. Open a Terminal window, `cd` into the notebook directory and type

    jupyter notebook

A browser window opens and shows the directory contents. Click on the notebook and play with it!

### Uninstalling prerequisites
 
To uninstall Jupyter open a Terminal window and type

    pip3 uninstall jupyter

To uninstall Python 3 type

    brew uninstall python

To uninstall homebrew type

    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"
