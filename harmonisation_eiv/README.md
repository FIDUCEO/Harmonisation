# FIDUCEO Harmonisation Software
Software for solving large, non-linear errors-in-variables regression problems

## Dependencies

Before installing with pip the `Cython` and `numpy` modules must be installed to ensure the cython code can compiled during the install. Install these with pip as:

`$ pip install Cython`

`$ pip install numpy`

All other Python-based dependencies should be installed automatically.

## Installation

Clone the Harmonisation repository

`$ git clone https://github.com/FIDUCEO/Harmonisation.git`

Change to src/main directory

`cd src/main`

Install with pip

`$ pip install harmonisation/`

NB: First make sure pip is upgraded to the latest version `pip install --upgrade pip`

## Usage

After install the following command-line scripts are available:

* `full_eiv_solver`
* `harm_plot`
* `harm_comp_plot`

For more details try:

`$ <script_name> --help`