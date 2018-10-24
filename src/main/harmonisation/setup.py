from setuptools import find_packages, setup

from os import path
from io import open
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "harmonisation.core.matchup.matchupProcessing.sample2ind.cython_lib_sample2ind",
        ["harmonisation/core/matchup/matchupProcessing/sample2ind/cython_lib_sample2ind.pyx"],
    ),
    Extension(
        "harmonisation.core.matchup.matchupToolbox.harmonisation.harmonisationProcessing.harmonisation_eiv.cython_lib_GN_algo",
        ["harmonisation/core/matchup/matchupToolbox/harmonisation/harmonisationProcessing/harmonisation_eiv/cython_lib_GN_algo.pyx"],
    )
]
# Get package __version__.
# Same effect as "from s import __version__",
# but avoids importing the module which may not be installed yet:
version = None
here = path.abspath(path.dirname(__file__))
with open('harmonisation/version.py') as f:
    exec(f.read())
    __version__ = version

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='harmonisation',
      version=__version__,
      description='FIDUCEO sensor series harmonisation software',
      author='Sam Hunt',
      long_description=long_description,
      author_email='sam.hunt@npl.co.uk',
      url='http://www.fiduceo.eu',
      keywords="FIDUCEO FCDR harmonisation metrology climate",
      packages=find_packages(exclude=['contrib', 'docs', 'tests', 'lsf']),
      ext_modules=cythonize(extensions),
      include_dirs=[np.get_include()],
      package_data={"": ["*.dat"]},
      install_requires=['numpy>=1.11.0', 'netCDF4>=1.1.0', 'scipy>=0.19', 'matplotlib>=2.2.2'],
      entry_points={
          'console_scripts': [
              'full_eiv_solver =  harmonisation.tools.full_eiv_solver:main',
              'harm_plot =  harmonisation.tools.harm_plot:main',
              'harm_comp_plot =  harmonisation.tools.harm_comp_plot:main',
          ],
})