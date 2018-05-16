from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(ext_modules=cythonize("cython_lib_PC_algo.pyx"), include_path=[numpy.get_include()])
