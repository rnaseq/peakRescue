from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'readToGeneAssignmentWithCython',
  ext_modules = cythonize("readToGeneAssignmentWithCython.pyx"),
)
