from __future__ import absolute_import, division, print_function

import os
import sys

import setuptools
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize


#from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Build import cythonize

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
#print("FILE_DIR ==", FILE_DIR)
TEST_DIR = os.path.dirname(FILE_DIR)
#print("TEST_DIR ==", TEST_DIR)
PACKAGE_DIR = os.path.dirname(TEST_DIR)
#print("PACKAGE_DIR ==", PACKAGE_DIR)
SRC_DIR = os.path.dirname(PACKAGE_DIR)


sys.path.append(SRC_DIR)   # for all the stuff in the 'mmgroup' package
from mmgroup.tests import test_hadamard
from mmgroup.tests.tempdir import temp_dir_name
assert sys.path.pop() == SRC_DIR

C_subdir = temp_dir_name(test_hadamard)


ext_modules=[
    Extension("test_hadamard",
        sources=[
           "test_hadamard.pyx",
           "test_hadamard_op.c",
        ],
        #libraries=["m"], # Unix-like specific
        include_dirs = [ "." ],
    )
]

setup(
  name = "test_hadamard",
  ext_modules =  cythonize(ext_modules)
)

