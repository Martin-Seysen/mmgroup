"""Configutation for development"""

import sys
import os

# Bit length of an integer used in the representation of the monster 
# This should be 32 for a 32-bit and 64 for a 64-bit system.
INT_BITS = 64

# Representations of the monster are built modulo the following numbers
PRIMES = [3,7,15,127]

# Extra compile arguments for .c files for Cython to be used in setup.py
EXTRA_COMPILE_ARGS = [
    "-m64", "-Ofast", "-flto", "-march=native", "-funroll-loops"
]




# directory names
DEV_DIR = os.path.dirname(os.path.realpath(__file__))
PACKAGE_DIR = os.path.dirname(DEV_DIR)
SRC_DIR = os.path.dirname(PACKAGE_DIR)

C_DIR = os.path.join(DEV_DIR, "c_files")
DOC_DIR = os.path.join(DEV_DIR, "c_doc")









