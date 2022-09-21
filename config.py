"""Configutation for development"""

import sys
import os
from collections import defaultdict

###############################################################################
# Set NATIVE = True for a version optimized for your PC
# Set NATIVE = False for a generic version
NATIVE = False

###############################################################################
# Bit length of an integer used in the representation of the monster 
# This should be 32 for a 32-bit and 64 for a 64-bit system.
INT_BITS = 64

###############################################################################
# Representations of the monster are built modulo the following numbers
PRIMES = [3, 7, 15, 31, 127, 255]



###############################################################################
# Extra compile arguments for .c files to be used in setup.py
EXTRA_COMPILE_ARGS = defaultdict(list)

# Update the EXTRA_COMPILE_ARGS dictionary with entries
#     compiler : <list of options> .
# Here 'compiler' is a string. For a list of available compilers, run:
#     python setup.py build_ext --help-compiler

mingw32_compile_args =  ["-m64", "-Ofast", "-flto", "-funroll-loops"]
unix_compile_args = ["-Ofast", "-flto", "-funroll-loops"]
msvc_compile_args = ["/O2"]

if NATIVE:
     mingw32_compile_args.append("-march=native") 
     unix_compile_args.append("-march=native") 
     msvc_compile_args.append("/arch:AVX2")

EXTRA_COMPILE_ARGS.update({
    'mingw32' : mingw32_compile_args,
    'unix' : unix_compile_args,
    'msvc': msvc_compile_args,
})


###############################################################################
# Extra link arguments for .c files to be used in setup.py
EXTRA_LINK_ARGS = defaultdict(list)
# Similiar to updating EXTRA_COMPILE_ARGS, we may update 
# the EXTRA_LINK_ARGS dictionary with a list of pairs
#     compiler : <list of options> .


###############################################################################
# Directories

ROOT_DIR = "."
SRC_DIR = os.path.join(ROOT_DIR, "src")
PACKAGE_DIR = os.path.join(SRC_DIR, "mmgroup")
DEV_DIR = os.path.join(PACKAGE_DIR, "dev")

REAL_SRC_DIR = os.path.realpath(SRC_DIR)


C_DIR = os.path.join(DEV_DIR, "c_files")
PXD_DIR = os.path.join(DEV_DIR, "pxd_files")











