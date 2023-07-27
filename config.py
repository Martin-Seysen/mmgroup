"""Configutation for development"""

import sys
import os
from collections import defaultdict

###############################################################################
# Set NATIVE = True for a version optimized for your PC
# Set NATIVE = False for a generic version
NATIVE = False


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












