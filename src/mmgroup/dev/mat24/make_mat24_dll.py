from __future__ import absolute_import, division, print_function
#from __future__ import  unicode_literals

import sys
import os

#print(__file__)
FILE_DIR = os.path.dirname(os.path.realpath(__file__))
#print("FILE_DIR ==", FILE_DIR)
DEV_DIR = os.path.dirname(FILE_DIR)
#print("DEV_DIR ==", DEV_DIR)
PACKAGE_DIR = os.path.dirname(DEV_DIR)
#print("PACKAGE_DIR ==", PACKAGE_DIR)
SRC_DIR = os.path.dirname(PACKAGE_DIR)
#print("SCR_DIR ==", SRC_DIR)

sys.path.append(SRC_DIR)   # for all the stuff in the 'mmgroup' package
from mmgroup.dev import config
from mmgroup.dev.config import C_DIR , DOC_DIR, EXTRA_COMPILE_ARGS
from mmgroup.dev.generate_c.make_dll import make_dll
assert sys.path.pop() == SRC_DIR


make_dll(
    name = "mat24", 
    sources=[
        os.path.join(C_DIR, "mat24_functions.c"),
    ],
    libraries = [], 
    include_dirs = [".",C_DIR],
    lib_path = C_DIR,
    dll_path = PACKAGE_DIR,   
    extra_args = EXTRA_COMPILE_ARGS + ["-DMAT24_DLL_EXPORTS"]
    )
 



