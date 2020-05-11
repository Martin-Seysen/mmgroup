from __future__ import absolute_import, division, print_function
#from __future__ import  unicode_literals

import os
import sys

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

sys.path.append("..")
import config
from config import C_DIR, EXTRA_COMPILE_ARGS

from codegen import  PYX_FILE_NAME, list_all_c_files

C_sources =  list_all_c_files()



ext_C_files = ["mat24_functions.c", 
               "mat24_xi_functions.c",
              ]

ext_sources = [os.path.join(C_DIR, f) for f in ext_C_files]

sources = list(map(str, [PYX_FILE_NAME] + C_sources + ext_sources))


ext_modules=[
    Extension("mm_op",
        sources = sources,
        #libraries=["m"] # Unix-like specific
        include_dirs = ["./", C_DIR ],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
    )
]

setup(
  name = "mm_op",
  ext_modules = cythonize(ext_modules)
)