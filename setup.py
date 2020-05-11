from __future__ import absolute_import, division, print_function

import sys
import os
import subprocess
import numpy as np
from glob import glob
import shutil

import setuptools
from setuptools import setup, find_packages
from distutils.errors import *
from build_ext_steps import Extension, CustomBuildStep, SharedExtension
from build_ext_steps import BuildExtCmd




import config
from config import EXTRA_COMPILE_ARGS, EXTRA_LINK_ARGS
from config import ROOT_DIR, SRC_DIR, PACKAGE_DIR, DEV_DIR
from config import REAL_SRC_DIR
from config import C_DIR, DOC_DIR,  PXD_DIR
from config import PRIMES



####################################################################
# Delete files
####################################################################


# The following files are before building the extension
# if the command line option -f or --force has been set
ext_delete = [
    os.path.join("C_DIR", "*.*"),
    os.path.join("DOC_DIR", "*.*"),
    os.path.join("PXD_DIR", "*.*"),
    os.path.join(PACKAGE_DIR, "*.dll"), 
    os.path.join(PACKAGE_DIR, "*.pyd"),
]


def force_delete():
    """Delete some files before command 'build_ext'"""
    if not "-f" in sys.argv and not  "--force" in sys.argv:
        return
    for file_pattern in ext_delete:
        for file in glob(file_pattern):
            try:
                #print(file)
                os.remove(file)
            except:
                pass
   


####################################################################
# create directories
####################################################################

def make_dir(*args):
    """Create subdirectory if it does not exist

    The path is given by the arguments
    """
    directory = os.path.realpath(os.path.join(*args))
    if not os.path.exists(directory):
        os.makedirs(directory)
    fname = os.path.join(directory, "readme.txt")
    with open(fname, "wt") as f:
        f.write(
"""The files in this directory have been created automatically
or copied from some other place.
So it is safe to delete all files in this directory.
"""
        )   

####################################################################
# extend path
####################################################################

def extend_path():
    sys.path.append(REAL_SRC_DIR)

####################################################################
# Check if we are in a 'readthedocs' environment
####################################################################


on_readthedocs = os.environ.get('READTHEDOCS') == 'True'



####################################################################
# Desription of the list 'mat24_presteps'.
#
# This is a list of programs to be run before executing the 'build_ext' 
# command. Each entry of list 'custom_presteps' is a list which we call 
# a program lists. A program list ia a list of strings corresponding to 
# a program to be executed with:
#     subprocess.call(program_list) . 
# The first entry of a program list is the name of the program to be 
# executed; here sys.executable means the current python version. 
# Subsequents entries correspond to command line arguments.
#
# If the first entry in that list is not a string then it is 
# interpreted as a function to be called with the arguments
# given by the subsequent entries of that list.
####################################################################


pyx_sources = [
    os.path.join(DEV_DIR, "mat24", "mat24fast.pyx"),
    os.path.join(DEV_DIR, "mat24_xi", "mat24_xi.pyx"),
    os.path.join(DEV_DIR, "mm_basics", "mm_basics.pyx"),
]

def copy_pyx_sources():
    for filename in pyx_sources:
        shutil.copy(filename, PXD_DIR)

mat24_presteps = CustomBuildStep("Starting code generation",
  [make_dir, "src", "mmgroup", "dev", "c_files"],
  [make_dir, "src", "mmgroup", "dev", "c_doc"],
  [make_dir, "src", "mmgroup", "dev", "pxd_files"],
  [force_delete],
  [copy_pyx_sources],
  [extend_path],
  [sys.executable, "codegen_mat24.py"],
)




mat24_shared = SharedExtension(
    name = "mmgroup.mmgroup_mat24", 
    sources=[
        os.path.join(C_DIR, "mat24_functions.c"),
    ],
    libraries = [], 
    include_dirs = [PACKAGE_DIR, C_DIR],
    library_dirs = [PACKAGE_DIR, C_DIR],
    extra_compile_args = EXTRA_COMPILE_ARGS,
    implib_dir = C_DIR,
    define_macros = [ ("MAT24_DLL_EXPORTS", None)],
)


mat24_extension = Extension("mmgroup.mat24",
        sources=[
            os.path.join(PXD_DIR, "mat24fast.pyx"),
        ],
        #libraries=["m"] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [PACKAGE_DIR, C_DIR ],
        libraries = ["libmmgroup_mat24"], 
        #runtime_library_dirs = ["."],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
        extra_link_args = EXTRA_LINK_ARGS, 
)


mat24_xi_extension = Extension("mmgroup.mat24_xi",
        sources=[
            os.path.join(PXD_DIR, "mat24_xi.pyx"),
            os.path.join(C_DIR, "mat24_xi_functions.c"),
        ],
        #libraries=["m"] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [PACKAGE_DIR, C_DIR ],
        libraries = ["libmmgroup_mat24"], 
        #runtime_library_dirs = ["."],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
        extra_link_args = EXTRA_LINK_ARGS, 
    )

mm_presteps =  CustomBuildStep("Code generation for modules mm and mm_op",
  [sys.executable, "codegen_mm.py"],
  [sys.executable, "codegen_mm_op.py"],
)


mm_shared =  SharedExtension(
    name = "mmgroup.mmgroup_mm_basics", 
    sources=[ os.path.join(C_DIR, f) for f in 
        [ "mm_aux.c", "mm_group_n.c", "mm_random.c", "mm_sparse.c",
          "mm_tables.c","mm_tables_xi.c",
        ]
    ],    
    libraries = ["libmmgroup_mat24"], 
    include_dirs = [PACKAGE_DIR, C_DIR],
    library_dirs = [PACKAGE_DIR, C_DIR],
    extra_compile_args = EXTRA_COMPILE_ARGS,
    implib_dir = C_DIR,
    define_macros = [ ("MM_BASICS_DLL_EXPORTS", None)],
)

mm_extension = Extension("mmgroup.mm",
    sources=[
            os.path.join(PXD_DIR, "mm_basics.pyx"),
    ],
    #libraries=["m"] # Unix-like specific
    include_dirs = [ C_DIR ],
    library_dirs = [ PACKAGE_DIR, C_DIR ],
    libraries = ["libmmgroup_mm_basics", "libmmgroup_mat24"], 
            # for openmp add "libgomp" 
    #runtime_library_dirs = ["."],
    extra_compile_args = EXTRA_COMPILE_ARGS, 
            # for openmp add "-fopenmp" 
    extra_link_args = EXTRA_LINK_ARGS, 
            # for openmp add "-fopenmp" 
)



ext_modules = [
    mat24_presteps,
    mat24_shared, 
    mat24_extension,
    mat24_xi_extension,
    mm_presteps,
    mm_shared, 
    mm_extension, 
]





PYX_SOURCE_P = "mm_op{P}.pyx"

C_SOURCES_P = [
    "mm{P}_op_pi",
    "mm{P}_op_misc",
    "mm{P}_op_xy",
    "mm{P}_op_t",
    "mm{P}_op_xi",
    "mm{P}_op_word",
]

def list_source_files(p):
    sources = [os.path.join(PXD_DIR, PYX_SOURCE_P.format(P = p))]
    for f in C_SOURCES_P:
         sources.append(os.path.join(C_DIR, f.format(P = p) + ".c"))
    return sources

    
for p in PRIMES:
    ext_modules.append(
        Extension("mmgroup.mm%d" % p,
            sources = list_source_files(p),
            #libraries=["m"] # Unix-like specific
            include_dirs = [ C_DIR ], 
            library_dirs = [PACKAGE_DIR, C_DIR ],
            libraries = ["libmmgroup_mm_basics", "libmmgroup_mat24"], 
                # for openmp add "libgomp" 
            #runtime_library_dirs = ["."],
            extra_compile_args = EXTRA_COMPILE_ARGS, 
                # for openmp add "-fopenmp" 
            extra_link_args = EXTRA_LINK_ARGS, 
                # for openmp add "-fopenmp" 
        )
    )




test_step = CustomBuildStep("import_all",
  [sys.executable, "import_all.py"],
  ["pytest",  "src/mmgroup/", "-v", "-s", "-m", "build"],
)


ext_modules.append(test_step)


if on_readthedocs:
    ext_modules = [ ]

   

setup(
    name = 'mmgroup',    
    version = '0.0.1',    
    license='BSD-2-Clause',
    description='Construction of the sporadic simple monster group.',
    long_description='yet to be done',
    author='Martin Seysen',
    author_email='m.seysen@gmx.de',
    url='yet unknown',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 1 - Planning',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        #'Operating System :: Unix',
        #'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        #Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    project_urls={
       # 'Changelog': 'yet unknown',
       # 'Issue Tracker': 'yet unknown',
    },
    keywords=[
        'sporadic group', 'monster group', 'finite simple group'
    ],
    python_requires='>=3.6',
    install_requires=[
         'numpy', 'scipy', 
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    setup_requires=[
        'numpy', 'scipy', 'pytest-runner', 'Cython',
        # 'sphinx',  'sphinxcontrib-bibtex',
    ],
    tests_require=[
        'pytest',
    ],
    cmdclass={
        'build_ext': BuildExtCmd,
    },
    ext_modules = ext_modules,
    include_dirs=[np.get_include()],  # This gets all the required Numpy core files
)



