from __future__ import absolute_import, division, print_function


####################################################################
# History
####################################################################

VERSION = '0.0.7' # 2021-12-01. Bugfix in version generation
# VERSION = '0.0.6' # 2021-12-01. Group operation accelerated
# VERSION = '0.0.5' # 2021-08-02. Word shortening in monster implemented
# VERSION = '0.0.4' # 2020-06-15. MSVC compiler is now supported
# VERSION = '0.0.3' # 2020-06-10. bugfixes in code generator
# VERSION = '0.0.2' # 2020-06-04. Order oracle added; bugfixes
# VERSION = '0.0.1' # 2020-05-20. First releae

# Version history must also be updated in the API Reference,
# section 'Version history'.

####################################################################
# Imports
####################################################################


import sys
import os
import re
import subprocess
import numpy as np
from glob import glob
import shutil

import setuptools
from setuptools import setup, find_packages
from collections import defaultdict
from distutils.errors import *
from build_ext_steps import Extension, CustomBuildStep, SharedExtension
from build_ext_steps import BuildExtCmd




import config
from config import EXTRA_COMPILE_ARGS, EXTRA_LINK_ARGS
from config import ROOT_DIR, SRC_DIR, PACKAGE_DIR, DEV_DIR
from config import REAL_SRC_DIR
from config import C_DIR, PXD_DIR
from config import PRIMES



####################################################################
# Delete files
####################################################################


# The following files are before building the extension
# if the command line option -f or --force has been set
ext_delete = [
    os.path.join("C_DIR", "*.*"),
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


codegen_args = ["mockup"] if on_readthedocs else []

####################################################################
# Set path for shared libraries in linux
####################################################################

if not on_readthedocs and os.name == "posix":    
    old_ld_path = os.getenv("LD_LIBRARY_PATH")
    old_ld_path = old_ld_path + ";" if old_ld_path else ""
    new_LD_LIBRARY_PATH = os.path.abspath(PACKAGE_DIR)
    os.environ["LD_LIBRARY_PATH"] =  old_ld_path + new_LD_LIBRARY_PATH 


####################################################################
# Add extensions and shared libraries to package data
####################################################################


if os.name in ["nt"]:
    extension_wildcards =  ["*.pyd", "*.dll"]     
elif os.name in ["posix"]:
    extension_wildcards =  ["*.so"]  
else:   
    extension_wildcards =  []  


package_data = {
        # If any package contains *.txt or *.rst files, include them:
        "mmgroup": extension_wildcards
}




####################################################################
# Desription of the list 'general_presteps'.
#
# This is a list of programs to be run before executing the 'build_ext' 
# command. Each entry of list 'custom_presteps' is a list which we call 
# a program list. A program list ia a list of strings corresponding to 
# a program to be executed with:
#     subprocess.check_call(program_list) . 
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
    os.path.join(DEV_DIR, "generators", "generators.pyx"),
    os.path.join(DEV_DIR, "mm_basics", "mm_basics.pyx"),
    os.path.join(DEV_DIR, "clifford12", "clifford12.pyx"),
]

def copy_pyx_sources():
    for filename in pyx_sources:
        shutil.copy(filename, PXD_DIR)

    

general_presteps = CustomBuildStep("Starting code generation",
  [make_dir, "src", "mmgroup", "dev", "c_files"],
  [make_dir, "src", "mmgroup", "dev", "pxd_files"],
  [force_delete],
  [copy_pyx_sources],
  [extend_path],
)



####################################################################
# List of libraries to be used an a code generation stage.
#
# This list contains the shared libraries to be used in a code
# generation stage. We have to divide the code generation process 
# into stages, since a library built in a certain stage may be 
# for generating the code used in a subsequent stage.
#
# For windows we list the import libraries instead.
####################################################################


if on_readthedocs:
    shared_libs_before_stage1 = []
    shared_libs_stage1 = shared_libs_stage2 = []
elif os.name in ["nt"]:
    shared_libs_before_stage1 = [
        "libmmgroup_mat24", 
    ]
    shared_libs_stage1 = shared_libs_before_stage1 + [
        "libmmgroup_clifford12",
    ]
    shared_libs_stage2 = shared_libs_stage1 + [
        "libmmgroup_mm_basics",
   ]
elif os.name in ["posix"]:
    shared_libs_before_stage1 = [
        "mmgroup_mat24", 
    ]
    shared_libs_stage1 = shared_libs_before_stage1 + [
        "mmgroup_clifford12",
    ]
    shared_libs_stage2 = shared_libs_stage1 + [
        "mmgroup_mm_basics",
    ]
else:
    raise DistutilsPlatformError(
        "I don't know how to build to the shared libraries "
        "in the '%s' operating system" % os.name
    )


####################################################################
# Building the extenstions at stage 1
####################################################################



mat24_presteps = CustomBuildStep("Generating code for extension 'mat24'",
  [sys.executable, "codegen_mat24.py"],
  [sys.executable, "codegen_clifford12.py"],
)

mat24_shared = SharedExtension(
    name = "mmgroup.mmgroup_mat24", 
    sources=[
        os.path.join(C_DIR, "mat24_functions.c"),
        os.path.join(C_DIR, "gen_xi_functions.c"),
        os.path.join(C_DIR, "mm_group_n.c"),
        os.path.join(C_DIR, "gen_leech.c"),
        os.path.join(C_DIR, "gen_leech3.c"),
        os.path.join(C_DIR, "gen_leech_reduce.c"),
        os.path.join(C_DIR, "gen_random.c"),
    ],
    libraries = [], 
    include_dirs = [PACKAGE_DIR, C_DIR],
    library_dirs = [PACKAGE_DIR, C_DIR],
    extra_compile_args = EXTRA_COMPILE_ARGS,
    implib_dir = C_DIR,
    define_macros = [ ("MAT24_DLL_EXPORTS", None)],
)


clifford12_shared = SharedExtension(
    name = "mmgroup.mmgroup_clifford12", 
    sources = [
        os.path.join(C_DIR, "bitmatrix64.c"),
        os.path.join(C_DIR, "qstate12.c"),
        os.path.join(C_DIR, "qmatrix12.c"),
        os.path.join(C_DIR, "xsp2co1.c"),
        os.path.join(C_DIR, "leech3matrix.c"),
        os.path.join(C_DIR, "xsp2co1_elem.c"),
        os.path.join(C_DIR, "involutions.c"),
    ],
    include_dirs = [PACKAGE_DIR, C_DIR],
    library_dirs = [PACKAGE_DIR, C_DIR],
    libraries = shared_libs_before_stage1, 
    extra_compile_args = EXTRA_COMPILE_ARGS,
    implib_dir = C_DIR,
    define_macros = [ ("CLIFFORD12_DLL_EXPORTS", None)],
)

mat24_extension = Extension("mmgroup.mat24",
        sources=[
            os.path.join(PXD_DIR, "mat24fast.pyx"),
        ],
        #libraries=["m"] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [PACKAGE_DIR, C_DIR ],
        libraries = shared_libs_stage1, 
        #runtime_library_dirs = ["."],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
        extra_link_args = EXTRA_LINK_ARGS, 
)

generators_extension = Extension("mmgroup.generators",
        sources=[
            os.path.join(PXD_DIR, "generators.pyx"),
        ],
        #libraries=["m"] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [PACKAGE_DIR, C_DIR ],
        libraries = shared_libs_stage1, 
        #runtime_library_dirs = ["."],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
        extra_link_args = EXTRA_LINK_ARGS, 
)


clifford12_extension =  Extension("mmgroup.clifford12",
        sources=[
            os.path.join(PXD_DIR, "clifford12.pyx"),
        ],
        #libraries=["m"] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [PACKAGE_DIR, C_DIR ],
        libraries = shared_libs_stage1, 
        #runtime_library_dirs = ["."],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
        extra_link_args = EXTRA_LINK_ARGS, 
)


####################################################################
# Building the extenstions at stage 2
####################################################################



mm_presteps =  CustomBuildStep("Code generation for modules mm and mm_op",
  [sys.executable, "codegen_mm.py"] + codegen_args,
  [sys.executable, "codegen_mm_op.py"] + codegen_args,
)


mm_shared =  SharedExtension(
    name = "mmgroup.mmgroup_mm_basics", 
    sources=[ os.path.join(C_DIR, f) for f in 
        [ "mm_aux.c",  "mm_group_word.c",
          "mm_tables.c","mm_tables_xi.c", "mm_crt.c",
        ]
    ],    
    libraries = shared_libs_stage1, 
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
    libraries = shared_libs_stage2, 
            # for openmp add "libgomp" 
    #runtime_library_dirs = ["."],
    extra_compile_args = EXTRA_COMPILE_ARGS, 
            # for openmp add "-fopenmp" 
    extra_link_args = EXTRA_LINK_ARGS, 
            # for openmp add "-fopenmp" 
)



ext_modules = [
    general_presteps,
    mat24_presteps,
    mat24_shared,
    clifford12_shared, 
    mat24_extension,
    generators_extension,
    clifford12_extension,
    mm_presteps,
    mm_shared, 
    mm_extension, 
]



####################################################################
# Adding the extension for the final stage.
#
# Here we build the representation of the monster modulo
# all small primes P in the set PRIMES. 
####################################################################


PYX_SOURCE_P = "mm_op{P}.pyx"

C_SOURCES_P = [
    "mm{P}_op_pi",
    "mm{P}_op_misc",
    "mm{P}_op_xy",
    "mm{P}_op_t",
    "mm{P}_op_xi",
    "mm{P}_op_word",
]


# Additions for the list C_SOURCES_P for spcific values p
C_SOURCES_SPECIFIC_P = defaultdict(list)
C_SOURCES_SPECIFIC_P.update( {
    15: [ "mm15_op_order", "mm15_op_axes_find", 
          "mm15_op_axes", "mm15_op_eval_A",
        ],
} )



def list_source_files(p):
    sources = [os.path.join(PXD_DIR, PYX_SOURCE_P.format(P = p))]
    for f in C_SOURCES_P + C_SOURCES_SPECIFIC_P[p]:
         sources.append(os.path.join(C_DIR, f.format(P = p) + ".c"))
    return sources

    
for p in PRIMES:
    ext_modules.append(
        Extension("mmgroup.mm%d" % p,
            sources = list_source_files(p),
            #libraries=["m"] # Unix-like specific
            include_dirs = [ C_DIR ] , 
            library_dirs = [PACKAGE_DIR, C_DIR ],
            libraries = shared_libs_stage2, 
                # for openmp add "libgomp" 
            #runtime_library_dirs = ["."],
            extra_compile_args = EXTRA_COMPILE_ARGS, 
                # for openmp add "-fopenmp" 
            extra_link_args = EXTRA_LINK_ARGS, 
                # for openmp add "-fopenmp" 
        )
    )





####################################################################
# After building the externals we add a tiny little test step.
####################################################################



test_step = CustomBuildStep("import_all",
  [sys.executable, "import_all.py"],
  ["pytest",  "src/mmgroup/", "-v", "-s", "-m", "build"],
)


ext_modules.append(test_step)



####################################################################
# Don't build any externals when building the documentation.
####################################################################


if on_readthedocs:
    ext_modules = [ 
        general_presteps,
        mat24_presteps,
        mm_presteps,
    ]

   

def read(fname):
    """Return the text in the file with name 'fname'""" 
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


####################################################################
# The main setup program.
####################################################################


### Prelimiary!!!
"""
ext_modules = [
    mat24_presteps,
    clifford12_shared, 
    clifford12_extension,
]
"""

setup(
    name = 'mmgroup',    
    version = VERSION,    
    license='BSD-2-Clause',
    description='Implementation of the sporadic simple monster group.',
    long_description=read('README.rst'),
    author='Martin Seysen',
    author_email='m.seysen@gmx.de',
    url='https://github.com/Martin-Seysen/mmgroup',
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
        #'Programming Language :: Python :: 3.4',
        #'Programming Language :: Python :: 3.5',
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
         'numpy',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    setup_requires=[
        'numpy', 'scipy', 'pytest-runner', 'cython',
        # 'sphinx',  'sphinxcontrib-bibtex',
    ],
    tests_require=[
        'pytest', 'scipy', 
    ],
    cmdclass={
        'build_ext': BuildExtCmd,
    },
    ext_modules = ext_modules,
    package_data = package_data,
    include_dirs=[np.get_include()],  # This gets all the required Numpy core files
)


####################################################################
# In linux, use 'auditwheel' afterwards, to build the final wheel.
#
# The current version is still erroneous.
####################################################################



def build_posix_wheel():   
    assert  os.name == "posix" 
    if not on_readthedocs and "bdist_wheel" in sys.argv:
        PROJECT_NAME = r"mmgroup"
        SUFFIX_MATCH = r"[-0-9A-Za-z._]+linux[-0-9A-Za-z._]+\.whl"
        DIST_DIR = "dist"
        w_match = re.compile(PROJECT_NAME + SUFFIX_MATCH)
        wheels = [s for s in os.listdir(DIST_DIR) if w_match.match(s)]
        for wheel in wheels:
            wheel_path = os.path.join(DIST_DIR, wheel)
            args = ["auditwheel", "-v", "repair", wheel_path]
            print(" ".join(args))
            subprocess.check_call(args)



if os.name == "posix":
    ## build_posix_wheel()  # This does not work
    pass


