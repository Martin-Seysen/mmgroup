####################################################################
# History
####################################################################

VERSION = '1.0.2' # 2024-02-27. Changes in build system regarding shared libraries
# VERSION = '1.0.1' # 2024-02-23. Support for exporting C functions from libraries
# VERSION = '1.0.0' # 2024-01-23. Documentation updated
# VERSION = '0.0.14' # 2024-01-18. Demonstration code for reduction in the Monster added
# VERSION = '0.0.13' # 2023-10-27. Supports numbering of elements of the Monster.
                     #             Adjusted for Python 3.12 
# VERSION = '0.0.12' # 2023-01-09. Support for the group Y_555 and the Bimonster
# VERSION = '0.0.11' # 2022-10-19. Bugfixes and macOS version added
# VERSION = '0.0.10' # 2022-10-11. Support for cibuildwheel added
# VERSION = '0.0.9' # 2022-09-01. Performance improved
# VERSION = '0.0.8' # 2022-07-12. Performance improved
# VERSION = '0.0.7' # 2021-12-01. Bugfix in version generation
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
import time
import subprocess
import numpy as np
from glob import glob

import setuptools
from setuptools import setup, find_namespace_packages
from collections import defaultdict



######################################################################
# Directories and inports relative to these driectories
######################################################################

ROOT_DIR = os.path.realpath(os.path.dirname(__file__))
SRC_DIR =  os.path.realpath(os.path.join(ROOT_DIR, 'src'))
PACKAGE_DIR = os.path.join(SRC_DIR, 'mmgroup')
DEV_DIR = os.path.join(PACKAGE_DIR, 'dev')
C_DIR = os.path.join(DEV_DIR, 'c_files')
LIB_DIR = os.path.join(DEV_DIR, 'lib_files')
PXD_DIR = os.path.join(DEV_DIR, 'pxd_files')
HEADERS_DIR = os.path.join(DEV_DIR, 'headers')


sys.path.append(ROOT_DIR)
sys.path.append(SRC_DIR)



####################################################################
# Print patform and command line arguments (if desired)
####################################################################


def print_platform():
    import platform
    print('Build running on platform')
    uname = platform.uname()
    for key in [
        'system', 'node', 'release', 'version', 'machine'
        ]:
        print(" " + key + ":", getattr(uname, key, None))
    import numpy as np
    print("numpy version:", np.version.version)
    print("")


def print_commandline_args():
    print('Command line arguments of setup.py (in mmgroup project):')
    for arg in sys.argv:
        print(' ' + arg)
    print('Current working directory os.path.getcwd():')
    print(' ' + os.getcwd())
    print('Absolute path of file setup.py:')
    print(' ' + os.path.abspath(__file__))
    print('')

print_platform()
print_commandline_args()

####################################################################
# Global options
####################################################################

STAGE = 1
STATIC_LIB = False
NPROCESSES = 16
COMPILER = None
CFLAGS = None
MOCKUP = False
VERBOSE = False
# Parse a global option '--stage=i', '--compiler=c', and set variable 
# ``STAGE`` to the integer value i if such an option is present.
for i, s in enumerate(sys.argv[1:]):
    if s.startswith('--n='):
        NPROCESSES = int(s[4:])
        sys.argv[i+1] = None
    if s.startswith('--stage='):
        STAGE = int(s[8:])
        sys.argv[i+1] = None
    elif s.startswith('--compiler='):
        COMPILER = s[11:]
        sys.argv[i+1] = None
    elif s.startswith('--cflags='):
        CFLAGS = s[9:]
        sys.argv[i+1] = None
    elif s.startswith('--static'):
        STATIC_LIB = True
        sys.argv[i+1] = None
    elif s == '--mockup':
        MOCKUP = True
        sys.argv[i+1] = None
    elif s == '-v':
        VERBOSE = True
        sys.argv[i+1] = None
    elif s[:1].isalpha:
        break
while None in sys.argv: 
    sys.argv.remove(None)

if COMPILER and COMPILER not in ['unix','msvc', 'mingw32']:
    raise ValueError("Unknown compiler '%s'" % COMPILER)


####################################################################
# Check if we are in a 'readthedocs' environment
####################################################################

on_readthedocs = MOCKUP or os.environ.get('READTHEDOCS') == 'True'


####################################################################
# Remove old stuff and import build tools
####################################################################


# Remove old shared libraries before(!) anybody else can grab them!!!

CLEANUP_ARGS = defaultdict(list)
CLEANUP_ARGS.update({
    1: ['-pcxs'],
    2: ['--rm-ext', 'mm_op', 'mm_reduce'],
    3: ['--rm-ext', 'mm_reduce'],
})

subprocess.check_call([sys.executable, 'cleanup.py',
   '--check-uninstalled'] +  CLEANUP_ARGS[STAGE])


####################################################################
# import build tools
####################################################################



from mmgroup.generate_c.build_ext_steps import Extension, CustomBuildStep
from mmgroup.generate_c.build_ext_steps import BuildExtCmd
from mmgroup.generate_c.build_shared import shared_lib_name
from mmgroup.generate_c.copy_shared import copy_shared_libs
from config import EXTRA_COMPILE_ARGS, EXTRA_LINK_ARGS


####################################################################
# Add extensions and shared libraries to package data
####################################################################

header_wildcards = ['*.h']
extension_wildcards =  []
if os.name in ['nt']:
    extension_wildcards =  ['*.pyd', '*.dll']
    header_wildcards = ['*.h', '*.lib']
#elif os.name in ['posix']:
#    extension_wildcards =  ['*.so']


package_data = {
        # If any package contains files as given above, include them:
        'mmgroup': extension_wildcards,
        'mmgroup.dev.headers': header_wildcards,
}





####################################################################
# Initialize list of external build operations
####################################################################



ext_modules = []


####################################################################
# We have to divide the code generation process 
# into stages, since a library built in a certain stage may be 
# for generating the code used in a subsequent stage.
####################################################################

DIR_DICT = {
   'SRC_DIR' : SRC_DIR,
   'C_DIR' : C_DIR,
   'LIB_DIR' : LIB_DIR,
   'DEV_DIR' : DEV_DIR,
   'PXD_DIR' : PXD_DIR,
   'PACKAGE_DIR': PACKAGE_DIR,
   'STATIC_LIB' : int(STATIC_LIB),
   'NPROCESSES' : str(NPROCESSES),
   'HEADERS_DIR' : HEADERS_DIR,
}

DIR_DICT['MOCKUP'] = '--mockup\n' if on_readthedocs else ''
DIR_DICT['COMPILER'] = '--compiler %s\n' % COMPILER if COMPILER else ''
DIR_DICT['CFLAGS'] = '--cflags=' + CFLAGS + "\n" if CFLAGS else ''
DIR_DICT['VERBOSE'] = '-v\n' if VERBOSE else ""

GENERATE_START = '''
 {VERBOSE}
 {MOCKUP}
 --py-path {SRC_DIR}
 --out-dir {C_DIR}
 --out-pxd-dir {PXD_DIR}
 --library-path  {PACKAGE_DIR}
'''.format(**DIR_DICT)

SHARED_START = '''
    {COMPILER}
    {CFLAGS}
    {MOCKUP}
    --source-dir {C_DIR}
    --include-path {PACKAGE_DIR} {LIB_DIR}  
    --library-path {PACKAGE_DIR} {LIB_DIR}
    --library-dir {LIB_DIR}
    --shared-dir {PACKAGE_DIR}
    --rpath $ORIGIN
    --define
    --static {STATIC_LIB}
    --n {NPROCESSES}
'''.format(**DIR_DICT)



####################################################################
# Building the extenstions at stage 1
####################################################################

DIR_DICT["DLL_NAME"] = "None" if STATIC_LIB else "MAT24"


MAT24_SOURCES = '''
   mat24_tables.c
   mat24_functions.c
   mat24_random.c
   check_endianess.c
'''

MAT24_GENERATE = GENERATE_START + '''
 --dll {DLL_NAME}
 --source-path {SRC_DIR}/mmgroup/dev/mat24
 --tables mmgroup.dev.mat24.mat24_ref 
 --sources mat24_functions.h
 --sources
'''.format(**DIR_DICT) + MAT24_SOURCES + '''
 --pxd  mat24_functions.pxd
 --pyx  mat24fast.pyx
'''

GENERATORS_SOURCES = '''
   gen_xi_functions.c mm_group_n.c gen_leech.c 
   gen_leech_type.c gen_leech3.c gen_leech_reduce.c
   gen_leech_reduce_22.c
   gen_leech_reduce_n.c gen_random.c
   gen_union_find.c
   gen_leech_reduce_mod3.c
'''


GENERATORS_GENERATE = GENERATE_START + '''
 --dll {DLL_NAME}
 --source-path {SRC_DIR}/mmgroup/dev/generators
 --tables mmgroup.dev.generators.gen_cocode_short
          mmgroup.dev.generators.gen_leech_reduce_n 
          mmgroup.dev.generators.gen_xi_ref
          mmgroup.dev.generators.gen_leech_22
          mmgroup.dev.generators.gen_leech_reduce_mod3
 --sources mmgroup_generators.h
 --sources 
'''.format(**DIR_DICT) + GENERATORS_SOURCES + '''
 --pxd  generators.pxd
 --pxi
 --pyx generators.pyx
'''


CLIFFORD12_SOURCES = '''
  qstate12.c qstate12io.c qmatrix12.c
  bitmatrix64.c uint_sort.c xsp2co1.c
  xsp2co1_word.c 
  leech2matrix.c leech3matrix.c xsp2co1_elem.c
  involutions.c xsp2co1_traces.c
  xsp2co1_map.c
'''


CLIFFORD12_GENERATE = GENERATE_START + '''
 --dll {DLL_NAME}
 --source-path {SRC_DIR}/mmgroup/dev/clifford12
 --tables mmgroup.dev.clifford12.bit64_tables
 --sources clifford12.h
 --sources  
'''.format(**DIR_DICT) + CLIFFORD12_SOURCES + '''
 --pxd  clifford12.pxd
 --pxi 
 --pyx  clifford12.pyx
'''


STAGE1_SOURCES = MAT24_SOURCES + GENERATORS_SOURCES + CLIFFORD12_SOURCES
STAGE1_LIBS = ['mmgroup_mat24']
STAGE1_BUILD_LIBS = shared_lib_name(STAGE1_LIBS, 'build', 
                   static=STATIC_LIB)
STAGE1_BUILD_EXT = shared_lib_name(STAGE1_LIBS, 'build_ext', 
                   static=STATIC_LIB)


MAT24_SHARED = SHARED_START + '''
    --name mmgroup_mat24 
    --sources 
''' + STAGE1_SOURCES + '''
    --libraries 
''' + " ".join([])


mat24_steps = CustomBuildStep(
   'Build libraries for extensions mat24, generate. and clifford12',
   [sys.executable, 'generate_code.py'] + MAT24_GENERATE.split(),
   [sys.executable, 'generate_code.py'] + GENERATORS_GENERATE.split(),
   [sys.executable, 'generate_code.py'] + CLIFFORD12_GENERATE.split(),
   [sys.executable, 'build_shared.py'] + MAT24_SHARED.split(),
)

mat24_extension = Extension('mmgroup.mat24',
        sources=[
            os.path.join(PXD_DIR, 'mat24fast.pyx'),
        ],
        #libraries=['m'] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [PACKAGE_DIR, LIB_DIR ],
        libraries = STAGE1_BUILD_EXT, 
        #runtime_library_dirs = ['.'],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
        extra_link_args = EXTRA_LINK_ARGS, 
)

generators_extension = Extension('mmgroup.generators',
        sources=[
            os.path.join(PXD_DIR, 'generators.pyx'),
        ],
        #libraries=['m'] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [PACKAGE_DIR, LIB_DIR ],
        libraries = STAGE1_BUILD_EXT, 
        #runtime_library_dirs = ['.'],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
        extra_link_args = EXTRA_LINK_ARGS, 
)

clifford12_extension =  Extension('mmgroup.clifford12',
        sources=[
            os.path.join(PXD_DIR, 'clifford12.pyx'),
        ],
        #libraries=['m'] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [PACKAGE_DIR, LIB_DIR ],
        libraries = STAGE1_BUILD_EXT, 
        #runtime_library_dirs = ['.'],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
        extra_link_args = EXTRA_LINK_ARGS, 
)



MAKE_ENDIANESS_ARGS = os.path.join(C_DIR,'mmgroup_endianess.h') + """
 --library-path {PACKAGE_DIR} {LIB_DIR}
 {MOCKUP}
 -v
""".format(**DIR_DICT)

if STAGE < 2:
    ext_modules += [
        mat24_steps,
    ]
    if not on_readthedocs:
        ext_modules += [
            mat24_extension,
            generators_extension,
            clifford12_extension,
        ]
    ext_modules +=  [CustomBuildStep(
        'Generate header describing the endianess of the machine',
        [sys.executable, 'make_endianess_header.py'
            ] + MAKE_ENDIANESS_ARGS.split()
    )]


####################################################################
# Building the extensions at stage 2
####################################################################

DIR_DICT["DLL_NAME"] = "None" if STATIC_LIB else "MM_OP"


MM_SOURCES = '''
    mm_index.c mm_aux.c mm_tables.c mm_group_word.c
    mm_tables_xi.c mm_crt.c
'''

MM_GENERATE = GENERATE_START + '''
 --dll {DLL_NAME}
 --source-path {SRC_DIR}/mmgroup/dev/mm_basics
               {SRC_DIR}/mmgroup/dev/mm_op
 --tables mmgroup.dev.mm_basics.mm_basics
          mmgroup.dev.mm_basics.mm_tables_xi
          mmgroup.dev.mm_basics.mm_aux
          mmgroup.dev.mm_basics.mm_tables
          mmgroup.dev.mm_basics.mm_crt
 --sources mm_basics.h
 --sources
'''.format(**DIR_DICT) + MM_SOURCES + '''
 --pxd mm_basics.pxd
 --pxi
'''


MM_OP_SUB_SOURCES = ''

for p in [3, 7, 15, 31, 127, 255]:
   MM_OP_SUB_SOURCES += '''
      mm{p}_op_misc.c
      mm{p}_op_pi.c                
      mm{p}_op_xy.c
      mm{p}_op_t.c
      mm{p}_op_xi.c
      mm{p}_op_word.c
      mm{p}_op_scalprod.c
      '''.format(p=p)

for p in [3, 15]:
   MM_OP_SUB_SOURCES += '''
      mm{p}_op_rank_A.c
      mm{p}_op_eval_A.c
      '''.format(p=p)

for p in [15]:
   MM_OP_SUB_SOURCES += '''
      mm{p}_op_eval_X.c
      '''.format(p=p)

MM_OP_SUB_GENERATE = GENERATE_START + '''
 --dll {DLL_NAME}
 --source-path {SRC_DIR}/mmgroup/dev/mm_op
 --subst mm(?P<p>[0-9]+)_op mm_op
 --tables mmgroup.dev.mm_op.mm_op
          mmgroup.dev.mm_op.mm_op_xi
          mmgroup.dev.mm_op.mm_op_pi
          mmgroup.dev.mm_op.mm_op_xy
          mmgroup.dev.hadamard.hadamard_t
          mmgroup.dev.hadamard.hadamard_xi
          mmgroup.dev.mm_op.mm_op_scalprod
 --sources mm_op_sub.h
 --sources
'''.format(**DIR_DICT) + MM_OP_SUB_SOURCES



MM_OP_P_SOURCES = '''
    mm_op_p_vector.c mm_op_p_axis.c
'''

MM_OP_P_GENERATE = GENERATE_START + '''
 --dll {DLL_NAME} -v
 --source-path {SRC_DIR}/mmgroup/dev/mm_op
 --set C_DIR={C_DIR}
 --tables mmgroup.dev.mm_basics.mm_basics
          mmgroup.dev.mm_op.dispatch_p
 --sources mm_op_p.h
 --sources  
'''.format(**DIR_DICT) + MM_OP_P_SOURCES + '''
 --pxd  mm_op_p.pxd
 --pxi
 --pyx  mm_op_p.pyx
'''


STAGE2_SOURCES = MM_SOURCES + MM_OP_SUB_SOURCES + MM_OP_P_SOURCES
STAGE2_LIBS = STAGE1_LIBS + ['mmgroup_mm_op']
STAGE2_BUILD_LIBS = shared_lib_name(STAGE2_LIBS, 'build', 
                    static=STATIC_LIB)
STAGE2_BUILD_EXT = shared_lib_name(STAGE2_LIBS, 'build_ext', 
                   static=STATIC_LIB)


MM_OP_SHARED = SHARED_START + '''
    --name mmgroup_mm_op 
    --sources 
''' + STAGE2_SOURCES + '''
    --libraries 
''' + " ".join(STAGE1_BUILD_LIBS)

    

mm_op_steps =  CustomBuildStep(
  'Build shared library and mm_op',
   [sys.executable, 'generate_code.py'] + MM_GENERATE.split(),
   [sys.executable, 'generate_code.py'] + MM_OP_SUB_GENERATE.split(),
   [sys.executable, 'generate_code.py'] + MM_OP_P_GENERATE.split(),
   [sys.executable, 'build_shared.py'] + MM_OP_SHARED.split(),
)

mm_op_extension = Extension('mmgroup.mm_op',
    sources=[
            os.path.join(PXD_DIR, 'mm_op_p.pyx'),
    ],
    #libraries=['m'] # Unix-like specific
    include_dirs = [ C_DIR ],
    library_dirs = [ PACKAGE_DIR, LIB_DIR ],
    libraries = STAGE2_BUILD_EXT, 
            # for openmp add 'libgomp' 
    #runtime_library_dirs = ['.'],
    extra_compile_args = EXTRA_COMPILE_ARGS, 
            # for openmp add '-fopenmp' 
    extra_link_args = EXTRA_LINK_ARGS, 
            # for openmp add '-fopenmp' 
)

mm_op_poststeps =  CustomBuildStep(
   'Create substituts for legacy extensions',
   [sys.executable, 'make_legacy_extensions.py', '--out-dir',
       os.path.join(SRC_DIR, 'mmgroup')
   ] 
)


if STAGE < 3:
    ext_modules += [
        mm_op_steps,
    ]
    if not on_readthedocs:
        ext_modules += [
            mm_op_extension,
            mm_op_poststeps, 
        ]


####################################################################
# Building the extensions at stage 3
####################################################################

DIR_DICT["DLL_NAME"] = "None" if STATIC_LIB else "MM_REDUCE"


MM_REDUCE_SOURCES = '''
   mm_order_vector.c
   mm_order.c
   mm_compress.c
   mm_reduce.c
   mm_suborbit.c
   mm_shorten.c
   mm_vector_v1_mod3.c
'''

MM_REDUCE_GENERATE = GENERATE_START + '''
 --dll {DLL_NAME} -v
 --source-path {SRC_DIR}/mmgroup/dev/mm_reduce
 --set p=15
 --tables mmgroup.dev.mm_op.mm_op
          mmgroup.dev.mm_reduce.order_vector_tables
          mmgroup.dev.mm_reduce.vector_v1_mod3
 --sources mm_reduce.h
 --sources  
'''.format(**DIR_DICT) + MM_REDUCE_SOURCES + '''
 --pxd   mm_reduce.pxd
 --pxi
 --pyx   mm_reduce.pyx
'''



STAGE3_SOURCES = MM_REDUCE_SOURCES
STAGE3_LIBS = STAGE2_LIBS + ['mmgroup_mm_reduce']
STAGE3_BUILD_LIBS = shared_lib_name(STAGE3_LIBS, 'build', 
                   static=STATIC_LIB)
STAGE3_BUILD_EXT = shared_lib_name(STAGE3_LIBS, 'build_ext', 
                   static=STATIC_LIB)

MM_REDUCE = SHARED_START + '''
    --name mmgroup_mm_reduce
    --sources 
''' + STAGE3_SOURCES + '''
    --libraries 
''' + " ".join(STAGE2_BUILD_LIBS)


mm_reduce_steps =  CustomBuildStep(
   'Building shared library mm_reduce',
   [sys.executable, 'generate_code.py'] + MM_REDUCE_GENERATE.split(),
   [sys.executable, 'build_shared.py'] + MM_REDUCE.split(),
)

mm_reduce_extension = Extension('mmgroup.mm_reduce',
    sources=[
            os.path.join(PXD_DIR, 'mm_reduce.pyx'),
    ],
    #libraries=['m'] # Unix-like specific
    include_dirs = [ C_DIR ],
    library_dirs = [ PACKAGE_DIR, LIB_DIR],
    libraries = STAGE3_BUILD_EXT, 
            # for openmp add 'libgomp' 
    #runtime_library_dirs = ['.'],
    extra_compile_args = EXTRA_COMPILE_ARGS, 
            # for openmp add '-fopenmp' 
    extra_link_args = EXTRA_LINK_ARGS, 
            # for openmp add '-fopenmp' 
)


if STAGE < 4:
    ext_modules += [
        mm_reduce_steps,
    ]
    if not on_readthedocs:
        ext_modules += [
            mm_reduce_extension,
        ]

####################################################################
# Copy generated header files
####################################################################

COPY_HEADERS = """
  --header-path {C_DIR}
  --header-dir {HEADERS_DIR}
  --lib-path {LIB_DIR}
  --lib-dir {HEADERS_DIR}
""".format(**DIR_DICT)


if not STATIC_LIB and not on_readthedocs:
    copy_headers_step = CustomBuildStep(
       'Copy header files',
      [sys.executable, 'shared_headers.py'] + COPY_HEADERS.split(),
      [sys.executable, 'copy_shared.py', r'${build_lib}'],
    )
    ext_modules.append(copy_headers_step)


####################################################################
# Don't build any externals when building the documentation.
####################################################################


def read(fname):
    '''Return the text in the file with name 'fname' ''' 
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

####################################################################
# The main setup program.
####################################################################

if os.name ==  'posix': 
   EXCLUDE = ['*.dll', '*.pyd', '*.*.dll', '*.*.pyd'] 
elif os.name ==  'nt': 
   EXCLUDE = ['*.so', '*.*.so'] 
else:
   EXCLUDE = [] 


setup(
    name = 'mmgroup',    
    version = VERSION,    
    license='BSD-2-Clause',
    description='Implementation of the sporadic simple monster group.',
    long_description=read('README.rst'),
    author='Martin Seysen',
    author_email='m.seysen@gmx.de',
    url='https://github.com/Martin-Seysen/mmgroup',
    packages=find_namespace_packages(
        where = 'src',
        exclude = EXCLUDE
    ),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=False,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 1 - Planning',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        #'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        #'Programming Language :: Python :: 3.6',
        #'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
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
       'Issue Tracker': 'https://github.com/Martin-Seysen/mmgroup/issues',
    },
    keywords=[
        'sporadic group', 'monster group', 'finite simple group'
    ],
    python_requires='>=3.6',
    install_requires=[
         'numpy', 'regex',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    setup_requires=[
        'numpy', 'pytest-runner', 'cython', 'regex',
        # 'sphinx',  'sphinxcontrib-bibtex',
    ],
    tests_require=[
        'pytest', 'numpy', 'regex', 'pytest-xdist'
    ],
    cmdclass={
        'build_ext': BuildExtCmd,
    },
    ext_modules = ext_modules,
    package_data = package_data,
    include_dirs=[np.get_include()],  # This gets all the required Numpy core files
)



