from __future__ import absolute_import, division, print_function

import sys
import os
import subprocess
import numpy as np
from glob import glob

import setuptools
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Distutils import build_ext




ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
SRC_DIR = os.path.join(ROOT_DIR, "src", "mmgroup")
DEV_DIR = os.path.join(SRC_DIR, "dev")
C_DIR = os.path.join(DEV_DIR, "c_files")

sys.path.append(SRC_DIR)
from dev import config
from dev.config import EXTRA_COMPILE_ARGS, PRIMES
assert sys.path.pop() == SRC_DIR



# Desription of the list 'custom_presteps'.
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
custom_presteps = [
 ##  [sys.executable, "setup_mat24.py", "build_ext", "--inplace"],
  [sys.executable, os.path.join(DEV_DIR, "mm_basics", "codegen.py")],
  [sys.executable, os.path.join(DEV_DIR, "mm_basics", 
                                      "make_mm_basics_dll.py")],
  [sys.executable, os.path.join(DEV_DIR, "mm_op", "codegen.py")],
]



build_ext_delete = [
     os.path.join(DEV_DIR, "mm_basics", "mm*.c"),
     os.path.join(DEV_DIR, "mm_basics", "mm*.o"),
     os.path.join(DEV_DIR, "mm_basics", "mm*.def"),
     os.path.join(DEV_DIR, "mm_op", "mm*.c"),
     os.path.join(DEV_DIR, "mm_op", "mm*.o"),
     os.path.join(DEV_DIR, "mm_op", "mm*.def"),
]




class GenCodeCommand(setuptools.Command):
    """Custom 'gen_code' command to run a fixed sequence of programs.

    These programs will be run automatically as a customized pre-
    processing step before executing the 'build_ext' command; they 
    will generate the C code to be processed by the 'build_ext' command.

    You may also execute the 'gen_code' setup command for just running
    all these preprocessing steps without the subsequent 'build_ext'
    step.
    
    The programs to be run are given by the global list with name
    'custom_presteps' in this file. We make no attempt to pass the 
    list of these programs to the setup() function as an argument.

    The basic idea for this class is taken from:
    https://jichu4n.com/posts/how-to-add-custom-build-steps-and-commands-to-setuppy/
    """

    description = "run code generating steps for 'build_ext' command"
    user_options = [
      # The format is (long option, short option, description).
    ]

    def run(self):
        for prestep_command in custom_presteps:
            if subprocess.call(prestep_command) != 0:
                sys.exit(-1)
        print("Code generation done")

    def initialize_options(self):
        """Set default values for options."""
        # Each user option must be listed here with their default value.
        pass

    def finalize_options(self):
        """Post-process options."""
        pass

class  BuildExtCommand(build_ext):
    """Customize the 'build_ext' command.

    We will always automatically run the 'gen_code' custom command 
    before the 'build_ext' command. The purprose of the 'gen_code' 
    custom command is to generate the C code used by the 'build_ext' 
    command for building an extension.

    Class GenCodeCommand implements the 'gen_code' custom command.

    Befor and after building the extension, all files listed in the 
    global list 'build_ext_delete' will bedeleted. Wildcards in that
    list are expanded using glob.golob()
    """

    def run(self):
        self.custom_delete_intermediate_files()
        self.run_command('gen_code')
        super(BuildExtCommand, self).run()
        self.custom_delete_intermediate_files()
 
    def custom_delete_intermediate_files(self):
        for file_pattern in build_ext_delete:
            for file in glob(file_pattern):
                try:
                    os.remove(file)
                except:
                    pass

ext_modules=[
    Extension("mmgroup.mm",
        sources=[
            os.path.join(DEV_DIR, "mm_basics", "mm_basics.pyx"),
        ],
        #libraries=["m"] # Unix-like specific
        include_dirs = [ C_DIR ],
        library_dirs = [C_DIR ],
        libraries = ["libmm_basics", "libmat24"], 
        #runtime_library_dirs = ["."],
        extra_compile_args = EXTRA_COMPILE_ARGS, 
    )
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
    sources = [os.path.join(DEV_DIR, "mm_op", PYX_SOURCE_P.format(P = p))]
    for f in C_SOURCES_P:
         sources.append(os.path.join(C_DIR, f.format(P = p) + ".c"))
    return sources

    
for p in PRIMES:
    ext_modules.append(
        Extension("mmgroup.mm%d" % p,
            sources = list_source_files(p),
            #libraries=["m"] # Unix-like specific
            include_dirs = [ C_DIR ], 
             library_dirs = [C_DIR ],
            libraries = ["libmm_basics", "libmat24"], 
            #runtime_library_dirs = ["."],
            extra_compile_args = EXTRA_COMPILE_ARGS, 
        )
    )

#TODO: organize the package following
#https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure


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
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    python_requires='>=3.4',
    install_requires=[
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    setup_requires=[
        'pytest-runner',
    ],
    cmdclass={
        'gen_code': GenCodeCommand,
        'build_ext': BuildExtCommand,
    },
    ext_modules = ext_modules,
    include_dirs=[np.get_include()],  # This gets all the required Numpy core files


)



