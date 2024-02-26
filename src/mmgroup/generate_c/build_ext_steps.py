r"""Module ``build_ext_steps`` provides a customized version of 
the 'build_ext' command for setup.py.

It should be placed in the same directory as the module ``setup.py``.

Distributing python packages
............................

The standard toolkit for distributing python packages is the 
``setuptools`` package. Here the user types::

   python setup.py build_ext

at the console for building the extensions to the python package, which 
are typically written in a language like C or C++ for the sake of speed. 
We may use e.g. the ``Cython`` package to write python wrappers for the 
functions written in C or C++. The ``setuptools`` package supports the
integration of ``Cython`` extensions.

The ``setup.py`` script describes a list of extensions, where each
extension is an instance of class ``Extension`` which is provided  
by  ``setuptools``. Then it calls the ``setup`` function which builds
all these extensions:

.. code-block:: python

  from setuptools import setup
  from setuptools.extension import Extension   

  ext_modules = [
      Extension(
          ...  # description of first extension
      ),
      Extension(
          ...  # description of second extension
      ),
      ...
  ]

  setup(
      ..., # Standard arguments for the setup function
      ext_modules = ext_modules,   # List of extensions to be built
  )

We assume that the reader is familiar with the standard python setup 
process. For background, we refer to


https://setuptools.readthedocs.io/en/latest/


Added functionality for building python packages
...................................................


This ``build_ext_steps`` module supports a new paradigm for building a 
python extension:

  * A python program ``make_stage1.py`` creates a C program 
    ``stage1.c``.

  * We create a python extension ``stage1.so`` (or ``stage1.pyd`` 
    in Windows) that makes the functionality of ``stage1.c``
    available in python.

  * A python program ``make_stage2.py`` creates a C program 
    ``stage2.c``. Here ``make_stage2.py`` may import ``stage1.so`` 
    (or ``stage1.pyd``).

  * We create a python extension  that makes the functionality of 
    ``stage2.c``  available in python.

  * etc.

This paradigm is not supported by the ``setuptools`` package.

Using class ``BuildExtCmd`` for the new paradigm
................................................

For using the new building paradigm we have to replace the standard 
class ``build_ext`` by the class ``build_ext_steps.BuildExtCmd``. 

.. code-block:: python

  from setuptools import setup
  from build_ext_steps import Extension   
  from build_ext_steps import BuildExtCmd

  ext_modules = [
      # description of extension as above
  ]

  setup(
      ..., # Standard arguments for the setup function
      ext_modules = ext_modules,   # List of extensions to be built
      cmdclass={
         'build_ext': BuildExtCmd, # replace class for build_ext
      },
  )

This change has a few consequences:

  * It is guaranteed that the extension are build in the given order

  * Extensions are always build in place (option ``build_ext --inplace``)
    (The current version does not support building the extension in
    a special build directory.)

  * The building of all extensions is now forced 
    (option ``build_ext --f``), regardless of any time stamps.


Apart from these changes, an extension is created in the same way 
as with ``setuptools``.

For a documentation of the ``Extension`` class in the ``setuptools``
package, see

https://docs.python.org/3/distutils/apiref.html?highlight=extension#distutils.core.Extension


Inserting user-defined functions into the build process
.......................................................


Module ``build_ext_steps`` provides a class ``CustomBuildStep``
for adding user-defined functions to the build process.

In the list ``ext_modules`` of extensions, instances of class 
``CustomBuildStep`` may be mixed with instances of class ``Extension``,
Class ``CustomBuildStep`` models an arbitrary  sequence of functions to 
be executed.

The constructor for that class takes a string 'name' describing the
action of these functions, followed by an arbitrary number of lists, 
where each list describes a function to be executed.

Here the first entry should be a string. Then a subprocess with that
name is called. Subsequent entries in the list are arguments passed to
the subprocess.

If such an argument contains the string '${build_lib}' then that string
is replaced by the name of the root directory of the path where the
extionsion is to be built. If the 'build_ext' command of  ``setup.py``
is invoked with the '--inplace' option then the string '${build_lib}'
is replaced by 'null'.

Such a subprocess may be e.g. a step that generates C code to be
used for building a subsequent python extension.    

Its recommended to use the string ``sys.executable`` (provided
by the ``sys`` package) instead of the string ``'python'`` for 
starting a python subprocess.

The following functionality is deprecated:

If the first entry of a list as descibed above is a python function
then that function is called.  Subsequent entries in the list are
passed as arguments to the function.


Using a shared library in an extension
......................................

The user may have to perform some os-specific steps for making the 
library available for python extension. Details are out of the scope 
of this documentation.

Using an extension in a subsequent build step
.............................................

Once a python extension has been built, it can also be used in
a subsequent step of the build process, e.g. for calculating large
arrays of constants for C programs. Details are out of the scope 
of this documentation..

"""


import sys
import os
import re
import shutil
import subprocess
from subprocess import CalledProcessError
from glob import glob

import setuptools
from setuptools.extension import Extension as _Extension
from Cython.Distutils import build_ext as _build_ext

from mmgroup.generate_c.build_shared import shared_lib_name



_compiler = None


# Class ``Extension`` is equal to class ``Cython.Distutils.build_ext``.
# This may change in future versions.
Extension = _Extension

class CustomBuildStep(_Extension):
    """Model a custom build step as described in the header of this module 
    """
    def __init__(self, name = "Customer_step", *args, **kwds):
        # Set attribute name (for displaying the corresponding build step)
        self.name = name
        if not isinstance(name, str):
            raise AssertionError("'name' must be a string")
        # Set the list of functions
        self.function_list = list(args) or []
        # Set the dict of keywords
        self.keyword_dict = kwds
        # Mark this as a nonstandard instance of class Extension
        self.is_non_standard = True
        # Set all other attributes of class Extension to default values
        self.sources = []
        self.include_dirs = []
        self.define_macros =  []
        self.undef_macros = []
        self.library_dirs =  []
        self.libraries =  []
        self.runtime_library_dirs =  []
        self.extra_objects = []
        self.extra_compile_args = []
        self.extra_link_args =  []
        self.export_symbols =  []
        self.swig_opts =  []
        self.depends =  []
        self.language = None
        self.optional = None
          



class BuildExtCmdObj:
    """Placeholder for argument in class CustomBuildStep

    When class BuildExtCmd is used as a substitute for
    class setuptools.command.build_ext then there is subclass
    CustomBuildStep of class Extension that may simply execute a 
    python function. This object ``BuildExtCmdObj`` may be passed
    to such a python function as a placeholder for the
    object of class BuildExtCmd that has invoked that python
    function.
    """
    pass



def _get_default_compiler():
    """Simplified version of distutils.ccompiler.get_default_compiler()

    Note that distutils will be removed in python 3.12. So this function
    does what we need from distutils.ccompiler.get_default_compiler().
    """
    if re.match('nt', os.name):
        return 'msvc'
    return 'unix'



class  BuildExtCmd(_build_ext):
    """Substitute for class setuptools.command.build_ext 

    The ``run`` function of this class treats an instance  of class
    Extension in the usual way, i.e. it builds a python extension. 

    Instances of class ``CustomBuildStep`` are treated as in the 
    description of this module.
    """
    WARN_PY = "Executing a python function with class BuildExtCmd is deprecated!"
    user_options = _build_ext.user_options + [
        ("nprocesses=", None, 
        "Number of prcesses used by class ParallelSteps"
        )
    ]

    def initialize_options(self):
        super(BuildExtCmd, self).initialize_options()
        self.nprocesses = 1

    def finalize_options(self):
        super(BuildExtCmd, self).finalize_options()
        self.nprocesses = max(1, int(self.nprocesses))
   

    def run(self):
        """Implementation of the run() function.

        See class distutils.cmd.Command for background.
        """
        # At this point, 'self.compiler' is None or a string,  but this 
        # will change while running the corresponding distutils command
        self.compiler_name = self.compiler
        if self.compiler_name is None:
            # Then try to get the default compiler name
            # The following line will no longer work in python 3.12. 
            # As a future remedy, we might try:
            # from setuptools._distutils.ccompiler import get_default_compiler
            try:
                from ddistutils.ccompiler import get_default_compiler
                self.compiler_name = get_default_compiler()
            except (ModuleNotFoundError, ImportError, NameError):
                self.compiler_name = _get_default_compiler()

        # Now process all extensions in the given order
        for ext in self.extensions[:]:
            if isinstance(ext, CustomBuildStep):
                # Then we execute a sequence of subprocesses or functions 
                print("\nExecuting custom build step '%s'" %  ext.name)
                sys.stdout.flush()
                sys.stderr.flush()
                for args in ext.function_list:
                    if isinstance(args[0], str):
                        # Then execute a subprocess
                        args1 = self._convert_subprocesss_args(args)
                        print(" ".join(map(str, args)))
                        try:
                            subprocess.check_call(args1)
                        except CalledProcessError:
                            err = "\nSubprocess '%s' has failed!"
                            print(err % str(args[0]))
                            raise             
                    else:
                        # Then execute a python function
                        warning.warn(self.WARN_PY, UserWarning)
                        f, f_args = args[0], args[1:]
                        for i, arg in enumerate(f_args):
                            ## Substitute BuildExtCmdObj by self
                            if arg == BuildExtCmdObj:
                                f_args[i] = self
                        f(*f_args) 
                sys.stdout.flush()
                sys.stderr.flush()
            elif isinstance(ext, str):
                # A string is taken as a comment
                pass 
            else:
                # A standard Extension is build in the same way as in the
                # base class. We rely on the fact that setuptools always 
                # builds the python extension in the build directory,
                # even if the 'inplace' option has not been set.

                # Evaluate extra_compile_args and extra_link_args to a string
                ext.extra_compile_args = self.eval_extra_args(
                    ext.extra_compile_args)
                ext.extra_link_args = self.eval_extra_args(
                    ext.extra_link_args)
                # Save the old value of self.compiler because the self.run()
                # method appears to change it so that it self.run() method
                # cannot be applied twice.
                compiler = self.compiler
                extensions = self.extensions
                self.extensions = [ext] # Build one extension only
                inplace = self.inplace
                self.inplace = True     # Always build inplace
                force = self.force      
                self.force = True       # Always force building
                # Run the corresponding mathod of the base class
                super(BuildExtCmd, self).run() 
                # Restore old attributes
                self.extensions = extensions
                self.inplace = inplace
                self.force = True
                self.compiler = compiler


    def _convert_subprocesss_args(self, args):
        def subs(s):
            s_subst = 'null' if self.inplace else str(self.build_lib)
            f_subst = lambda x: str(s_subst)
            if r'${build_lib}' in s:
                print(s_subst)
                return re.sub(r'\${build_lib}', f_subst, s)
            return s
        return [subs(arg) for arg in args]

    def eval_extra_args(self, extra_args):
        """Evaluate the dictionary ``extra_args`` to a string.

        Here ``extra_args`` is one of the keyword arguments 
        ``extra_compile_args`` or ``extra_link_args`` passed
        to the constructor of class Extension. The  dictionary
        has entries:

            'compiler' : <List if arguments>

        If ``self.compiler_name`` matches a key 'compiler'  then  
        the list ``extra_args[self.compiler_name]`` is returned. 
        If ``self.compiler_name`` is a string then that string is
        returned.
        
        For a list of compilers, run

            ``python setup.py build_ext --help-compiler`` .

        """
        if isinstance(extra_args, dict):
            return extra_args[self.compiler_name]
        return extra_args
 



    def get_build_directory(self):
        """Returns name of directory for the output to be built
        """ 
        print("build_dir is", self.build_lib)
        return self.build_lib

    def get_package_dir(self):
        """Return the package directory of the build command"""
        build_py = self.get_finalized_command('build_py')
        package_dir = os.path.abspath(build_py.get_package_dir('.'))
        #print(package_dir)
        return package_dir

###########################################################################################



