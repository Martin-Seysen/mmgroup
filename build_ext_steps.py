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


A new paradigm for building python packages
...........................................


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

  * A keyword argument ``extra_compile_args`` and ``extra_link_args`` 
    for  an instance of class ``Extension`` may be a dictionary

       'compiler' : <List if arguments>

    instead of a list of arguments. Here 'compiler' is a string
    describing a compiler. For a list of compilers, run

      ``python setup.py build_ext --help-compiler`` .


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

Here the first entry of each list is either a string or a callable
python function. If the first entry is a string then a subprocess
with that name is called. Otherwise the corresponding python function
is executed. Subsequent entries in the list are arguments given to
the subprocess or to the function.

Such a subprocess may be e.g. a step that generates C code to be
used for building a subsequent python extension.    

Its recommended to use the string ``sys.executable`` (provided
by the ``sys`` package) instead of the string ``'python'`` for 
starting a python subprocess.


Building shared libraries
.........................

Module ``build_ext_steps`` provides another class ``SharedExtension``
for building  a shared library (or a DLL in Windows).

In the list ``ext_modules`` of extensions, instances of class 
``SharedExtension`` may be mixed with instances of other  classes.

Arguments for the constructor of class ``SharedExtension`` are the
same as for class ``Extension``. Especially, the following keyword 
arguments are recognized:


  * ``name``, ``sources``, ``include_dirs``, ``libraries``,
    ``define_macros``, ``undef_macros``

Here ``name`` is a Python dotted  name without any extension as in
class ``Extension``. The command uses the C compiler to build
the shared library and stores it at the location given by ``name``
in the same way as a python extension is built. The appropriate
extension (e.g. '.so' for unix and '.dll' for Windows) is 
automatically appended to the file name of the shared library.

The user should provide an additional keyword argument ``implib_dir`` 
specifying a directory where to store the import library for a Windows 
DLL. In Windows, the import library for ``foo.dll`` has the name 
``libfoo.lib``. If a program uses a Windows DLL then it should be 
linked to that import library. In unix operating systems there is no 
concept similar to an import library.

Caution!!

The current class ``SharedExtension`` supports Windows DLLs containing
C programs (no C++) compiled with the ``mingw32`` compiler only.


Using a shared library in an extension
......................................

The reason for building a shared library is that several python
extensions may use the same shared library.

The way how a shared library (or a Windows DLL) is linked to a program
using that library depends on the operating system.

The user may have to perform some os-specific steps for making the 
library available for python extension. Therefore he may read
the variable ``os.name`` (which has value ``'nt'`` for Windows and
``'posix'`` for unix) and use class ``CustomBuildStep`` for performing
the appropriate steps.

For Windows one has to add the directory of the import library
(discussed in the previous section) to the search path for the library
by specifying that directory in the ``'library_dirs'`` keyword 
argument for class Extension. Also, one has to specify the name
of the import library in the  ``'libraries'`` keyword argument for 
class Extension. This is sufficient if the python extension and
the DLLs used by that extension are in the same directory.

More involved cases of shared libraries and the procedures required
for other operating systems are out of the scope of this documentation.

Using an extension in a subsequent build step
.............................................

Once a python extension has been built, it can also be used in
a subsequent step of the build process, e.g. for calculating large
arrays of constants for C programs.

This approach works well on a Windows system, but it might not work
on other operating systems. Here it is a good idea to write a
pure-python substitute for any C extension to be used in a subsequent
build step. This may slow down the build process considerably. But it 
is better to have a slow build process than no build process at all.


Technical remarks about shared libraries
........................................

The functionality for building a Windows DLL with the ``migw32`` 
compiler in coded in function ``make_dll_win32_gcc``. One could have 
used the functionality of class  ``distutils.ccompiler`` instead, 
but the author has decided not to dive any deeper into the source 
code of the ``distutils`` package. It may be worth using class  
``distutils.ccompiler`` for porting the functionality of classes 
``BuildExtCmd`` and  ``SharedExtension`` to other compilers and 
operating systems.

The ``mingw32`` compiler is not the standard compiler for python on
Windows. C++ files should be compiled with the standard compiler
(which is ``msvc`` for Windows) in order to avoid trouble with 
name mangling.
"""

import sys
import os
import subprocess

from glob import glob

import setuptools
from setuptools.extension import Extension as _Extension


from Cython.Distutils import build_ext as _build_ext
from distutils.file_util import copy_file
from distutils.errors import *


_compiler = None


# Class ``Extension`` is equal to class ``Cython.Distutils.build_ext``.
# This may change in future versions.
Extension = _Extension

class CustomBuildStep(_Extension):
    """Model a custom build step as described in the header of this module 
    """
    def __init__(self, name = "Customer_step", *args):
        # Set attribute name (for displaying the corresponding build step)
        self.name = name
        if not isinstance(name, str):
            raise AssertionError("'name' must be a string")
        # Set the list of functions
        self.function_list = list(args) or []
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
          


class SharedExtension(_Extension):
    """Model a shared library as described in the header of this module 
    """
    def __init__(self,  *args, **kwds):
        self.implib_dir = None
        # If keyword argument "implib_dir" is given, set self.implib_dir
        # to the corresponding value
        if "implib_dir" in kwds:
            self.implib_dir = kwds["implib_dir"]
            del kwds["implib_dir"]
        # Treat all other arguments as in the base class 'Extension'.
        super(SharedExtension, self).__init__(*args, **kwds)
        self.is_non_standard = True






class  BuildExtCmd(_build_ext):
    """Substitute for class setuptools.extension import Extension

    The ``run`` function of this class treats an instance  of class
    Extension in the usual way, i.e. it builds a python extension. 

    Instances of class ``CustomBuildStep`` or ``SharedExtension``
    are treated as in the description of this module.
    """
    def run(self):
        """Implementation of the run() function.

        See class ditutils.cmd.Command for background.
        """
        # At this point, 'self.compiler' is None or a string,  but this 
        # will change while running the corresponding distutils command
        self.compiler_name = self.compiler
        if self.compiler_name is None:
            # Then try to get the default compiler name
            from distutils.ccompiler import get_default_compiler
            self.compiler_name = get_default_compiler()

        # Now process all extensions in the given order
        for ext in self.extensions[:]:
            if isinstance(ext, CustomBuildStep):
                # Then we execute a sequence of subprocesses or a functions 
                print("\nExecuting custom build step '%s'" %  ext.name)
                for args in ext.function_list:
                    if isinstance(args[0], str):
                        # Then execute a subprocess
                        print(" ".join(map(str, args)))
                        if subprocess.call(list(map(str, args))) != 0:
                            err = "Subprocess %s failed"
                            raise ValueError(err %  args[0])              
                    else:
                        # Then execute a python function
                        args[0](*(args[1:])) 
            elif isinstance(ext, SharedExtension):
                # Then we build a shared librry using method
                # self.build_shared_extension()

                # Evaluate extra_compile_args and extra_link_args to a string
                ext.extra_compile_args = self.eval_extra_args(
                    ext.extra_compile_args)
                ext.extra_link_args = self.eval_extra_args(
                    ext.extra_link_args)
                print("\nBuilding shared library '%s'" %  ext.name)
                self.build_shared_extension(ext)
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
                self.inplace = True     # Always bild inplace
                force = self.force      
                self.force = True       # Always force building
                # Run the corresponding mathod of the base class
                super(BuildExtCmd, self).run() 
                # Restore old attributes
                self.extensions = extensions
                self.inplace = inplace
                self.force = True
                self.compiler = compiler


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
 
    def build_shared_extension(self, ext):
        """Build a shared library

        The shared libarary is described by the instance ``ext``
        of class SharedExtension.
        """ 
        compiler = self.compiler_name
        if os.name == "posix":
            if compiler == 'unix':
                make_so_posix_gcc(self, ext)
            else:
                raise DistutilsPlatformError(
                  "I don't know how to build a posix shared library "
                  "with the '%s' compiler" % compiler)
        elif os.name == "nt":
            if compiler == 'mingw32':
                make_dll_nt_mingw32(self, ext)
            else:
                raise DistutilsPlatformError(
                  "I don't know how to build a Windows DLL "
                  "with the '%s' compiler" % compiler)
        else:
            raise DistutilsPlatformError(
                  "I don't know how to build a shared library "
                  "on platform '%s'" % os.name)

    @classmethod
    def get_shared_ext_filename(cls, ext_name):
        """Convert a python module to the name of a shared library

        Here ``ext_name`` is the   Python dotted  name without any 
        extension as in class ``Extension``.

        The function returns the full path name of the corresponding
        shared library.
        """
        if os.name == "posix":
            suffix = ".so"
        elif os.name == "nt":
            suffix = ".dll"
        else:
            raise DistutilsPlatformError(
                  "I don't know the suffix for a shared library "
                  "on platform '%s'" % os.name)
        ext_path = ext_name.split('.')
        return os.path.join(*ext_path) + suffix

    def get_outputs(self):
        """This extends the corresponding method of class Extension

        It returns a list of output file names. Apart from the names 
        of the created python extensions, that list also contains 
        the names of the created shared libraries.
        """
        extensions = self.extensions
        self.extension = [ x for  x in extensions
            if not getattr(x, "is_non_standard", None)
        ]
        self.extensions = extensions
        standard_outputs = super(BuildExtCmd, self).get_outputs()
        return standard_outputs + self.get_non_standard_outputs()

    def get_non_standard_outputs(self):
        """Auxilary method for self. get_outputs()

        Returns list of the shared libraries (to be) built.
        """
        outputs = []
        for ext in self.extensions[:]:
            if isinstance(ext, SharedExtension):
                outputs.append(self.get_shared_ext_filename(ext.name))
        return outputs

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







def make_dll_nt_mingw32(cmd, ext):
    """Create a Windows DLL with the mingw compiler"""
    compile_args = [ "-c", "-Wall", "-DMS_WIN64"] + ext.extra_compile_args 
    for ipath in ext.include_dirs:
        compile_args.append("-I" + os.path.realpath(ipath))
    link_args = ["-shared",  "-Wall", "-DMS_WIN64"] + ext.extra_link_args 
    objects = []

    # add define macros
    for name, value in ext.define_macros:
        if value is None:
            compile_args.append("-D" + name)
        else:
            compile_args.append("-D" + name + "=" + value)

    # add undef macros
    for name in ext.undef_macros:
        compile_args.append("-U" + name)

    # Compile sources and add objects to list 'objects'
    for source in ext.sources:
        args = compile_args[:]
        args.append(os.path.realpath(source))
        args.append("-o")
        obj = os.path.realpath(os.path.splitext(source)[0] + ".o")
        args.append(obj)
        objects.append(obj)
        print("gcc " + " ".join(args))
        subprocess.call("gcc " + " ".join(args)) 

    # Link
    largs = link_args[:] + objects 
    for inc_dir in ext.library_dirs:
        # same search path for include files an libraries
        largs.append("-L" + inc_dir)
    for library in ext.libraries: 
        largs.append("-l" + library )
    dll_name = cmd.get_shared_ext_filename(ext.name)  
    dll_path =   os.path.join(cmd.get_package_dir(), dll_name)
    largs +=  ["-o",  dll_path ]
    if not ext.implib_dir:
        raise DistutilsSetupError(
                  "Class SharedExtension requires keyword "
                  "argument 'implib_dir' for building Windows DLL" )
    pure_dll_name = ext.name.split(".")[-1]
    def_path =   os.path.join(ext.implib_dir, "%s.def" % pure_dll_name)
    largs += [ "-Wl,--output-def," + def_path ]
    implib_path = os.path.join(ext.implib_dir, "lib%s.lib" % pure_dll_name)
    largs += [ "-Wl,--out-implib," + implib_path ]
    print("gcc " + " ".join(largs))
    subprocess.call("gcc " + " ".join(largs)) 

    print("removing object files...")
    for obj in objects:
        os.remove(obj)
    print(dll_path + "\nhas been generated.\n")

    #Todo: copy dll to build path if requested
    


def make_so_posix_gcc(cmd, ext):
    """Create a posix shared library with gcc"""
    compile_args = [ "-c", "-Wall"] + ext.extra_compile_args 
    for ipath in ext.include_dirs:
        compile_args.append("-I" + os.path.realpath(ipath))
    link_args = ["-shared",  "-Wall"] + ext.extra_link_args 
    objects = []

    # add define macros
    for name, value in ext.define_macros:
        if value is None:
            compile_args.append("-D" + name)
        else:
            compile_args.append("-D" + name + "=" + value)

    # add undef macros
    for name in ext.undef_macros:
        compile_args.append("-U" + name)

    # Compile sources and add objects to list 'objects'
    for source in ext.sources:
        args = compile_args[:] + ["-fPIC"]
        args.append(os.path.realpath(source))
        args.append("-o")
        obj = os.path.realpath(os.path.splitext(source)[0] + ".o")
        args.append(obj)
        objects.append(obj)
        print("gcc " + " ".join(args))
        subprocess.call(["cc"] + args) 

    # Link
    largs = link_args[:] + objects 
    for inc_dir in ext.library_dirs:
        # same search path for include files an libraries
        largs.append("-L" + inc_dir)
    for library in ext.libraries: 
        largs.append("-l" + library )
    dll_name = cmd.get_shared_ext_filename(ext.name).split("/")
    dll_name[-1] =  "lib" + dll_name[-1]  
    dll_path =   os.path.join(cmd.get_package_dir(), *dll_name)
    largs +=  ["-o",  dll_path ]
    subprocess.call(["cc"] + largs) 

    print("removing object files...")
    for obj in objects:
        os.remove(obj)
    print(dll_path + "\nhas been generated.\n")

    #Todo: copy dll to build path if requested
    







