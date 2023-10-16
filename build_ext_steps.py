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

Keyword arguments in the constructor are passed to function
``subprocess.check_call`` when executing a subprocess and 
ignored when executing a python function.

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

The user should provide an additional keyword argument ``lib_dir`` 
specifying a directory where to store any static libraries, object
files or import library for a Windows DLL. In Windows, the import 
library for ``foo.dll`` has the name ``libfoo.lib``. If a program 
uses a Windows DLL then it should be linked to that import library. 
In unix operating systems there is no concept similar to an 
import library.

Caution!!

The current class ``SharedExtension`` supports Windows DLLs containing
C programs (no C++) compiled with the ``mingw32`` or the ``msvc``
compiler only.


Using a shared library in an extension
......................................

The reason for building a shared library is that several python
extensions may use the same shared library.

The way how a shared library (or a Windows DLL) is linked to a program
using that library depends on the operating system.

After creation, an object of class  ``SharedExtension`` contains an
attribute ``lib_name``. Attribute ``lib_name`` contains name of the 
library (or of the import library in Windows) that has to be linked
to a program using the shared library.

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
import re
import subprocess
from subprocess import CalledProcessError

from glob import glob

import setuptools
from setuptools.extension import Extension as _Extension



def shared_lib_name(name, mode, static=False, os_name=None, pymod=False, dir=True):
    """Adjust name of a (shared) library

    Here parameter ``name`` is a name of a (shared) library to be
    adjusted to various situations in the build process. This name
    may be prefixed by some path information, but the name itself
    should not contain any suffixes, e.g. '.a', '.so', '.lib', '.dll'.
    The part of parameter ``name`` indicating the library name
    should not contain any prefix strings, e.g. 'lib'.

    The function adjusts the name of the library as configured by the
    other parameters of the function and returns adjusted name as a
    string. If ``name`` is a list of strings then the function returns
    the list of adjusted strings.

    Parameter ``mode`` describes the part of the build process
    where the adjusted name is to be used:

    mode = 'load':      The file name name of the (shared) library
                        including all prefexes and extensions.

    mode = 'build':     The name of the library to be linked to
                        a program using that library. In Unix-like
                        systems this is the (shared) library itsalf.
                        In case of a Windows DLL this is an import
                        library generated together with the DLL.

    mode = 'build_ext'  Similar to mode 'build'; to be used in
                        ``setuptools`` instead of mode 'build'.

    Parameter ``static`` should be True in case of a static library
    and False (default) in case of a shared library.

    Parameter ``os_name`` defaults to the string returned by function
    ``os.name`` in the ``os`` package. It may be set to model the
    behaviour of a different os.

    If parameter ``pymod`` is True then ``name`` is parsed as a python
    module separated by '.' characters. Default is False.

    If parameter ``dir`` is False then any information in parameter
    ``name`` refering to a directory is removed. Default is True.
    """
    if not os_name:
        os_name = os.name
    if isinstance(name, list):
        return [shared_lib_name(x, mode, static, os_name, pymod, dir)
                 for x in name]
    if pymod:
        name = re.sub('\.', '/', name)
    path = os.path.split(name)
    if len(path) == 0:
        raise ValueError("No library name specified")
    path, name = list(path[:-1]), path[-1]
    if not dir:
        path = []
    platform_found = new_name = False
    if os.name == "posix":
        platform_found = True
        if mode == "build_ext":
            new_name = name
        elif mode in ["build", "load"]:
            suffix = ".a" if static else ".so"
            new_name =  "lib" + name + suffix
    elif os.name == "nt":
        platform_found = True
        if static:
            if mode == "build_ext":
                new_name = name
            elif mode in ["build", "load"]:
                new_name = name + ".lib"
        else:
           if mode == "build_ext":
               new_name = "lib" + name
           if mode == "build":
               new_name =  "lib" + name + ".lib"
           if mode == "load":
               new_name =  "lib" + name + ".dll"

    if not platform_found:
        err = "Don't know the suffix for a %slibrary on platform %s"
        sh = '' if static else 'shared '
        raise ValueError(err % (sh, os_name))
    if not new_name:
        err = "Illegal mode %d in function shared_lib_name"
        raise ValueError(err % mode)
    new_name = os.path.normpath(os.path.join(*(path + [new_name])))
    #print('%s' % new_name)
    return new_name


from Cython.Distutils import build_ext as _build_ext

from parallel_processes import SimpleProcessWorkers

class DistutilsError (Exception):
    """The root of all Distutils exceptions.

    Specifically, I do not know how and when exactly 'distutils'
    will eventually become deprecated.

    Since the autor is not willing to dig any deeper into the details
    of setuptools/distutils, I simply copy the relevant exceptions
    here. 
    """
    pass

class DistutilsPlatformError (DistutilsError):
    """We don't know how to do something on the current platform (but
    we do know how to do it on some platform) -- eg. trying to compile
    C files on a platform not supported by a CCompiler subclass."""
    pass


class DistutilsSetupError (DistutilsError):
    """For errors that can be definitely blamed on the setup script,
    such as invalid keyword arguments to 'setup()'."""
    pass



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
          


class SharedExtension(_Extension):
    """Model a shared library as described in the header of this module 
    """
    def __init__(self,  *args, **kwds):
        self.lib_dir = None
        self.static_lib = False
        # If keyword argument "lib_dir" is given, set self.lib_dir
        # to the corresponding value
        if "lib_dir" in kwds:
            self.lib_dir = kwds["lib_dir"]
            del kwds["lib_dir"]
        if "static_lib" in kwds:
            self.static_lib = bool(kwds["static_lib"])
            del kwds["static_lib"]
        # Treat all other arguments as in the base class 'Extension'.
        super(SharedExtension, self).__init__(*args, **kwds)
        self.is_non_standard = True
        self.lib_name = shared_lib_name(self.name, 'build_ext',
               static=self.static_lib,  pymod=True, dir=False)

    def get_lib_dir(self):    
        if not self.lib_dir:
            raise DistutilsSetupError(
                  "Class SharedExtension requires keyword "
                  "argument 'lib_dir' for building libary" )
        return self.lib_dir


    def linker_library_options(self, compiler):
        if compiler in ['mingw32', 'unix']:
            inc, lib = "-L", "-l"
        elif compiler in ['msvc']:
            inc, lib = "/LIBPATH:", "/DEFAULTLIB:"
        else:
            raise DistutilsSetupError(
                  "Don't know how to deal with %s compiler" % compiler )
        largs = []
        for inc_dir in self.library_dirs:
            # same search path for include files an libraries
            largs.append(inc + inc_dir)
        for library in self.libraries: 
            largs.append(lib + library )
        return largs





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

    Instances of class ``CustomBuildStep`` or ``SharedExtension``
    are treated as in the description of this module.
    """
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
                for args in ext.function_list:
                    if isinstance(args[0], str):
                        # Then execute a subprocess
                        print(" ".join(map(str, args)))
                        try:
                            subprocess.check_call(list(map(str, args)), 
                                **ext.keyword_dict)
                        except CalledProcessError:
                            err = "\nSubprocess '%s' has failed!"
                            print(err % str(args[0]))
                            raise             
                    else:
                        # Then execute a python function
                        f, f_args = args[0], args[1:]
                        for i, arg in enumerate(f_args):
                            ## Substitute BuildExtCmdObj by self
                            if arg == BuildExtCmdObj:
                                f_args[i] = self
                        f(*f_args) 
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
        os.makedirs(ext.get_lib_dir(), exist_ok=True)
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
            elif compiler == 'msvc':
                make_dll_nt_msvc(self, ext)
            else:
                raise DistutilsPlatformError(
                  "I don't know how to build a Windows DLL "
                  "with the '%s' compiler" % compiler)
        else:
            raise DistutilsPlatformError(
                  "I don't know how to build a shared library "
                  "on platform '%s'" % os.name)


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
            if isinstance(ext, SharedExtension) and not ext.static_lib:
                #outputs.append(self.get_shared_ext_filename(ext.name))
                name = shared_lib_name(self.name, 'load', pymod=True)
                outputs.append(name)

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





def make_dll_nt_msvc(cmd, ext):
    """Create a Windows DLL with the mingw compiler"""
    compile_args = [ "/c",  "/W4", "/DMS_WIN64"] + ext.extra_compile_args 
    for ipath in ext.include_dirs:
        compile_args.append("/I " + os.path.realpath(ipath))
    objects = []
    lib_dir = ext.get_lib_dir()

    # add define macros
    for name, value in ext.define_macros:
        if value is None:
            compile_args.append("-D" + name)
        else:
            compile_args.append("-D" + name + "#" + value)

    # add undef macros
    for name in ext.undef_macros:
        compile_args.append("-U" + name)

    # Compile sources and add objects to list 'objects'
    arglist = []
    for source in ext.sources:
        args = ["cl"] + compile_args[:]
        args.append(os.path.realpath(source))
        objname = os.path.splitext(os.path.split(source)[-1])[0]
        obj = os.path.realpath(os.path.join(lib_dir, objname + ".obj"))
        args.append('/Fo%s' % obj)
        objects.append(obj)
        arglist.append(args)
    workers = SimpleProcessWorkers(cmd.nprocesses)
    workers.run(arglist)

    # Link
    if ext.static_lib:
        lib_name = shared_lib_name(ext.name, 'load',
           static=True, pymod=True, dir=False)
        dll_path =   os.path.join(lib_dir, lib_name)
        lcmd =  ["lib"] + objects + ["/OUT:" + dll_path ]
    else:
        lcmd = ["link", "/DLL"] + ext.extra_link_args + objects
        lcmd += ext.linker_library_options('msvc')
        dll_name = shared_lib_name(ext.name, 'load', pymod=True)
        dll_path =   os.path.join(cmd.get_package_dir(), dll_name)
        lcmd +=  ["/OUT:" + dll_path ]
        implib_name = shared_lib_name(ext.name, 'build', pymod=True, dir=False)
        implib_path = os.path.join(lib_dir, implib_name)
        lcmd += [ "/IMPLIB:%s" % os.path.abspath(implib_path) ]
    print(" ".join(lcmd))
    subprocess.check_call(" ".join(lcmd)) 
    print(dll_path + "\nhas been generated.\n")

    


def make_dll_nt_mingw32(cmd, ext):
    """Create a Windows DLL with the mingw compiler"""
    compile_args = [ "-c", "-Wall", "-DMS_WIN64"] + ext.extra_compile_args 
    for ipath in ext.include_dirs:
        compile_args.append("-I" + os.path.realpath(ipath))
    objects = []
    lib_dir = ext.get_lib_dir()

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
        objname = os.path.splitext(os.path.split(source)[-1])[0]
        obj = os.path.realpath(os.path.join(lib_dir, objname + ".o"))
        args.append(obj)
        objects.append(obj)
        print("gcc " + " ".join(args))
        subprocess.check_call("gcc " + " ".join(args)) 

    # Link
    largs = objects + ext.linker_library_options('mingw32')
    if ext.static_lib:
        raise DistutilsSetupError("Cannot build static library")
    else:
        dll_name = shared_lib_name(ext.name, 'load', pymod=True)
        dll_path =   os.path.join(cmd.get_package_dir(), dll_name)
        largs +=  ["-o",  dll_path ]
        implib_name = shared_lib_name(ext.name, 'build', pymod=True, dir=False)
        implib_path = os.path.join(lib_dir, implib_name)
        largs += [ "-Wl,--out-implib," + implib_path ]
        print("gcc " + " ".join(largs))
        subprocess.check_call("gcc " + " ".join(largs)) 
        print(dll_path + "\nhas been generated.\n")

    


def make_so_posix_gcc(cmd, ext):
    """Create a posix shared library with gcc"""
    compile_args = [ "-c", "-Wall"] + ext.extra_compile_args 
    for ipath in ext.include_dirs:
        compile_args.append("-I" + os.path.realpath(ipath))
    objects = []
    lib_dir = ext.get_lib_dir()

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
        objname = os.path.splitext(os.path.split(source)[-1])[0]
        obj = os.path.realpath(os.path.join(lib_dir, objname + ".o"))
        args.append(obj)
        objects.append(obj)
        print("cc " + " ".join(args))
        subprocess.check_call(["cc"] + args) 

    # Link
    if ext.static_lib:
        lib_name = shared_lib_name(ext.name, 'load',
           static=True, pymod=True, dir=False)
        dll_path =   os.path.join(lib_dir, lib_name)
        lcmd =  ["ar", "rcs", dll_path ] + objects
    else:
        lcmd = ["cc", "-shared",  "-Wall"] + ext.extra_link_args
        lcmd += objects + ext.linker_library_options('unix')
        dll_name = shared_lib_name(ext.name, 'load', pymod=True)
        dll_path =   os.path.join(cmd.get_package_dir(), dll_name)
        lcmd +=  ["-o",  dll_path ]
    print(" ".join(lcmd))
    subprocess.check_call(lcmd) 
    print(dll_path + "\nhas been generated.\n")

    







