The standard build process for python packages is based on the
``setuptools`` package. We use ``setuptools`` to build the ``mmgroup`` 
package. We have added the following functionality 
to the standard build process:

  * The user may call arbitrary functions or processes, e.g.
    for the automatic generation of C code.

  * The user may add a shared library that may be used by
    several extensions.

Module ``build_ext_steps`` in the root directory contains the
python classes required for that extension. The following table 
lists some important modules in the root directory.


  .. table:: Files in the root directory used by the build process
    :widths: 30 70

    ========================= ============================================
    File name                 Purpose
    ========================= ============================================
    ``build_ext_steps.py``    Extension for the ``setuptools`` package
    ``build_shared.py``       Tool for building a shared library
    ``cleanup.py``            Clean up intermediate files
    ``config.py``             A configuration file of the project
    ``generate_code.py``      Tool for generating C code
    ``MANIFEST.in``           List of files to be added to the source
                              distribution
    ``pyproject.toml``        Main configuration file for the project
    ``README.rst``            Main documentation file for GitHub
    ``setup.py``              Main file for the build process with
                              ``setuptools``
    ========================= ============================================





The tool build_shared.py
------------------------

The tool ``build_shared.py`` is a simple stand-alone tool for
generating a shared library with a standard C compiler. It supports
Windows, Linux and macOS. This job could also be done with a tool 
like ``CMake`` or ``Meson``; but for the time being it has been 
easier to extract this functionality from an older tool with a 
similar functionality.

Here are the options of the tool ``build_shared.py``:

.. code-block:: text

  -h, --help            show this help message and exit
  --name MAME           Set NAME of library to be generated. NAME should
                        have no extension and no prefix such as 'lib'.
  --sources [SOURCES ...]
                        Add list SOURCES of to list of source files. Each
                        source file should be a .c file with extension '.c'.
  --source-dir DIR      Set root directory DIR for all source files.
  --include-path [PATHS ...]
                        Set list of PATHS for finding .h files to be used.
  --library-path [PATHS ...]
                        Set list of PATHS for finding libraries to be used.
  --define [VARS ...]   Define variables given by the list VARS for the
                        compiler. Each variable must be a string VAR or
                        VAR=VALUE.
  --undef [VARS ...]    Undefine variables given by the list VARS for the
                        compiler.
  --libraries [LIBS ...]
                        Search libraries with names in the list given by
                        LIBS. This corresponds to the gcc option '-l'.
  --library-dir DIR     Set directory DIR for storing object files and
                        static libraries.
  --shared-dir DIR      Set directory DIR for storing shared libraries.
  --static [STATIC]     Create static instead of shared library. Optional
                        argument STATIC may be 0 (=shared) or 1 (=static).
  --n N                 Use N parallel processes for compiling.
  --compiler C          Specify name of default compiler. C must be 'unix',
                        'msvc', or 'mingw32'.
  --cflags FLAGS        Add extra arguments FLAGS for the compiler. E.g. 
                        '--cflags=-c,foo=bar' adds arguments '-c' and
                        'foo=bar'.
  --display MODE        display (full) name of generated library and exit.
                        MODE = 'load' displays the full file name of the
                        library. MODE = 'build' displays the name of the
                        library to be linked to a program using that
                        library.




Creating a new version
----------------------

The ultimate goal of building a new version is to upload a new python 
version to  the **pypi** server. 

At present we will upload a source distribution and python 3 wheels
for 64-bit Windows, Linux and macOS, say, for the latest two or three 
python versions.

Before creating a new version, (at least) the following test should 
be executed in a shell in Windows 64, and also in some Linux version:

 .. code-block:: text

    pytest src/mmgroup/ -v -s -m "very_slow"

Version numbering 
.................

We assume that we want to create Version 0.0.8 at date 2022-07-12,
with the short version description 'Performance improved'.

You should update the version number in file **setup.py** by 
writing e.g:

 .. code-block:: python

    VERSION = '0.0.8' # 2022-07-12. Performance improved

into that file.
You should also comment out older version descriptions in that file.
In section **Version history** of file **docs/source/api.rst**
you should add the following line to the version history:

 .. code-block:: text

    | Version 0.0.8, 2022-07-12. Performance improved

Then you should upload the new version with these changes to the
``master`` branch in the github repository:

https://github.com/Martin-Seysen/mmgroup

After uploading, you should create a new release in the github repository.
Therefore, click **Create a new release** in  main page of the github
repository. Here you should write the tag **v0.0.8** into the field
**Choose a tag**. The **Target** of the release should be **master**,
referring to the master branch in git. You should enter the title 
**mmgroup v0.0.8** into the field **Release title**.
We recommend to enter (at least)
the version description (which is 'Performance improved' in our case)
into the field **Describe this release**. Finally, you should click the 
button **Publish release** in the github window.


Generating the wheels 
.....................


This subsection describes how to create wheels manually. This process
has now been automated to some extent by using ``GitHub`` actions
that trigger the python tool ``cibuildwheel``. So the reader may skip
this section.

Here you must generate a wheel for each python version, and also for
each operating system that you want to support. Here we assume that
Anaconda is used for creating wheels for Windows 64 for various 
python versions. An Anaconda environment e.g for python 3.9 is
created by typing the following command

 .. code-block:: text

    conda create --name python39 python = 3.9

Then we may switch to that python version by typing:

 .. code-block:: text

    conda activate python39

Environments for other python versions are created similarly.
One has to install all required python packages for each version.
For uploading a version to **pypi** we also have to install **twine**
with 

 .. code-block:: text

    pip install twine

In each version to be supported we have to type:

 .. code-block:: text

    python setup.py build_ext bdist_wheel
    python setup.py sdist

Here the first line creates the wheel for the selected python version.
Before doing so in Windows, you must install a C compiler for cython
as described in section :ref:`install_from_source_label`.

The second line creates a source distribution; this must be done
only once. The wheels and source distributions are stored in 
subdirectory **dist**. The wheel for mmgroup version 0.0.8 for python
3.9 has a name similar to ``mmgroup-0.0.9-cp37-cp37m-win_amd64.whl``;
and the source distribution has a name like  ``mmgroup-0.0.8.zip``.


Uploading the version to **pypi** 
..................................

You may upload the new version with the following commend:


 .. code-block:: text

    twine upload twine upload dist/*

This uploads all files from subdirectory **dist** to **pypi**.
So you'd better cleanup that directory before uploading.

.. warning::
   Uploading with **twine** is irreversible. If your uploaded version
   is buggy, you will have to create a new version!


Description of module ``build_ext_steps.py``
--------------------------------------------


.. automodule:: build_ext_steps


