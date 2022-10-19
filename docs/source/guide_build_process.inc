The standard build process for python packages is based on the
``setuptools`` package. We use ``setuptools`` to build the ``mmgroup`` 
package,  especially, for building extensions coded in the C language 
with the ``Cython`` package. We have added the following functionality 
to the standard build process:

  * The user may call arbitrary functions or processes, e.g.
    for the automatic generation of C code.

  * The user may generate a shared library that may be used by
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
    ``codegen_clifford12.py`` Generates C code for extension
                              ``clifford12``
    ``codegen_mat24.py``      Generates C code for extensions ``mat24``
                              and ``generators``
    ``codegen_mm.py``         Generates C code for extension ``mm``
    ``codegen_mm_op.py``      Generates C code for extensions ``mm<p>``,
                              where ``p`` is a modulus of a representation
                              of the monster
    ``config.py``             Main configuration file of the project
    ``import_all.py``         Imports all extensions (after building them)
    ``MANIFEST.in``           List of files to be added to the source
                              distribution
    ``pytest.ini``            Main configuration file for testing 
                              with ``pytest``
    ``README.rst``            Main documentation file for GitHub
    ``setup.py``              Main file for the build process with
                              ``setuptools``
    ========================= ============================================





Porting the project
-------------------

Porting a larger project to a different operating system or compiler
toolchain requires a thorough knowledge of the target system. 
Here the ``setuptools`` package comes handy, since it supports
several operating systems and compilers. 

However, we use some specific extensions for ``setuptools`` coded
in file ``build_ext_steps.py``; and the build process for the 
``mmgroup`` project may take a long time.

There is a much simpler demonstration project (with name ``miniproject``)
using exactly the same extensions for its build process with 
``setuptools``, see: 

https://github.com/Martin-Seysen/mmgroup_miniproject .

The files ``build_ext_steps.py`` are exactly the same in both projects.
For porting the ``mmgroup`` project we recommend to port the 
``miniproject`` first. This will almost surely require modifications
in file ``build_ext_steps.py``. These modifications can then be 
copied to the  ``mmgroup`` project. 



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


