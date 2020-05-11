The standard build process for python packages is based on the
``setuptools`` package. We use ``setuptools`` to build the ``mmgroup`` 
package,  especially, for building extensions coded in the C language 
with the ``Cython`` package. We have added the following functionality 
to the standard build process:

  * The user may call arbitrary functions or processes, e.g.
    for the atomatic generation of C code.

  * The user may generate a shared library that may be used by
    several extensions.

Module ``build_ext_steps`` in the root directory contains the
python classes required for that extension. The following table 
lists some important modules in the root directory.


  .. table:: Files in the root directory used by the build process
    :widths: 30 70

    ======================= =============================================
    File name               Purpose
    ======================= =============================================
    ``build_ext_steps.py``  Extension for the ``setuptools`` package
    ``codegen_mat24.py``    Generates C code for extensions ``mat24``
                            and ``mat24_xi``
    ``codegen_mm.py``       Generates C code for extension ``mm``
    ``codegen_mm_op.py``    Generates C code for extensions ``mm<p>``,
                            where ``p`` is a modulus of a representation
                            of the monster
    ``config.py``           Main configuration file of the project
    ``import_all.py``       Imports all extensions (after building them)
    ``MANIFEST.in``         List of files to be added to the source
                            distribution
    ``pytest.ini``          Main configuration file for testing 
                            with ``pytest``
    ``README.rst``          Main documentation file for GitHub
    ``setup.py``            Main file for the build process with
                            ``setuptools``
    ======================= =============================================





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

Description of module ``build_ext_steps.py``
--------------------------------------------


.. automodule:: build_ext_steps


