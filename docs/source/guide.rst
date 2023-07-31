

================================
The mmgroup guide for developers
================================


Introduction
============

This *guide for developers* explains some selected topics of the
``mmgroup`` project which are relevant for developers who want to make 
changes or bugfixes. It also explains some mathematical aspects of the 
implementation which are not covered by :cite:`Seysen20`.

The *guide* is far from complete and we hope that we can provide 
more information in future versions of it.

Directory structure
===================

.. include:: guide_directory.inc


Some mathematical aspects of the implementation
===============================================

.. include:: guide_mathematical.inc


.. include:: guide_mathematical_2.inc



.. _install_from_source_label:

Installation from a source distribution
=======================================


.. include:: installation_from_source.inc



.. _build-process-label:


The build process
=================


Warning! At present the build process is under construction.

We plan to switch to the ``Meson`` build system. Therefore a much more
declarative style is required for describing the build operations.
Thus the description here is pretty much outdated!


.. include:: guide_build_process.inc



.. _code-generation-label:

Code generation
===============

Warning! At present the code generation process is under construction.

We plan to switch to the ``Meson`` build system. Therefore a much more
declarative style is required for describing the code generation.
Also, a strict separation between input and output directories is
required.
Thus the description here is pretty much outdated!


The main code generation tool *generate_code.py*
------------------------------------------------

An invocation of the main code generation script *generate_code.py*
generates a set of .c files and also a single .h file from a set of
source files. Here each .c file is generated from a single source file
which has extension .ske. Prototypes for the functions in the .c files
may also be generated form these source files; they are copied into
the genrated .h file. The user may also specify a source file with an
extension .h. Such a .h file usually contains global declarations;
and it is also copied into the generated header file.
 
The user may also specifiy a set of python modules to be used for
code generation with the ``--tables`` option. Any such module must
have a class ``Tables`` containing information how to enter tables
and code snippets into the generated C code.

TODO: Explain basic structure of class ``Tables`` !!!

We use the *Cython* language to integrate the generated .c file into
our python project. In a Cython project, a .pxd file must be generated
from a header file. Here the .pxd file contains essentially the same
information as the header file. The ``--pxd`` option generates a .pxd
file from the generated header file automatically. Here the user
may specify a source file (with extension .pxd), which will also be 
copied into the generated .pxd file.  

The user may also create a .pxi file from the .pxd file with the
``--pxi`` option . Such a .pxi file will contains wrappers for the 
C functions declared in the .pxi file. These wrappers can be used 
directly from python by including the .pxi file into a .pyx file
that defines a Cython extension. Note that this automatic wrappping
mechanism works for rather simple prototypes only.

With the ``--pyx`` option, the user may simply copy a .pyx file from
the source path to a destination directory.

The Meson build system requires strict separation between the input
and the output directory. The user may specify a path where to look 
for the input files (with extensions .ske, .h, .pxd, and .pyx).
The user may specify a directory where to store the generated
.c and .h files; and also a directory where to store the generated
.pxd, .pyi, and .pyx files. 

In the following table we list the options available in the
*generate_code.py* tool.


.. code-block:: text

  -h, --help            show this help message and exit
  --sources [SOURCE ...]
                        List of SOURCE files to be generated. For each
                        SOURCE with extension '.c' a '.c' file is
                        generated from a file with the same name and
                        extension '.ske'. A SOURCE with extension '.h' is
                        copied into the common header file. Each SOURCE
                        is seached in the path set by parameter '--
                        source-path'. Output is written to the directory
                        set by parameter '--out-dir'.
  --out-header HEADER   Set name of output header file to HEADER. By
                        default we take the first file with extension .h
                        in the list given in the argument '--sources'.
  --pxd PXD             Set input '.pxd' file PXD for generating '.pxd'
                        file with same name from that input file and from
                        the generated header.
  --pxi                 Generate '.pxi' file from '.pxd' file if this
                        option is set.
  --pyx PYX             Copy input '.pyx' file PYX from source path to
                        output directory.
  --tables [TABLES ...]
                        Add list TABLES of table classes to the tables to
                        be used for code generation.
  --set VAR=VALUE [VAR=VALUE ...]
                        Set variable VAR to value VALUE. When generating
                        code with subsequent '--sources' options then the
                        table classes will set VAR=VALUE.
  --subst [PATTERN SUBSTITUTION ...]
                        Map the name of a '.c' or '.h' file to be
                        generated to the name of a file used as a source
                        for the generation process. Example: "--subst
                        mm(?P<p>[0-9]+)_op mm_op" maps e.g.'mm3_op' to
                        'mm_op'. We substitute PATTERN in the name of a
                        generated file by SUBSTITUTION for obtaining the
                        name of that source. The part "(?P<p>[0-9]+)"
                        creates a variable 'p' that takes the decimal
                        string following the intial letters 'mm' in the
                        file name. Then variable 'p' will be passed to
                        the table classes used for code generation.
                        PATTERN must be given in python regular
                        expression syntax.
  --source-path [PATHS ...]
                        Set list of PATHS for finding source files to be
                        used for generation process.
  --py-path [PATHS ...]
                        Set list of PATHS for finding python scripts
  --out-dir DIR         Store output files with extensions .c and .h in
                        directory DIR.
  --out-pxd-dir DIR     Store output files with extensions .pxd, .pxi,
                        and .pyx in directory DIR.
  --dll DLL_NAME        Generate code for exporting C functions to a DLL
                        or to a shared library with name DLL_NAME.
                        Parameter DLL_NAME must be the same for all
                        generated C files to be placed into the same DLL.
                        This parameter is not used for any other
                        purposes.
  --nogil               Optional, declare functions from .pxi file as
                        'nogil' when set.
  --mockup              Use tables and directives for Sphinx mockup if
                        present
  -v, --verbose         Verbose operation



Generation of code snippets and tables
---------------------------------------

In this section we describe the most important functions and classes
used for the automatic generation of C code.
  

.. automodule:: mmgroup.generate_c.make_c_tables_doc


Classes and functions provided by the code generator
----------------------------------------------------
 

.. autoclass:: mmgroup.generate_c.TableGenerator
   :members: generate


.. autoclass:: mmgroup.generate_c.UserDirective

.. autoclass:: mmgroup.generate_c.UserFormat

          
.. autofunction::   mmgroup.generate_c.c_snippet  

.. autofunction::   mmgroup.generate_c.format_item

.. autofunction::   mmgroup.generate_c.make_doc

.. autofunction::   mmgroup.generate_c.make_table

.. autofunction::   mmgroup.generate_c.prepend_blanks
   
.. autofunction::   mmgroup.generate_c.pxd_to_pxi
  
.. autofunction::   mmgroup.generate_c.generate_pxd

    

How the code generator is used
------------------------------

.. include:: guide_use_code_generation.inc



Description of some typical bugs
================================

.. include:: bugs.inc

