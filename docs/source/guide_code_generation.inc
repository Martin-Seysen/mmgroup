

The main code generation tool *generate_code.py*
------------------------------------------------

An invocation of the main code generation script *generate_code.py*
generates a set of .c files and also a single .h file from a set of
source files. Here each .c file is generated from a single source file
which has extension .ske. Prototypes for the functions in the .c files
may also be generated form these source files; they are copied into
the generated .h file. The user may also specify a source file with an
extension .h. Such a .h file usually contains global declarations;
and it is also copied into the generated .h file.
 
The user may also specifiy a set of python modules to be used for
code generation with the ``--tables`` option. Any such module must
have a class ``Tables`` containing information how to enter tables
and code snippets into the generated C code. The structure of such a
code-generating module is described in the next section.

We use the *Cython* language to integrate the generated .c file into
our python project. In a Cython project, a .pxd file should be 
generated from a header file. Here the .pxd file contains essentially
the same information as the header file. The ``--pxd`` option generates
a .pxd file from the generated header file automatically. Here the user
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


Python modules used for code generation
----------------------------------------

In the main code generation tool  *generate_code.py* the user may specify
a list of python modules to be used for code generation. Any such module
must have a class ``Tables`` that contains two dictionaries ``tables`` and
``directives`` to be used for code generation. Both of these dictionaries
map identifiers to objects implementing tables or directives to be used in 
the code generation process. Here the identifiers should be valid python 
identifiers beginning with an alphabetic character. If several code 
generation modules are specified then the corresponding dictionaries in
the ``Tables`` classes of  these modules are merged into a single
dictionary.

Dictionary ``tables`` maps identifiers to values. If we have
e.g. ``tables[P] = 3`` then the expression ``%{P}`` in a source file
with extension .ske evaluates to the string ``3`` in the generated .c
file. More elaborated use cases for the ``tables`` dictionary are given 
in the following sections. A value in the ``tables`` dictionary may also 
be a list of integers; then there ia a way to to initialize an array of
integers in the generated .c files with these integer values.   

Dictionary ``directives`` maps identifiers to directives. Invoking
such a directive in a source file causes the code generator to put
a code snippet into the generated .c file. Here a typical use case is the
multiplication of a part of a vector of the representation of the Monster
group with a fixed matrix.

A class ``Tables`` as described above should accept arbitrary keyword
arguments in the constructor. (A standard way to do this is to
provide a parameter ``**kwds`` in the constructor.) The 
*generate_code.py* tool may pass such keyword arguments to
the constructor of a class ``Tables`` either with the ``--set``
option or with the ``--subst`` option. Using the ``--set`` option
is straightforward. Using ``--subst`` is a bit more involved.
The following paragraph provides an example for using the ``--subst``
option.

The source file ``mm_op_pi.ske`` supports certain oparations of the 
Monster group on vectors modulo several small odd numbers ``p``. The
module ``mm_op_pi.py`` provides a class ``Tables`` that supports
the generation of code snippets to be used in these operations. 
From the source file ``mm_op_word.ske`` we want to generate files
``mm3_op_pi.c``, ``mm7_op_pi.c``, ... etc. for the operation modulo
3, 7, ..., etc. In the corresponding code generation process there
is an option 

``--subst mm(?P<p>[0-9]+)_op mm_op``

that primarily describes a subtitution of strings that maps both, 
``mm3_op_pi``, and  ``mm7_op_pi`` to ``mm_op_pi``. So this
substitution maps the name of a generated .c file to the name of a
.ske file (without extension) from which that .c file is to be
generated. The first argument ``mm(?P<p>[0-9]+)_op`` of the 
``--subst`` option is the pattern to be substituted given in the
standard python regular expression syntax. Note that the part
``(?P<p>[0-9]+)`` of that argument also creates a
dictionary ``{'p' : '3'}`` mapping identifier 'p' to the string
'3', when we substitute the name ``mm3_op_pi`` using this pattern.
When such a dictionary is created, the code generation tool
passes the keyword argument ``p = '3'`` to the constructor of
every class ``Tables`` used for the generation of file
``mm3_op_pi.c``. This way the the ``Tables`` classes can be
instructed to generate code snippets for operation modulo 3.




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

    

Support for Sphinx *mockup* feature when generating documentation
--------------------------------------------------------------------

We use the *Sphinx* package on the *Read the Docs* server for
generating documentation. On that server we do not want to compile
.c files. However, we want the .c files to be present for extracting
documentation from them, using the *Doxygen* tool. Here it suffices
if the prototypes and comments relevant for documentation are
present in the .c files. 

This means that we have to perform the code generation process also
on the documentation server. Note that some code generation procedures
in the ``Tables`` classes in the python modules used for code generation
rely on compiled Cython extensions. But Cython extensions will not
be compiled on the documentation server. To avoid using compiled code,
in each python module to
be used for code generation we may add a ``MockupTables`` class to
the ``Tables`` class already present. Such a ``MockupTables`` class
may generate the parts of a .c file relevant for documentation only,
without using any compiled Cython extensions. 

We may set the ``--mockup`` option in the code generation tool
*generate_code.py*. Then any ``Tables`` class in a python module to
be used for code generation is replaced by the ``MockupTables``
class in that module, if present.
