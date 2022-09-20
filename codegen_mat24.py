r"""Generation of  C code dealing with the Mathieu group Mat24

Generating the ``mmgroup.mat24`` extension
..........................................

Function ``mat24_make_c_code()`` generates C code for basic computations
in the Golay code, its cocode, and the Mathieu group Mat24. It also
generates code for computations in the Parker loop and in its 
automorphism group. 

The generated C modules are used in the python extension 
``mmgroup.mat24``. The functions used by that extension are contained 
in a shared library with name ``mmgroup_mat24.dll``. The reason for
creating such a shared library is that that these functions are also
called by C functions written for other python extensions.

We use the C code generation mechanism in class 
``generate_c.TableGenerator``. Here a .c file  and a .h file is
crreated from file ``mat24_functions.ske`` in subdirectory
``src/mmgroup/dev/mat24``. The .ske file is like a .c file, but 
augmented with some code generation statements for entering tables 
and automatically generated code into the .c file to be generated. 
This .ske file may also have statements for automatically 
generating a .h file declaring the exported functions.  

We create an instance ``tg`` of class ``TableGenerator`` for 
generating the .c files. The table generator ``tg`` takes two
dictionaries ``tables`` and ``directives`` as arguments. These
dictionaries provide user-defined tables and directives for the
code generator. Class  ``Mat24`` in module
``mmgroup.dev.mat24.mat24_ref`` has methods ``tables()`` 
and ``directives()`` creating the required tables and directives.


Generating the ``mmgroup.generators`` extension
.............................................

Function ``generators_make_c_code`` generates C code for computing 
the monomial part of the operation of the elements :math`\xi` and
:math`\xi^2`of the monster group. These C functions are used for
computing (rather large) tables required for the implmentation of
the functions that compute the operation :math`\xi` and :math`\xi^2` 
on a representation of the monster.

The generation of the ``mmgroup.generators`` extension is similar
to the generation of the ``mmgroup.mat24`` extension. Here the
list of .ske file is given in the list GENERATORS_C_FILES.
For each file in that list a C file is created. 

A common header with name given by H_GENERATORS_NAME is created 
from all these .ske files, prependend by the header files in
the list GENERATORS_H_FILES. A .pxd file with name
PXD_GENERATORS_NAME is created from that header file. That .pxd
file will also contain the declarations in the string
PXD_DECLARATIONS.

All input files are read fom the directory SKE_DIR.


 
Location of the output files
............................

The location of the generated output files is controlled by certain
variables in module config.py. Each of these variables specifies the
name of a directory.

Files with extension .c, .h go to the directory ``C_DIR``. Files with
extension .pxd, .pxi, .pyx go to the directory ``PXD_DIR``. 
"""



from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import types
import sys
import re
import os
import subprocess
import shutil



from config import SRC_DIR, DEV_DIR,  C_DIR, PXD_DIR
from config import REAL_SRC_DIR 
sys.path.append(REAL_SRC_DIR)

from mmgroup.dev.mat24.mat24_ref import Mat24
from mmgroup.dev.generators.gen_xi_ref import GenXi
from mmgroup.dev.generators.gen_leech_reduce_n import GenLeechReduceY
from mmgroup.dev.generators.gen_cocode_short import ShortCocodeTables
from mmgroup.generate_c import TableGenerator, make_doc
from mmgroup.generate_c import pxd_to_pyx


########################################################################
# Generate mat24_functions.c
########################################################################

pxd_declarations = """
from libc.stdint cimport uint32_t, uint16_t, uint8_t

"""



def mat24_make_c_code():
    """Create .c and .h file with the functionality of class Mat24

    The input of this function is the file MAT24_C_FILE.ske that 
    contains a (much faster) C version of the functions in class
    mmgroup.dev.mat24.mat24_ref.Mat24.

    The functions in the .ske file make use of the tables that have 
    been generated for this module and also of some functions for
    generating C code automatically. An example where generating
    code automatically makes sense is the matrix multiplication with
    a constant bit matrix.
     
    The code generating process is described in class TableGenerator
    in module make_c_tables.  
    """ 
    print("Creating C source from file mat24_functions.ske\n")
    MAT24_C_FILE = "mat24_functions"
    SKE_DIR = os.path.join(DEV_DIR, "mat24")
    # The following two tables can't easily be computed earlier
    Mat24.tables["Mat24_doc_basis"] = Mat24.str_basis()
    generator = TableGenerator(Mat24.tables, Mat24.directives)
    f = os.path.join(SKE_DIR, MAT24_C_FILE)
    path_ = os.path.join(C_DIR, MAT24_C_FILE)  
    #print("pwd", os.getcwd())
    #print(os.path.realpath(path_ + ".c"))
    generator.generate(f + ".ske", path_ + ".c", path_ + ".h")
    ## generator.export_tables(file_name = "mat24_export.py")

    generator.generate_pxd(
        os.path.join(PXD_DIR, MAT24_C_FILE + ".pxd"), 
        MAT24_C_FILE + ".h",  
        pxd_declarations   
    )
    print("C files for extension mat24 have been created" )
    


########################################################################
# Generate c files for module 'generators'
########################################################################

SKE_DIR =  os.path.join(DEV_DIR, "generators")

GENERATORS_C_FILES = [
    "gen_xi_functions",
    "mm_group_n",
    "gen_leech",
    "gen_leech_type",
    "gen_leech3",
    "gen_leech_reduce",
    "gen_leech_reduce_n",
    "gen_random",
]


GENERATORS_H_START = """
// %%GEN h
#ifndef MMGROUP_GENERATORS_H
#define MMGROUP_GENERATORS_H
// %%GEN c

"""

GENERATORS_H_END = """
// %%GEN h
#endif // ifndef MMGROUP_GENERATORS_H
// %%GEN c
"""

GENERATORS_H_FILES = [
    GENERATORS_H_START,
    "mmgroup_generators.h",
]

GENERATORS_TABLE_CLASSES = [
     GenXi, GenLeechReduceY, ShortCocodeTables
]

H_GENERATORS_NAME = "mmgroup_generators.h"
PXD_GENERATORS_NAME = "generators.pxd"
PXI_GENERATORS_NAME = "generators.pxi"

PXD_DECLARATIONS = """

from libc.stdint cimport uint64_t, uint32_t, uint16_t, uint8_t
from libc.stdint cimport int64_t, int32_t

"""


def generators_make_c_code():
    """Create .c and .h file with the functionality of class Mat24Xi

    """
    print("Creating C sources for the 'generators' extension\n")

    # Setp table and directives for code generation
    GenXi.tables["GenXi_doc"] = GenXi # can't do this earlier
    tables = {}
    directives = {}
    for table_class in GENERATORS_TABLE_CLASSES:
        table_instance = table_class()
        tables.update(table_instance.tables)
        directives.update(table_instance.directives)
    print(tables.keys())
    tg = TableGenerator(tables, directives)

    # Generate c files
    all_ske_files = [os.path.join(SKE_DIR, name) 
        for name in GENERATORS_H_FILES]
    for name in GENERATORS_C_FILES:
        ske_file = name + ".ske"
        ske_path = os.path.join(SKE_DIR, ske_file)
        c_file = name + ".c"
        c_path = os.path.join(C_DIR, c_file)
        print("Creating %s from %s" % (c_file, ske_file))
        tg.generate(ske_path, c_path)
        all_ske_files.append(ske_path)

    # generate .h file
    all_ske_files.append(GENERATORS_H_END)
    h_file =  H_GENERATORS_NAME
    h_path =  os.path.join(C_DIR, h_file)
    pxd_file =  PXD_GENERATORS_NAME
    print("Creating %s from previous .ske files" % h_file)
    tg.generate(all_ske_files, None, h_path)

    # generate .pxd file
    tg.generate_pxd(
        os.path.join(PXD_DIR, PXD_GENERATORS_NAME), 
        h_file, 
        PXD_DECLARATIONS
    )
    print("C files for extension 'generators' have been created" )
        
    # generate .pxi file
    def pxi_comment(text, f):
        print("\n" + "#"*70 + "\n### %s\n" % text + "#"*70 + "\n\n",
        file=f
    )
    f_pxi = open(os.path.join(PXD_DIR, PXI_GENERATORS_NAME), "wt")
    pxi_comment(
        "Wrappers for C functions from file %s" % PXD_GENERATORS_NAME, 
        f_pxi
    )
    print(PXD_DECLARATIONS, file = f_pxi)
    pxi_content = pxd_to_pyx(
        os.path.join(PXD_DIR, PXD_GENERATORS_NAME),
        os.path.split(PXD_GENERATORS_NAME)[0],
        select = True
    )
    print(pxi_content, file = f_pxi)
    f_pxi.close()

 
########################################################################
# Main program
########################################################################

if __name__ == "__main__":
    mat24_make_c_code()    
    generators_make_c_code()    
    