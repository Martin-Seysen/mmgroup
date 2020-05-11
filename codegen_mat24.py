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


Generating the ``mmgroup.mat24_xi`` extension
.............................................

Function ``mat24_xi_make_c_code`` generates C code for computing the
monomial part of the operation of the elements :math`\xi` and
:math`\xi^2`of the monster group. These C functions are used for
computing (rather large) tables required for the implmentation of
the functions that compute the operation :math`\xi` and :math`\xi^2` 
on a representation of the monster.

The generation of the ``mmgroup.mat24_xi`` extension is similar
to the generation of the ``mmgroup.mat24`` extension. Here the
.ske file is file ``mat24_xi_functions.ske`` in subdirectory
``src/mmgroup/dev/mat24_xi``. Class  ``Mat24Xi`` in module
``mmgroup.dev.mat24_xi.mat24_xi_ref`` has methods ``tables()`` and
``directives()`` for creating the required tables and directives.
The functions used by the ``mmgroup.mat24_xi` extension are 
statically linked; so they are part of that entension. 
 
Location of the output files
............................

The location of the generated output files is controlled by certain
variables in module config.py. Each of these variables specifies the
name of a directory.

Files with extension .c, .h go to the directory ``C_DIR``. Files with
extension .pxd, .pxi, .pyx go to the directory ``PXD_DIR``. 
Documentaion files with extension .txt go to the directory ``DOC_DIR``.
"""



from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import types
import sys
import re
import os
import subprocess
import shutil



from config import SRC_DIR, DEV_DIR,  C_DIR, DOC_DIR, PXD_DIR
from config import REAL_SRC_DIR 
sys.path.append(REAL_SRC_DIR)

from mmgroup.dev.mat24.mat24_ref import Mat24
from mmgroup.dev.mat24_xi.mat24_xi_ref import Mat24Xi
from mmgroup.generate_c import TableGenerator, make_doc



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
    # The follwing two tables can't easily be computed earlier
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
    doc_path_ = os.path.join(DOC_DIR, MAT24_C_FILE)  
    make_doc(path_ + ".c", doc_path_ + "_c_doc.txt")
    print("C files for extension mat24 have been created" )
    


def mat24_xi_make_c_code():
    """Create .c and .h file with the functionality of class Mat24Xi

    """
    print("Creating C source from file mat24_xi_functions.ske\n")
    MAT24_XI_C_FILE = "mat24_xi_functions"
    SKE_DIR = os.path.join(DEV_DIR, "mat24_xi")
    Mat24Xi.tables["Mat24Xi_doc"] = Mat24Xi # can't do this earlier
    generator = TableGenerator(Mat24Xi.tables, Mat24Xi.functions)
    f = os.path.join(SKE_DIR, MAT24_XI_C_FILE)
    path_ = os.path.join(C_DIR, MAT24_XI_C_FILE)  
    generator.generate(f + ".ske", path_ + ".c", path_ + ".h")
    #generator.export_tables(file_name = "mat24_xi_export.py")

    generator.generate_pxd(
        os.path.join(PXD_DIR, MAT24_XI_C_FILE + ".pxd"), 
        MAT24_XI_C_FILE + ".h",  
        pxd_declarations   
    )
    doc_path_ = os.path.join(DOC_DIR, MAT24_XI_C_FILE)  
    make_doc(path_ + ".c", doc_path_ + "_c_doc.txt")
    print("C files for extension mat24_xi have been created" )


if __name__ == "__main__":
    mat24_make_c_code()    
    mat24_xi_make_c_code()    
    