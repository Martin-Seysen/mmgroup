from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import types
import sys
import re
import os
import subprocess
import shutil


FILE_DIR = os.path.dirname(os.path.realpath(__file__))
DEV_DIR = os.path.dirname(FILE_DIR)
PACKAGE_DIR = os.path.dirname(DEV_DIR)
SRC_DIR = os.path.dirname(PACKAGE_DIR)

sys.path.append(SRC_DIR)
from mmgroup.dev.config import C_DIR , DOC_DIR
from mmgroup.dev.generate_c.make_c_tables import TableGenerator
from mmgroup.dev.generate_c.make_c_tables import make_doc
from mmgroup.dev.mat24.mat24_ref import Mat24
assert sys.path.pop() == SRC_DIR

os.chdir(FILE_DIR)


MAT24_C_FILE = "mat24_functions"



pxd_declarations = """
from libc.stdint cimport uint32_t, uint16_t, uint8_t

"""

def mat24_make_c_code():
    """Create .c and .h file with the functionality of class Mat24

    The input of this function is the file MAT24_C_FILE.ske that 
    contains a (hopefully much faster) C version of the functions
    in this class.

    The functions in the .ske file make use of the tables that have 
    been generated for this module and also of some functions for
    generating C code automatically. An example where generating
    code automatically makes sense is the matrix multiplication with
    a constant bit matrix.
     
    The code generating process is described in class TableGenerator
    in module make_c_tables.  
    """ 
    # The follwing two tables can't easily be computed earlier
    Mat24.tables["Mat24_doc"] = Mat24.__doc__ 
    Mat24.tables["Mat24_doc_basis"] = Mat24.str_basis()
    generator = TableGenerator(Mat24.tables, Mat24.directives)
    f = MAT24_C_FILE
    path_ = os.path.join(C_DIR, MAT24_C_FILE)  
    generator.generate(f + ".ske", path_ + ".c", path_ + ".h")
    ## generator.export_tables(file_name = "mat24_export.py")

    generator.generate_pxd(MAT24_C_FILE + ".pxd", 
           MAT24_C_FILE + ".h",  pxd_declarations   )
    doc_path_ = os.path.join(DOC_DIR, MAT24_C_FILE)  
    make_doc(path_ + ".c", doc_path_ + "_c_doc.txt")
    print("Golay code tables of %d bytes length have been created" % ( 
                generator.table_size() ) )
    


if __name__ == "__main__":
    print("Creating C source from file mat24_functions.ske\n")
    mat24_make_c_code()
    