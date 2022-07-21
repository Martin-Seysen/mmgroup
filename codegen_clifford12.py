r"""Generation of  C code dealing with the Mathieu group Mat24

Generating the ``mmgroup.clifford12`` extension
..........................................

Function ``clifford12_make_c_code()`` generates C code for basic 
computations in the Clifford group of up to 12 qubits.

The generated C modules are used in the python extension 
``mmgroup.clifford``. The functions used by that extension are 
contained in a shared library with name ``mmgroup_clifford12.dll``. 
The reason for creating such a shared library is that that these 
functions are also called by C functions written for other python 
extensions.

We use the C code generation mechanism in class 
``generate_c.TableGenerator``. Here a .c file  and a .h file is
crreated from files ``*.ske`` in subdirectory
``src/mmgroup/dev/clifford12``. The .ske file is like a .c file, but 
augmented with some code generation statements for entering tables 
and automatically generated code into the .c file to be generated. 
This .ske file may also have statements for automatically 
generating a .h file declaring the exported functions.  


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

import mmgroup
from mmgroup.generate_c import TableGenerator, make_doc
from mmgroup.generate_c import pxd_to_pyx
from mmgroup.dev.clifford12.bit64_tables import Bit64Tables



VERBOSE = 0

##########################################################################
# Generating .c files
##########################################################################




SKE_DIR = os.path.join(DEV_DIR, "clifford12")

SIMPLE_SKE_FILES = ["qstate12", "qstate12io", "qmatrix12"]
PXI_SKE_FILES = [
    "bitmatrix64",  "uint_sort", "xsp2co1", 
    "leech3matrix", "xsp2co1_elem",
    "involutions", "xsp2co1_traces",
]
SKE_FILES = SIMPLE_SKE_FILES + PXI_SKE_FILES


##########################################################################
# Generating .h and .pxd files files
##########################################################################



H_FILE_NAME = "clifford12.h"      # name of generated .h file


PXD_FILE_NAME = "clifford12.pxd"  # name of generated .pxd file


H_FILE_BEGIN = """

// %%EXPORT_KWD CLIFFORD12_API

// %%GEN h

/** @file clifford12.h
 File ``clifford.h`` is the header file for shared library 
 ``mmgroup_clifford12``. This comprises the C modules
 in the lists SIMPLE_SKE_FILES and PXI_SKE_FILES.
*/

#ifndef CLIFFORD12_H
#define CLIFFORD12_H

/// @cond DO_NOT_DOCUMENT 

#include <stdint.h>

#define CLIFFORD12_DLL  // We want a DLL!!


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define CLIFFORD12_HELPER_DLL_IMPORT __declspec(dllimport)
  #define CLIFFORD12_HELPER_DLL_EXPORT __declspec(dllexport)
  #define CLIFFORD12_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define CLIFFORD12_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define CLIFFORD12_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define CLIFFORD12_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define CLIFFORD12_HELPER_DLL_IMPORT
    #define CLIFFORD12_HELPER_DLL_EXPORT
    #define CLIFFORD12_HELPER_DLL_LOCAL
  #endif
#endif


// Now we use the generic helper definitions above to define CLIFFORD12_API 
// and CLIFFORD12_LOCAL.
// CLIFFORD12_API is used for the public API symbols. It either DLL imports 
// or DLL exports (or does nothing for static build). 
// CLIFFORD12_LOCAL is used for non-api symbols.

#ifdef CLIFFORD12_DLL // defined if CLIFFORD12 is compiled as a DLL
  #ifdef CLIFFORD12_DLL_EXPORTS // defined if we are building the CLIFFORD12 DLL 
                           // (instead of using it)
    #define CLIFFORD12_API CLIFFORD12_HELPER_DLL_EXPORT
  #else
    #define CLIFFORD12_API CLIFFORD12_HELPER_DLL_IMPORT
  #endif // CLIFFORD12_DLL_EXPORTS
  #define CLIFFORD12_LOCAL CLIFFORD12_HELPER_DLL_LOCAL
#else // CLIFFORD12_DLL is not defined: this means CLIFFORD12 is a static lib.
  #define CLIFFORD12_API
  #define CLIFFORD12_LOCAL
#endif // CLIFFORD12_DLL


/// @endcond 

"""


H_FILE_END = """
// %%GEN h

#ifdef __cplusplus
}
#endif
#endif  // #ifndef CLIFFORD12_H
"""



### Declarations for the generated .pxd files


PXD_DECLARATIONS = """

from libc.stdint cimport uint64_t, uint32_t, uint16_t, uint8_t
from libc.stdint cimport int64_t, int32_t, int8_t


cdef extern from "clifford12.h":
    enum: QSTATE12_MAXCOLS
    enum: QSTATE12_MAXROWS
    enum: QSTATE12_UNDEF_ROW
    ctypedef struct qstate12_type:
        uint32_t maxrows
        uint32_t nrows
        uint32_t ncols
        int32_t  factor
        uint32_t shape1
        uint32_t reduced
        uint64_t *data 


""".format()




##########################################################################
# Entering tables and automatically generated code
##########################################################################


TABLE_CLASSES = [
    Bit64Tables,
]



##########################################################################
# Generating the .pyx file
##########################################################################


##########################################################################
# Generating the .pyx file
##########################################################################



PXI_FILE_NAME = "clifford12.pxi"



##########################################################################
# Functions of this module for generating .c, .h and .pxi files
##########################################################################



def list_c_files():
    """Return list of names of c files generated by make_clifford12()"""
    def c_file(name):
         return os.path.join(C_DIR, name + ".c")
    return [c_file(name) for name in SKE_FILES]




def make_clifford12():
    """Generate clifford12 .c files from  .ske files 

    The relevant .ske files are listed in variable
    SKE_FILES. One corresponding .h and .pxd file
    is also generated. The .c files and the .h file are written 
    to directory given by C_DIR, the .pxd file is written to the 
    current directory.

    Tables and directives for automatically generated C code
    are taken from the classes listed in BASIC_TABLE_CLASSES.

    Return pair of lists, one of the .c files and one of the.
    .pxd files generated
    """
    tables = {}
    directives = {}
    global generated_tables
    for table_class in TABLE_CLASSES:
        table_instance = table_class()
        tables.update(table_instance.tables)
        directives.update(table_instance.directives)
    # print("Clifford 12 functions:\n",  directives.keys())
    tg = TableGenerator(tables, directives, verbose = VERBOSE)
    # first generate C files
    c_files = []
    all_ske_files = [H_FILE_BEGIN]
    for name in SKE_FILES:
        ske_file = name + ".ske"
        ske_path = os.path.join(SKE_DIR, ske_file)
        c_file = name + ".c"
        c_path = os.path.join(C_DIR, c_file)
        print("Creating %s from %s" % (c_file, ske_file))
        tg.generate(ske_path, c_path)
        all_ske_files.append(ske_path)
        c_files.append(c_path)
    all_ske_files.append (H_FILE_END)

    # generate .h and .pxd file
    h_file =  H_FILE_NAME
    h_path =  os.path.join(C_DIR, h_file)
    pxd_file =  PXD_FILE_NAME
    print("Creating %s from previous .ske files" % h_file)
    tg.generate(all_ske_files, None, h_path)
    tg.generate_pxd(
        os.path.join(PXD_DIR, pxd_file), 
        h_file, 
        PXD_DECLARATIONS
    )
    return c_files,  [ pxd_file ]


##########################################################################
# The main function for generating code
##########################################################################


def generate_files():
    """The main function of this module for generating code

    This function generates all reqired .c files and also the required
    headers.

    It also generates the .pxi file with name given by PYX_FILE_NAME and
    stores it in the current directory. Cython will use the .pxi file to
    build a wrapper for all generated C functions. 
    """
    def pxi_comment(text, f):
        print("\n" + "#"*70 + "\n### %s\n" % text + "#"*70 + "\n\n",file=f)
    c_files,  pxd_files =  make_clifford12() 
    f_pxi = open(os.path.join(PXD_DIR, PXI_FILE_NAME), "wt")
    print(PXD_DECLARATIONS, file = f_pxi)
    
    for pxd_f in pxd_files:
        pxi_comment("Wrappers for C functions from file %s" % pxd_f, f_pxi)
        pxi_content = pxd_to_pyx(
            os.path.join(PXD_DIR, pxd_f),
            os.path.split(pxd_f)[0],
            select = True
        )
        print(pxi_content, file = f_pxi)
     
    f_pxi.close()
    return c_files



##########################################################################
# Main program
##########################################################################


if __name__ == "__main__":
    generate_files()






