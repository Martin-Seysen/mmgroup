from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import types
import sys
import re
import os
import subprocess


FILE_DIR = os.path.dirname(os.path.realpath(__file__))
DEV_DIR = os.path.dirname(FILE_DIR)
PACKAGE_DIR = os.path.dirname(DEV_DIR)
SRC_DIR = os.path.dirname(PACKAGE_DIR)

sys.path.append(SRC_DIR)
from mmgroup.dev.config import C_DIR , DOC_DIR
from mmgroup.dev.generate_c.make_c_tables import TableGenerator
from mmgroup.dev.generate_c.make_c_tables import make_doc
from mmgroup.dev.mat24.mat24_ref import Mat24
from mmgroup.dev.mat24_xi.mat24_xi_ref import Mat24Xi
assert sys.path.pop() == SRC_DIR

os.chdir(FILE_DIR)


MAT24_XI_C_FILE = "mat24_xi_functions"

pxd_declarations = """
from libc.stdint cimport uint32_t, uint16_t, uint8_t

"""



def mat24_xi_make_c_code():
    """Create .c and .h file with the functionality of class Mat24Xi

    will be moved to directory ma24_xi!!!
    """
    Mat24Xi.tables["Mat24Xi_doc"] = Mat24Xi # can't do this earlier
    generator = TableGenerator(Mat24Xi.tables, Mat24Xi.functions)
    f = MAT24_XI_C_FILE
    path_ = os.path.join(C_DIR, MAT24_XI_C_FILE)  
    generator.generate(f + ".ske", path_ + ".c", path_ + ".h")
    #generator.export_tables(file_name = "mat24_xi_export.py")

    generator.generate_pxd(MAT24_XI_C_FILE + ".pxd", 
           MAT24_XI_C_FILE + ".h",  pxd_declarations   )
    doc_path_ = os.path.join(DOC_DIR, MAT24_XI_C_FILE)  
    make_doc(path_ + ".c", doc_path_ + "_c_doc.txt")
    print("Mat24Xi tables of %d bytes length have been created" % ( 
                generator.table_size() ) )




def cythonize(del_intermediate_files = True):
    """Build the python extension Mat24Fast

    This is one of the few places where OS-specific code is required
    """
    res = subprocess.call("python setup.py build_ext --inplace -c mingw32")
    if res:
        print("\nError in:\npython setup.py build_ext --inplace -c mingw32")
        print("Process terminated with errorlevel", res)
        raise RuntimeError( "Could not build the Goly24Sub extension" )
    intermediate_files = [
              "./mat24fast.c",
    ]
    return
    if del_intermediate_files:
        subprocess.call("python setup.py clean")
        for f in intermediate_files:
            print( "removing" , f )
            os.remove(f) 

def test_fast_version(gcf):
    """Tiny regression test of fast version mat24fast.Mat24.

    Python dll mat24fast.pyd must exist.
    This dll is created e.g. by method Mat24.fast_version()

    gcf must be equal to mat24fast.Mat24.
    """
    from random import randint
    for i in range(10):
        v = randint(0,0xffffff)
        c1 = Mat24.vect_to_vintern(v)
        c2 = gcf.vect_to_vintern(v)
        assert c1 == c2
    #print( "Test of Cython Math24 code version passed" )


def fast_version(del_intermediate_files = True, rebuild = True): 
    """Create mat24fast.pyd, return class mat24fast.Mat24.

    That class is a fast C/Cython version of class Mat24. 
    It is built on the fly if not present and Cython and Mingw32 
    is available.

    Inermediate files of the build process are deleted if
    del_intermediate_files is True. The fast version 
    mat24fast.pyd is always built if 'rebuild' is True.

    """
    try:
        #assert not rebuild
        #import mat24fast
        1/0
    except:
        print("Creating C source from file mat24_xi_functions.ske\n")
        #mat24_make_c_code()
        mat24_xi_make_c_code()
        #print( "Build the python extension mat24fast" )
        #cythonize(del_intermediate_files)
        #print( "mat24fast extension created successfully\n" )
        #import mat24fast 
        #print("mat24fast extension imported")
    #fast = mat24fast.Mat24
    #print( "The fast Mat24.pyd module is now available" )
    #test_fast_version(fast) # make a tiny regeression test
    #return fast


if __name__ == "__main__":
    fast_version(del_intermediate_files = False, rebuild = True)
    