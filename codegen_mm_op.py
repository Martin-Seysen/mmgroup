r"""Generation of C code for the monster representations modulo ``p``

This script generates the C code for computing in the monster 
representations modulo several small odd numbers ``p`` of the form
``p = 2**l - 1``, ``1 < k <= 8``. The list of integers ``p``
for which code will be genrated is given by the variable
``PRIMES`` in file ``config.py``. The generated C modules are used 
in the python extension ``mmgroup.mm{P}``. Here ``{P}`` is to be
replaced by the actaul modulus ``p``, so the extension for ``p = 3``
has the name ``mmgroup.mm3``. 

Therefore several .c, .h, .pyd files are created, and also one .pyx
file is generated for each modulus ``p``. Then the Cython package  
is used for integrating all C functions into a python extension. 
Here we create one extension for each modulus ``p``. The action of 
each extension is controlled by a .pyx file, which also includes 
the generated .pyd files.

Function ``make_c_h_pxd(p)`` in this module generates the .c, .h
and .pxd files, and function  ``generate_pyx(p, ...)`` generates
the .pyx file for modulus ``p``. 
    
Variable INT_BITS in file config.py should be set to 32 or 64 for
32-bit or 64-bit target systems. There might be performace issues
if INT_BITS is not chosen properly. 


Generating .c files
-------------------

We use the C code generation mechanism in class 
``generate_c.TableGenerator``. There .c and .h files are generated
from files with extension .ske. These .ske files al like .c files, 
but augmented with some code generation statements for entering 
(usually rather large) tables and automatically generated code into 
the .c file to be generated. These .ske files may also have statements 
for automatically generating .h files declaring the exported funtions.  
    
Variable C_FILE_SKELETONS is a list of all .ske files. We will 
create  one  .c file from each .ske file in that list. The variable
C_FILE_SKELETONS is a list of strings where the substring
``'{P}'`` occurs in each string. From the string e.g.
``"mm{P}_op_pi"`` we derive the name ``"mm_op_pi.ske"`` of the
.ske file and the names ``"mm3_op_pi.c"``, ``"mm15_op_pi.c"``
of the C files generated for modulus ``3`` and ``5``.

Entering tables and automatically generated code
................................................

We create an instance ``tg`` of class ``TableGenerator`` for 
generating the .c files. For each modulus ``p`` a specific instance 
of class  ``TableGenerator`` is created. A table generator ``tg`` 
takes two dictionaries ``tables`` and ``directives`` as arguments. 
These dictionaries provide user-defined tables and directives for 
the code generator. There are several table-providing classes that
provide dictionaries ``tables`` and ``directives`` for the 
code genertor ``tg``. Variable ``TABLE_CLASSES`` contains 
the list of table-providing classes. For each of these
table-providing classes an instance is created, one instance for 
each modulus ``p``.  

Each table-providing class has methods ``tables()`` and 
``directives()`` that return the corresponding dictionaries.
The union of the dictionaries  ``tables()`` or ``directives()`` 
of all table-providing classes is passed as the constructor
of the instance ``tg`` as argument ``tables()`` or 
``directives()``, respectively.

A simple example for a table-generating class is class 
``Lsbit24Function`` in module ``mmgroups.dev.mat24.mat24aux``.



Generating the .h and the .pxd file
...................................

For each modulus ``p`` we create a single .h file containing the 
prototypes for the functions in all .c files generated for that
modulus. The name of that .h file is given by the variable 
``H_BASICS_NAME``. In case  ``H_FILE_NAME = "mm{P}_op.h"`` the
.h file ``"mm3_op.h"`` will be created for modulus ``p = 3``.

The content of a .h file is also  generated by the corresponding
instance ``tg`` of class  ``TableGenerator`` described above. 

Each .h file has a prefix given by the string ``H_FILE_BEGIN``
and suffix given by the string ``H_FILE_END``. 

For each .h file the Cython package requires a .pxd file for
generating a python extension. The .pxd file contains essentially 
the same information as the .h file. Method ``generate_pxd()``
of the code generator ``tg`` generates the .pxd file. The name
of the .pxd file is given by ``PXD_FILE_NAME``. The generated
.pxd file is prefixed with the content of the string
``PXD_DECLARATIONS``. In case ``PXD_FILE_NAME = "mm_op{P}.pxd"`` 
the .pxd file ``"mm3_op.pxd"`` will be created for  the
modulus ``p = 3``.


Generating the .pyx file
........................

The functions exported by the python extension ``mmgroup.mm{P}`` 
are just simple wrappers for the C functions generated by the code
generator. These wrappers must be coded in the Cython language.
They are simple enough so that they can be generated automatically
by function  ``pxd_to_pyx()`` in module ``mmgroup.generate_c``.
Function ``pxd_to_pyx()`` takes the .pxd files for one modulus 
``p`` (see last section) as input and it creates the .pyx file
for that modulus with name given by ``PXI_FILE_NAME``. in case
``PYX_FILE_NAME = "mm_op{P}.pyx"`` the file ``"mm_op3.pyx"``
for ``p = 3``.


Augmenting the .pyx file
........................

The .pyx file is augmented by some constants for its specific modulus
``p``. The list of constants is given by ``MM_BASICS_CONSTANTS``.
For each modulus ``p`` we create an instance of class ``MM_Basics`` in
module ``mmgroup.dev.mm_basics.mm_basics``. Here we pass the modulus
``p`` to the constructor of that class. Then the instance of class
``MM_Basics`` contains the required constants as attributes, and we
simply define these contants with its appropriate values in the .pyx
file. So contants with the same name may have different values in
extensions corresponding to different  moduli ``p``. 


 
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


import sys
import os
import re
import collections
import warnings
from numbers import Integral

from config import SRC_DIR, DEV_DIR,  C_DIR, DOC_DIR, PXD_DIR
from config import REAL_SRC_DIR 
sys.path.append(REAL_SRC_DIR)

from config import INT_BITS, PRIMES
SKE_DIR = os.path.join(DEV_DIR, "mm_op")



from mmgroup.generate_c import TableGenerator, make_doc, format_item
from mmgroup.generate_c import pxd_to_pyx, pxd_to_function_list

from mmgroup.dev.mm_basics import mm_aux, mm_tables, mm_basics
from mmgroup.dev.mm_basics import mm_tables_xi, mm_random

from mmgroup.dev.hadamard import hadamard_t 
from mmgroup.dev.hadamard import hadamard_xi 

from mmgroup.dev.mm_op import mm_op, mm_op_pi, mm_op_xy, mm_op_xi






VERBOSE = 0

##########################################################################
# Generating .c files
##########################################################################




H_FILE_SKELETONS = [
    "mm{P}_op_defines",
]

C_FILE_SKELETONS = [
    "mm{P}_op_pi",
    "mm{P}_op_misc",
    "mm{P}_op_xy",
    "mm{P}_op_t",
    "mm{P}_op_xi",
    "mm{P}_op_word",
]

##########################################################################
# Generating .h and .pxd files files
##########################################################################


H_FILE_NAME = "mm_op{P}.h"


PXD_FILE_NAME = "mm_op{P}.pxd"




### Generate the headers for specific p

H_FILE_BEGIN = """
// %%GEN h
// This header has been created automatically, do not edit!

#ifndef MM_OP{P}_H
#define MM_OP{P}_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mm_basics.h"

"""


H_FILE_END  = """
// %%GEN h
#ifdef __cplusplus
}
#endif
#endif  // #ifndef MM_OP{P}_H

"""


### Declarations for the generated .pxd files


PYX_DECLARATIONS = """
# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from cython.parallel cimport parallel
cimport openmp

from libc.stdint cimport uint32_t, uint16_t, uint8_t, int32_t
from libc.stdint cimport uint16_t
from libc.stdint cimport uint{INT_BITS}_t as uint_mmv_t


"""


##########################################################################
# Entering contants from class MM_Basics 
##########################################################################

# List of constants added to classes for a specific p from MM_Basics(p)
MM_BASICS_CONSTANTS = [
   "MMV_ENTRIES", "INT_BITS", "LOG_INT_BITS",
   "P", "FIELD_BITS", "LOG_FIELD_BITS", 
   "INT_FIELDS", "LOG_INT_FIELDS", "P_BITS", "MMV_INTS", 
]   

##########################################################################
# Entering tables and automatically generated code
##########################################################################




TABLE_CLASSES = [
    mm_op.MM_Op,
    mm_op_pi.Perm24_Benes,
    mm_op_pi.SmallPerm64,
    mm_op_pi.ScalarProd2048,    
    mm_op_xy.Perm64_xy, 
    hadamard_t.HadamardOpT64,
    hadamard_t.HadamardOpT3,
    hadamard_xi.HadamardOpXi64,
    hadamard_xi.HadamardOpXi16,
    mm_op_xi.MonomialOp_xi,
]



    
##########################################################################
# Generating the .pxd file
##########################################################################



PYX_FILE_NAME = "mm_op{P}.pyx"




##########################################################################
# Generating the .c, .h and .pyx file
##########################################################################



def make_c_h_pxd(p):
    """Generate .c files for a specific p from  .ske files 

    The relevant .ske files are listed in variable
    C_FILE_SKELETONS. One corresponding .h and .pxd file
    is also generated. The .c files and the .h file are written 
    to directory given by C_DIR, the .pxd file is written to the 
    current directory.

    Tables and directives for automatically generated C code
    are taken from the classes listed in TABLE_CLASSES.

    Return pair of lists, one of the .c files and one of the.
    .pxd files generated
    """
    tables = {}
    directives = {}
    global generated_tables
    for table_class in TABLE_CLASSES:
        table_instance = table_class(p)
        tables.update(table_instance.tables)
        directives.update(table_instance.directives)
    tg = TableGenerator(tables, directives, verbose = VERBOSE)
    # first generate C files
    c_files = []
    all_ske_files = [H_FILE_BEGIN]
    for name in H_FILE_SKELETONS:
        ske_file = name.format(P = "") + ".ske"
        ske_path = os.path.join(SKE_DIR, ske_file)
        all_ske_files.append(ske_path)
    for name in C_FILE_SKELETONS:
        ske_file = name.format(P = "") + ".ske"
        ske_path = os.path.join(SKE_DIR, ske_file)
        c_file = name.format(P = p) + ".c"
        c_path= os.path.join(C_DIR, c_file)
        print("Creating %s from %s" % (c_file, ske_file))
        tg.generate(ske_path, c_path)
        all_ske_files.append(ske_path)
        c_files.append(c_path)
    all_ske_files.append (H_FILE_END)

    # generate .h and .pxd file
    h_file =  H_FILE_NAME.format(P = p)
    h_path =  os.path.join(C_DIR, h_file)
    pxd_file =  PXD_FILE_NAME.format(P = p)
    print("Creating %s from previous .ske files" % h_file)
    tg.generate(all_ske_files, None, h_path)
    pyx_declarations = PYX_DECLARATIONS.format(
        INT_BITS = INT_BITS, P = p
    )
    tg.generate_pxd(
       os.path.join(PXD_DIR, pxd_file), 
       h_file, pyx_declarations, nogil = True
    )
    return c_files,  [ pxd_file ]







##########################################################################
# Generating the .pyx file
##########################################################################




m_translate_pxd = re.compile("mm_op(\d+)_(\w+)")



def translate_pxd(name):
    m = m_translate_pxd.match(name)
    if (m):
        p, method_name = m.groups()
        return  "op_"  + method_name 
    else:
        return name



def generate_pyx(p, pxd_files):
    """The main function of this module for generating code

    This function generates the .pyx file with name given by PYX_FILE_NAME 
    and stores it in the current directory. Cython will use the .pyx file 
    to build a wrapper for all generated C functions. 
    """
    def pyx_comment(text, f):
        print("\n" + "#"*70 + "\n### %s\n" % text + "#"*70 + "\n\n",file=f)

    f_pyx = open(os.path.join(PXD_DIR, PYX_FILE_NAME.format(P = p)), "wt")
    declarations = PYX_DECLARATIONS.format(INT_BITS = INT_BITS, P = 0)
    print(declarations, file = f_pyx)

    for pxd_f in pxd_files:
        pyx_comment("Wrappers for C functions from file %s" % pxd_f, f_pyx)
        pyx_code = pxd_to_pyx(
            os.path.join(PXD_DIR, pxd_f), 
            None, translate_pxd, nogil=True
        )
        print(pyx_code, file = f_pyx)

    pyx_comment("Constants", f_pyx)
    basis = mm_op.MM_Op(p)
    for name in MM_BASICS_CONSTANTS:
        value = getattr(basis, name)
        declaration = "%s = %d\n" % (name, value) 
        print(declaration,  file = f_pyx)

    f_pyx.close()




##########################################################################
# Exported functions for generating code
##########################################################################






if __name__ == "__main__":
    for p in PRIMES:
        c_files,  pxd_files = make_c_h_pxd(p)
        generate_pyx(p, pxd_files)
