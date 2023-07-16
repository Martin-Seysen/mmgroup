
from __future__ import absolute_import, division, print_function
#from __future__ import  unicode_literals



import sys
import os
import itertools
from collections import defaultdict, Counter
from numbers import Integral
import numpy as np
import subprocess
import shutil

import pytest
 

from mmgroup.dev.mm_basics.mm_basics import INT_BITS
from mmgroup.dev.hadamard.hadamard_codegen import HadamardMatrixCode
from mmgroup.dev.hadamard.hadamard_t import HadamardOpT64
from mmgroup.dev.hadamard.hadamard_t import HadamardOpT3
from mmgroup.dev.hadamard.hadamard_xi import HadamardOpXi64
from mmgroup.dev.hadamard.hadamard_xi import HadamardOpXi16

from mmgroup.tests import test_hadamard
from mmgroup.tests.tempdir import make_temp_dir, kill_temp_dir
from mmgroup.tests.tempdir import temp_dir_name


FILE_NAME = "test_hadamard_op"
C_FILE_NAME = FILE_NAME + ".c"
H_FILE_NAME = FILE_NAME + ".h"
PXD_FILE_NAME = FILE_NAME + ".pxd"
PXD_MODULE_NAME =  "test_hadamard_op"
PYX_FILE_NAME =    "test_hadamard.pyx"

################################################################
# Defining the names of the test functions
################################################################

SIZES = [2, 4, 6]
PRIMES = [3, 7, 15,127, 255]

def hadamard_function_name(p, n_qubits):
    return "mod%d_hadamard%d_test" % (p, 1 << n_qubits)

def xch_parity_function_name(p):
    return "mod%d_xch_parity64_test" % p

def op_t64_function_name(p):
    return "mod%d_op_t64_test" % p
  
def op_t3_function_name(p):
    return "mod%d_op_t3_test" % p

def op_xi64_function_name(p):
    return "mod%d_op_xi64_test" % p

def op_xi16_function_name(p):
    return "mod%d_op_xi16_test" % p

def iter_function_names():
    for q in SIZES:
        for p in PRIMES:
             yield hadamard_function_name(p, q)
    for p in PRIMES:
         yield  xch_parity_function_name(p)
    for p in PRIMES:
         yield  op_t64_function_name(p)
    for p in PRIMES:
         yield  op_t3_function_name(p)
    for p in PRIMES:
         yield  op_xi64_function_name(p)
    for p in PRIMES:
         yield  op_xi16_function_name(p)



################################################################
# Coding test routines in C file
################################################################


def code_hadamard_test(p, n_qubits):
    array_name = "p_data"
    function_name = hadamard_function_name(p, n_qubits)
    s =  """

void %s(uint_mmv_t *%s)
{
""" % (function_name, array_name) 

    cc =  HadamardMatrixCode(p, n_qubits)
    cc.load_vector_direct(array_name)
    cc.hadamard_op()
    cc.mul_pwr2(-n_qubits >> 1)
    cc.store_vector_direct(array_name)
    cc.comment_statistics()
    return  s + cc.generate() +  "}\n"


def code_xch_parity_test(p):
    array_name = "p_data"
    function_name = xch_parity_function_name(p)
    s =  """

void %s(uint_mmv_t *%s)
{
""" % (function_name, array_name) 

    cc = HadamardOpT64(p, 6)
    cc.load_vector_direct(array_name)
    cc.swap_parities()
    cc.store_vector_direct(array_name)
    cc.comment_statistics()
    return  s + cc.generate() +  "}\n"



def code_op_t64_test(p):
    array_name = "p_data"
    function_name = op_t64_function_name(p)
    cc = HadamardOpT64(p, 6)
    s =  """

void %s(uint_mmv_t *%s)
{
uint_mmv_t exp1 = -(%s[%d] & 1ULL);
""" % (function_name, array_name, array_name, cc.V64_INTS) 

    s += cc.make_code(array_name, "exp1", array_name)
    return  s +   "}\n"



def code_op_t3_test(p):
    cc = HadamardOpT3(p)
    d = {
        "array_name" : "p_data",
        "function_name" : op_t3_function_name(p),
        "a_len" : 24 * cc.V24_INTS,
    }
    s =  """

void {function_name}(uint_mmv_t *{array_name})
""".format(**d)  # 1
    s += "{\n"
    s += """
uint_mmv_t exp1 = -({array_name}[3] & 1ULL);
uint_mmv_t a[3 * {a_len}];
a[0] = {array_name}[0];
a[{a_len}] = {array_name}[1];
a[2*{a_len}] = {array_name}[2];
""".format(**d)  # 2
    s += cc.make_code("a", "exp1", "a")
    s += """
{array_name}[0] = a[0];
{array_name}[1] = a[{a_len}];
{array_name}[2] = a[2*{a_len}];
""".format(**d) 
    return  s +   "}\n"






def code_op_xi64_test(p):
    array_name = "p_data"
    function_name = op_xi64_function_name(p)
    cc = HadamardOpXi64(p)
    data = ",\n".join(map(cc.hex, cc.exp_table() + [0] * 8))
    s =  """

void %s(uint_mmv_t *%s)
{
static uint_mmv_t mask0[16] = {
   %s
}; 
int index = (%s[%d]) %% 3;
uint_mmv_t *p_mask = mask0 + (index == 2 ? 8 : index);
""" % (function_name, array_name, data, array_name, 16*cc.V24_INTS) 
    #s += 'printf("OFS = %d\\n", (int)(index));\n' 
    s += cc.make_code(array_name, "p_mask", array_name)
    return  s +   "}\n"




def code_op_xi16_test(p):
    array_name = "p_data"
    function_name = op_xi16_function_name(p)
    cc = HadamardOpXi16(p)
    s =  """

void %s(uint_mmv_t *%s)
{
uint_mmv_t exp1 = -(%s[%d] & 1ULL);
""" % (function_name, array_name, array_name, 4*cc.V24_INTS) 

    s += cc.make_code(array_name, "exp1", array_name)
    return  s +   "}\n"


################################################################
# Creating the .c, .h, .pyd, and .pyx file
################################################################


def generate_test_c_file():
    with open(C_FILE_NAME, "wt") as f:
        f.write("""
#include <stdint.h>
typedef uint%d_t uint_mmv_t;
#include "%s"

"""  % (INT_BITS, H_FILE_NAME))
        for q in SIZES:
            for p in PRIMES:
                 f.write(code_hadamard_test(p, q))
        for p in PRIMES:
            f.write(code_xch_parity_test(p))
        for p in PRIMES:
            f.write(code_op_t64_test(p))
        for p in PRIMES:
            f.write(code_op_t3_test(p))
        for p in PRIMES:
            f.write(code_op_xi64_test(p))
        for p in PRIMES:
            f.write(code_op_xi16_test(p))

def generate_test_h_file():
    with open(H_FILE_NAME, "wt") as f:
        f.write("""
typedef uint%d_t uint_mmv_t;

"""  % INT_BITS)
        for f_name in iter_function_names():
            f.write("void %s(uint_mmv_t*);\n" % f_name)





def generate_test_pxd_file():
    with open(PXD_FILE_NAME, "wt") as f:
        f.write("""# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from libc.stdint cimport uint%d_t as uint_mmv_t

cdef extern from "%s":

"""  % (INT_BITS, "test_hadamard_op.h" ))
        for f_name in iter_function_names():
            f.write("    void %s(uint_mmv_t*)\n" % f_name)



def generate_test_pyx_file():
    with open(PYX_FILE_NAME, "wt") as f:
        f.write("""# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from libc.stdint cimport uint%d_t as uint_mmv_t

from mmgroup.tests.test_hadamard.temp cimport %s


"""  % (INT_BITS, PXD_MODULE_NAME ))
        for f_name in iter_function_names():
            f.write("""
def %s(v):
    cdef uint_mmv_t[::1] pv = v
    %s.%s(&pv[0])

""" % (f_name, PXD_MODULE_NAME, f_name))





################################################################
# Creating the python library with cython
################################################################




def cythonize():
        """Build the python extension Hadamard64Sub

        This is one of the few places where OS-specific code is required
        """
        if os.name in ["nt"]:
            compilers = ["msvc", ""]
        else:
            compilers = [""]
        for compiler in compilers:
            try:
                command = [sys.executable, 
                   os.path.join("..", "setup_test_hadamard.py"),
                     "build_ext", "--inplace",
                ]
                if compiler: 
                    command.append("-c" + compiler)
                if subprocess.call(command) == 0:
                    return
            except FileNotFoundError:
                pass
        raise ValueError("setup has failed")

def cleanup():
    for fname in [
             C_FILE_NAME, H_FILE_NAME, PXD_FILE_NAME, PYX_FILE_NAME,
             "test_hadamard.c",
        ]: 
        try:
            os.remove(fname)
        except:
            pass
    try:
        shutil.rmtree("build")
    except:
        pass
        

################################################################
# Main program
################################################################


def build():
    old_dir = os.getcwd()
    temp_dir = make_temp_dir(test_hadamard)
    os.chdir(temp_dir)
    generate_test_c_file()
    generate_test_h_file()
    generate_test_pxd_file()
    generate_test_pyx_file()
    cythonize()
    cleanup()
    os.chdir(old_dir)


def kill():
    kill_temp_dir(temp_dir_name(test_hadamard))

 