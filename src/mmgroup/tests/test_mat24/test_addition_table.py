from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint
import types
import sys
import re
import os
from operator import __or__, __xor__
from functools import reduce
import pytest


from mmgroup.bitfunctions import bitlen, bitweight, hibit, v2


        
                        
 
 
from mmgroup.dev.mat24.mat24_ref import Mat24
from mmgroup.dev.mat24.make_addition_table import make_addition_operations
from mmgroup.dev.mat24.make_addition_table import check_addition_operations
from mmgroup.dev.mat24.make_addition_table import display_addition_operations
from mmgroup.dev.mat24.make_addition_table import BitMatrixMulFix


def table_testcases():
    yield Mat24.basis[13:], False, 12
    yield [0], True, 4
    yield [0], False, 4
    for i in range(1, 10):
         mask =  (1 << i) - 1
         for j in range(1, 10):
             m = [randint(0,  mask) for k in range(j)]
             yield m, randint(0,1), randint(4,6)


@pytest.mark.mat24
def test_make_additon_table(verbose =  0):
    print("\nTesting module make_addition_table")
    for i, (m, singleton, granularity)  in enumerate(table_testcases()):
        if verbose:
            nrows, ncols = len(m), bitlen(reduce(__or__, m, 0))
            print("\nTest case %d, a %d time %d matrix, singleton = %s" %
                (i, nrows, ncols, bool(singleton)))
            print("m =", [hex(x) for x in m])
        ops = make_addition_operations(m, singleton, granularity)
        if verbose:
            display_addition_operation(ops) 
        check_addition_operations(m, ops, singleton)

    """
    data = [ 6, 10, 12, 15]
    d, stage = make_addition_tree(data, granularity = 4)  
  
    d, stage = make_addition_tree(b, granularity = 4) 
    if verbose: 
       show_addition_tree(d)


    bm = BitMatrixMulFix(granularity = 4)
    bm.set_fixed_matrix(b)
    bm.selftest()
    bm.set_names("uint32_t", "c", "a")
    if verbose:
        print( bm.make_c_program() )
        print( "===============================================" )
        print( bm.n_registers , "registers needed" )
        print( len(bm.ops), "operations needed" )
    print("Test passed")
    """
    
    
     
    