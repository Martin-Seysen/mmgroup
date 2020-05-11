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
from mmgroup.dev.mat24.make_addition_table import make_addition_tree
from mmgroup.dev.mat24.make_addition_table import show_addition_tree
from mmgroup.dev.mat24.make_addition_table import BitMatrixAddChain


@pytest.mark.mat24
def test_make_additon_table(verbose =  0):

    print("\nTesting module make_addition_table")
    b = Mat24.basis[13:]

    data = [ 6, 10, 12, 15]
    d, stage = make_addition_tree(data, granularity = 4)  
  
    d, stage = make_addition_tree(b, granularity = 4) 
    if verbose: 
       show_addition_tree(d)


    bm = BitMatrixAddChain(b, granularity = 4)
    bm.selftest()
    bm.set_names("uint32_t", "c", "a")
    if verbose:
        print( bm.make_c_program() )
        print( "===============================================" )
        print( bm.n_registers , "registers needed" )
        print( len(bm.ops), "operations needed" )
    print("Test passed")

    
    
     
    