"""Generating tables for operations on the rep of the monster

This includes also some functions that use these table independent
of characteristic p
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals




import os
import sys
import numpy as np


from mmgroup.bitfunctions import bitweight, bitparity

from mmgroup.dev.mm_basics.mm_basics import MM_Const


def make_perm64_table():
    """Make table for conversion from bit vector to suboctad

    Let o be an octad. Let {{i},{j}}, 0 <= i, j < 8, i != j be the 
    even cocode word which is a subset of octad  o,  such that 
    precisely the i-th and  the j-th bit  of octad  o  is set, as 
    in module mm_aux.c. Then  entry  8 * i + j  in the returned
    table is the number of the subocatad  {{i},{j}} in the standard 
    numbering given in module mm_aux.c.
    """
    step = [0, 1, 2, 4, 8, 16, 32, 63]
    return [step[i] ^ step[j] for j in range(8) for i in range(8)] 
  

def make_perm64_weights():
    """Return a 64-bit integer containing halved suboctad weights

    Bit i of the returned integer is set if the suboctade with
    number i has bit weight 2 (mod 4) and cleared if that
    suboctad has bit weight 0 (mod 4).    
    """
    return sum((((bitweight(i) + 1) >> 1) & 1) << i for i in range(64))


def make_suboctad_table():
    """Make table for conversion from suboctad to bit vector

    Let o be an octad. Let b_0,...,b_7 be the elements of octad 
    o in natural order. Let x = sum 2**i * x_i, i = 0,...,7.
    Then sum (x_i * b_i) is a suboctad as defined in function
    mat24_suboctad_to_cocode() in module mat24_functions.c.  
    In that function we also define a natural numbering of the
    suboctads of an octad. Entry x of following table is the
    number of the suboctad sum (x_i * b_i) for 0 <= x < 128.
    Note that bit b_7 can be computed from bits b_0,...,b_6
    as the even parity bit.
    """
    return [((i >> 1) ^ -bitparity(i)) & 0x3f for i in range(128)]



class MM_OctadTable:    
    perm64_table = make_perm64_table()
    perm64_weights = make_perm64_weights() 
    suboctad_table = make_suboctad_table()

    tables = {
        "MM_PERM64_TABLE": perm64_table,
        "MM_SUBOCTAD_TABLE": suboctad_table,
    }

    directives = {}



class Tables(MM_OctadTable):
    def __init__(self, *args, **kwds):
        pass



    
    
 


