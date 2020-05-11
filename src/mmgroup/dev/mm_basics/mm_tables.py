"""Generating tables for operations on the rep of the monster

This includes also somw functions that use these table independent
of characteristic p
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals




import os
import sys
import numpy as np


from mmgroup import mat24 
from mmgroup.bitfunctions import bitweight, bitparity

from mmgroup.dev.mm_basics.mm_basics import MM_Const


def make_octads_table():
    octads = [mat24.octad_to_vect(i) for i in range(759)] 
    bl = [mat24.vect_to_bit_list(v1)[0] for v1 in octads]
    return np.array(sum(bl, []), dtype = np.uint8)


def make_perm64_table():
    step = [0, 1, 2, 4, 8, 16, 32, 63]
    return [step[i] ^ step[j] for j in range(8) for i in range(8)] 
  
def make_inv_table():
    return [
       sum((mat24.gcode_weight(32*i + j) & 1) << j for j in range(32))
       for i in range(64)
    ]

def make_perm64_weights():
    return sum((((bitweight(i) + 1) >> 1) & 1) << i for i in range(64))


def make_suboctad_table():
    return [((i >> 1) ^ -bitparity(i)) & 0x3f for i in range(128)]


#234567890123456789012345678901234567890123456789012345678901234567890  

class MM_OctadTable(MM_Const):    
    octad_elements = make_octads_table()
    perm64_table = make_perm64_table()
    inv_table = make_inv_table()
    perm64_weights = make_perm64_weights() 
    suboctad_table = make_suboctad_table()

    new_tables = {
        "MM_INV_2048_TABLE": inv_table,
        "MM_OctadElementsTable" :  octad_elements,
        "MM_PERM64_TABLE": perm64_table,
        "MM_SUBOCTAD_TABLE": suboctad_table,
    }

    def __init__(self):
        super(MM_Const, self).__init__()
        self.make_tables(self.new_tables, {})



