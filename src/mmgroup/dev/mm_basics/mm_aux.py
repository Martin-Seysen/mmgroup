"""Accessing vectors of the representation of the monster

We deal with the rational representation [Seysen20]  of the monster 
group MM, which is  based on the Conway's construction [Con85] of
the monster, modulo various small integers p = 2**n-1, 2 <= n <= 8.
Here the integer p is called the modulus. 

All representations R_p of MM are 196884-dimensional matrix 
representations over the integers modulo p. Here the matrices
operate on vectors which are represented as arrays of integers.
The matrix operations on these vectors are performed by C functions.

For vectors there is an external representation as an array of 
integers of type uint8_t and there is also an internal representation. 
The internal representation is fairly close to the description of the 
basis vectors of the representation in [Seysen20] , so that matrix 
operations can take advantage of the internal structure of that 
repesentation.

The purpose of this module is access to the internal representation 
of such a vector and conversion between the external and the 
internal representation.  

There are versions optimized for 32-bit and for 64-bit computers. 
This version is optimized for {INT_BITS} bits. It supports 
moduli p = {P_LIST}.


References
----------
see file mmgroup.pdf
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
import collections
"""TODO: write module documentation

"""


import warnings
from numbers import Integral
import numpy as np


from mmgroup.dev.mm_basics import mm_basics   


from mmgroup.dev.mm_basics.mm_basics import MM_Basics, MM_Const
from mmgroup.generate_c import c_snippet, TableGenerator
from mmgroup.generate_c import format_item

from mmgroup.bitfunctions import bitlen, hibit

from mmgroup.dev.mm_basics.mm_doc import __doc__ as MM_DOC
from mmgroup.dev.mm_basics.mm_basics import __doc__ as MM_BASICS_DOC


INT_BITS =  MM_Basics.INT_BITS 




#2345678901234567890123456789012345678901234567890123456789012345678901

class MM_IO24:
    """Generate code to access vectors of R_p in internal representation

    This class provides code generation functions for fast data
    transfer between a small array in a vector in R_p and an
    array of bytes. 

    We use the code generation mechanism in class  
    make_c_tabels.TableGenerator. The generated functions take
    an input parameter 'log_field_bits' coresponding to
    LOG_FIELD_BITS. This indicates the a bit field containing an
    entry contains 2**log_field_bits.
    """
    directives = {}

    @staticmethod
    def iter_ext_ABC():
        """Yield triples (tag, i, j) for external rep for tags A, B, C

        Triples are given in the order of the external representation
        for tags A, B and C.
        tag = 0, 1, 2 means 'A', 'B', 'C', respectively
        """
        for i in range(24):
            yield 0, i, i
        for tag in range(3):
            for i in range(24):
                for j in range(i):
                     yield tag, i, j

    @classmethod
    def abc_table(cls):
        """Table for expanding entries for tags 'A', 'B', 'C'.

        Entry k0 of the external representation of the monster is
        mapped to location k1 in the internal representation with
        k1 = (Table[k0] & 0x7ff) + k0 - 24. Entry k0 is also copied 
        to location k1 - 31 * (Table[k0] >> 11) of the internal rep.

        See comments 'Internal representation of a vector' and
        'External representation of a vector' for internal and
        external representation of a vector.
        """
        a = np.zeros(24 + 3 * 276, dtype = np.uint16)
        for k0, (tag, i, j) in enumerate(cls.iter_ext_ABC()):
            k1 = (24 * tag + i) * 32 + j
            k1_entry = k1 - k0 + 24
            d_t  = i - j
            assert 0 <= k1_entry <  2048
            a[k0] = k1_entry  + (d_t << 11)
        return a
               


    @classmethod
    def reduce_table(cls):
        """Masks for reducing the fields of an uint_mmv_t modulo p

        Table entries 2*i-2 and 2*i-1 refer to modulus p = 2**i-1 for 
        i = 2,...,8. In the standard case, when i is not a power of 
        two, entries 2*i-2 and  2*i-1 have the following values:

        Index  Table entry
        2*i-2: A mask containg the value 1 in each field. 
        2*i-1: A mask containg the value 2**i-1 in each field. 

        If i is a power of two then table entry 2*i-1 is mask containing 
        the  value 2**(i/2)-1 in each field.

        The reason for that special case is that there is no space for 
        a carry bit between two adjacent fields if i is a power of two.
        In that case we need a special trick for detecting the value 
        2**i-i in a field.        
        """ 
        tbl = []
        for p_bits in range(2, 9):
            field_bits = p_bits
            while (field_bits & (field_bits - 1)): 
                field_bits += 1 
            mask_1 = sum(1 << i for i in range(0, INT_BITS, field_bits))
            q = p_bits >> (field_bits == p_bits)
            mask_q = ((1 << q) - 1) * mask_1
            tbl += [mask_1, mask_q]
        return tbl      
                    
    def __init__(self):
        """For code generation with class make_c_tables.TableGenerator

        """
        self.tables = {
            "MM_DOC": MM_DOC,
            "MM_BASICS_DOC": MM_BASICS_DOC,
            "MM_AUX": sys.modules[__name__],
            "MM_AUX_IO24": self,
            "MM_AUX_TBL_ABC" : self.abc_table(), 
            "MM_AUX_TBL_REDUCE": self.reduce_table(),
        }
        

class Tables(MM_IO24):
    def __init__(self, *args, **kwds):
        super().__init__()
        pass





