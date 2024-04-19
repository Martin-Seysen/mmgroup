from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
import re
import warnings

from mmgroup import mat24
from mmgroup.dev.mm_op.mm_op import MM_Op


###########################################################################
# Table for operation xy on the 64 entries of for tag T
###########################################################################

class Perm64_xy(MM_Op):
    """Yet to be documented!!!!!!!!!!!!!!!!!!!!!

    """
    weights = sum(mat24.suboctad_weight(j) << j for j in range(64))
    mask = (1 << 64) - 1
    hi_list = [ 0, mask, weights, weights ^ mask ]
    directives = {}

    def __init__(self, **kwds):
        """Initialise for calulations with small integers modulo p

        p+1 must be a power of two. Calculations modulo p are described 
        in more detail in the base classes of this class.
        """
        p = int(kwds.get('p', 3))
        super(Perm64_xy, self).__init__(p = p)
        self.t_hi = self.hi_table()
        self.t_lo = self.lo_table()
        self.tables.update(self.make_tables())

    def table_value(self, index):
        i0 = index & 0x3f
        v = sum(mat24.suboctad_scalar_prod(j, i0) << j for j in range(64))
        return self.hi_list[ (index >> 6) & 3 ] ^ v

    def table_part(self, index):
        tbl = self.table_value(index)
        p, i_fields = self.P, self.INT_FIELDS
        return [self.smask(p, tbl >> i) for i in range(0, 64, i_fields)]

    def hi_table(self):
        return sum((self.table_part(i) for i in range(0, 256, 16)), [])
        
    def lo_table(self):
        return sum((self.table_part(i) for i in range(16)), [])

    def make_tables(self):
        return {
            "TABLE_PERM64_XY_LOW" : self.t_lo,
            "TABLE_PERM64_XY_HIGH" : self.t_hi,
        }
 


Tables = Perm64_xy

class MockupTables:
    tables, directives = {}, {}
    def __init__(self, **kwds):
        pass
