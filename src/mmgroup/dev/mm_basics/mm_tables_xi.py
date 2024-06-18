r"""Generate tables for operation xi on the rep 196884x of the monster

Class ``MM_TablesXi`` in this module contains the arrays PERM_TABLES,
SIGN_TABLES, OFFSETS, SHAPES, as secribed in the
**mmgroup guide for developpers**, section
**Some mathematical aspects of the implementation**, subsection
**Implementing generators of the Monster group**, subsubsection
**Monomial operation of the generators \xi^e**.
We adopt the terminology from that subsection.

These arrays are computed by the methods of class ``Pre_MM_TablesXi``
in this module. Her we use teh implementation of the operation of
generator :math:`\xi^{\pm ^}` in module ``mmgroup.generators``.
"""

#234567890123456789012345678901234567890123456789012345678901234567890

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import os
import sys
import numpy as np
from itertools import product

from mmgroup.generate_c import UserDirective, UserFormat
from mmgroup.generate_c import ConstUserFormat


try:
    # Try importing the fast C function
    from mmgroup import generators as gen 
except (ImportError, ModuleNotFoundError):
    # Use the slow python function if the C function is not available
    print("\nUsing slow Python functions for table generation!!!\n")
    from mmgroup.dev.generators.gen_xi_ref import GenXi
    gen = GenXi




def make_table_bc_symmetric(table):
    def b(i, j):
        return 32 * i + j
    def c(i, j):
        return 32 * (i + 24) + j
    for i in range(24):
        for j in range(i):
            table[b(j,i)] = table[b(i,j)] 
            table[c(j,i)] = table[c(i,j)] 
        table[b(i,i)] = b(i,i)
        table[c(i,i)] = c(i,i)


def cut24(table):
    table = table.reshape(-1, 32)
    #assert sum(table[:,24:].reshape(-1)) == 0, (table.shape, table[:,24:])
    table = table[:,:24]
    table = table.reshape(-1)
    return table



def check_table(table, blocks, row_length):
    assert row_length in (24, 32)
    if row_length == 24:
        table = cut24(np.copy(table))
    length = len(table)
    image_length = (int(max(table & 0x7fff)) + 31) & -32
    blocklen, r0 = divmod(length, blocks)
    image_blocklen, r1 = divmod(image_length, blocks)
    assert r0 == r1 == 0
    for i, start in enumerate(range(0, length, blocklen)):
         part = table[start:start+blocklen] & 0x7fff
         t_min, t_max = i*image_blocklen, (i+1)*image_blocklen
         assert t_min <= min(part) <= max(part) < t_max



def table_immage_box(u_box, u_exp):
    return gen.gen_xi_op_xi_short(u_box << 16, u_exp) >> 16


class Pre_MM_TablesXi: 
    def __init__(self):
        """
        BC, T0 = 24*32, 72*32 + 15*64
        T1, X0  = 72*32 + 375*64,  72*32 + 759*64
        X1 = X0 + 1024*32
        """
        BC, T0, T1, X0, X1 = 1, 2, 3, 4, 5
        _OFS_X0 = 72*32 + 759*64
        OFFSETS = [None, 24*32, 72*32 + 15*64, 72*32 + 375*64, 
                   _OFS_X0, _OFS_X0 + 1024*32]
        self.TAG_NAMES = {BC:"BC", T0:"T0", T1:"T1", X0:"X0", X1:"X1"} 

        self.BOX_SHAPES = {
            BC: (1, 78, 32), T0: (45, 16, 32), T1: (64, 12, 32),
            X0: (64, 16, 24), X1: (64, 16, 24)
        }

        self.MAP_XI = [
            [[BC, BC], [BC, BC]],
            [[T0, T0], [T0, T0]],
            [[T1, X0], [T1, X1]],
            [[X0, X1], [X1, X0]],
            [[X1, T1], [X0, T1]]
        ]   

        for n in range(5):
            _m = self.MAP_XI[n]
            for j in range(2):
                assert self.BOX_SHAPES[_m[0][j]] == self.BOX_SHAPES[_m[1][j]]
            for exp1 in range(2):
                box = _m[exp1][0]
                img = gen.gen_xi_op_xi_short(box << 16, exp1 + 1) >> 16
                assert _m[exp1][1] == img

        self.PERM_TABLES = [[0,0], [0,0], [0,0], [0,0], [0,0]]
        self.SIGN_TABLES = [[0,0], [0,0], [0,0], [0,0], [0,0]]
        self.TABLE_BOX_NAMES = {} # items are pairs of source and destination box
        self.OFFSETS = np.zeros((5,2,2), dtype = np.uint32).tolist()
        self.SHAPES =  np.zeros((5,2, 3), dtype = np.uint32).tolist()

        for n in range(5):
            box = self.MAP_XI[n][0][0]
            img = self.MAP_XI[n][0][1]
            self.SHAPES[n][0] = self.BOX_SHAPES[box]
            self.SHAPES[n][1] = self.BOX_SHAPES[img]
            #print(n, box, img, self.SHAPES[n])
        

        for n in range(5):
            for exp1 in range(2):
                box = self.MAP_XI[n][exp1][0]
                img = self.MAP_XI[n][exp1][1]
                assert self.SHAPES[n][0] == self.BOX_SHAPES[box]
                assert self.SHAPES[n][1] == self.BOX_SHAPES[img]
                table = gen.make_table(box, exp1 + 1)
                shape =  self.SHAPES[n][0]
                img_shape = self.SHAPES[n][1]
                #print(shape, img_shape)
                assert shape[0] == img_shape[0], (shape[0], img_shape[0]) 
                assert len(table) == shape[0] * shape[1] * 32
                img_len = img_shape[0] * img_shape[1] * 32
                check_table(table, shape[0], shape[2])
                inv_table = gen.invert_table(table, shape[2], img_len)
                if box == BC:
                    make_table_bc_symmetric(inv_table)
                t_perm, t_sign = gen.split_table(inv_table, shape[1]*32)
                del table
                del inv_table
                if img_shape[2] == 24:
                    t_perm = cut24(t_perm)
                self.PERM_TABLES[n][exp1] = t_perm
                self.SIGN_TABLES[n][exp1] = t_sign
                self.TABLE_BOX_NAMES[n, exp1] = [self.TAG_NAMES[box],
                                                  self.TAG_NAMES[img]]
                self.OFFSETS[n][exp1][0] = OFFSETS[box]
                self.OFFSETS[n][exp1][1] = OFFSETS[img]




class MM_TablesXi:
    done_ = False
    directives = {}
    tables = {}
    def __init__(self):
        cls = self.__class__
        if cls.done_:
            return
        Pre = Pre_MM_TablesXi()

        cls.PERM_TABLES = Pre.PERM_TABLES
        cls.SIGN_TABLES = Pre.SIGN_TABLES
        cls.TABLE_BOX_NAMES = Pre.TABLE_BOX_NAMES
        cls.OFFSETS = Pre.OFFSETS
        cls.SHAPES = Pre.SHAPES

        cls.tables["MM_TABLE_PERM_XI"] = cls.PERM_TABLES
        cls.tables["MM_TABLE_SIGN_XI"] = cls.SIGN_TABLES
        cls.tables["MM_TABLE_OFFSETS_XI"] = cls.OFFSETS
        cls.tables["MM_TABLE_XI_COMMENT"] = UserFormat(
             cls.comment_table_mapping, "ii")
        cls.done_ = True


    @classmethod
    def comment_table_mapping(cls, i, j):
         s = "table for xi**%d: map tag %s to tag %s"
         src, dest = cls.TABLE_BOX_NAMES[i, j]
         return s % (j + 1, src, dest)

    @classmethod
    def display_config(cls):
        cls.__init__()
        print("source:")
        print(" shapes:", cls.SHAPES)
        print(" start: ", cls.OFFSETS)
        print(" boxes: ", cls.TABLE_BOX_NAMES)


class Mockup_MM_TablesXi:
    directives = {}
    tables = {}


class Tables:
     mockup_tables = Mockup_MM_TablesXi.tables
     mockup_directives = Mockup_MM_TablesXi.directives

     def __init__(self, *args, **kwds):
         self._table_class = None

     def _load_tables(self):
         if self._table_class is None:
             self._table_class =  MM_TablesXi()
         return self._table_class

     @property
     def tables(self):
         return self._load_tables().tables

     @property
     def directives(self):
         return self._load_tables().directives


