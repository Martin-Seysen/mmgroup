from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import collections
import re
import warnings
from numbers import Integral
import numpy as np


from mmgroup.bitfunctions import isinteger, v2, hibit

from mmgroup.dev.generate_c.make_c_tables import c_snippet, TableGenerator
from mmgroup.dev.generate_c.make_c_tables import make_table
from mmgroup.dev.generate_c.generate_functions import UserDirective
from mmgroup.dev.generate_c.generate_functions import UserFormat

from mmgroup.dev.mm_basics.mm_tables_xi import MM_TablesXi
from mmgroup.dev.mm_basics.mm_basics import MM_Basics, MM_Op

class MonomialOpTableInfo:
    def __init__(self, shape, start, op_diff):
        self.SHAPE = shape
        self.START = start 
        self.OP_DIFF = op_diff


class MonomialOp_t(MM_Op):
    def __init__(self, p):
        super(MonomialOp_t, self).__init__(p)
        self.make_table_info()

    def load_bytes(self, source, bytes, n_bytes):
        s = ""
        n_registers = 8 // self.FIELD_BITS
        registers =  ", ".join(["r%d" % i for i in range(n_registers)])
        s = "uint_mmv_t %s;\n"% registers
        LOAD_MASK = self.hex(self.smask(self.P, -1, 8))
        for i, start in enumerate(range(0, n_bytes, self.INT_FIELDS)):
            s += "r0 =  (%s)[%d];\n" % (source, i)
            for j in range(1, n_registers):
                s += "r%d = (r0 >> %d) & %s;\n" % (
                   j, j * self.FIELD_BITS, LOAD_MASK)
            s += "r0 &= %s;\n" % LOAD_MASK
            for j in range(self.INT_FIELDS):
                if start + j < n_bytes:
                    s += "(%s)[%d] = r%d >> %d;\n" % (bytes,  start + j, 
                        j % n_registers, 8 * (j // n_registers))
        return s;


    def store_bytes(self, bytes, perm, sign, dest, n_bytes):
        s = "uint_mmv_t r0, r1;\n"

        d = ["((uint_mmv_t)(%s[%s[%d]]) << %d)" % (bytes, perm, i, 
                  self.FIELD_BITS * (i % self.INT_FIELDS)) 
                  for i in range(n_bytes)
            ]
        for i, start in enumerate(range(0, n_bytes, self.INT_FIELDS)):
            s += "r0 = " 
            s += "\n  + ".join(d[start:start+self.INT_FIELDS])
            s += ";\n" 
            s += "r1 = (%s)[0] >> %s;\n" % (sign, start)
            s += self.gen_mmv_uint_spread("r1", "r1")
            s += "(%s)[%d] = r0 ^ r1;\n" % (dest, i)



        #s +=("// i=%d\n" % i)
        for j in range(i+1, 32 // self.INT_FIELDS):
            s += "%s[%d] = 0;\n" % (dest,j)
        return s


    def comment(self,i):
        src_t = MM_TablesXi.SOURCE_TAGS[i]
        dest_t = MM_TablesXi.DEST_TAGS[i]
        if (src_t[0] == src_t[1] and dest_t[0] == dest_t[1]):
            return "// Map tag %s to tag %s."% (src_t[0], dest_t[0])
        s = "// Map tag %s to tag %s if e = 1\n" % (src_t[0], dest_t[0])
        s += "// Map tag %s to tag %s if e = 2" % (src_t[1], dest_t[1])
        return s


    def make_table_info(self):
        self.table_info = []
        sh = self.LOG_INT_FIELDS
        for i in range(1,6):
            c = self.comment(i)
            info_source =  MonomialOpTableInfo(
                MM_TablesXi.SOURCE_SHAPES[i],
                MM_TablesXi.SOURCE_START_1[i] >> sh, 
                MM_TablesXi.SOURCE_OP_DIFF[i]
            )
            info_dest =  MonomialOpTableInfo(
                MM_TablesXi.DEST_SHAPES[i],
                MM_TablesXi.DEST_START_1[i] >> sh, 
                MM_TablesXi.DEST_OP_DIFF[i]
            )
            self.table_info.append((i-1, info_source, info_dest, c))    
        self.table_diff = MM_TablesXi.MAX_ABS_START_DIFF >> sh

    def tables(self):
        return {
            "OP_XI_TABLE_INFO": self.table_info, 
            "OP_XI_TABLE_DIFF": self.table_diff
        }
     
    def directives(self):
        return {
            "OP_XI_LOAD": UserDirective(self.load_bytes, "ssi"),
            "OP_XI_STORE": UserDirective(self.store_bytes, "ssssi"),
        }

                
#234567890123456789012345678901234567890123456789012345678901234567890

if __name__ == "__main__":
    print("#include <stdint.h>")
    print("typedef uint%d_t uint_mmv_t;\n\n" % MM_Basics.INT_BITS)
    primes = [3, 7, 15, 255]
    for p in primes:
        for n in [24, 32]:
            mot = MonomialOp_t(p)
            s = """void f_%d_%d(uint_mmv_t *src, uint32_t *signs, uint8_t *b)
{"""
            print(s % (p, n)) 
            print(mot.load_bytes("src", "signs", "b", n))
            print("}\n")

    for p in primes:
        for n in [24, 32]:
            mot = MonomialOp_t(p)
            s = """void g_%d_%d(uint8_t *b, uint16_t *perm, uint_mmv_t *dest)
{"""
            print(s % (p, n)) 
            print(mot.store_bytes("b", "perm", "dest", n))
            print("}\n")
         