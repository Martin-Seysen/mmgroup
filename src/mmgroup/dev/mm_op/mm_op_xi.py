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

from mmgroup.generate_c import c_snippet, TableGenerator, make_table
from mmgroup.generate_c import UserDirective, UserFormat


from mmgroup.dev.mm_basics.mm_tables_xi import MM_TablesXi
from mmgroup.dev.mm_op.mm_op import MM_Op


class MonomialOpTableInfo:
    def __init__(self, shape, start, op_diff):
        self.SHAPE = shape
        self.START = start 
        self.OP_DIFF = op_diff


class MonomialOp_xi_uint8_t(MM_Op):
    """Standard implementation for class MonomialOp_xi

    This is the standard implementation for the monomial cases
    of operation xi.
      
    In this case the main loop is:

    uint8_t b[2496], *p_b;
    for (i = 0; i < {i_src.SHAPE[0]}; ++i) {
        p_b = b;
        
        for (j = 0; j < {i_src.SHAPE[1]}; ++j) {
           // %%OP_XI_LOAD p_src, p_b, {i_src.SHAPE[2]}
           p_src += {V24_INTS};
           p_b += 32;
        }
        
        for (j = 0; j < {i_dest.SHAPE[1]}; ++j) {
           // %%OP_XI_STORE b, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           p_dest += {V24_INTS};
           p_perm += {i_dest.SHAPE[2]};
           p_sign += 1;
        }        
    }

    """

    def __init__(self, p):
        super(MonomialOp_xi_uint8_t, self).__init__(p)
        MM_TablesXi()
        self.make_table_info()
        # make tables and directives for code generation
        self.tables.update(self.make_tables())
        self.directives.update(self.make_directives())

    def load_bytes(self, source, bytes, n_bytes, type_=""):
        s = ""
        n_registers = 8 // self.FIELD_BITS
        registers =  ", ".join(["r%d" % i for i in range(n_registers)])
        s = "uint_mmv_t %s;\n"% registers
        LOAD_MASK = self.hex(self.smask(self.P, -1, 8))
        type_ = "(%s)" % type_ if len(type_) else ""
        for i, start in enumerate(range(0, n_bytes, self.INT_FIELDS)):
            s += "r0 =  (%s)[%d];\n" % (source, i)
            for j in range(1, n_registers):
                s += "r%d = (r0 >> %d) & %s;\n" % (
                   j, j * self.FIELD_BITS, LOAD_MASK)
            s += "r0 &= %s;\n" % LOAD_MASK
            for j in range(self.INT_FIELDS):
                if start + j < n_bytes:
                    s += "(%s)[%d] = %s(r%d >> %d);\n" % (bytes,  start + j, 
                        type_, j % n_registers, 8 * (j // n_registers))
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

    def make_tables(self):
        return {
            "OP_XI_TABLE_INFO": self.table_info, 
            "OP_XI_TABLE_DIFF": self.table_diff
        }
     
    def make_directives(self):
        return {
            "OP_XI_LOAD": UserDirective(self.load_bytes, "ssis"),
            "OP_XI_STORE": UserDirective(self.store_bytes, "ssssi"),
        }



class MonomialOp_xi_uint_fast8_t(MonomialOp_xi_uint8_t):
    """This is a variant of class MonomialOp_xi_uint8_t

    Here the main loop is as in class ``MonomialOp_xi_uint8_t``.
    But the variable ``b`` may also be integer type.
 
    It might save a *partial register stall* when storing data to a
    variable of type ``uint_fast8_t`` instead of type ``uint8_t``
    """
    n_regs = {1:1, 2:2, 4:2, 8:2, 16:4, 32:4, 64:8}

    def __init__(self, p):
        super(MonomialOp_xi_uint_fast8_t, self).__init__(p)
        self.n_registers = self.n_regs[self.INT_FIELDS]

    def load_bytes(self, source, bytes, n_bytes):
        s = ""
        n_registers = self.n_registers
        registers =  ", ".join(["r%d" % i for i in range(n_registers)])
        s = "uint_mmv_t %s;\n"% registers
        WIDTH = n_registers * self.FIELD_BITS
        LOAD_MASK = self.hex(self.smask(self.P, -1, WIDTH))
        for i, start in enumerate(range(0, n_bytes, self.INT_FIELDS)):
            s += "r0 =  (%s)[%d];\n" % (source, i)
            for j in range(1, n_registers):
                s += "r%d = (r0 >> %d) & %s;\n" % (
                   j, j * self.FIELD_BITS, LOAD_MASK)
            s += "r0 &= %s;\n" % LOAD_MASK
            for j in range(self.INT_FIELDS):
                if start + j < n_bytes:
                    s += "(%s)[%d] = r%d >> %d;\n" % (bytes,  start + j, 
                        j % n_registers, WIDTH * (j // n_registers))
        return s;

    def store_bytes(self, bytes, perm, sign, dest, n_bytes):
        s = "uint_mmv_t r0, r1;\n"

        d = ["((uint_mmv_t)(%s[%s[%d]]) << %d)" % (bytes, perm, i, 
                  self.FIELD_BITS * (i % self.INT_FIELDS)) 
                  for i in range(n_bytes)
            ]
        n_registers = self.n_registers

        for i, start in enumerate(range(0, n_bytes, self.INT_FIELDS)):
            reg = "r0 = ("
            d1 = d[start:start+self.INT_FIELDS]
            M = self.smask(self.P, range(n_registers))
            while len(d1):
                hd, d1 = d1[:n_registers], d1[n_registers:]
                s += reg + "\n  + ".join(hd) + ") & %s;\n" % self.hex(M)
                reg = "r0 += ("
                M <<= n_registers * self.FIELD_BITS

            s += "r1 = (%s)[0] >> %s;\n" % (sign, start)
            s += self.gen_mmv_uint_spread("r1", "r1")
            s += "(%s)[%d] = r0 ^ r1;\n" % (dest, i)

        #s +=("// i=%d\n" % i)
        for j in range(i+1, 32 // self.INT_FIELDS):
            s += "%s[%d] = 0;\n" % (dest,j)
        return s



class MonomialOp_xi_alternative(MM_Op):
    """Alternative implementation to class MonomialOp_xi
 
    In this case the main loop is:

    for (i = 0; i < {i_dest.SHAPE[0]}; ++i) {
        for (j = 0; j < {i_dest.SHAPE[1]}; ++j) {
           // %%OP_XI p_src, p_perm, p_sign, p_dest, {i_dest.SHAPE[2]}
           p_dest += {V24_INTS};
           p_perm += {i_dest.SHAPE[2]};
           p_sign += 1;
        }
        p_src += {int:i_src.SHAPE[1] * V24_INTS};
    }

    Although here fewer data are being moved, this leads to a slower
    implementation on a 64 bit 80x86 CPU.
    """
    def __init__(self, p):
        super(MonomialOp_xi, self).__init__(p)
        self.make_table_info()
        self.tables.update(self.make_tables())
        self.directives.update(self.make_directives())

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

    def comment(self,i):
        src_t = MM_TablesXi.SOURCE_TAGS[i]
        dest_t = MM_TablesXi.DEST_TAGS[i]
        if (src_t[0] == src_t[1] and dest_t[0] == dest_t[1]):
            return "// Map tag %s to tag %s."% (src_t[0], dest_t[0])
        s = "// Map tag %s to tag %s if e = 1\n" % (src_t[0], dest_t[0])
        s += "// Map tag %s to tag %s if e = 2" % (src_t[1], dest_t[1])
        return s


    def op(self, src, perm, sign, dest, length):
        assert length in (24, 32)
        INT_FIELDS = self.INT_FIELDS
        s = "uint_mmv_t r0, rd, rp;\n"
        for i in range(0, length, INT_FIELDS): 
            op = ""
            for j in range(0, min(INT_FIELDS, length - i)):
                s += "rp = (%s)[%d];\n" % (perm, i + j)
                s += "rd = (%s)[rp >> %d] >> ((rp & %d) << %d);\n" % (
                      src, self.LOG_INT_FIELDS, INT_FIELDS - 1, 
                      self.LOG_FIELD_BITS
                )
                s += "r0 %s= ((rd & %d) << %s);\n" % (
                     op, self.p, j * self.FIELD_BITS
                )
                op = "+"
            s += "rd = (%s)[0] >> %s;\n" % (sign, i)
            s += self.gen_mmv_uint_spread("rd", "rd")
            s += "(%s)[%d] = r0 ^ rd;\n" % (dest, i // INT_FIELDS)
        for i in range(1 + (length - 1) // INT_FIELDS, 32 // INT_FIELDS):
            s += "(%s)[%d] = 0;\n"   % (dest, i)
        return s    

    def make_tables(self):
        return {
            "OP_XI_TABLE_INFO": self.table_info, 
            "OP_XI_TABLE_DIFF": self.table_diff
        }
     
    def make_directives(self):
        return {
            "OP_XI": UserDirective(self.op, "ssssi"),
        }


 

MonomialOp_xi =  MonomialOp_xi_uint8_t

class Mockup_MonomialOp_xi:
    def __init__(srlf, *args):
        pass
    tables =  {
        "OP_XI_TABLE_INFO": [], 
        "OP_XI_TABLE_DIFF": 0
    }
    def op(self, *args):
        return "\n"
    directives = {
        "OP_XI": UserDirective(op, "ssssi")
    }
        
        