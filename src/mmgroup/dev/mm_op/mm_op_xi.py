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
    r"""Standard implementation for class MonomialOp_xi

    This is the standard implementation for the monomial cases
    of operation xi.
      
    In this case the main block for computing the monomial part of
    ``v_in * xi**(e1 + 1)`` in ``v_out`` is:

    .. code-block:: c

      uint_mmv_t *p_src, *p_dest;
      uint_fast32_t j, k;
      uint8_t b[2496], *p_b;
      mm_sub_table_xi_type *p_tables;
      uint16_t *p_perm;
      uint32_t *p_sign;
      uint32_t L =  %{LOG_INT_FIELDS};

      // %%FOR* i  in range(5)
      // %%WITH* SRC_SHAPE = MM_TABLE_SHAPES_XI[i][0]
      // %%WITH* DEST_SHAPE = MM_TABLE_SHAPES_XI[i][1]
      p_src = v_in + (MM_SUB_OFFSET_TABLE_XI[%{i}][e1][0] >> L);
      p_dest = v_out + (MM_SUB_OFFSET_TABLE_XI[%{i}][e1][1] >> L);
      p_tables =  &MM_SUB_TABLE_XI[%{i}][e1];
      p_sign = p_tables->p_sign;
      p_perm = p_tables->p_perm;

      for (j = 0; j < %{SRC_SHAPE[0]}; ++j) {
          p_b = b;
          for (k = 0; k < %{SRC_SHAPE[1]}; ++k) {
             // %%OP_XI_LOAD p_src, p_b, SRC_SHAPE[2], uint8_t
             p_src += %{V24_INTS};
             p_b += 32;
          }
        
          for (k = 0; k < %{DEST_SHAPE[1]}; ++k) {
             // %%OP_XI_STORE b, p_perm, p_sign, p_dest, DEST_SHAPE[2]
             p_dest += %{V24_INTS};
             p_perm += %{DEST_SHAPE[2]};
             p_sign += 1;
          }
      }
      // %%END WITH
      // %%END WITH
      // %%END FOR

    """

    def __init__(self, p):
        super(MonomialOp_xi_uint8_t, self).__init__(p = p)
        self.tables.update(self.make_tables())
        sub = MM_TablesXi()
        self.tables["MM_TABLE_SHAPES_XI"] = sub.SHAPES 
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

    def comment(self, i):
        """This does not work!!!"""
        boxes = [MM_TablesXi.TABLE_BOX_NAMES[i,j] for j in (0,1)]
        if boxes[0] == boxes[1]:
            args = boxes[0][0], boxes[0][1]
            return "// Map box %s to box %s\n" % args
        else:
            fmt =  "// Map box %s to box %s if exponent is %d\n"
            s = ""
            for exp1, (src, dest) in enumerate(boxes):
                s += fmt % (src, dest, exp1 + 1)
            return s

    def make_tables(self):
        return {}
     
    def make_directives(self):
        return {
            "OP_XI_LOAD": UserDirective(self.load_bytes, "ssis"),
            "OP_XI_STORE": UserDirective(self.store_bytes, "ssssi"),
            "OP_XI_COMMENT": UserDirective(self.comment, "i"),
        }



class MonomialOp_xi_uint_fast8_t(MonomialOp_xi_uint8_t):
    """This is a variant of class MonomialOp_xi_uint8_t

    Here the main loop is as in class ``MonomialOp_xi_uint8_t``.
    But the variable ``b`` may also be integer type.
 
    It might save a *partial register stall* when storing data to a
    variable of type ``uint_fast8_t`` instead of type ``uint8_t``.

    Warning: This class has not been tested and will probably fail!!
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

        for j in range(i+1, 32 // self.INT_FIELDS):
            s += "%s[%d] = 0;\n" % (dest,j)
        return s




class Tables(MonomialOp_xi_uint8_t):
    def __init__(self, *args, **kwds):
        p = int(kwds.get('p', 3))
        super().__init__(p)
     

class MockupTables:
    tables =  {}
    directives = {}
    def __init__(self, *args, **kwds):
        pass
        
        