"""The module contains class ``BitMatrixMulTransp``.

Class ``BitMatrixMulTransp`` contains a directive that generates
code to multiply a bit vector with the transposed matrix of a
fixed bit matrix. Depending on the size of the matrix and the
bit length of the underlying integer type, several bit vectors
can be multiplied with the same matrix simultaneously.

The C function ``mat24_autpl_set_qform`` uses that matrix 
multiplication, see section :ref:`implement-autpl-label`
in the *Guide for developers* for background.
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint
import types
import sys
import re
import os
from operator import __or__, __xor__

from mmgroup.generate_c import UserDirective


class BitMatrixMulTransp(object):
    """Multiplication of bit vector v by the transpose of bit matrix m.

    Bit matrix m and bit vector v have fixed size. 
    Vector v is lv bits long, where lv must be a power of two. Matrix m 
    has lm rows and lv columns, with lm <= lv. Bitvectors and matrices 
    are stored according to the conventions in bitfunctions.py.

    The purpose of this class is to generate code for computing
    v = v * transposed(m). There must be a C integer type uint_fast<x> 
    that supports unsigned integers of x bit length with  x = lv.
    input/output variable v should be an unsigned in of bit length
    at least max(lv, lm). 

    Actually, we can do several multplications  v_i * transposed(m_i), 
    i = 0,...,n-1 in with a single call to this function, provided that
    an integer is at least n * lv bits long. We support n 
    simultaneous multiplications for n being a power of two.
    Then the i-th component of v_i and of a row of m_i is stored in
    bits  i*lv, ..., i*lv + lv-1  of the corresponding integer.

    """ 

    

    def set_matrix(self, lv, lm, n=1):
        """Set lv, lm, n  as in the description of the class"""
        assert 0 <= lm <= lv <= lv * n <= 64, (lm, n, lv)
        assert (lv * n) & (lv * n - 1) == 0  
        self.lv, self.lm, self.n = lv, lm, n
        self.nsteps = 0
        while 1 << self.nsteps < lv:
            self.nsteps += 1
        assert  1 << self.nsteps == lv
        self.masks = {}
        self.nops = 0
        self.vtype = "uint_fast%d_t" % (max(8, n * lv))
        self.v = self.m = self.temp = "??"
       

    def setnames(self, v, m, temp = "_t"):
        """Set names of bit vector v and bit matrix m

        Optionally, a prefix for temporary variables may be set.
        """
        self.v, self.m, self.temp =  v, m, temp
        self.maxvar = 0

    def isnull(self, start):
        """Check if a row of bit matrix m i supposed to be zero"""
        return not (start < self.lm)


    def var(self, i):
        """Return name of the i-th temporary variable"""
        self.maxvar = max(i+1, self.maxvar)
        return self.temp + str(i)

    def mk_mask(self, st, high):
        """Return fixed mask for stage st as a C integer literal.

        In case high=0 the mask has bits  k*sh, ..., k*sh + sh - 1  
        set for  sh = 1 << st  and  k = 0,2,4,6,... . 

        The mask for high=1 is obtained by shifting the corresponding
        mask for high=0 left by  1 << st  bits.
        """ 
        try:
            return self.masks(sh, high)
        except:
            sh = 1 << st
            m = (1 << sh) - 1
            nsteps = self.n << self.nsteps
            m = sum(m << i for i in range(0, nsteps, 2 * sh))
            if high:
                m <<= sh

            m = hex(m)
            if m[-1] in "lL": m = m[:-1]
            m = m + "U"
            if nsteps > 16: m = m + "L"
            if nsteps > 32: m = m + "L"
            self.masks[(sh, high)] = m
            return m
    
    


    def gen_start(self, i, tmp):
        """Genrate code to process row i of matrix m at stage 0.

        Result is stored in temporary variable with number tmp.
        This amounts to setting "tmp" = m[i] & v. If m[i] is supposed
        to be zero, an empty string is returned. 
        """
        if self.isnull(i):   
            return ""
        c = "{tmp} =  {m}[{i}] & {v};\n".format( tmp = self.var(tmp),
             m = self.m, i = str(i), v = self.v )
        self.n_ops += 1
        return c

    def gen_step(self, i=0, st=None, tmp=0, result=None):
       """Genreate code to process m[i],...,m[i+2**st-1] at stage st.

       Number i must be a multiple of sh, where sh = 2**st. The
       result of this stage is stored in the temporary variable with 
       number tmp. The content of that variable after stage  st  is:

       bit[k*sh + j] 
         =  sum(v[k*sh + l] & m[i + j, k*sh + l] for l in range(sh)),

       for 0 <= j < sh.

       The function recursively calls gen_step(i, st-1) and
       gen_step(i + sh/2, st-1), and computes the result in stage
       st  from these two results in stage st - 1.

       So in case st=0 the result is just v & m[i], which is the
       content of the variable returned by of member function 
       self.gen_start(i).

       If  st  is suffciently high, then the content of bit j of
       the  variable  returned by this function with i=0 is 

          sum(v[l] & m[j,l] for l in range(sh)),

       which is just bit j of the requested matrix product.

       Temproary variable tmp is also used in recursive calculations.
       If 'result' is set, the final result is stored in variable
       with name 'result'.

       Though the code generation process is recursive, the generated 
       code is very efficient and, obvioulsy, not recursive.
       """
       if st is None: 
            st, self.n_ops, self.maxvar = self.nsteps, 0, 0
       if result is None:
           result = self.var(tmp)
       if st == 0:
           return self.gen_start(i, tmp)
       st -= 1
       sh = 1 << st
       c1 = self.gen_step(i, st, tmp) 
       tmp2 = tmp + bool(c1)
       c2 = self.gen_step(i + sh, st, tmp2)
       if c2:
           ct2 = "(({tmp} << {sh}) ^ {tmp}) & {mask}".format(
            tmp = self.var(tmp2), sh = str(sh), mask = self.mk_mask(st,1)
           ) 
           self.n_ops += 3 
       if c1:
           ct1 = "(({tmp} >> {sh}) ^ {tmp}) & {mask}".format(
            tmp = self.var(tmp), sh = str(sh), mask = self.mk_mask(st,0) 
           ) 
           self.n_ops += 3 
           if c2:
               c = "%s = (%s)\n    | (%s);\n" % (result, ct1, ct2)
               self.n_ops += 1
           else:
               c = "%s = %s;\n" % (result, ct1)
       else:
           if c2:
               c = "%s = %s;\n" % (result, ct2)
           else:
                return  ""
       return c1 + c2 + c

    def define_temps(self):
        """Generate code for defining temporary variablea"""
        return "%s %s;\n" % (self.vtype, 
             ", ".join([self.var(i) for i in range(self.maxvar)]) )

    def generate(self):
        """Generate code for doing the requested matrix multiplication"""
        self.n_ops = 0
        code = self.gen_step(result=self.v)
        ops = "// %d operations\n" % self.n_ops
        s = self.comment() + "{\n" + self.define_temps() + code + ops + "}\n"
        return s

    def n_operations(self):
        """Return number of Boolean operations required"""
        self.gen_step()
        return self.n_ops

    def parity(self, v):
        """Return parity of least 1 << self.nsteps bits of v"""
        for i in range(self.nsteps):
             v ^= v >> (1 << i)
        return v & 1

    def compute(self, v, m):
        """Compute the requested matrix product of v and transpose(m)"""
        res = 0
        rows = m[:self.lm]
        for i in range(self.n):
            sh = i * self.lv
            v1 = sum(self.parity((v&r) >> sh) << j for j,r in enumerate(rows))
            res  |= v1 << sh 
        return res

    def comment(self):
        s = """// Compute v = v * transpose(m), v a bit vector, m a bit matrix,
// v coded as an integer with {0} valid bits, m a {1} x {0} bit matrix
// coded as an array of integers.
""".format(self.lv, self.lm) 
        return s
        


    def generate_c_mul_transp(self, v, m, lv, lm, n=""):
        n = max(1, int(n) if n else 0)
        self.set_matrix(int(lv), int(lm), n)
        self.setnames(str(v), str(m)) 
        c = self.generate()
        #print(b.n_ops, "operations for generate_c_mul_transp", lv, lm, n)
        return c
    

    def tables(self):
        """This class generates no tables"""
        return {}

    def directives(self):
        """Code generation directive has name "BITVMULTRANSP". """
        return {
           "BITVMULTRANSP" : UserDirective(self.generate_c_mul_transp, "ss"),
        }

