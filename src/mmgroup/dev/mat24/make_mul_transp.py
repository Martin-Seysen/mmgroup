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

from mmgroup.generate_c import UserDirective, UserFormat


class BitMatrixMulTransp(object):
    """Generate code for multipling a bit vector with a bit matrix

    The generated code multiplies several bit vectors ``v_i`` by the 
    transposed of some corresponding bit matrices ``m_i``.

    Bit matrices ``m_i`` and bit vectors ``v_i`` have fixed size. 
    Vector ``v_i`` is ``lv`` bits long, where ``lv`` must be a power 
    of  two  or three times a power of two. Matrix ``m_i`` has ``lm`` 
    rows and ``lv`` columns, with ``lm <= lv``. We will actually 
    multiply  ``n`` vectors ``v_i`` by ``n`` matrices ``m_i``. Bit 
    ``v_i[k]`` is stored in bit ``i *lv + k`` of the integer ``v`` 
    and bit ``m_i[j, k]`` is stored in bit ``i *lv + k`` of the 
    integer``m[j]`` in the integer array ``m``. ``n * lv  <= 64`` 
    must hold. The result is stored in the variable ``v`` in the 
    same format as the input ``v``.

    The sizes ``lv, lm, n`` must be set with method ``set_matrix``.
    Variable names ``m`` and ``v`` must be set with method 
    ``set_names``.
     
    ``v`` should be an integer and ``m`` should be an array of 
    integers of type ``uint_fast64``. In case ``n * lv <= 32`` 
    that type may also be ``uint_fast32``.  
    """ 

    def set_matrix(self, lv, lm, n=1):
        """Set lv, lm, n  as in the description of the class"""
        assert 0 < lm <= lv <= lv * n <= 64, (lm, n, lv)
        self.lv, self.lm, self.n = lv, lm, n
        self.wordlength = 32 if n * lv <= 32 else 64
        self.int_fmt = "0x%xUL" if self.wordlength == 32 else "0x%xULL"
        self.int_type = "uint_fast%d_t" % self.wordlength
        self.masks = {}
        self.all_mask = (1 << (n * lv)) - 1
        self.nops = 0
        self.temp_name = "t_"
        self.n_temps = self.max_temps = 0

    def set_names(self, v, m, temp = "t"):
        """Set names of bit vector ``v`` and bit matrix ``m``

        Optionally, a prefix for temporary variables may be set.
        """
        self.v, self.m, self.temp =  v, m, temp
        self.n_temps = 0
        self.nops = 0

    def mask(self, space, width, pos):
        """Return a certain bit mask as a C literal

        In the returned bit mask we set the bits at positions

        ``k * space + pos * width + j,   0 <= j < width``,

        for all ``k >= 0`` such that the bit position is less than
        ``n * lv``, with ``n, lv`` as given by method ``set_matrix``.

        ``space`` must be divisible by ``lv`` and ``width`` must be
        divisible by ``space``.
 
        Since the same mask is used frequently, we use a dictionary
        for storing the masks already created. 
        """
        try:
            return self.masks[(space, width, pos)]
        except KeyError:
            assert self.lv % space == space % width == 0
            m = self.all_mask // ((1 << space) - 1) 
            m *= (1 << width) - 1
            m = m << (pos * width)
            s = self.masks[(space, width, pos)] = self.int_fmt % m
            return s
 
    def var(self, i):
        """Return the name of the i-th temporary variable"""
        return  self.temp_name + str(i)

    def new_var(self):
        """Return the name of the next unused temporary variable

        The returned variable is marked as used.
        """
        v = self.var(self.n_temps)
        self.n_temps += 1
        self.max_temps = max(self.n_temps, self.max_temps)
        return v

    def release_vars(self, n_vars = 1):
        """Mark temporary variables as unused

        We mark the ``n_vars`` last recently used variables as unused.
        ``n_vars`` defaults to 1.
        """
        self.n_temps -= n_vars
   
    def gen_start(self, i):
        """Generate code to process ``m[i]`` 

        The function returns a pair ``(var, code)``. Here ``var``
        is a new temporary variable containing the result and ``code``
        is the C code generated for computing that result.

        Bit ``[k]`` of the result ``var`` is equal to``v[k] & m[i,k] 
        for all relevant indices ``k``.  ``var`` is None if the result 
        is always zero. This happens if ``i >= self.lm``.
        """
        if i >= self.lm:   
            return None, ""
        tmp = self.new_var()
        c = "{tmp} =  {m}[{i}] & {v};\n".format( tmp = tmp,
             m = self.m, i = str(i), v = self.v )
        self.n_ops += 1
        return tmp, c

    def shl(self, var, shift):
        """Return the C expression ``(var << shift)``

        ``var`` must be the name of a variable and ``shift`` must
        be an integer corresponding to a fixed shift factor.  In case 
        ``sh < 0`` we generate the expression ``(var >> -shift)``.
        """
        if shift == 0:
            return var
        elif shift > 0:
            self.n_ops += 1
            return "(%s << %d)" % (var, shift)
        else:
            self.n_ops += 1
            return "(%s >> %d)" % (var, -shift)
     

    def merge(self, var, width, blocks, pos):
        """Generate code for XORing blocks of variable ``var``

        For the integer variable ``var`` the generated code computes
        the XOR sum ``x`` sum of the values ``var >> (i * width)``,
        ``0 <= i < blocks``. Then ``x`` is masked with the mask
        ``self.mask(blocks * width, width, 0)`` and shifted
        left by ``pos * blocks`` bits. The (masked and shifted) 
        result ``x`` is stored in variable ``var``.
        """
        assert blocks > 0
        if not var or blocks == 1:
            return ""
        parts = [self.shl(var, width * (pos-i)) for i in range(blocks)]
        self.n_ops += blocks
        mask = self.mask(width * blocks, width, pos)
        return "%s = (%s) & %s;\n" % (var, " ^ ".join(parts), mask)
  

    def gen_step(self, st, width):
        """Generate code to process ``m[st],...,m[st + width - 1]`` 

        The function returns a pair ``(var, code)``. Here ``var``
        is a new temporary variable containing the result and ``code``
        is the C code generated for computing that result.

        Bit ``[k*width + j]`` of the result ``var`` is equal to

        ``sum(v[l] & m[st+j, l] for l in range(k*width, (k+1)*width))``

        for ``0 <= j < width``. Here ``width`` must be a factor of
        the bit length of ``v``  and it should be a power of two
        or three times a power of two.

        The function recursively divides the blocks of length 
        ``width`` into two or three block of equal length
        and processes the parts of block.
       
        So in case ``width = 1`` the result is just ``v & m[st]``, 
        which is the content of the variable returned by of member 
        function ``self.gen_start(st).

        In case ``width = lv, st = 0`` the content of bit ``[i*lv + j]``
        of the the  variable ``var``  returned by this function with
        is eqaual to

        ``sum(v[l] & m[j, l] for l in range(i*lv, (i+1)*lv)``

        which is just bit ``j`` of the requested matrix product
        of the bit vector ``v_i`` with the transposed  matrix of 
        the bit matrix ``m_i``. 

        The function returns the pair ``(None,"")`` if the result
        is equal to zero. This happens in case ``st >= self.lm``

        Though the code generation process is recursive, the generated 
        code is very efficient and, obvioulsy, not recursive.
        """
        if st >= self.lm:   
            return None, ""
        if width == 1:
            return self.gen_start(st)
        if width & (width - 1) == 0:
            blocks = 2
        elif width % 3 == 0:
            blocks = 3
        else:
            err = "Illegal vector length in class BitMatrixMulTransp"
            raise ValueError(err)
        width = width // blocks
        var, c = self.gen_step(st, width) 
        m = self.merge(var, width, blocks, 0)
        code = c + m
        for j in range(1, blocks):
            var1, c = self.gen_step(st + j * width, width)
            if var1:
                m = self.merge(var1, width, blocks, j) 
                code += c + m
                code += "%s |= %s;\n" % (var, var1)
                self.n_ops += 1
                self.release_vars(1)
        return var, code  
            
    def define_temps(self):
        """Generate code for defining temporary variables"""
        return "%s %s;\n" % (self.int_type, 
             ", ".join([self.var(i) for i in range(self.max_temps)]))

    def generate(self):
        """Generate code for doing the requested matrix multiplication"""
        self.n_ops = 0
        var, code = self.gen_step(0, self.lv)
        out = "// %d operations\n%s = %s;\n" % (self.n_ops, self.v, var)
        s = self.comment() + "{\n" + self.define_temps() + code + out + "}\n"
        return s

    def comment(self):
        s = """// Compute ``v_i = v_i * transpose(m_i), i = 0,...,{n1}``, ``v_i`` 
// a bit vector of length {lv}, ``m_i`` a {lm} times {lv}  bit 
// matrix. ``v_i`` is coded in the integer ``{name_v}``, ``m_i`` 
// is coded in the array  ``{name_m}`` of integers. Bits ``v_i[k]`` 
// and ``m_i[j][k]`` are coded in bits ``{lv}*i + k`` of the 
// variables  ``{name_v}`` and ``{name_m}[j]``, ``j = 0,...,{lm1}``, 
// respectively. 
""".format(n1 = self.n-1, lv = self.lv, lm = self.lm, lv1 = self.lv - 1, 
        lm1 = self.lm - 1, name_v =  self.v, name_m = self.m ) 
        return s
        
    def n_operations(self):
        """Return number of Boolean operations required"""
        self.gen_step()
        return self.n_ops

    def generate_c_mul_transp(self, v, m, lv, lm, n=1):
        """Generate code for multiplying bit vector with matrix

        ``v`` is the name of an integer encoding ``n`` bit vectors
        of length ``lv``. ``m`` is the name of an array of integers 
        encoding ``n`` bit matrices of shape ``lm`` times ``lv``.

        More details of the computation are given in the description
        of the class. The computedbit vectors are encoded in ``v``.
        """
        n = max(1, int(n) if n else 0)
        self.set_matrix(int(lv), int(lm), n)
        self.set_names(str(v), str(m)) 
        c = self.generate()
        #print(b.n_ops, "operations for generate_c_mul_transp", lv, lm, n)
        return c

    def parity(self, v):
        """Return bit parity of the unsigned 64-bit integer ``v``"""
        v ^= v >> 32; v ^= v >> 16; v ^= v >> 8; v ^= v >> 4
        return (0x6996 >> (v & 0x0f)) & 1

    def compute(self, v, m):
        """Compute requested matrix product of ``v`` and ``transpose(m)``

        Here ``v`` is an integer and ``m`` is a list of integers. The
        matrix product given by ``lv, lm, n`` in method ``set_matrix``
        (as described in the header of the class) is computed and
        returned as an integer.
        """
        res = 0
        rows = m[:self.lm]
        for i in range(self.n):
            sh = i * self.lv
            v1 = sum(self.parity((v & r) >> sh) << j 
                for j, r in enumerate(rows))
            res  |= v1 << sh 
        return res
    
    def extract(self, var, start, length, dest = 0):
        """Return C expression containing some bits of ``var``
        
        Here ``var`` is an integer variabbe. The function extracts
        bits ``start, ..., start + length - 1`` variable and returns
        an integer expression where these bits are copied to bits
        ``start, ..., start + length - 1``. The other bits  of the  
        returned expression are cleared.
        """
        start, length = int(start), int(length)
        dest = int(dest) if dest else 0
        v = self.shl(var, dest - start)
        mask = ((1 << length) - 1) << dest
        fmt = "0x%xUL" if dest + length <= 32 else "0x%xULL"
        return "(%s & %s)" % (v, fmt % mask)

    def tables(self):
        """This class generates no tables"""
        return {
           "BITV_EXTRACT" : UserFormat(self.extract, "sii")
        }

    def directives(self):
        """Code generation directive has name "BITVMULTRANSP". """
        return {
           "BITVMULTRANSP" : UserDirective(self.generate_c_mul_transp, "ss"),
        }

