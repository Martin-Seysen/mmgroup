from __future__ import absolute_import, division, print_function
#from __future__ import  unicode_literals


import sys
import os

from mmgroup.bitfunctions import bitparity, bitweight
from mmgroup.dev.hadamard.hadamard_codegen import HadamardMatrixCode
from mmgroup.generate_c import UserDirective, UserFormat


class HadamardOpT64(HadamardMatrixCode):
    """Code for multiplying vector by a modified Hadamard matrix

    This class supports the operation of the triality element t**e,
    e = 1,2  on a block v with tag T in the representation 196884x 
    of the monster group as described in[Seys19]. We have 
    t**3 == 1. Here t is a 64 times 64 matrix acting on v by right
    multiplication We have t = M * D and t**2 = D * M for a
    diagonal matrix D and a Hadamard-linke matrix M as described 
    below. 

    For D we have D[i,i] = (-1)**((bitweight(i) + 1) >> 1).

    M is the parity-adjusted 2**n times 2**n Hadamard matrix 
    for an even number n as defined in the header of module
    define_matrices.py. This class supports the generation of
    C code for multiplying a vector of integers modulo p
    by a parity-adjusted Hadamard matrix M.

    Notation and implementation of vectors is as in the base 
    class HadamardMatrixCode, which supports right-multiplication 
    of a vector by a Hadamard  matrix.

    To 'upgrade' multiplication by a Hadamard matrix to
    multiplication by a parity-adjusted Hadamard matrix M, we 
    need to code a permutation of the entries of a vector v 
    that exchanges v[i] with v[63 - i] if i has odd parity 
    and that fixes v[i] if i has even parity.
    """
    def __init__(self, p, log_vlength = 6, verbose = 0):
        """Create an instance of the class.

        Parameters are as in the base class HadamardMatrixCode,
        with the restriction that log_vlength must be even.
        """
        args = (p, log_vlength, verbose)
        super(HadamardOpT64, self).__init__(*args)
        assert log_vlength & 1== 0
        self.make_parity_masks()
        self.make_diagonal_mask()
        self.directives.update(self.make_directives())

    def make_parity_masks(self):
        """Create parity masks

        self.PARITY_MASK[0] and self.PARITY_MASK[1] select the 
        entries v[j] of a vector stored in a single variable, 
        where index j has even and odd parity, respectively.
        """
        parity = sum((bitparity(i) << i for i in range(self.INT_FIELDS)))
        mask0 = self.smask(self.P, ~parity)
        mask1 = self.smask(self.P, parity)
        self.PARITY_MASK = [mask0, mask1]

    def make_diagonal_mask(self):
        """Create masks for multiplication with diagonal matrix

        For the diagonal matrix D we have 

           D[i,i] = (-1)**((bitweight(i) + 1) >> 1),   0 <= i < 64,
        
        where bitweight(i) is the bit weight of the binary number i.

        For j, k let i = j * self.INT_FIELDS + k. Then we set
        bit field k in self.DIAG_MASK[j] to 0 if D[i,i] = 1 and to
        self.P if D[i,i] = -1.        
        """
        diag = sum((bool((bitweight(i)+1) & 2) << i for i in range(64)))
        self.DIAG_MASK = []
        for i in range(0, 64, self.INT_FIELDS):
            self.DIAG_MASK.append(self.smask(self.P, diag >> i)) 

    def merge_parities(self, dest, var0, var1, parity = 0):
        """Merge entries of dest from var0 and var1

        Here dest, var0, and var1 are vectors stored in a
        single variable.

        If the integer 'parity' has even bit parity then the
        entry var0[i] copied to dest[i] if i has even parity
        and the entry var1[i] is copied to dest[i] if i has
        odd parity. If the integer 'parity' has odd parity, 
        the roles of var0 and var1 are echanged.
        """
        mask0, mask1 = self.PARITY_MASK
        if bitparity(parity) == 1:
            mask0, mask1 = mask1, mask0
        dest.assign((var0 & mask0) | (var1 & mask1))

    def swap_fields(self, var):
        """Swap entries of vector stored in variable 'var'.

        Here var is a variable storing a vector 
        var[0], ..., var[m - 1], for m = self.INT_FIELD.
        The function exchanges var[i] with var[m - i - 1].
        """
        for i0 in range(self.LOG_INT_FIELDS - 1, -1, -1):
            i = 1 << i0
            msk = self.B0_MASK[i]
            sh = i << self.LOG_FIELD_BITS
            if sh << 1 == self.INT_BITS:
                var.assign((var << sh) | (var >> sh))
            else:
                var.assign(((var & msk) << sh) | ((var >> sh) & msk))

    def swap_odd_parities_log_nvars_even(self):
        """Auxiliary method for method swap_parities.

        It handles the case when self.LOG_INT_FIELDS is odd.
        """
        l = len(self.vars)
        t0, t1 = self.vars.temp(0), self.vars.temp(1)
        for i in range(0, l >> 1, 2):
            v0, v1 = self.vars[i], self.vars[i + 1]
            v2, v3 = self.vars[l - i - 2], self.vars[l - i - 1]
            self.merge_parities(t0, v0, v1, i ^ 1)
            self.swap_fields(t0)
            self.merge_parities(t1, v2, v3, i)
            self.swap_fields(t1)
            self.merge_parities(v0, v0, t1, i)
            self.merge_parities(v1, v1, t1, i ^ 1)
            self.merge_parities(v2, v2, t0, i ^ 1)
            self.merge_parities(v3, v3, t0, i)
            
            
    def swap_odd_parities_log_nvars_odd(self):
        """Auxiliary method for method swap_parities.

        It handles the case when self.LOG_INT_FIELDS is even.
        """
        l = len(self.vars)
        t0 = self.vars.temp(0)
        for i in range(0, l >> 1, 1):
            v0, v1 = self.vars[i], self.vars[l - i - 1]
            self.merge_parities(t0, v0, v1, i ^ 1)
            self.swap_fields(t0)
            self.merge_parities(v0, v0, t0, i)
            self.merge_parities(v1, v1, t0, i ^ 1)
            
    def swap_parities(self):
        """Do a certain permutaition on the entries of the vector

        The function operates in a vector v[0],...,v[2**n-1] 
        stored in variables as described in the base class
        HadamardMatrixCode. It exchanges v[i] with v[2**n-i-1]
        if i has odd parity and it fixes v[i] if i has odd
        parity.
        
        Here n is the parameter log_vlength given in the 
        constructor. The function is defined for even n only.
        """
        if self.LOG_VLEN & 1:
            s = "Size of Hadamard matrix must be a power of 4"
            raise ValueError(s)
        if len(self.vars) < 2:
            s = "Size of Hadamard matrix too small for p = %d"
            raise NotImplementedError(s % self.P)
        odd =  self.LOG_INT_FIELDS & 1
        self.comment(
"""Exchange component i with component %d-i if i 
has odd parity; fix it if i has even parity."""
           % ((1 << self.LOG_VLEN) - 1)
        ) 
        if odd:
            self.swap_odd_parities_log_nvars_odd()
        else: 
            self.swap_odd_parities_log_nvars_even()


    def load_vector_mul_diagonal(self, array_name, mask):
        """Load a vector v form an (external) array A of integers.

        Here 'mask' is an integer of type uint_mmv_t wich must be
        0 or (uint_mmv_t)(-1). If mask == -1 then the vector is 
        multiplied with the dagonal matrix D described in method
        make_diagonal_mask().

        'array_name' must a variable of type uint_mmv_t 
        referring to array A as in method load_vector_direct(). 
        """
        self.comment(
"""Loading vector v from array %s; multiply v
with diagonal matrix if %s == -1."""
            % (array_name, mask)
        )
        for i,v in enumerate(self.vars):
            s = "%s = %s[%d] ^ ((%s) & %s);\n" % (v, 
                array_name, i, mask, self.hex(self.DIAG_MASK[i])) 
            self.add(s, 1, 2)
        self.comment_vector()  

    def store_vector_mul_diagonal(self, mask, array_name):
        """Store a vector v to an (external) array A of integers.

        Here 'mask' is an integer of type uint_mmv_t wich must be
        0 or (uint_mmv_t)(-1). If mask == -1 then the vector is 
        multiplied with the diagonal matrix D described in method
        make_diagonal_mask().

        'array_name' must a variable of type uint_mmv_t 
        referring to array A.
        """ 
        self.comment(
"""Storing vector v to array %s; multiply v
with diagonal matrix if %s == -1."""
            % (array_name, mask)
        )
        s = ""
        for i,v in enumerate(self.vars):
            s = "%s[%d] = %s ^ ((%s) & %s);\n" % (array_name, i, 
                v, mask, self.hex(self.DIAG_MASK[i])) 
            self.add(s, 1, 2)


    def make_code(self, source, mask, dest):
        """Apply triality operation on vector.

        Right multiply the vector of integers mod self.P stored in 
        'src' by t**e, where t is the 64 times 64 triality matrix 
        operating on blocks of the rep 196884x with tag. 

        For e = 1, parameter 'mask' must be 0; for e = 2, parameter
        'mask' must bei (uint_mmv_t)(-1). The result is stored
        in 'dest'.

        'source' and 'dest' must be pointers of type uint_mmv_t*,
        'mask' must be a variable of type uint_mmv_t.
        """
        self.reset_vars()
        self.comment(
"""Multiply the vector of integers mod {p} stored
in ({src}) by t**e, where t is the 64 times 64 
triality matrix and e = 1 if {mask} = 0, e = 2 if
{mask} = (uint_mmv_t)(-1). The result is stored
in ({dest}).
""".format(p = self.P, src = source, mask = mask, dest = dest)
        )
        self.load_vector_mul_diagonal(source, mask)
        self.complement_variable(mask)
        self.swap_parities() 
        self.hadamard_op(shift = -3)
        #self.mul_pwr2(-3)
        self.store_vector_mul_diagonal(mask, dest)
        self.complement_variable(mask)
        self.comment_statistics()
        return self.generate()


    def make_directives(self):
        return {
          "MUL_MATRIX_T64" : UserDirective(self.make_code, "sss"),
        }
        



class HadamardOpT3(HadamardMatrixCode):
    """Apply triality element to tags A, B, C

    Yet to be documented!!!!
    """
        
    def __init__(self, p, verbose = 0):
        """Create an instance of the class.

        Parameters are as in the base class HadamardMatrixCode,
        with the restriction that log_vlength must be even.
        """
        super(HadamardOpT3, self).__init__(p, 0, verbose)
        self.set_vars()
        self.directives.update(self.make_directives())

    def set_vars(self):
        self.vlen = 3 << self.NO_CARRY
        self.vars.resize(self.vlen)
        

    def reset_vars(self):
        super(HadamardOpT3, self).reset_vars()
        self.set_vars()


    def main_op(self):
        """yet to be documented!!!"""
        if self.NO_CARRY:
            for i in range(3):
                self.expand_hadamard_extern(self.vars[i], self.vars[i+3])
            self.external_butterfly(self.vars[4], self.vars[5])
            self.mul_var_pwr2(self.vars[4], -1)
            self.mul_var_pwr2(self.vars[5], -1)
            self.external_butterfly(self.vars[3], self.vars[5])
        self.external_butterfly(self.vars[1], self.vars[2])
        self.mul_var_pwr2(self.vars[1], -1)
        self.mul_var_pwr2(self.vars[2], -1)
        self.external_butterfly(self.vars[0], self.vars[2])
        if self.NO_CARRY:
            for i in range(3):
                self.compress_hadamard_extern(self.vars[i], self.vars[i+3])
        self.vars.xch(0,1)


    def load_vector_mul_diagonal(self, array_name, mask):
        """Load vector v from tags A, B, C of rep 196884x.

        Here 'mask' is an integer of type uint_mmv_t wich must be
        0 or (uint_mmv_t)(-1). If mask == -1 then the vector v is 
        multiplied with the dagonal matrix D described in method
        make_diagonal_mask().

        'array_name' must a variable of type uint_mmv_t 
        referring to array A as in method load_vector_direct(). 
        More specifically, A contains a set of componets of
        a vector for the rep 196884x with tag A.
        The function also loads the corresponding sets of 
        components of the vector with tags B and C.
        """
        self.comment(
"""Loading vector from rep 196884x with tags A,B,C
to v[0...2]. Here %s refers to the tag A part. 
Negate v[2] if %s == -1."""
            % (array_name, mask)
        )
        s = "%s = (%s)[0];\n" % (self.vars[0], array_name)
        s += "%s = (%s)[%d];\n" % (self.vars[1], array_name,
              24 * self.V24_INTS )
        s += "%s = (%s)[%d] ^ ((%s) & %s);\n" % (self.vars[2], 
                array_name, 48 * self.V24_INTS, mask, 
                self.hex(self.smask(self.P)) ) 
        self.add(s, 3, 2)
        self.comment_vector()  

    def store_vector_mul_diagonal(self, mask, array_name):
        """Store vector v to tags A, B, C of rep 196884x.

        Here 'mask' is an integer of type uint_mmv_t wich must be
        0 or (uint_mmv_t)(-1). If mask == -1 then the vector v is 
        multiplied with the dagonal matrix D described in method
        make_diagonal_mask().

        'array_name' must a variable of type uint_mmv_t 
        referring to array A as in method load_vector_direct(). 
        More specifically, A contains a set of componets of
        a vector for the rep 196884x with tag A.
        The function also loads the corresponding sets of 
        components of the vector with tags B and C.
        """
        self.comment(
"""Store vector v[0...2] to rep 196884x with 
tags A,B,C. Here %s refers to the tag A part. 
Negate v[2] if %s == -1."""
            % (array_name, mask)
        )
        s = "(%s)[0] = %s;\n" % (array_name, self.vars[0])
        s += "(%s)[%d] = %s;\n" % (array_name, 
              24 * self.V24_INTS, self.vars[1] )
        s += "(%s)[%d]  = %s ^ ((%s) & %s);\n" % ( 
                array_name, 48 * self.V24_INTS, self.vars[2],
                mask,  self.hex(self.smask(self.P)) ) 
        self.add(s, 3, 2)


    def make_code(self, source, mask, dest):
        """Apply triality operation on vector.

        Right multiply the vector of integers mod self.P stored in 
        'src' by t**e, where t is the 64 times 64 triality matrix 
        operating on blocks of the rep 196884x with tag. 

        For e = 1, parameter 'mask' must be 0; for e = 2, parameter
        'mask' must bei (uint_mmv_t)(-1). The result is stored
        in 'dest'.

        'source' and 'dest' must be pointers of type uint_mmv_t*,
        'mask' must be a variable of type uint_mmv_t.
        """
        self.reset_vars()
        self.comment(
"""Multiply the vector of integers mod {p} stored in
({src}) by t**e, where t is the 3 times 3 triality
matrix [[0, 2,  -2], [1, 1, 1], [1,  -1, -1]] / 2.
and e = 1 if {mask} = 0, e = 2 if {mask} = 
(uint_mmv_t)(-1). The result is stored in ({dest}).

{src} and {dest} are pointers of type *uint_mmv_t.
Components with tags A, B, C referred by ({src}) 
are processed, one integer of type uint_mmv_t
for each tag.

""".format(p = self.P, src = source, mask = mask, dest = dest)
        )
        self.load_vector_mul_diagonal(source, mask)
        self.complement_variable(mask)
        self.main_op()
        self.store_vector_mul_diagonal(mask, dest)
        self.complement_variable(mask)
        self.comment_statistics()
        return self.generate()

    def make_directives(self):
        return {
          "MUL_MATRIX_T3" : UserDirective(self.make_code, "sss"),
        }
        
########################################################################################


class HadamardOpT3A(HadamardMatrixCode):
    """Apply triality element to tags A, B, C

    This is a simplified version of class HadamardOpT3. It
    computes the A part of the vector (multipled by a triality
    element) only.

    Yet to be documented!!!!
    """
        
    def __init__(self, p, verbose = 0):
        """Create an instance of the class.

        Parameters are as in the base class HadamardMatrixCode,
        with the restriction that log_vlength must be even.
        """
        super(HadamardOpT3A, self).__init__(p, 0, verbose)
        self.set_vars()
        self.directives.update(self.make_directives())

    def set_vars(self):
        self.vlen = 2 << self.NO_CARRY
        self.vars.resize(self.vlen)
        

    def reset_vars(self):
        super(HadamardOpT3A, self).reset_vars()
        self.set_vars()


    def external_add(self, v1, v2):
        """Compute  v1 += v2.
   
        Perform addition of the two C variables v1 and v2.  For all 
        components v1[i], v2[i] of of v1, v2
        (stored in v1 and v2) we put

            v1[i] = v1[i] + v2[i], 

        Each result is reduced modulo p. This is a simplified version
        of method external_butterfly().
        """
        if self.FAST_MOD3:
            self.add_mod3(v1, v2, v1)
        else:
            v1.assign(v1 + v2)
            self.reduce_butterfly(v1)

    def main_op(self):
        """yet to be documented!!!"""
        if self.NO_CARRY:
            for i in range(2):
                self.expand_hadamard_extern(self.vars[i], self.vars[i+2])
            self.external_add(self.vars[2], self.vars[3])
        self.external_add(self.vars[0], self.vars[1])
        if self.NO_CARRY:
            self.compress_hadamard_extern(self.vars[0], self.vars[2])
        self.mul_var_pwr2(self.vars[0], -1)




    def load_vector_mul_diagonal(self, array_name, mask):
        """Load vector v from tags A, B, C of rep 196884x.

        Here 'mask' is an integer of type uint_mmv_t wich must be
        0 or (uint_mmv_t)(-1). If mask == -1 then the vector v is 
        multiplied with the dagonal matrix D described in method
        make_diagonal_mask().

        'array_name' must a variable of type uint_mmv_t 
        referring to array A as in method load_vector_direct(). 
        More specifically, A contains a set of componets of
        a vector for the rep 196884x with tag A.
        The function also loads the corresponding sets of 
        components of the vector with tags B and C.
        """
        self.comment(
"""Loading vector from rep 196884x with tags A,B,C
to v[0...2]. Here %s refers to the tag A part. 
Negate v[2] if %s == -1."""
            % (array_name, mask)
        )
        s = "%s = (%s)[%d];\n" % (self.vars[0], array_name,
              24 * self.V24_INTS )
        s += "%s = (%s)[%d] ^ ((%s) & %s);\n" % (self.vars[1], 
                array_name, 48 * self.V24_INTS, mask, 
                self.hex(self.smask(self.P)) ) 
        self.add(s, 2, 2)
        self.comment_vector()  

    def store_vector(self, array_name):
        """Store vector v to tag A, B, C of  of rep 196884x.

        'array_name' must a variable of type uint_mmv_t 
        referring to array A as in method load_vector_direct(). 
        More specifically, A contains a set of componets of
        a vector for the rep 196884x with tag A.
        The function also loads the corresponding sets of 
        components of the vector with tags B and C.
        """
        self.comment(
"""Store vector v[0] to rep 196884x with 
tags A. Here %s refers to the tag A part. """
            % (array_name)
        )
        s = "(%s)[0] = %s;\n" % (array_name, self.vars[0])
        self.add(s, 1)


    def make_code(self, source, mask, dest):
        """Apply triality operation on vector.

        Right multiply the vector of integers mod self.P stored in 
        'src' by t**e, where t is the 64 times 64 triality matrix 
        operating on blocks of the rep 196884x with tag. 

        For e = 1, parameter 'mask' must be 0; for e = 2, parameter
        'mask' must bei (uint_mmv_t)(-1). The part with tag 'A' of
        the result result is stored in 'dest'.

        'source' and 'dest' must be pointers of type uint_mmv_t*,
        'mask' must be a variable of type uint_mmv_t.
        """
        self.reset_vars()
        self.comment(
"""Put dest_A =  (src_B + mask * src_C) / 2   (mod {p})

Here src_B and src_C are the part of a vector of integers 
mod {p} stored in ({src}) with tag B and C, and dest_A is 
the part of a vector of integers mod {p} stored in ({dest}),  
with tag A. Here {mask} must be 0 or -1.

This means that the function computes the part with tag A of
the vector ({dest}) = ({src}) * t**e, where e = 1 - mask.

{src} and {dest} are pointers of type *uint_mmv_t.
Components with tags B, C referred by ({src}) 
are processed, one integer of type uint_mmv_t
for each tag.

""".format(p = self.P, src = source, mask = mask, dest = dest)
        )
        self.load_vector_mul_diagonal(source, mask)
        self.main_op()
        self.store_vector(dest)
        self.comment_statistics()
        return self.generate()

    def make_directives(self):
        return {
          "MUL_MATRIX_T3A" : UserDirective(self.make_code, "sss"),
        }
        

########################################################################################
# Summarizing the table classes given above
########################################################################################


class Tables:
    def __init__(self, **kwds):
        p = kwds.get('p', 3)
        self.tables = {}
        self.directives = {}
        table_classes = [HadamardOpT64(p), HadamardOpT3(p),
                         HadamardOpT3A(p)]
        for t in table_classes:
            self.tables.update(t.tables)
            self.directives.update(t.directives)

class MockupTables:
    tables = {}
    directives = {}
    def __init__(self, **kwds):
        pass







if __name__ == "__main__":
    for p in (3, 7, 127):
        cg = HadamardOpT64(p)
        cg.swap_parities()
        cg.comment_statistics()
        print("\n", cg.generate())
