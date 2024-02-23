from __future__ import absolute_import, division, print_function
#from __future__ import  unicode_literals


import sys
import os
from numbers import Integral

from mmgroup.bitfunctions import bitparity, bitweight
from mmgroup.dev.hadamard.hadamard_codegen import HadamardMatrixCode
from mmgroup.dev.hadamard.hadamard_codegen import C_UintVarPool
from mmgroup.generate_c import UserDirective, UserFormat



class C_UintVarArray(C_UintVarPool):
    """Model a (usually short) local array of integer variables

    The entries of the array can be treated in the same way as the
    variable of a pool, where the pool is an instance of class
    C_UintVarPool.
    """
    def __init__(self, c_type = "uint_mmv_t", name = "a[%d]", length = 0):
        """Create an array of integer variables of type 'ctype'.

        The array will contain 'length'  variables. These
        variables are instances of class C_UintVar and they are
        obtained as self[0], ..., self[length - 1].

        The names of the variables are name % 0, name % 1, name % 2, ..
        where name is given by parameter 'name'. So for obtaining
        a[0], a[1], a[2], ... parameter 'name' should be "a[%d]"
        """
        super(C_UintVarArray, self).__init__(c_type, name, length)
        self.array_name = self.name[:self.name.index("[")]
             

    def declare(self):
        """Declare the array in C

        The function returns a string containing that declaration.
        """
        if len(self) == 0:
            return ""
        name = self.array_name
        return "%s %s[%d];\n" % (self.c_type, name, len(self))

    def load_pool(self, start,  pool):
        """Load a 'pool' from the array of variables
 
        Here 'pool' is an instance of class C_UintVarPool
        modelling a pool of integer vriables. The function loads
        all 'official' variables of the pool from the array, 
        starting at the index 'start'.
        """
        if isinstance(start, Integral):
            indices = map(str, range(start, start + len(pool)))
        else:
            indices = ["(%s) + %d" % (start, i) for i in range(len(pool))]
        name = self.array_name
        for i, index in enumerate(indices):
            pool[i].assign("%s[%s]" %  (name, index)) 

    def store_pool(self, pool, start):
        """Store a 'pool' to the array of variables
 
        Here 'pool' is an instance of class C_UintVarPool
        modelling a pool of integer vriables. The function stores
        all 'official' variables of the pool to the array, starting 
        at the index 'start'.
        """
        if isinstance(start, Integral):
            indices = map(str, range(start, start + len(pool)))
        else:
            indices = ["(%s) + %d" % (start, i) for i in range(len(pool))]
        name = self.array_name
        for i, index in enumerate(indices):
            s = "%s[%s] = %s;\n" %  (name, index, pool[i])
            self.context.add(s, 1)

    def _bad(self, *args, **kwds):
        raise NotImplementedError

    temp = _bad
    xch = _bad


class HadamardOpXi64(HadamardMatrixCode):
    """Code for multiplying vector by a modified Hadamard matrix

    This class supports the operation of the  element l**e,
    e = 1,2  on a block of a vector v of length 64 with tag 'Y' 
    or 'Z' in  representation  196884x  of the monster group as 
    described in [Seys19].

    The group element l is called xi in [Seys19], but we cannot 
    use greek letters as tags in python or C.  We have xi**3 == 1. 
    Here a block of xi is a 64 times 64 matrix acting on a block 
    of v by right multiplication. For xi acting on a components
    of v with tags 'Y' and 'Z', we can partition xi  into blocks 
    of identical 64  times 64  matrices xi64, where each block 
    xi64  has shape xi64 = DD4 * M * DD16. Similarly, we can 
    partition  xi**2  into identical blocks xi64**2, with 
    xi64**2 = DD16 * M * DD4. Here  DD16 and DD4 are certain 
    diagonal matrices, and M is a Hadamard-like matrix as 
    described below. 

    We have DD16 = kron(D16, U4), DD4 = kron(U16, D4) with 
    D16 the diagonal matrix 
    ( -1, -1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1)
    and D4 the diagonal matrix (-1, 1, 1, 1), and U<n> the
    n times n unit matrix. kron() is the kronecker product
    as given by the numpy function numpy.kron().

    M is equal to kron(M4, M2), where M<n> is the parity-adjusted 
    2**n times 2**n Hadamard matrix  for an even number n as 
    defined in the header of module define_matrices.py. This class 
    supports the generation of C code for multiplying a vector v
    of integers modulo p by a matrix xi64 or xi64**2.

    Notation and implementation of vectors is as in the base 
    class HadamardMatrixCode, which supports right-multiplication 
    of a vector by a Hadamard  matrix.

    More precisely, M acts on each components of a vector labeled
    (U, i0 + i, j0 + j), i = 0,...,15, j = 0,..3, with column
    4 * i + j of matrix M acting on (U, i0 + i, j0 + j). Here the
    tag U must be 'Y' or 'Z', i0 must be divisible by 16, and j0
    must be divisible by 4.

    See  [Seys19], section 9 and module define_matrix.py for
    background.

    The code generated by method self.make_code() multiplies 16
    consecutive rows (U, i0 + i, j), i = 0,...,15, j = 0,..23,
    with 6 identical matrix blocks, where each of these blocks
    is either xi64 or x64**2. The code block generated by that
    method is rather large. On the other hand, this is one of the
    most time-critical operations in the monster group. Here 
    dividig a row (U, i, j), j = 0,...,23 into blocks of four
    entries is also quite time critical, requiring a fair amount 
    of optimization.
    """
    _PRE_MDIAG16 = [
            -1, -1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1
    ]
    MDIAG16 = [ int(i < 0) for i in _PRE_MDIAG16 ]

    def __init__(self, p, verbose = 0):
        """Create an instance of the class.

        Parameters are as in the base class HadamardMatrixCode.
        p is the modulus, which must be 2**k - 1, 2 <= k <= 8.
        """
        super(HadamardOpXi64, self).__init__(p, 0, verbose)
        # Next set sizes for Hadamard-like matrices
        self.set_sizes()
        # Next set the default position for method expand_hadamard    
        self.free_cy_pos = -1  
        # Next set positions where Hadamard operation muse  be done
        self.hadamard_operations = 3 + (15 << self.LOG_INT_FIELDS)
        # make tables and directives for code generation
        self.tables.update(self.make_tables())
        self.directives.update(self.make_directives())


    def set_sizes(self):
        array_size, log_vars = 0, 4
        if self.V24_INTS_USED == 1:
            array_size, log_vars = 24, 3
        elif self.V24_INTS_USED == 2:
            cy = self.NO_CARRY
            array_size, log_vars = 24 << cy, 3 + cy
        self.reset_vars(self.LOG_INT_FIELDS + log_vars)
        self.array = C_UintVarArray(length = array_size) 
        self.array.set_context(self)   

    def load_var(self, source, i, ofs, p_mask, dest, tmp = None):
        if tmp is None:
            tmp = self.vars.temp()
        i1 = 15 - i if bitparity(i) else i
        index = (i1 << self.LOG_V24_INTS) + ofs
        mindex = 2 * self.MDIAG16[i1] 
        dest.assign("%s[%d] ^  %s[%s]" % (source, index, p_mask, mindex))
        sh = self.FIELD_BITS
        swap_mask = self.smask(self.P, range(1, 24, 4))
        tmp.assign((dest ^ (dest >> sh)) & swap_mask)
        dest.assign_xor(tmp | (tmp << sh))


    def store_var(self, source, p_mask, dest, i, ofs, tmp = None):
        if tmp is None:
            tmp = self.vars.temp()
        index = (i << self.LOG_V24_INTS) + ofs
        mindex = 2 * self.MDIAG16[i] + 4
        s = ""
        if self.V24_INTS_USED < 3 and ofs == self.V24_INTS_USED - 1:
            mask1 = self.hex(self.smask(self.P, range(24 % self.INT_FIELDS)))
            s += "%s[%d] = (%s ^  %s[%s]) & %s;\n" % (
                dest, index, source, p_mask, mindex, mask1)
            self.n_operations += 1
        else:
            s += "%s[%d] = %s ^  %s[%s];\n" % (
                dest, index, source, p_mask, mindex)
        self.add(s, 1, 1)

    def load_all_v24_1(self, source, p_mask):
        t = self.vars.temp(0)
        for i in range(16):
            v0 = self.vars[i]
            self.load_var(source, i, 0, p_mask, v0, t)
        self.hadamard_operations = 3 + (15 << 5)

    def store_all_v24_1(self, p_mask, dest):
        t = self.vars.temp(0)
        for i in range(16):
            v0 = self.vars[i]
            #self.mul_var_pwr2(v0, -3) 
            self.store_var(v0, p_mask, dest, i, 0, t)   


    def load_all_v24_1_no_cy(self, source, p_mask):
        t, t1 = self.vars.temp(0), self.vars.temp(1)
        sh1 = 8 * self.FIELD_BITS
        sh2 = 16 * self.FIELD_BITS
        m0, m1, m2 = [self.smask(self.P, range(8*i, 8*i+8)) for i in range(3)]
        for i in range(8):
            v0, v1, v2 = self.vars[i], self.array[i+8], self.array[i+16]
            self.load_var(source, 2*i, 0, p_mask, v0, t1)
            self.load_var(source, 2*i+1, 0, p_mask, t, t1)
            v1.assign(((v0 >> sh1) & m0) | ((t << sh1) & m2))
            v2.assign(((v0 >> sh2) & m0) | ((t) & m2))
            v0.assign(((v0) & m0) | ((t << sh2) & m2))
        self.free_cy_pos = 3
        self.hadamard_operations = 3 + (15 << 4)
            
    def store_all_v24_1_no_cy(self, p_mask, dest):
        t, t1 = self.vars.temp(0), self.vars.temp(1)
        sh1 = 8 * self.FIELD_BITS
        sh2 = 16 * self.FIELD_BITS
        m0, m1, m2 = [self.smask(self.P, range(8*i, 8*i+8)) for i in range(3)]
        for i in range(8):
            v0, v1, v2 = self.array[i], self.array[i+8], self.vars[i]
            t.assign((v0 & m0) | ((v1 << sh1) & m1) | (v2 << sh2))
            #self.mul_var_pwr2(t, -3)           
            self.store_var(t, p_mask, dest, i+i, 0, t1)   
            v0.assign(((v0 >> sh2) & m0) | ((v1 >> sh1) & m1) | (v2 & m2))
            #self.mul_var_pwr2(v0, -3)           
            self.store_var(v0, p_mask, dest, i+i+1, 0, t1)


    def load_all_v24_2(self, source, p_mask):
        t, t1 = self.vars.temp(0), self.vars.temp(1)
        sh1 = 8 * self.FIELD_BITS
        m0, m1 = [self.smask(self.P, range(8*i, 8*i+8)) for i in range(2)]
        for i in range(8):
            v0, v1, v2 = self.vars[i], self.array[i+8], self.array[i+16]
            self.load_var(source, 2*i, 0, p_mask, v0, t1)
            self.load_var(source, 2*i+1, 0, p_mask, t, t1)
            v1.assign(((v0 >> sh1) & m0) | ((t) & m1))
            v0.assign(((v0) & m0) | ((t << sh1) & m1))
            self.load_var(source, 2*i, 1, p_mask, v2, t1)
            self.load_var(source, 2*i+1, 1, p_mask, t, t1)
            v2.assign(((v2) & m0) | ((t << sh1) & m1))
        self.hadamard_operations = 3 + (15 << 3)
            

    def store_all_v24_2(self, p_mask, dest):
        t, t1 = self.vars.temp(0), self.vars.temp(1)
        sh1 = 8 * self.FIELD_BITS
        m0, m1 = [self.smask(self.P, range(8*i, 8*i+8)) for i in range(2)]
        for i in range(8):
            v0, v1, v2 = self.array[i], self.array[i+8], self.vars[i]
            t.assign((v0 & m0) | ((v1 << sh1) & m1)) 
            #self.mul_var_pwr2(t, -3)           
            self.store_var(t, p_mask, dest, 2*i, 0, t1)   
            v0.assign(((v0 >> sh1) & m0) | ((v1) & m1))
            #self.mul_var_pwr2(v0, -3)           
            self.store_var(v0, p_mask, dest, 2*i+1, 0, t1)
            #self.mul_var_pwr2(v2, -3)           
            t.assign((v2) & m0)
            self.store_var(t, p_mask, dest, 2*i, 1, t1)
            v2.assign((v2 >> sh1) & m0)
            self.store_var(v2, p_mask, dest, 2*i+1, 1, t1)


    def load_all_v24_2_no_cy(self, source, p_mask):
        t = self.vars.temp(0)
        sh1 = 8 * self.FIELD_BITS
        m0, m1 = [self.smask(self.P, range(8*i, 8*i+8)) for i in range(2)]
        for i in range(16):
            v0, v1, v2 = self.vars[i], self.array[i+16], self.array[i+32]
            self.load_var(source, i, 0, p_mask, v0, t)
            v1.assign((v0 >> sh1) & m0)
            v0.assign((v0) & m0)
            self.load_var(source, i, 1, p_mask, v2, t)
        self.hadamard_operations = 3 + (15 << 4)
        self.free_cy_pos = 3

    def store_all_v24_2_no_cy(self, p_mask, dest):
        t = self.vars.temp(0)
        sh1 = 8 * self.FIELD_BITS
        m0, m1 = [self.smask(self.P, range(8*i, 8*i+8)) for i in range(2)]
        for i in range(16):
            v0, v1, v2 = self.array[i], self.array[i+16], self.vars[i]
            v0.assign((v0 & m0) | ((v1 << sh1) & m1))
            #v2.assign(v2)
            #self.mul_var_pwr2(v0, -3)           
            #self.mul_var_pwr2(v2, -3)           
            self.store_var(v0, p_mask, dest, i, 0, t)
            self.store_var(v2, p_mask, dest, i, 1, t)

    def small_load_all(self, source, p_mask):
        if self.V24_INTS_USED == 1:
            if self.NO_CARRY:
                self.load_all_v24_1_no_cy(source, p_mask)
            else:
                raise NotImplmentedError("V24_INTS_USED=1")
        elif self.V24_INTS_USED == 2:
            if self.NO_CARRY:
                self.load_all_v24_2_no_cy(source, p_mask)
            else:
                self.load_all_v24_2(source, p_mask)
        else:
            raise NotImplmentedError("V24_INTS_USED > 2")
  

    def small_store_all(self, p_mask, dest):
        if self.V24_INTS_USED == 1:
            if self.NO_CARRY:
                self.store_all_v24_1_no_cy(p_mask, dest)
            else:
                raise NotImplmentedError("V24_INTS_USED=1")
        elif self.V24_INTS_USED == 2:
            if self.NO_CARRY:
                self.store_all_v24_2_no_cy(p_mask, dest)
            else:
                self.store_all_v24_2(p_mask, dest)
        else:
            raise NotImplmentedError("V24_INTS_USED > 2")
         

    def add_mask(self, p_mask, value):
        d = {1: "++", -1: "--"}
        self +=  "%s%s;\n" %(d[value], p_mask)
        
    def label(self, index):
        return "l_mmv%d_op_l64_%d"  % (self.P, index)


    def very_small_op(self, source, p_mask, dest):
        self.reset_vars(4+5)
        self.load_all_v24_1(source, p_mask)
        self.comment_statistics()
        self.hadamard_op(self.hadamard_operations, shift = -3)
        self.comment_statistics()
        self.store_all_v24_1(p_mask, dest)
        self += "%s += %d;\n" % (source, 16 * self.V24_INTS)
        if dest != source:
            self += "%s += %d;\n" % (dest, 16 * self.V24_INTS)
        self.comment_statistics()


        
    def small_op(self, source, p_mask, dest):
        self.small_load_all(source, p_mask)
        l1, l2 = self.label(1), self.label(2)
        self.comment_statistics()
        a = self.array
        self += "i = 0;\ngoto %s;\n" % l2
        self += l1 + ":\n"
        a.store_pool(self.vars, "i")
        self += "i += %d;\n" % len(self.vars)
        a.load_pool("i",  self.vars)
        self += l2 + ":\n"
        self.expand_hadamard(self.free_cy_pos)   
        self.hadamard_op(self.hadamard_operations, shift = -3, 
            compress = True)
        #self.compress_hadamard()   
        self.comment_statistics()
        i_end = 2 * len(a) // 3
        self += "if (i < %d) goto %s;\n" % (i_end, l1)
        self.small_store_all(p_mask, dest)
        self += "%s += %d;\n" % (source, 16 * self.V24_INTS)
        if dest != source:
            self += "%s += %d;\n" % (dest, 16 * self.V24_INTS)
        self.comment_statistics()

    def large_op_on_int(self, source, p_mask, dest):
        t = self.vars.temp(0)
        for i in range(16):
            self.load_var(source, i, 0, p_mask, self.vars[i], t)    
        self.hadamard_op(self.hadamard_operations, shift = -3)
        for i in range(16):
            #self.mul_var_pwr2(self.vars[i], -3)           
            self.store_var(self.vars[i], p_mask, dest, i, 0, t)
        self += "%s++;\n" % source
        if dest != source:
            self += "%s++;\n" % dest


    def large_op(self, source, p_mask, dest):
        USED = self.V24_INTS_USED
        self += "for (i = 0; i < %d; ++i) {\n" % USED
        self.large_op_on_int(source, p_mask, dest)
        self += "}\n"
        for i in range(16):
            for j in range(self.V24_INTS - USED):
                k = i * self.V24_INTS + j
                self += "%s[%d] = 0;\n" % (dest, k)
        k = 16 * self.V24_INTS - USED
        self += "%s += %d;\n" % (source, k)
        if dest != source:
            self += "%s += %d;\n" % (dest, k)


    def additional_declare(self):
        self += self.array.declare()
        self += ("uint_fast32_t i;\n\n")

    def make_code(self, source, p_mask, dest):
        """Apply operation on vector.

        YET TO BE DOCUMENTED!!!

        Right multiply the vector of integers mod self.P stored in 
        'source' by t**e, where t is the 64 times 64 triality matrix 
        operating on blocks of the rep 196884x with tag. 


        Let A_T be the table created by function  self.exp_table(). 
        'p_mask' must be a pointer of type uint_mmv_t* with value
        p_mask = &A_T[e-1] for e = 1, 2.
         
        'source' and 'dest' must be pointers of type uint_mmv_t*. 
        """
        self.reset_vars()
        if not self.FAST_MOD3:
            self.additional_declare()
        self.comment(
"""TODO: write comment!!!
""".format(p = self.P, source = source, p_mask = p_mask, dest = dest)
        )
        if self.V24_INTS_USED < 3:
            if self.FAST_MOD3:
                self.very_small_op(source, p_mask, dest)
            else:
                self.small_op(source, p_mask, dest)
        else: 
            self.large_op(source, p_mask, dest)  
        return self.generate()

    def exp_table(self):
        d4 = self.smask(self.P, range(0, 24, 4))
        d16 = self.smask(self.P, (1 << 24) - 1)
        m = [0] * 8
        m[0:8:2] = [d4, d4, 0, d16] 
        m[1:8:2] = [0, d16, d4, d4] 
        return m


    def make_tables(self):
        return {
            "TABLE_MUL_MATRIX_XI64" : self.exp_table()
        }


    def make_directives(self):
        return {
            "MUL_MATRIX_XI64" : UserDirective(self.make_code, "sss"),
        }
        


class HadamardOpXi16(HadamardMatrixCode):

    def __init__(self, p, verbose = 0):
        """Create an instance of the class.

        Parameters are as in the base class HadamardMatrixCode.
        p is the modulus, which must be 2**k - 1, 2 <= k <= 8.
        """
        super(HadamardOpXi16, self).__init__(p, 0, verbose)
        # Next set sizes for Hadamard-like matrices
        self.reset_vars(self.LOG_INT_FIELDS + 2)
        # Next set positions where Hadamard operation muse  be done
        self.hadamard_operations = 3 + (3 << self.LOG_INT_FIELDS)
        # set masks for diagonal operation
        self.MASKS = [self.smask(self.P, range(0,24,4))]
        self.MASKS.append(~self.MASKS[0] & self.smask(self.P, range(24)))
        # make tables and directives for code generation
        #self.tables.update(self.make_tables())
        self.directives.update(self.make_directives())
        

    def load_var(self, source, i, mask, dest):
        tmp = self.vars.temp()
        i1 = [0, 2, 1, 3][i]
        index = (i1 << self.LOG_V24_INTS)
        m = self.hex(self.MASKS[i1 == 0])
        dest.assign("%s[%d] ^  (%s & %s)" % (source, index, m, mask))
        sh = self.FIELD_BITS
        swap_mask = self.smask(self.P, range(1, 24, 4))
        tmp.assign((dest ^ (dest >> sh)) & swap_mask)
        dest.assign_xor(tmp | (tmp << sh))


    def store_var(self, source, mask, dest, i):
        index = (i << self.LOG_V24_INTS)
        m = self.hex(self.MASKS[i == 0])
        s =  "%s[%d] = %s ^ (%s & %s);\n" % (
            dest, index, source, mask, m)
        self.add(s, 1, 2)



    def op(self, source, mask, dest):
        self.complement_variable(mask)
        for i in range(4):
            self.load_var(source, i, mask, self.vars[i])
        self.hadamard_op(self.hadamard_operations)
        self.complement_variable(mask)
        for i in range(4):
            self.mul_var_pwr2(self.vars[i], -2)           
            self.store_var(self.vars[i], mask, dest, i)  
        self.comment_statistics()
        self += "%s++;\n" % source
        if dest != source:
            self += "%s++;\n" % dest


    def make_code(self, source, mask, dest):
        """Apply operation on vector.

        YET TO BE DOCUMENTED!!!

        Right multiply the vector of integers mod self.P stored in 
        'source' by t**e, where t is the 64 times 64 triality matrix 
        operating on blocks of the rep 196884x with tag. 

        For e = 1, parameter 'mask' must be 0; for e = 2, parameter
        'mask' must bei (uint_mmv_t)(-1). The result is stored
        in 'dest'.

        'source' and 'dest' must be pointers of type uint_mmv_t*,
        'mask' must be a variable of type uint_mmv_t.
        """
        self.reset_vars()
        USED = self.V24_INTS_USED
        if USED > 1:
            self += "uint_fast32_t i;\n"
        self.comment(
"""TODO: write comment!!!
""".format(p = self.P, source = source, mask = mask, dest = dest)
        )
        if USED > 1:
            self += "for (i = 0; i < %d; ++i) {\n" % USED
        self.op(source, mask, dest)
        if USED > 1:
            self += "}\n" 
        if USED == 2:
            for i in range(-1, 7, 2):
                mask1 = self.hex(self.smask(self.P, range(8)))
                self.add("%s[%s] &= %s;\n" % (dest, i, mask1), 1, 1)
        elif USED >= 3:
            for i in range(4):
                for j in range(self.V24_INTS - USED):
                    k = i * self.V24_INTS + j
                    self += "%s[%d] = 0;\n" % (dest, k)
        k = 4 * self.V24_INTS - USED
        self += "%s += %d;\n" % (source, k)
        if dest != source:
            self += "%s += %d;\n" % (dest, k)
        return self.generate()

    def make_directives(self):
        return {
            "MUL_MATRIX_XI16" : UserDirective(self.make_code, "sss"),
        }
  



########################################################################################
# Summarizing the table classes given above
########################################################################################


class Tables:
    def __init__(self, **kwds):
        p = kwds.get('p', 3)
        self.tables = {}
        self.directives = {}
        table_classes = [HadamardOpXi64(p), HadamardOpXi16(p)]
        for t in table_classes:
            self.tables.update(t.tables)
            self.directives.update(t.directives)

class MockupTables:
    tables = {}
    directives = {}
    def __init__(self, **kwds):
        pass







if __name__ == "__main__":
    for p in (3,7, 15, 127):
        print( """typedef unsigned int uint_fast32_t;
typedef unsigned long long  uint_mmv_t;
void mm_f%d(uint_mmv_t *source, uint_mmv_t mask, uint_mmv_t *dest)
""" % p )
        cg = HadamardOpXi16(p)
        print("\n", cg.make_code("source", "mask", "dest"))
