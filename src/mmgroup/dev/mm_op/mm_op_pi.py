from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import collections
import re
import warnings
from numbers import Integral
from array import array
from random import randint



from mmgroup.dev.generate_c.make_c_tables import c_snippet, TableGenerator
from mmgroup.dev.generate_c.make_c_tables import make_table
from mmgroup.dev.generate_c.generate_functions import UserDirective
from mmgroup.dev.generate_c.generate_functions import UserFormat



from mmgroup import mat24

from mmgroup.dev.mm_basics.mm_basics import MM_Basics, MM_Op
from mmgroup.dev.mm_basics.mm_tables import MM_OctadTable


###########################################################################
# Standard permutations of the 24 entries of a small array 
###########################################################################

class Perm24_Standard(MM_Op):
    """Permute 24 small integers   deprecated!!!!!!!!!!!!!!!!!!!!!

    A modulus p for the small integers is given in the constructor of the 
    class.  An variable of type uint_mmv_t stores self.INT_FIELDS small
    integers, with each of them stored in a field of self.FIELD_BITS bits.
    Sizes self.INT_FIELDS and self.FIELD_BITS are calculated by a base
    class of this class depending on modulus p. The 24 small integers
    to be permuted are stored in a contiguous array of integers of
    type uint_mmv_t.

    Given a permutation P of 24 integers 0,...,23 as an array of 8-bit 
    integers, method self.prepare(P, M) generates code  for computing
    array M of 48 8-bit integers from array P. That code is optimized
    for computing a permutation of 24 small integers modulo p. Given a
    source permutation and a fast permutation array M, method
    self.expression() essentially returns a list of code snippets for 
    computing all components of the the array of type unit_mmv_t[] 
    representing the permuted array 'source'. Using the code generator
    in module make_c_tables you may generate C code for storing the
    permuted array 'source' in an array 'dest' with

    // %%FOR i, expr in PERM24_EXPR(perm_fast, source)
    dest[{i}] = {expr};
    // %%END FOR

    Method self.expression() also allows to negate all small integers 
    before permuing them. 

    The recommened way to use this class in the code generator is:

    // The following array contains the permutation to be performed
    uint8_t permutation[24]

    // The following array will store an version of that permutation 
    // optimized for fast operation on a vector of 24 small integers: 
    uint8_t perm_fast[48];

    // Prepare permutation for fast operation
    // %%PERM24_PREPARE "permutation", "perm_fast"

    // The following code permutes the small array 'source' and stores
    // the result to 'dest'. The result is negated if bit 0 of 'sign'
    // is set.
    uint_mmv_t *soure, *dest, sign;

    // %%FOR i, expr in PERM24_EXPR("source", "perm_fast", "sign")
    dest[{i}] = {expr};
    // %%END FOR
    """
    def __init__(self, p):
        """Initialise for calulations with small integers modulo p

        p+1 must be a power of two. Calculations modulo p are described 
        in more detail in the base classes of this class.
        """
        super(Perm24_Standard, self).__init__(p)

    def prepare(self, perm, perm_fast):
        """Generates code for the precomputation needed for a permutation

        Parameter 'perm' is a permutation coded as an array of 24
        integers with values 0,...,23. Parameter 'perm_fast' is an 
        array of 48 integers of at least 8 bit computed by the function. 
        Array 'perm_fast' encodes the inverse permutation 'inv' of 'perm'.
        Here entries 2*i and 2*i+1 encode the result of the inverse
        permutation 'inv' of 'perm' applied to i as follows:

        If 'source' is an array of type uint_mmv_t[] that holds the source
        data and 'dest' is an array of type uint_mmv_t[] that shall hold 
        the result, the small integer in field of source[2*i] starting at 
        bit [2*i+1] is mapped to the i-th small integer in dest.   
  
        The function generates a C block in curly brackets.  It uses the
        identifiers 'i', 'j' inside the block.
        """
        return  self.snippet("""{
    uint_fast8_t i, j;
    for(i = 0; i < 24; ++i) {
        j = {perm}[i];
        {result}[j+j] = (i >> {LOG_INT_FIELDS});
        {result}[j+j+1] = (i & {hex:INT_FIELDS-1}) << {LOG_FIELD_BITS};
    }
}"""    , perm = perm, result = perm_fast)


    def sign_term(self, sign, i):
        """Return sign term contributed by 'sign' for destination[i]"""
        if i == 0:
            return self.snippet(
                "({0} = (-({0} & {hex:1})) & {smask:P, range(24)})",
                sign)
        if (i == 1 and self.V24_INTS_USED == 2):
            return self.snippet( "({0} & {smask:P, range(8)})", sign)
        return sign


    def expressions(self, source, perm_fast, sign = None):
        """Generate code for a permutation of 24 small integers

        Given a vector of 24 small integers coded in the bit fields
        of an array 'source' of type uint_mmv_t[], and an array 
        'perm_fast' encoding an optimized version of a fast
        permutation of 24 entries as returned by method perpare(),
        this method essentially returns a list expressions that make 
        up the permuted vector 'source'.

        If you want to store the permuted vector 'source' in an array
        array 'dest' of the same type as the array containing 'source' 
        you may simply code:

            dest[0] = expr[0]; dest[1] = expr[1]; ...

        where expr[0], expr[1], ... is the list of expressions
        returned by this method.

        Since this method is designed for use in the code generator,
        it actually returns a lsit of pairs

            [ (0, expr[0]), (1, expr[1]), ...,  (N, expr[N]) ],
        
        where N = self.V24_INTS_USED - 1.

        If an addititional optional argument 'sign' is passed to
        the method, this must be an lvalue of type uint_mmv_t.
        Then all 24 small integers are negated before being permuted 
        if bit 0 of 'sign' is set. 
        """
        def inv_perm_term(j):
            """Return contribution of inverse_permuation[j]""" 
            return self.snippet(
               "(((({0})[{1}[{int:2*j}]] >> {1}[{int:2*j+1}]) & {P})" + (
               " << {int:j*FIELD_BITS % INT_BITS})"), 
               source, perm_fast, j = j
            )
        results = []
        INT_FIELDS = self.INT_FIELDS
        for i in range(self.V24_INTS_USED):
            terms = [ inv_perm_term(j) for j in
                range(INT_FIELDS * i, min(24, INT_FIELDS * (i+1))) ]
            result = "\n  | ".join(terms)
            if sign:
                result = "(%s)\n  ^ %s" % (result, self.sign_term(sign,i))
            results.append((str(i), result))
        return results            


    def cond_neg(self, sign):
        """Generate code for conditional negation of 24 small integers

        This is equivalent to

            self.expressions(source, perm_fast, sign)

        with  perm_fast  equal  to the identity permutation and source
        equal to the zero vector. 
        """
        results = []
        for i in range(self.V24_INTS_USED):
            result = self.sign_term(sign,i)
            results.append((str(i), result))
        return results            



    def tables(self):
        return {}
        return {
            "PERM24_EXPR": UserFormat(self.expressions, "sss"),
            "PERM24_COND_NEG": UserFormat(self.cond_neg, "ss"),
        }
 
    def directives(self):
        return {}
        return {
           "PERM24_PREPARE" : UserDirective(self.prepare, "ss"),
        }


###########################################################################
# Permutations of the 24 entries of a small array with a Benes network
###########################################################################


class Perm24_Benes(MM_Op):
    """Permute 24 small integers with a Benes network

    This is a modified version of class Perm24_Standard for permuting
    vectors of small integers.

    A modulus p for the small integers is given in the constructor of the 
    class.  An variable of type uint_mmv_t stores self.INT_FIELDS small
    integers, with each of them stored in a field of self.FIELD_BITS bits.
    Sizes self.INT_FIELDS and self.FIELD_BITS are calculated by a base
    class of this class depending on modulus p.

    Here a Benes network may be used in case self.INT_FIELDS >= 8 only,
    i.e. an integer of type uint_mmv_t muts be able to store at least
    8 small integers. The 24 small integers to be permuted will be 
    copied into a set of at most 3 local variables v0, v1, v2.

    The Benes network for performing the permutation can be computed by
    the C function mat24_perm_to_net(). The result of this function is
    a array B of type uint32_t[9] describing a Benes network for a 
    permutation. From array B we must compute another array M type
    uint_mmv_t[] and length at most 21, with the exact length given by
    self.n_benes_operations(). Array M is a more suitable mask for
    permuting fields of bitlength self.FIELD_BITS, see method
    self.iter_shift_op() for details.
 
    Method self.prepare(benes_net, result) generates code for computing
    array M from array B. Here benes_net is the name of a variable
    referring to the input array B and  result is the name of a variable
    referring to the output array M.

    The recommended way to use this class in the code generator is:

    // The following array contains the permutation to be performed
    uint8_t permutation[24];

    // The following code prepares the permutation operation
    uint32_t benes_net[9];     // Benes network
    // The following mask is used by the actual permutation code
    uint_mmv_t benes_mask[{PERM24_BENES_MASK_LEN}]; 
    // Compute Benes network
    mat24_perm_to_net(permutation, benes_net);
    // Prepare mask array from Benes network
    // %%PERM24_BENES_PREPARE "benes_net", "benes_mask"

    // The following code permutes the small array 'source' and stores
    // the result to 'dest'. The result is negated if bit 0 of 'sign'
    // is set.
    uint_mmv_t *soure, *dest, sign;

    // Declare temporary variables t0, t1,...
    // %%PERM24_BENES_DECLARE "t"
    // Load 'source' to temporary variables
    // %%PERM24_BENES_LOAD "source"
    // Permute and possibly negate data in temp. variables
    // %%PERM24_BENES_PERMUTE "benes_mask", "sign"
    // Store temporary variables to 'dest'
    // %%PERM24_BENES_STORE "dest"

    Some details of the internal operation of the benes network are
    described in method iter_shift_cases.
    """

    BENES_SHIFTS = [1,2,4,8,16,8,4,2,1]
    permute_comment = """// Permute the 24 small integers in '{vars}' 
// using the Benes network given by '{mask}'. All small   
// integers are negated if bit 0 of 'sign' is set.
"""

    def __init__(self, p, prefer_benes_net = True):
        """Initialise for calulations with small integers modulo p

        p+1 must be a power of two. Calculations modulo p are described 
        in more detail in the base classes of this class.

        If 'prefer_benes_net' then self.use_benes_net is set True
        whenever Benes networks may be used for a given p.
        Ohterwise self.use_benes_net is set True only if using 
        a Benes network is an obvious advantage compared to using
        the methods of class Perm24_Standard.
        """
        super(Perm24_Benes, self).__init__(p)
        self.prep_table = []
        self.use_benes_net = (
            self.V24_INTS_USED <=  2 + bool(prefer_benes_net) )
        if self.V24_INTS_USED <= 3:
            self.prep_table = self.make_prep_table()
            #print("Benes table:", list(map(hex, self.prep_table)))

    def not_supported(self, *args, **kwds):
        err = "Benes network not supported for p=%d on a %d-bit system"
        raise TypeError(err % (self.P, self.INT_BITS)) 

    def n_benes_operations(self):
        """Return No of computation steps required for Benes network"""
        return len(self.prep_table)

    def iter_shift_cases(self):
        """Yield pairs (layer, offset) to be used for shift opration

        Given a permutation  p1 of 24 entries, the C function 
        mat24_perm_to_net() computes a Benes network with 9 layers for 
        that permutation. Each layer is  represented as an integer which 
        encodes an array of 24 bits. If bit i is in layer l of the 
        Benes network is set, then entry i must be exchanged with entry
        i + self.BENES_SHIFTS[l] in layer i. The first layer is layer 0.

        In layers 3, 4 and 5 bits are set at positions i < 8 only. 

        In the context of a small permutation, the 24 entries of a
        small array are stored in self.V24_INTS_USED integers of type
        uint_mmv_t, with 1 <= self.V24_INTS_USED <= 3. The integer at 
        offset y contains entries 

            y * INT_FIELDS, ..., y * INT_FIELDS + INT_FIELDS - 1 .

        This function yields pairs (layer, offset) describing a step.
        Steps are given in the order in which the Benes network will be 
        processed. Here offset = y in a step means that the integer at 
        offset y must be processed together with the integer at offset 
        y + k  for k = floor(self.BENES_SHIFTS[layer] / INT_BITS).
        This means that up to INT_FIELDS entries of a layer are
        processed simultaneously in one step. 

        The order of the operations is chosen so that we obtain a 
        good locality.
        """
        if self.V24_INTS_USED == 1: 
            for i in range(9): yield i, 0
        elif self.V24_INTS_USED == 2: 
            for ofs in [0,1]:
                for i in range(4 - ofs): yield i, ofs
            yield 4, 0
            for ofs in [0,1]:
                for i in range(5 + ofs, 9): yield i, ofs
        elif self.V24_INTS_USED == 3: 
            for ofs in [0,1,2]:
                for i in range(3): yield i, ofs
            for i in range(3,6): yield i, 0
            for ofs in [0,1,2]:
                for i in range(6, 9): yield i, ofs
        else:
            raise ValueError("Illegal parameter for Benes network")



    def make_prep_table(self):
        """Return the table required by method prepare()

        This is a compressed version of the table returend by method
        iter_shift_cases() with each pair (layer, offset) stored in
        on byte. Here the lower four bits of that byte store the value 
        'layer' and the higher four bits of that byte store the value
        offset/16. In case self.V24_INTS_USED == 3 the higher four 
        bits store the value offset/8 instead.
        """ 
        d = self.V24_INTS_USED == 3
        return [  (offset << (self.LOG_INT_FIELDS + d)) + layer
              for layer, offset in self.iter_shift_cases()]  

    def prep_declare(self):
        """Delare local and static variables for self.prepare()"""
        s = """uint_mmv_t tmp; 
uint_fast8_t i;
"""
        if self.V24_INTS_USED > 1:
            s +=  self.snippet("""static uint8_t tbl[] = {
    // %%TABLE table_prepare_perm24, uint8_t
};
"""     , table_prepare_perm24 = self.prep_table)
        return s


    def prep_step(self, benes_net):
        """Generate code to compute entry i of output array of self.prepare(). 

        """
        if self.V24_INTS_USED == 1: 
            s = "tmp = ({net})[i];"
        elif self.V24_INTS_USED == 2: 
            s = "tmp = tbl[i]; tmp = ({net})[tmp & 15] >> (tmp & 0xf0);" 
        elif self.V24_INTS_USED == 3: 
            s = "tmp = tbl[i];\n"
            s += "tmp = ({net})[tmp & 15] >> ((tmp & 0xf0) >> 1);" 
        return s.format(net = benes_net)

    def prepare(self, benes_net, result):
        """Generates code for the precomputation needed for a permutation

        Given an array B describing a Benes network for a permutation 
        of 24 elements, this method genrates code for the necessary
        precomputation for performing such permutations. The result of
        the precomputation is an array M of type uint_mmv_t[] and length
        at most 21, with the exact length given by 
        self.n_benes_operations().

        The use of array M is described in method self.iter_shift_op(). 

        Parameter 'benes_net'  is the name of a variable referring to
        input array B and parameter 'result' is the name of a variable
        referring to the output array M.

        The function generates a C block in curly brackets. That block
        may contain a static array of up to 21 bytes. It uses the
        identifiers 'i', 'tmp' and 'tbl' inside the block.
        """
        s = """{
    {declare}
    for(i = 0; i < {nsteps}; ++i) {
        {step}
        // %%MMV_UINT_SPREAD tmp, tmp
        {result}[i] = tmp;
    }
}
"""
        return self.snippet(s, declare = self.prep_declare(), 
            nsteps = len(self.prep_table), 
            step = self.prep_step(benes_net),
            result = result)


    def iter_shift_op(self):
        """Yield instructions for the actual Benes network operation

        This function yields a sequence of triples

            (offset1, offset2, shift)

        For the permuation of an array v of type uint_mmv_t[] representing 
        the small integers in bit fields of length self.FIELD_BITS
        this means:

        Exchange bit  i  in v[offset1] bit  i + shift  in v[offset2]
        if bit i in the corresponding mask is set.

        Method self.prepare() generates code for computing the 
        corresponding array of masks for a Benes network prepresenting
        a permutation of 24 small integers.

        As usual in a Benes network, we have  shift == 0  if 
        offset1 != offset 2.
        """
        for (layer, offset) in self.iter_shift_cases():
            diff, sh = divmod(self.BENES_SHIFTS[layer], self.INT_FIELDS)
            yield offset, offset+diff, sh*self.FIELD_BITS

    def cond_negate(self, var, sign):
        """Generate code to negate  array of small integers

        Parameter 'var' describes the set of  of 24 small integers to be
        negated as in method benes_permutation(). All these integers
        are negated if bit 0 of variable 'sign' is set. 'sign' must be 
        of type uint_mmv_t. 
        """
        s = "{sign} = (-({sign} & {hex:1})) & {smask:P, range(24)};\n"
        for i in range(self.V24_INTS_USED):
            if i == 1 and self.V24_INTS_USED == 2:
                s += "{sign} &= {smask:P, range(8)};\n"
            s += "{var}%d ^= {sign};\n" %  i
        return self.snippet(s, var = var, sign = sign) 


    def permute(self, mask, sign, var = None):
        """Generate code for permutation with Benes network.

        The generated code permutes a sequence of 24 small integers. All
        these integers are negated if bit 0 of the variable of type
        uint_mmv_t with name 'sign' is set. Variable 'sign' is destroyed.

        Parameter 'mask' is an array M of type uint_mmv_t[] describing the 
        permutation. 'var' describes a set of at most three variables of
        type uint_mmv_t containing the sequence of 24 small integers to  
        be permuted. If var="v", these variables have names v0, v1, v2.
        There are precisely self.V24_INTS_USED such variables. By default, 
        'var' is the last recent argument passed to method declare()

        The description 'mask' of the permutation should be generated 
        using the code generated by method prepare().
        """
        if var is None: var = self.last_var
        vars = ", ".join(var + str(i) for i in range(self.V24_INTS_USED))
        vars = "(" + vars + ")"
        s = self.permute_comment.format(vars=vars, mask=mask, sign=sign)
        s += self.cond_negate(var, sign)
        for i, (ofs1, ofs2, sh) in enumerate(self.iter_shift_op()):
            if sh:
                assert ofs1 == ofs2
                s +=  """{t} = ({v}{ofs1} ^ ({v}{ofs1} >> {sh})) & {m}[{i}]; 
{v}{ofs1} ^=  {t} ^ ({t} << {sh});
""".format( t = sign, v = var, ofs1 = ofs1, sh = sh, m = mask, i = i )
            else: 
                s += """{t} = ({v}{ofs1} ^ {v}{ofs2}) & {m}[{i}]; 
{v}{ofs1} ^=  {t};  {v}{ofs2} ^=  {t};
""".format( t = sign, v = var, ofs1 = ofs1, ofs2 = ofs2, m = mask, i = i)
        s += "// Permutation of small integers done.\n"
        return s

    def declare(self, var):
        vars = ", ".join(var + str(i) for i in range(self.V24_INTS_USED)) 
        self.last_var = var
        return "uint_mmv_t %s;\n" % vars

    def load(self, source, var = None):
        if var is None: var = self.last_var
        s = ""
        for i in range(self.V24_INTS_USED):
            s += "{var}{i} = ({source})[{i}];\n".format(
                 source = source, var = var, i = i)
        return s;

    def store(self, destination, var = None):
        if var is None: var = self.last_var
        s = ""
        for i in range(self.V24_INTS_USED):
            s += "({destination})[{i}] = {var}{i} ;\n".format(
                 destination = destination, var = var, i = i)
        return s;


    def tables(self):
        return {
            "PERM24_USE_BENES_NET": self.use_benes_net,
            "PERM24_BENES_MASK_LEN": len(self.prep_table),
        }
 
    def directives(self):
        d = {
           "PERM24_BENES_PREPARE" : UserDirective(self.prepare, "ss"),
           "PERM24_BENES_DECLARE" : UserDirective(self.declare, "s"),
           "PERM24_BENES_LOAD" : UserDirective(self.load, "ss"),
           "PERM24_BENES_PERMUTE" : UserDirective(self.permute, "sss"),
           "PERM24_BENES_STORE" : UserDirective(self.store, "ss"),
        }
        if len(self.prep_table):
            return d
        else:
            d0 = {}
            for key_ in d: d0[key_] = self.not_supported
            return d0
        


###########################################################################
# Scalar products and inversions of a 2048 x 24 vector
###########################################################################


class ScalarProd2048(MM_Op):

    def __init__(self, p):
        super(ScalarProd2048, self).__init__(p)

    ######################################################################
    # Generate tables for scalar products
    ######################################################################

    def scalar_prod_table_entry(self, v1):
        v = mat24.gcode_to_vect(2*v1)
        return [ self.smask(self.P, v >> j)
                     for   j in range(0, 32, self.INT_FIELDS) ]

    def scalar_prod_table(self, dist, length, start=0):
        return sum( (self.scalar_prod_table_entry(v1)
           for v1 in range(start, start + length*dist, dist)), [] )

    def scalar_prod_2048_unroll(self, p_high, p_low, v):
        s = "uint_mmv_t v_t;\n"
        e1_3 = [self.scalar_prod_table_entry(i) for i in range(1,4)]
        for i in range(self.V24_INTS_USED):
            s += "{v}[{i}] ^= (v_t = {ph}[{i}] ^ {pl}[{i}]);\n".format( 
                 v = v, ph = p_high, pl = p_low, i = i)
            for (j, e) in tuple(zip(range(1,4), e1_3)):
                s += "{v}[{k}] ^= v_t ^ {e};\n".format( v = v,
                   k = j*self.V24_INTS + i, e = self.hex(e[i]) )
        return s



    def directives(self):
        return {
            "SCALAR_PROD_2048_UNROLL": 
                   UserDirective(self.scalar_prod_2048_unroll, "sss")
        }

    def tables(self):
        return {
            "MMV_TBL_SCALPROD_HIGH": self.scalar_prod_table(64, 32),
            "MMV_TBL_SCALPROD_LOW": self.scalar_prod_table(4, 16),
        }



###########################################################################
# Permutations of 64 entries in a row
###########################################################################


def v2_list(log_n):
    v2list = []
    for i in range(log_n):
        v2list = v2list + [i] + v2list
    return [None] + v2list
     

class SmallPerm64(MM_Op):
    N_REGS = 3
    v2 = v2_list(6)

    def __init__(self, p, shift = True):
        super(SmallPerm64, self).__init__(p)

    def generate_load_perm64(self, source, bytes, sign):
        """Load 64 enries of a vector 'source' to a byte array 'bytes'

        The vector of type uint_mmv_t[] is negated if bit 12 of the
        'sign' is set.
        """
        s = """
// Load 64 entries of vector %s to byte array %s.
// The byte array is negated if bit 12 of %s is set.
""" % (source, bytes, sign)
        n_registers = 8 // self.FIELD_BITS
        registers =  ", ".join(["r%d" % i for i in range(n_registers)])
        s = "uint_mmv_t %s;\n" % registers
        LOAD_MASK = self.hex(self.smask(self.P, -1, 8))
        SIGN_MASK = self.hex(self.smask(self.P, -1))
        s += "%s = (-((%s >> 12) & 1)) & %s;\n" % (sign, sign, SIGN_MASK)
        for i, start in enumerate(range(0, 64, self.INT_FIELDS)):
            s += "r0 =  (%s)[%d] ^ (%s);\n" % (source, i, sign)
            for j in range(1, n_registers):
                s += "r%d = (r0 >> %d) & %s;\n" % (
                   j, j * self.FIELD_BITS, LOAD_MASK)
            s += "r0 &= %s;\n" % LOAD_MASK
            for j in range(self.INT_FIELDS):
                s += "(%s)[%d] = r%d >> %d;\n" % (bytes,  start + j, 
                    j % n_registers, 8 * (j // n_registers))
        return s;

    def setup_store(self, op):
        N_REGS = max(self.N_REGS, self.LOG_INT_FIELDS)
        indices =  ["r%d" % i for i in range(N_REGS)]   
        registers = ["ri"] +  indices          
        indices += ["%s[%d]" % (op, i) for i in range(N_REGS,6)]
        s = "uint_mmv_t v;\n"
        s +=  "uint_fast8_t %s;\n" % ", ".join(registers);
        for i in range(N_REGS):
            s += "r%d = %s[%d];\n" % (i, op, i)
        s += "ri = %s;\n"  % indices[0]
        return s, indices


    def process_field(self, bytes, i, dest):
        s = "(uint_mmv_t)(%s[%s])" % (bytes, "ri" if i else "0")
        i0 = i &  (self.INT_FIELDS - 1)
        s = "v %s= %s << %d;\n" % ("+" if i0 else "",  
               s,  i0 << self.LOG_FIELD_BITS)
        if i0 == self.INT_FIELDS - 1:
            s += "%s[%d] = v ;\n" % (dest, i >> self.LOG_INT_FIELDS)
        return s     
        
    def modify_index(self, i, indices):
        if i < 2: return ""
        return "ri ^= %s;\n" % indices[self.v2[i]]

    def generate_store_perm64(self, bytes, dest, op):
        s, indices = self.setup_store(op)
        for i in range(64):
            s += self.modify_index(i, indices)
            s += self.process_field(bytes, i, dest)
        return s

    def generate_invert_perm64(self, src):
        s = ""
        wt = MM_OctadTable.perm64_weights
        for i in range(self.V64_INTS):
            s += self.snippet("{src}[{i}] ^= {smask:P,{mask}};\n",
                src = src, i = i,  mask = wt >> (i * self.INT_FIELDS))
        return s

    def directives(self):
        return {
           "LOAD_PERM64" :  UserDirective(self.generate_load_perm64, "sss"),
           "STORE_PERM64" :  UserDirective(self.generate_store_perm64, "sss"),
           "INVERT_PERM64" :  UserDirective(self.generate_invert_perm64, "s"),
        }

    def tables(self):
        return {
        }

    def _test(self):
        s = "#include <stdint.h>\ntypedef uint64_t uint_mmv_t;\n"
        s += "void f(int src[], int dest[], int op[])\n{\n"
        s += self.generate_load_perm64("src", "dest", "op") 
        s += "}\n"
        return s


######################################################################
# Test functions
######################################################################





if __name__ == "__main__":
    c = Perm24_Standard(7)
    print(c.prepare("perm", "perm_fast"))
    for i, x in c.expressions("source", "perm_fast", "sign"):
          print("dest[%s] =\n  %s;" % (i, x))
    print("1234567890"*7)
    print("\n//\n")

    c = Perm24_Benes(7)
    print(c.prepare("net", "result"))
    print(c.declare("v"))
    print(c.load("source"))
    print(c.permute("maask", "sign"))
    print(c.store("dest"))
    1/0


    with open("tmp.c", "wt") as f: f.write(SmallPerm64(3)._test())
    os.system("gcc -c tmp.c")
    os.system("del tmp.c tmp.o")
    1/0

    for p in [3, 7, 127]:
        print ("\ntest", p)
        a = SmallPerm24(p)
        a._test()





