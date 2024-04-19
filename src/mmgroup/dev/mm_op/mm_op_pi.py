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



from mmgroup.generate_c import c_snippet, TableGenerator, make_table
from mmgroup.generate_c import UserDirective, UserFormat



from mmgroup import mat24

from mmgroup.dev.mm_basics.mm_tables import MM_OctadTable
from mmgroup.dev.mm_op.mm_op import MM_Op




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
    uint_mmv_t benes_mask[%{PERM24_BENES_MASK_LEN}]; 
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
        super(Perm24_Benes, self).__init__(p = p)
        self.prep_table = []
        self.use_benes_net = (
            self.V24_INTS_USED <=  2 + bool(prefer_benes_net) )
        if self.V24_INTS_USED <= 3:
            self.prep_table = self.make_prep_table()
            #print("Benes table:", list(map(hex, self.prep_table)))
        self.tables.update(self.make_tables())
        self.directives.update(self.make_directives())

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
    // %%TABLE table_prepare_perm24, uint8
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
    %{declare}
    for(i = 0; i < %{nsteps}; ++i) {
        %{step}
        // %%MMV_UINT_SPREAD tmp, tmp
        %{result}[i] = tmp;
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
        s = "%{sign} = (0-(%{sign} & %{hex:1})) & %{smask:P, range(24)};\n"
        for i in range(self.V24_INTS_USED):
            if i == 1 and self.V24_INTS_USED == 2:
                s += "%{sign} &= %{smask:P, range(8)};\n"
            s += "%%{var}%d ^= %%{sign};\n" %  i
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

    def copy(self, var_src, var_dest):
        s = ""
        for i in range(self.V24_INTS_USED):
            s += "({dest})[{i}] = ({src})[{i}];\n".format(
                 src = var_src, dest = var_dest, i = i)
        for i in range(self.V24_INTS_USED, self.V24_INTS):
            s +=  "({dest})[{i}] = 0;\n".format(
                    dest = var_dest, i = i)
        return s;

    def store(self, destination, var = None):
        if var is None: var = self.last_var
        s = ""
        for i in range(self.V24_INTS_USED):
            s += "({destination})[{i}] = {var}{i} ;\n".format(
                 destination = destination, var = var, i = i)
        for i in range(self.V24_INTS_USED, self.V24_INTS):
            s +=  "({destination})[{i}] = 0;\n".format(
                    destination = destination, i = i)
        return s;

    def make_tables(self):
        return {
            "PERM24_USE_BENES_NET": self.use_benes_net,
            "PERM24_BENES_MASK_LEN": len(self.prep_table),
        }
 
    def make_directives(self):
        d = {
           "PERM24_BENES_PREPARE" : UserDirective(self.prepare, "ss"),
           "PERM24_BENES_DECLARE" : UserDirective(self.declare, "s"),
           "PERM24_BENES_LOAD" : UserDirective(self.load, "ss"),
           "PERM24_BENES_COPY" : UserDirective(self.copy, "ss"),
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
        super(ScalarProd2048, self).__init__(p = p)
        self.tables.update(self.make_tables())
        self.directives.update(self.make_directives())

    ######################################################################
    # Generate tables for scalar products
    ######################################################################

    def scalar_prod_table_entry(self, v1):
        v = mat24.gcode_to_vect(v1)
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



    def make_directives(self):
        return {
            "SCALAR_PROD_2048_UNROLL": 
                   UserDirective(self.scalar_prod_2048_unroll, "sss")
        }

    def make_tables(self):
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
        super(SmallPerm64, self).__init__(p = p)
        self.directives.update(self.make_directives())



    def generate_load_perm64(self, source, bytes, sign, type_ = ""):
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
        s += "%s = (0-((%s >> 12) & 1)) & %s;\n" % (sign, sign, SIGN_MASK)
        type_ = "(%s)" % type_ if len(type_) else ""
        for i, start in enumerate(range(0, 64, self.INT_FIELDS)):
            s += "r0 =  (%s)[%d] ^ (%s);\n" % (source, i, sign)
            for j in range(1, n_registers):
                s += "r%d = (r0 >> %d) & %s;\n" % (
                   j, j * self.FIELD_BITS, LOAD_MASK)
            s += "r0 &= %s;\n" % LOAD_MASK
            for j in range(self.INT_FIELDS):
                s += "(%s)[%d] = %s(r%d >> %d);\n" % (bytes,  start + j, 
                    type_, j % n_registers, 8 * (j // n_registers))
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
            s += self.snippet("%{src}[%{i}] ^= %{smask:P,{mask}};\n",
                src = src, i = i,  mask = wt >> (i * self.INT_FIELDS))
        return s


    def _test(self):
        s = "#include <stdint.h>\ntypedef uint64_t uint_mmv_t;\n"
        s += "void f(int src[], int dest[], int op[])\n{\n"
        s += self.generate_load_perm64("src", "dest", "op") 
        s += "}\n"
        return s


    def make_directives(self):
        return {
           "LOAD_PERM64" : UserDirective(self.generate_load_perm64, "ssss"),
           "STORE_PERM64" : UserDirective(self.generate_store_perm64, "sss"),
           "INVERT_PERM64" : UserDirective(self.generate_invert_perm64, "s"),
        }


######################################################################
# Summarizing the tables given above
######################################################################

class Tables:
    def __init__(self, **kwds):
        self.p = p = int(kwds.get('p', 3))
        self.tables = {}
        self.directives = {}
        table_classes = [Perm24_Benes(p), ScalarProd2048(p),
                         SmallPerm64(p)]
        for t in table_classes:
            self.tables.update(t.tables)
            self.directives.update(t.directives)

class MockupTables:
    tables = {}
    directives = {}
    def __init__(self, **kwds):
        pass

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





