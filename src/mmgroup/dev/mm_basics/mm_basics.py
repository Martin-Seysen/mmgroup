from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
import collections
import warnings
from numbers import Integral
import re
import ast
import functools
from functools import reduce, partial
from operator import __or__





from mmgroup.bitfunctions import bitlen, hibit, lmap
from mmgroup.dev import config
from mmgroup.dev.generate_c.make_c_tables import format_item, c_snippet
from mmgroup.dev.generate_c.generate_functions import UserDirective
from mmgroup.dev.generate_c.generate_functions import UserFormat

INT_BITS = config.INT_BITS


########################################################################
# Class MM_Basics and axuxluiar functions for that class
########################################################################

UINT_SUFFIX = "UL" if INT_BITS == 32 else "ULL"

def c_hex(x):
    """Convert the integer x to a useful unsigned C hex literal

    It is assumed that no integer type has more bits that
    an integer of type uint_mmv_t.
    """
    #print("c_hex", x, type(x))
    h = hex(x & ((1 << INT_BITS) - 1))
    while h[-1:] == "L": h = h[:-1] # for python 2
    return h + UINT_SUFFIX

def smask(value, fields, width):
    """Create a simple bit mask constant

    For 0 <= i and 0 <= j < 'width', the bit of that constant at 
    position i * 'width' + j is set if both, bit i of integer 'fields'
    and bit j of integer 'value' are set. A bit at 
    position >= INT_BITS is never set.
    
    Any of the arguments 'value' and 'fields' may either be an 
    integer or anything iterable that yields a list of integers. 
    Then this list is intepreted as a list of bit positions. A bit 
    of the argument is set if its position occurs in that list and 
    cleared otherwise.
    """
    try:
        value = reduce(__or__, [1 << i for i in value], 0)
    except TypeError:
        pass 
    try:
        fields = reduce(__or__, [1 << i for i in fields], 0)
    except TypeError:
        pass
    assert width > 0 
    value &= (1 << width) - 1
    result = sum(value << (width * i)
        for i in range(INT_BITS//width) if (fields >> i) & 1)
    return result & ((1 << INT_BITS) - 1)




class MM_Basics(object):
    """Representation of a vector of numbers modulo p in a C integer

    Such a vector is stored in an array of integers of type uint_mmv_t,
    which is uint32_t or uint64_t, depending on the constant INT_BITS
    in file config.py. 

    For the modulus p, the number p + 1 must be a power of two. So
    you should e.g. put p = 15 for calculating in characteristic 5.
    This greatly simplifies negation modulo p, which is just obtained
    as the bitwise complement of p.

    An integer of type uint_mmv_t (also called an 'int' for brevity)
    can store several fields, where each field stores one number 
    modulo p. The bit length of a field is always a power of two.
    So we have the following constants:

    P                 This is just p.
    INT_BITS          Bit length of an int; this is 32 of 64.   
    LOG_INT_BITS      Binary logarithm of INT_BITS.
    FIELD_BITS        Bit length of a field. This is a power of 2.
    LOG_FIELD_BITS    Binary logarithm of FIELD_BITS. 
    INT_FIELDS        Number of fields in an int. This a power 
                      of 2 and equal to INT_BITS / FIELD_BITS.
    LOG_INT_FIELDS    Binary logarithm of INT_FIELDS.    
    P_BITS            bit length of a number modulo p. We have    
                      p = 2 ** P_BITS - 1 and P_BITS <= FIELD_BITS.

    A vector af the 196884-dimensional representation 'MM' of the 
    monster group has a natural representation as a collection of 
    2-dimensional matrices, see e.g. [1]. Each of these matrices has 
    24 or 64 columns. Keeping track of these structures greatly
    simplifies the operation of the Monster group. Especially, we 
    align all rows of these matrices to the beginning of an int. 
    So a row of length 24 may or may not fill up all fields of an
    array of ints. To overcome this problem, we always reserve 32
    entries for a tow of length 24. We define more constants: 
       
    V24_INTS        No of ints required to store 24 fields. 
    V64_INTS        No of ints required to store 64 fields.
    MMV_INTS        Length of a vector of rep. MM in ints
    LOG_V24_INTS    binary logarithm of V24_INTS
    LOG_V64_INTS    binary logarithm of V64_INTS
    V24_INTS_USED   No of ints actually used by a vector of 24 fields

    Any instance of this class created for a modulus p contains
    these constants for the modulus p as attributes. The same
    names for these constants are also used for code generation,
    see method tables().

                      
    [1]  J.H. Conway, A simple construction for the Fischer-Griess
         monster group, Invent. Math. 79 (1985), 513-540.
    """
    MMV_ENTRIES = 759 * 64 + (3 * 24 + 3 * 2048) * 32
    INT_BITS = INT_BITS
    UINT_SUFFIX = "UL" if INT_BITS == 32 else "ULL"
    # must be 32 or 64 for uint32_t or uint64_t
    assert INT_BITS in (32, 64)
    LOG_INT_BITS = hibit(INT_BITS)
    table_dict = { 0 : { 
        "INT_BITS": INT_BITS,
        "LOG_INT_BITS": LOG_INT_BITS,
        "MMV_ENTRIES" : MMV_ENTRIES,
        "hex": UserFormat(hex, arg_fmt = c_hex),
        "smask": UserFormat(smask, fmt = c_hex),
        "directives": {}
    } }
    _sizes = {}
    p = 0  # p = 0 means that no modulus p has been specified

    def __init__(self, p = 0):
        """Create an instance of class MM_Basics for modulus p.

        """
        self.p = p
        self.const_format = "uint"+str(self.INT_BITS) 
        self.tables_ = {**self.sizes(p)}

    def __getattr__(self, attrib):
        return self.tables_[attrib]
 

    @classmethod
    def hex(cls, x):
        """Convert the integer x to a useful unsigned C hex literal

        It is assumed that no integer type has more bits that
        an integer of type uint_mmv_t.
        """
        return c_hex(x)




    @classmethod
    def sizes(cls, p):
        """Return dictionary of constants for a specific p.""" 
        try:
            return cls.table_dict[p]
        except:
            assert p > 1 and p & (p + 1) == 0, str(p)
            d = {}
            d.update(cls.table_dict[0])
            d["P"] = p
            d["P_BITS"] = P_BITS = bitlen(p)
            FIELD_BITS = P_BITS
            while (FIELD_BITS & (FIELD_BITS - 1)): 
                FIELD_BITS += 1 
            d["FIELD_BITS"] = FIELD_BITS
            d["LOG_FIELD_BITS"] = hibit(FIELD_BITS)
            d["INT_FIELDS"] = INT_FIELDS = cls.INT_BITS // FIELD_BITS
            d["LOG_INT_FIELDS"] = hibit(INT_FIELDS)
            V24_INTS = 32 // INT_FIELDS
            d["V24_INTS"] = V24_INTS
            d["LOG_V24_INTS"] = hibit(V24_INTS)
            d["V24_INTS_USED"] = V24_INTS - (V24_INTS >> 2)
            V64_INTS = 64  // INT_FIELDS
            d["V64_INTS"] = V64_INTS
            d["LOG_V64_INTS"] = hibit(V64_INTS)
            d["MMV_INTS"] = 3 * (2048 + 24) * V24_INTS + 759 * V64_INTS
            d["smask"] =  UserFormat(partial(cls.smask_p, p), fmt=c_hex)
            cls.table_dict[p] = d
            return d


    @classmethod
    def smask_p(cls, p, value = None, fields = -1, width = None):
        """Auxiliary function for method self.smask()

        self.smask is equivalent to partial(self.smask_p, self.p).
        """
        #print("This is smask_p")
        if value is None:
            value = p
        if width is None:
            width =  cls.table_dict[p]["FIELD_BITS"]
        return smask(value, fields, width)

    def smask(self, value = None, fields = -1, width = None):
        """Return smask(value, fields, width)

        'value' defaults to self.P, 'fields' defaults to -1, and 
        'width' defaults to self.FIELD_BITS.
        """
        return self.smask_p(self.p, value, fields, width)

    @staticmethod
    def mul_fix(data, i):
        """Generate expression <data>*i for a fixed unsigned int i"""
        assert i >= 0
        if i == 0:
            return "0"
        if i == 1:
            return "({0})".format(data)
        if i & (i-1) == 0:
            sh = hibit(i)
            return "(({0}) << {1})".format(data, sh)
        if i & (i+1) == 0:
            sh = hibit(i+1)
            return "((({0}) << {1}) - ({0}))".format(data, sh)
        if (i-1) & (i-2) == 0:
            sh = hibit(i-1)
            return "((({0}) << {1}) + ({0}))".format(data, sh)
        return "(({0}) * {1})".format(data, i)


    @staticmethod
    def shl_fix(data, i):
        """Generate expression <data> << i for a fixed signed int i

        Generate <data> >> (-i)  in case i < 0.
        """ 
        if i == 0:
            return "({0})".format(data)
        if i > 0:
            return "(({0}) << {1})".format(data, i)
        if i < 0:
            return "(({0}) >> {1})".format(data, -i)



    @classmethod
    def deref(cls, array, index = 0):
        """Return the string 'array[index]'

        'array' should be an identifer, i.e. a string and 'index'
        should be an integer. 

        If 'array' is "a + j" or "a - j" then "a[j + index]" or 
        "a[-j + index]" is returned.
        """
        try:
            a, i = cls.m_psum.match(array).groups() 
            while i and i[0] in " +":
                 i = i[1:]
            i = int(ast.literal_eval(i)) if i else 0
            return "%s[%d]" % (a, i + index)
        except:
            return "(%s)[%d]" % (array, index)  

    def make_tables(self, new_tables, new_directives = {}):
        self.tables_.update(new_tables)
        self.tables_["directives"].update(new_directives)

    def directives(self):
        return self.tables_["directives"]

    def tables(self):
        return self.tables_

    def snippet(self, source, *args, **kwds):
        tables = {**self.tables_, **kwds}
        return c_snippet(source, *args, **tables)
        pass


########################################################################
# Class MM_Const 
########################################################################


class MM_Const(MM_Basics):
    """Fast computation of MM_Basics constants for variable p

    The purpose of this program is to store the constants defined 
    in class MM_Basics, for all p of shape 2**k - 1 with
    2 <= k <= 8, so that they can be computed very efficiently for 
    any such p. The list of names of these constants is given in the 
    list  MM_Const.entries. We compute these constants and store them
    in the bit fields of a 32-bit integer for all p.

    We use the code generation mechanism in class in class 
    TableGenerator in module make_c_tables. A source file for
    generating and using these constants table may look as follows:

    // Store table of all constants in the array MMV_CONST_TABLE
    // %%USE_TABLE
    static const uint32_t MMV_CONST_TABLE[] = {
      // %%TABLE %MMV_CONST_TAB, uint32
    };
    
    uint32_t p, const_p, c;
    p = 7; // usually the variable p is an input of a function
    // The following directive MMV_LOAD_CONST loads the contstants
    // for the given p to variable const_p
    // %%MMV_LOAD_CONST p, const_p
    // Next assign constant FIELD_BITS to variable c for the given p.
    // Here we take the constants for that p from variable p_const.
    c = {MMV_CONST:FIELD_BITS,const_p}; 

    This class provides the directive MMV_CONST_TAB for generating
    all tables, the directive MMV_LOAD_CONST for storing the table
    of constants for a spcific p in an integer variable, and the
    function  MMV_CONST for extracting a specific contant from that
    variable. Here we follow the rules of class TableGenerator.

    On input p, the directive MMV_LOAD_CONST loads the contants
    to a variable of type uint32_t. Here bit fields of that
    variable correspond do the constants with names given by
    MM_Basics.entries. Start positions and lengths of the bit 
    fields are given by MM_Basics.pos and MM_Basics.lengths.
    Given p, the entry with index log(p+1) is taken from the
    array of constants created by the directive MMV_CONST_TAB.
    Note that p+1 is a power of 2, so that log(p+1) is an integer.
    We compute log(p+1) from p with a deBruijn sequence
    (given by the integer MM_Basics.deBruijnMult which represents
    a bit sequence). We also scramble the array of constants so 
    that so that access via a deBruijn sequence is most efficient. 
    """
    # Names of the constants to be stored in an integer for each p
    entries = [  
        "LOG_INT_FIELDS", "INT_FIELDS", "LOG_FIELD_BITS",
        "FIELD_BITS", "P_BITS",
    ] 

    primes = [(1 << k) - 1 for k in range(2,9)]  # list of all p
    lengths = [0] * len(entries)  # list of max bit length of constants
    for p in primes:
        data =  MM_Basics.sizes(p)
        for i, name in enumerate(entries):        
             lengths[i] |= data[name]
    lengths = lmap(bitlen, lengths)
    #print ("SUM MM_Const bitlengths", sum(lengths))

    assert sum(lengths) <= 32
    pos = {}                      # pos is a dictionary of shape
    start = 0                     #  {entry_name:(start_pos,max_value)}
    for i, length in enumerate(lengths):
        pos[ entries[i] ] = (start, (1 << length) - 1)  
        start += length
    deBruijnMult = 0xe8           # deBruijn sequence
    table = [0] * 8               # This will be the table of constants
    for p in primes:              # Assemble that table
        data =  MM_Basics.sizes(p)
        index = (((p + 1) * deBruijnMult) >> 8) & 7 # scramble indices of 
                                  # table access via deBruijn sequence
        value = 0                  
        for entry, (shift, mask) in pos.items(): # assemble entry for p
            assert 0 <= data[entry] <= mask
            value += data[entry] << shift
        table[index] = value      # store entry for p in table
    #print("pos", pos)
    #rint("tbl", lmap(hex,table))
    T_NAME = "MMV_CONST_TAB"      # python name of contant table
    F_NAME = "MMV_CONST"        # function name "MMV_CONST"
    LOAD_F_NAME = "MMV_LOAD_CONST" # directive name "MMV_LOAD_CONST"

    # Some more constants are definded as  dividend/INT_FIELDS
    # with dividends for contant names given by:
    dividend = {
        "V64_INTS": 64, 
        "V24_INTS": 32,
        "MMV_INTS" : MM_Basics.MMV_ENTRIES,
    }


    def __init__(self):
        super(MM_Const, self).__init__(0)
        new_tables = {
           self.T_NAME : self.table,
           "P_LIST" : config.PRIMES, 
           "MMV_CONST" : UserFormat(self.f, "ss"),
        }
        new_directives = {
          "MMV_LOAD_CONST" : 
              UserDirective(self.gen_get_const_table, "ss", "NAMES"),
        }
        self.make_tables(new_tables, new_directives)


    @classmethod
    def gen_get_const_table(cls, names, p,  const_p):
        """Code generating method for directive MMV_LOAD_CONST 

        According to rules of class TableGenerator the first argument 
        of this method may be an implicit dictionary that translates 
        the python name of the table of constants to the C name of that
        table. So names[cls.T_NAME] is the appropriate C name. 
       
        We calculate the index for the given input p and load the
        entry of that table corresponding to p (which is an integer)
        to the destination given by const_p.
        """
        s = "// Store constant table for {p} to {const_p}\n".format(
            const_p = const_p, p = p
        )    
        s += "{c} = {t}[(((({p}) + 1) * {mul}) >> 8) & 7];\n".format(
            c = const_p, p = p, t = names[cls.T_NAME],
            mul = cls.deBruijnMult
        )
        return s

    @classmethod
    def gen_get_const_expr(cls, const_name, const_p):
        """General code generating method for function MMV_CONST 

        It generates code that extracts the constant with name 
        const_name from the entry const_p for a specific p. 
        """
        sh, mask = cls.pos[const_name]
        s = "(({t}{sh}) & {mask})".format(
            t = const_p, mask = mask,
            sh = " >> " + str(sh) if sh else ""
        )
        return s


    @classmethod
    def gen_div_int_fields(cls, const_name, const_p):
        """Code generating for constants with names in self.dividend

        All these constants are defined as dividend/INT_FIELDS.
        The dividents for these contant names are taken from
        dictionary cls.dividend. The constant INT_FIELDS is calculated
        by method  gen_get_const_expr() of this class.
        """
        n_entries = cls.dividend[const_name]
        log_ = cls.gen_get_const_expr("LOG_INT_FIELDS", const_p)
        return "({0} >> {1})".format(n_entries, log_)
 
    @classmethod
    def f(cls, name, *args):
        """Code generating method for function MMV_CONST 

        It generates code that returns the constant with given
        'name'.
        """
        try:
             return cls.gen_div_int_fields(name, *args)
        except:
            return cls.gen_get_const_expr(name, *args)
      


########################################################################
# Class MM_Op
########################################################################

class MM_Op(MM_Basics):
    """Basic operations on integers of type uint_mmv_t.

    Here an integer of type uint_mmv_t is interpreted as an array of   
    INT_FIELDS bit fields, where each bit field is FIELD_BITS long 
    and stores an integer modulo p. Modulus p is given in the 
    constructor, and p + 1 must be a power of 2, see class 
    mm_basics.MM_Basics for details.

    Method of this class generate code for some basic operations
    on arrays of such bit fields.
    """
    m_psum = re.compile(r"\s*(\w+)\s*([+-].*)?$")

    def __init__(self, p):
        """Create instance for modulus p with p + 1 a power of two"""
        super(MM_Op, self).__init__(p)
        new_directives =  {
            "MMV_ROTL" : UserDirective(self.gen_mmv_rotl, "sis"),
            "MMV_XOR24": UserDirective(self.gen_mmv_xor24, "s"*100),
            "MMV_CNEG24": UserDirective(self.gen_mmv_cneg24, "ss"),
            "MMV_UINT_SPREAD": UserDirective(self.gen_mmv_uint_spread, "ss"),
        }
        new_tables = {}
        self.make_tables(new_tables, new_directives)

    def gen_mmv_rotl(self, src, count, dest):
        """Rotate 'src' left by count 'bits', store result to 'dest'

        This generatates code for computing 
        'dest' = 'src' * 2**k  (mod p)
        for any fixed integer k. 'src' and 'dest' of type uint_mmv_t.
        They are interpreted as arrays of bit fields. The operation 
        is done on all bit fields simultaneously.
        """
        count = int(count) % self.P_BITS
        if count == 0:
            return
        s = """// Put {dest} = {src} rotated left by {c} for all fields
{dest} = ((({src}) << {c}) & {smask: (-1 << {c}) & P})
       | ((({src}) >> {int:P_BITS-c}) & {smask:range(c)});
"""
        return self.snippet(s, src = src, dest = dest, c=count)


    def v24_len(self):
        res = self.V24_INTS
        return res if res < 4 else res - (res >> 2)

    def gen_mmv_xor24(self, *args):
        """Compute 'dest' = 'src1'  ^  'src2'   ^  ...

        Generarate code for computing 'dest' = 'src1' ^ 'src2' ^ ... .
        Here 'dest' is the last argument args[-1], and 
        'src1', 'src2', ... are the arguments in the tuple *args[:-2].
        All arguments are of  type uint_mmv_t[]. They are interpreted 
        as a small array of 24 integers modulo p. The operation is 
        done on all bit fields simultaneously. We need V24_INTS 
        integers of type uint_mmv_t to store 24 integers modulo p.
        """
        sources, dest  =  args[:-1], args[-1]
        s = """// Put {dest} = {0}, for arrays
// {dest}, {1}  of 24 mmv fields
""".format(" ^ ".join(sources), ", ".join(sources), dest=dest)
        
        for i in range(self.v24_len()):
            s += "{0} = {1};\n".format( self.deref(dest, i),
                " ^ ".join([self.deref(s, i) for s in sources]) )
        return s

    def gen_mmv_cneg24(self, sign, var):
        """Negate 'var' if 'sign' & 1 == 1
  
        Genrerate code for performing that operation on an integer
        'var' of type uint_mmv_t. 'var' is  interpreted as a small 
        array 24 integers modulo p. The operation is done on all bit
        fields (each representing an integer mod p) simultaneously. 
        We need V24_INTS integers of type uint_mmv_t to store
        24 integers modulo p.
        """
        #mask = self.simple_mask(nfields = 24) 
        #mask2 = self.simple_mask(nfields = 8)  
        s = """// if {sign} & 1 then negate ({var}), with 
// ({var}) an array of 24 mmv fields; destroying {sign}.
{sign} = -({sign} & 1) & {smask:-1,range(24)};
"""
        for j in range(self.v24_len()):
            if j == 1 and self.INT_FIELDS == 16:
                s += "{sign} &= {smask:-1,range(8)};\n"  
            s += "%s ^= {sign};\n" % self.deref(var, j)
        #return s.format(sign=sign, var=var, mask=mask, mask2=mask2)
        return self.snippet(s, sign=sign, var=var)



    def gen_mmv_uint_spread(self, src, dest, *args):
        """Spread bits of the integer 'src' to field array 'dest'

        'dest' is of type uint_mmv_t and interpreted an array of 
        bit fields. bit field i of 'dest' is set to -1 (mod p) if 
        bit i of 'src' is one and to 0 otherwise.

        An optional additional fixed integer  argument 'value' may 
        be given. Then  bit field i of 'dest' is set to  'value'
        if bit i of 'src' is one.
        """ 
        msize = 1 << (self.INT_BITS >> 1)

        value = self.P
        if len(args) and not args[0] is None:
            value = int(args[0])
            if value != self.P: value %= self.P
        if value == 0: 
            return "{dest} = 0;\n".format(dest = dest)
        s = """// Spread bits 0,...,{i} of {src} to the ({j}-bit long) fields
// of {dest}. A field of {dest} is set to {value} if its 
// corresponding bit in input {src} is one and to 0 otherwise.
""". format(i = self.INT_FIELDS - 1, value = hex(value), j = self.FIELD_BITS,
                      dest = dest, src= src)
        a, b = self.INT_FIELDS >> 1, self.INT_BITS
        while a:
            m0 = self.hex(self.smask(range(a), -1, b))
            m1 = self.hex(self.smask(range(a, 2*a), -1, b))
            b >>= 1
            s0 = "{dest} = ({src} & {m0}) + (({src} & {m1}) << {sh});\n"
            if len(m0 + m1) >= 17:
                s0 = "\n    + ".join(s0.split("+"))
            s += s0.format(dest=dest, src=src, m0=m0, m1=m1, sh=b-a)
            a >>= 1
            src = dest
        s += "%s = %s;\n"  % (dest, self.mul_fix(dest, value))
        s += "// Bit spreading done.\n"
        return s





########################################################################
# Test functions
########################################################################


if __name__ == "__main__":
    for p in (3,7):
        b =  MM_Basics(p)
        print(p, "\n", b.tables(), "\n" )
