r"""We describe the class ``MM_Basics`` and its subclass ``MM_Const``
 

Class ``MM_Basics`` in module ``mmgroup.dev.mm_basics.mm_basics`` defines 
several constants that describe the organization of vectors of integers
modulo :math:`p` used in the representation modulo :math:`\rho_p`. Here
several integers modulo  :math:`p` are stored in a single unsigned
integer of type ``uint_mmv_t``. Type  ``uint_mmv_t`` is equal to the
standard C integer type ``uint64_t`` on a 64-bit system. 

For a 32-bit system that type may be changed to ``uint32_t`` by adjusting
the variable ``INT_BITS`` in this module, but this has no longer been 
tested for several years.



    .. table:: Constants used for generating representations of the monster
      :widths: 25 75

      ==================== =====================================================
      Name                 Value
      ==================== =====================================================
      ``FIELD_BITS``       Number of bits used to store an integer modulo ``P``.
                           This is the smallest power of two greater than or
                           equal to ``P_BITS``.
      ``INT_BITS``         Number of bits available in an unsigned integer of
                           type ``uint_mmv_t``. This is equal to 64.
      ``INT_FIELDS``       Number of integers modulo ``P`` that is stored in a
                           variable of type ``uint_mmv_t``. ``INT_FIELDS`` is 
                           a power of two and equal to ``INT_BITS/FIELD_BITS``.
      ``LOG_FIELD_BITS``   Binary logarithm of the value ``FIELD_BITS``.
      ``LOG_INT_BITS``     Binary logarithm of the value ``INT_BITS``.
      ``LOG_INT_FIELDS``   Binary logarithm of the value ``INT_FIELDS``.
      ``LOG_V24_INTS``     Binary logarithm of the value ``V24_INTS``.
      ``LOG_V64_INTS``     Binary logarithm of the value ``V64_INTS``.
      ``MVV_ENTRIES``      Number of integers modulo ``P`` contained in a
                           vector of the representation of the monster,
                           including unused entries. 
      ``MVV_INTS``         Number of integers of type ``uint_mmv_t`` required
                           to store a vector of the representation of the
                           monster. This value depends on ``P``.
      ``P``                The modulus of the current representation being
                           generated. ``P`` is equal to ``2**P_BITS - 1``.
      ``P_BITS``           Bit length of modulus ``P``.
      ``V24_INTS``         Number of integers of type ``uint_mmv_t`` used up
                           to store a vector of ``24`` integers modulo ``P``.
                           ``V24_INTS`` is always a power of two.
      ``V64_INTS``         Number of integers of type ``uint_mmv_t`` used
                           to store a vector of ``64`` integers modulo ``P``.
                           ``V64_INTS`` is always a power of two.
      ``V64_INTS_USED``    Number of integers of type ``uint_mmv_t`` actually
                           used to store a vector of ``24`` integers modulo 
                           ``P``. This may be less than ``V24_INTS``.
      ==================== =====================================================

Class ``MM_Const`` in module ``mmgroup.dev.mm_basics.mm_basics`` is
a table-providing class for the C functions used by python extension
``mmgroup.mm``. Most C functions in that module take the modulus ``p`` as 
a parameter. They need fast access to the constants given in the table
above for all legal values of ``p``. Using the code generator with the 
tables and directives provided by class ``MM_Const``, we can easily 
obtain these values for any legal ``p`` as in the following example:

.. code-block:: c

   // Return a nonzero value if p is a bad modulus,  
   // i.e. not p = 2**k - 1 for some 2 <= k <= 8
   #define mm_aux_bad_p(p) (((p) & ((p)+1)) | (((p)-3) & ((0UL-256UL))))

   // Create a precomputed table containing the constants for modulus p
   // %%USE_TABLE
   static const uint32_t MMV_CONST_TABLE[] = {
   // %%TABLE MMV_CONST_TAB, uint32
   };

   int do_something_for_modulus_p(uint 32_t p, ...)
   {
       uint32_t modulus_data, p_bits, field_bits;
       // Reject any illegal modulus p
       if (mm_aux_bad_p(p)) return -1;
       // Load the constants for modulus p to variable modulus_data
       // %%MMV_LOAD_CONST  p, modulus_data;
       // Load value P_BITS for modulus P to p_bits (using modulus_data)
       p_bits = %{MMV_CONST:P_BITS,modulus_data};
       // Load value FIELD_BITS for modulus P to field_bits 
       field_bits = %{MMV_CONST:FIELD_BITS,modulus_data};
       // Now do the actual work of the function using these values
       ...
   }
 

"""
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
from mmgroup.generate_c import format_item, c_snippet
from mmgroup.generate_c import UserDirective, UserFormat

INT_BITS = 64  # Once and for all


########################################################################
# Class MM_Basics and auxiliary functions for that class
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

def smask_default(default_width, value, fields = -1, width = None):
    if width is None:
        width = default_width
    return smask(value, fields, width)


def shl_fix(data, i):
    """Generate expression <data> << i for a fixed signed int i

    Generate <data> >> (-i)  in case i < 0.
    """ 
    data = str(data)
    try:
       d = int(data)
       return c_hex(d >> -i if i < 0 else d << i) 
    except:
        if i == 0:
            return "({0})".format(data)
        if i > 0:
            return "(({0}) << {1})".format(data, i)
        if i < 0:
            return "(({0}) >> {1})".format(data, -i)
          



class MM_Basics(object):
    """The basic table-providing class for vectors in :math:`\rho_p` 

    This class provides the constants given in the table above.
    It is an abstract class. Subclasses of this class  may provide
    tables and directives for dealing with vectors in :math:`\rho_p`.
    """
    MMV_ENTRIES = 759 * 64 + (3 * 24 + 3 * 2048) * 32
    INT_BITS = INT_BITS
    UINT_SUFFIX = "UL" if INT_BITS == 32 else "ULL"
    # must be 32 or 64 for uint32_t or uint64_t
    assert INT_BITS in (32, 64)
    LOG_INT_BITS = hibit(INT_BITS)
    tables = { 
        "GENERATE_CODE" : True,
        "INT_BITS": INT_BITS,
        "LOG_INT_BITS": LOG_INT_BITS,
        "MMV_ENTRIES" : MMV_ENTRIES,
        "hex": UserFormat(hex, arg_fmt = c_hex),
        "shl": UserFormat(shl_fix, "si"),
    }
    table_dict = {}
    _sizes = {}
    const_format = "uint"+str(INT_BITS) 
    p = 0  # p = 0 means that no modulus p has been specified

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
            d = cls.tables.copy()
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
            partial_smask = partial(smask_default, FIELD_BITS)
            d["smask"] = UserFormat(partial_smask, fmt = c_hex)
            cls.table_dict[p] = d
            return d



    def snippet(self, source, *args, **kwds):
        """Generate a C code snippet

        This method calls function ``c_snippet`` in module
        ``mmgroup.generate_c`` to generate a C code snippet.
        Parameters are as in function ``c_snippet``. The method
        returns the generated C code snippet as a string.

        In addition, all tables and directives available in the given 
        instance of the table-providing class are also passed as
        arguments to function ``c_snippet``.
        """
        tables = {"directives" : self.directives, **self.tables, **kwds}
        return c_snippet(source, *args, **tables)





########################################################################
# Class MM_Const 
########################################################################


# Auxilary function for defining methods like MM_Const().P_BITS(p) etc.
_attr_from_table = lambda name, p : MM_Basics.sizes(p)[name]


class MM_Const(MM_Basics):
    """This is the basic table-providing class for module ``mmgroup.mm``

    The main purpose of this class is to provide the constants 
    defined in class ``MM_Basics, for a variable modulus ``p`` as
    shown in the example above.

    This class provides the directive ``MMV_CONST_TAB`` for generating
    all tables, and the directive ``MMV_LOAD_CONST`` for storing the 
    table of constants, for a specific modulus ``p``, in an integer 
    variable. The string formatting function ``MMV_CONST`` can be
    used for extracting a specific constant from that variable, as 
    indicated in the example above.

    Internally, we use a deBruijn sequence to translate the value ``p``
    to an index for the table generated via the directive 
    ``MMV_CONST_TAB``. 

    Constants not depending on the modulus ``p``, such as ``INT_BITS``, 
    ``LOG_INT_BITS``, and ``MMV_ENTRIES`` are available as attributes 
    of class ``MM_Const``. They can also be coded with the code 
    generator directly via string formatting, e.g.::

        uint_8_t  a[%{MMV_ENTRIES}];

    For a fixed modulus ``p`` the constants depending on ``p`` can also
    be coded with the code generator via string formatting, e.g.::

        a >>= %{P_BITS:3};

    That constant also availble in the form ``MM_Const().P_BITS(3)``.


    Class ``MM_Const`` provides a string-formatting function ``shl`` 
    which generates a shift expression  Here::

         %{shl:expression, i}

    generates an expression equivalent to ``((expression) << i)``. 
    Here ``i`` must be an integer. In case ``i < 0`` we generate 
    ``((expression) >> -i)`` instead. E.g. ``%{shl:'x',-3}``
    evaluates to ``x >> 3``.

    Class ``MM_Const`` provides another string-formatting function
    ``smask`` which generates an integer constant to be used as a 
    bit mask for integers of type ``uint_mmv_t``. Here::

         %{smask:value, fields, width}

    with integers ``value``, ``fields``, and ``width`` creates such a 
    bit mask. ``For 0 <= i`` and ``0 <= j < width``, the bit of that 
    mask at  position ``i * width + j`` is set if both, bit ``i`` of 
    the integer ``fields`` and bit ``j`` of the integer ``value`` are 
    set. A bit in the mask at  position ``INT_BITS`` or higher is 
    never set.
    
    Any of the arguments ``value`` and ``fields`` may either be an 
    integer or anything iterable that yields a list of integers. 
    Then this list is interpreted as a list of bit positions. A bit 
    of that argument is set if its position occurs in that list and 
    cleared otherwise.
  
    E.g. ``%{smask:3, [1,2], 4}`` evaluates to ``0x330``.
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

    

 
    def __init__(self, **kwds):
        new_tables = {
           self.T_NAME : self.table,
           "MMV_CONST" : UserFormat(self.f, "ss"),
           "smask":  UserFormat(smask, fmt = c_hex),
        }
        self.tables = self.tables.copy()
        self.tables.update(new_tables) 
        self.directives = {
          "MMV_LOAD_CONST" : 
              UserDirective(self.gen_get_const_table, "ss", "NAMES"),
        }
        for name in [
            "P", "P_BITS", "FIELD_BITS", "LOG_FIELD_BITS", 
            "INT_FIELDS", "LOG_INT_FIELDS", "V24_INTS", "LOG_V24_INTS",
            "V24_INTS_USED", "V64_INTS", "LOG_V64_INTS",
            ]:
            f = partial(_attr_from_table, name)
            setattr(self, name, f)
            self.tables[name] = UserFormat(f, "i")


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
      


class Mockup_MM_Const(MM_Const):
    def __init__(self, **kwds):
        super(Mockup_MM_Const, self).__init__()
        old_tables = self.tables
        self.tables = {}
        self.tables.update(old_tables)
        self.tables["GENERATE_CODE"] = False


Tables = MM_Const
MockupTables = Mockup_MM_Const

