r"""The module contains classes ``Lsbit24Function`` and ``MatrixToPerm``.

Class ``Lsbit24Function`` contains tables and directives for computing
the least significant bit of a ``24-bit`` integer using a
DeBruijn sequence. This may be overkill for such a simple function,
but it may also be considered as a didactic example of a 
table-providing class.

Class ``MatrixToPerm`` contains a directive that generates highly 
optimized code for converting a bit matrix acting on a Golay code 
word to a permutation in the Mathieu group ``Mat_24``. Here we 
assume the bit matrix actually encodes an element of the  
Mathieu group; otherwise garbage is returned. The Golay code
word must be given as an integer in ``gcode`` representation.
"""
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

# python 3 string formatting see:
# https://pyformat.info/


import types
import sys
import re
import os
from random import randint
from operator import __or__, __xor__, __and__
from functools import reduce

import numpy
from numpy import array, zeros, uint8, uint16, uint32


from mmgroup.bitfunctions import bw24, bits2list, iter_bitweight, v2
from mmgroup.bitfunctions import iter_ascending_bitweight, lmap



from mmgroup.generate_c import prepend_blanks, format_item
from mmgroup.generate_c import UserDirective, UserFormat



   


 


############################################################################
############################################################################
# Generating C code for bit functions
############################################################################
############################################################################


def lsbit24(x):
    """Return position of least significant bit of (x & 2**24)"""
    return v2(x | 0x1000000)


def _make_de_bruijn_sequence():
    """Create DeBruijn sequence for finding the lowest bit of an integer

    This is a modification of
    http://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightMultLookup
    which returns 24 if x & 0xffffff is zero.

    Here we just inspect the lower 24 bit of an integer and return 24
    if they are all zero.
    """
    lsbit_table = numpy.zeros(32, dtype = uint8)
    sequence = 0x077CB531
    for i in range(32):
        index =  ((sequence << i) >> 26) & 0x1f
        lsbit_table[index] = min(i,24)  
    assert lsbit_table[0] == 24
    return sequence,  lsbit_table


class Lsbit24Function(object):
    """Generate code for computing function ``lsbit24()``

    The class uses the ideas in function ``_make_de_bruijn_sequence()``
    for computing that  least significant bit of an integer.

    This documentation is intended to explain the cooperation of this
    class with the class ``TableGenerator`` in module
    ``mmgroup.generate_c``. Classes designed 
    for a similar cooperation may refer to this documentation.

    The class requires a table ``lsbit_table``, which is generated 
    by function ``_make_de_bruijn_sequence()``.

    This is clearly to much overhead fur such a tiny function, be we
    use this as an example for documenting an application of the code-
    generating class ``TableGenerator``.


    This class hat the following public member functions

    ================  ==================================================
    ``tables()``      Returns a dictionary containing the tables
                      required for generating the C code
    ``directives()``  Returns a dictionary containing directives 
                      required for generating the C code    
    ================  ==================================================

    This class is designed for cooperation with an instance of
    class ``TableGenerator``. A similar class should provode at least 
    the member functions ``tables()`` and ``directives()``.

 
    In order to create the table ``lsbit_table``, put e.g. the lines::

        static const uint8_t MAT24_LSBIT_TABLE[32] = { 
        // %%TABLE Mat24_lsbit_table, uint8
        };

    into the source file. This creates the table ``MAT24_LSBIT_TABLE``
    for  computing the lsbit quickly. 
 
    This class provides a function with name ``MAT24_LSBIT`` for 
    generating  C code that compute the lsbit of ``x``. After generating 
    the table given above, the expression ``{MAT24_LSBIT:x}`` in the input 
    file generates a very efficient expression that evaluates to the 
    lsbit of ``x``. 
    Note that we use standard python formatting syntax, with 
    ``MAT24_LSBIT`` refering to a member  function ``generate`` of class 
    ``Lsbit24Function``that generates the appropriate C code for 
    computing the lsbit of the 24-bit integer ``x``.
 
    The entry 

          "MAT24_LSBIT" : UserFormat(cls.generate, "s", prefix=True)
    
    in the attribute ``tables`` means the following:

    We create an instance of class ``mmgroup.generate_c.UserFormat``
    that behaves like member function ``generate`` of this class.
    Furthermore, this instance has a ``.format`` method so that
    the string ``'{MAT24_LSBIT:x}'`` in the source file for code 
    generation  evaluates to the string returned by applying member 
    function  ``Lsbit24Function.generate`` to the argument ``'x'``.
    Here function ``Lsbit24Function.generate`` returns a string with
    very efficient C code for computing the lsbit of variable ``x``.
    
    Note that python and C name spaces are strictly separated. So the
    table coded above has the name  ``Mat24_lsbit_table`` in python.
    In other words, the python object ``Mat24_lsbit_table`` contains the 
    data that make up the table. In the generated C file, the table 
    may have another name (which is ``MAT24_LSBIT_TABLE`` in the example
    above). The instance of the code generation class provides a 
    dictionary ``names`` that maps python names to C names. Dictionary 
    ``names`` must be passed to the member function ``generate`` as the 
    first parameter, so the member function ``generate`` can 
    infer the C name of the table from its python name 
    ``Mat24_lsbit_table``. The keyword argument ``prefix=True`` makes
    sure that dictionary ``names`` is prepended to the list of arguments
    for function ``Lsbit24Function.generate`` given by the user.

    The second argument ``"s"`` of ``UserFormat()`` means that the
    argument (which is ``x`` in the example ``{MAT24_LSBIT:x}``) is
    interpreted as the string ``"x"``. See documentation of class 
    ``UserFormat()`` for details.

    """
    CONST,  table = _make_de_bruijn_sequence()
    hCONST = format_item(CONST, "uint32")
   
    @classmethod
    def compute(cls, x):
        """Return index of least significant bit of an integer x

        If x & 0xffffff is zero then 24 is returned.        
        """
        x = int(x & 0xffffffff)
        return int(cls.table[(((x & -x) * cls.CONST) >> 26) & 0x1f])



    @classmethod
    def generate(cls, names, x, arg_is_pwr2 = False):
        """generate C code computing  lsbit24(x)

        Here 'x' is a string representing the input as C code.

        When calling this function from python, names is the
        dictionary for translating python names of tables into
        C names to be used by the generated code, see documention
        of this class.

        If 'arg_is_pwr2' then more efficient code is generated,
        assuming that x is a power of two.

        """
        if not arg_is_pwr2:
             x = "({0}) & (0-({0}))".format(x)
        return  """{1}[({2} *  \\
        ({0}) >> 26) & 0x1f]""".format(x, 
            names["Mat24_lsbit_table"], cls.hCONST)


    @classmethod
    def tables(cls):
        """Return dictionary of tables for class TableGenerator

        """ 
        return { 
            "Mat24_lsbit_table" :  cls.table[:], 
            "Mat24_lsbit_const" :  cls.CONST, 
            "MAT24_LSBIT" : UserFormat(cls.generate, "s", prefix=True),
        }

    @classmethod
    def directives(cls):
        """No directives here""" 
        return {} 
        


############################################################################
############################################################################
# Functions for analyzing basis vectors of the Golay code  
############################################################################
############################################################################
 
def split_golay_codevector(v, check=None):
    """Split Golay code vector into 'blackwhite' and 'colored' part
 
    The 'blackwhite' part has the same value (0 or 1) in the bottom
    three entries of each MOG column. 

    The 'colored' part has 0 or 2 entries with value 1 in the bottom
    three entries of each MOG column, and no entry with value 1 in
    the top row. 

    Returns pair of Golay code vectors ('blackwhite', 'colored ')
    with the sum of these vectors equal to v.

    If 'check' is set, it must behave lie an instance of class GolayCode.
    """
    color = [0, 6, 5, 3, 3, 5, 6, 0] 
    col = 0
    for i in range (1, 25, 4):
         col |= color[(v >> i) & 7] << i
    v ^= col
    if check:
        cn = check.vect_to_vintern(v) | check.vect_to_vintern(col) 
        assert cn & 0xfff == 0, "not a Golay code word"
    return v, col





############################################################################
############################################################################
# Functions for changing Colay code rep to permutation rep in Mat24
############################################################################
############################################################################


def generate_golay24_decode(table_names, s): 
        """Generate code for converting a Golay code word 's'.

        Generate C code for computing the 'vect' representation of a
        Golay code word 's' given in 'gcode' representation. See class
        Golay24.GolayCode for background.

        The generated code is and expression containing the input 's'
        and conversion tables from class Golay24.GolayCode.
        """ 
        return """({1}[({0} << 4) & 0xf0]
       ^ {2}[({0}  >> 4) & 0xff])
""".format(s, table_names["Mat24_dec_table1"],
                table_names["Mat24_dec_table2"])



class MatrixToPerm(object):
    """Support for conversion of Mat24 from matrix to permutation rep

    A basis of the Golay code is given with the constructor.
    The basis should contain 6 colored and 6 blackwhite basis vectors,
    exactly one of the blackwhite basis vectors should be odd.

    The basis is augmented by up to two colored vectors of weight 12,
    such that their XOR sum has also weight 12. Call these vectors
    c0 and c1. Then we calculate c0 & c1, c0 & ~c1, ~c0 & c1. These 3
    vectors, and the odd blackwhite vector and its complement form a 
    set of 5 vectors. Any singelton vector in the MOG an be obtained 
    as an intersection of one of these 5 vectors with a vector that 
    has ones in exactly one column of the MOG.

    Then the even blackwhite basis vectors are processed with 
    operators 'and' and 'not' so that for each MOG column there is a
    vector t[i]  that contains ones in precisely one column. 

    If an image of the standard basis of the Golay code is given,
    we can perform these XOR, AND and coplement operations with the 
    transformed basis vectores, yielding the transformed singletons.

    These transformed singletons are just the images of an element
    of the Mathieu group Mat24 considered ad a permutation. 

    Member functions tables(), directives() and bytelen() and also
    compute() and generate() follow the conventions in 
    class Lsbit24Function. 
    """
    def __init__(self, basis, verbose = 0):
        """Set basis for the Golay code and calculate tables
      
        separate (even and odd) blackwhite and colored basis vectors
        """ 
        self.basis = list( basis )
        self.verbose = verbose
        self.split_basis()   
        self.augment_colored_basis() 
        self.make_aux_table()
        self.test(verbose=0)
        #print (self.generate("b", "p"))


    def split_basis(self):
        """Split basis into blackwhite and colored vectors.

        self.blackwhite =  list of indices of even blackwhite basis vectors
        self.colored =     list of indices of colored basis vectors
        self.odd_v =       index of the (single) odd blackwhite basis vector
        """
        self.blackwhite = []
        self.colored = []
        self.odd_v = None
        for i, v in enumerate(self.basis):
            bw, col =  split_golay_codevector(v)
            if v == col:
                self.colored.append(i)
            if v == bw:
                if bw24(v & 0x111111) & 1:
                    self.odd_v = i
                else:
                    self.blackwhite.append(i)


    def and_basis(self, z):
        """Compute and operation (with complements) on basis vectors

        z is a list of pairs (index, sign). self.basis[index] is complemented
        if and only if sign is 1, and all these (possibly complemented)
        basis vectors are combined with the and operation.
        The result is returned. 
        """
        return  reduce(__and__, [int(self.basis[pos]) ^ -sign for (pos,sign) in z])
       

    def augment_colored_basis(self):
        """Augment self.basis by addding some vectors to the basis.
        
        The augmented vectors are XOR combinations of the basis vectors.
       
        The following tables are computed:

        self.basis_augment is a list of entries, where each entry represents
        an augmented basis vector as a list of the indices of basis vectors.
        The corresponding basis vectors are to be XORed to obtain the new
        basis vector.    

        self.col_aux_table is a list of 5 entries. Here entry i is the
        instruction how to compute the i-th auxiliary vector aux[i] from 
        the basis vectors with member function  and_basis.

        self.out_table is a table of 24 auxiliary vectors. For 
        obtaining the i-th singleteon we have to compute

            col[i >> 2] &  aux[self.out_table[i]] .
        
        Here col[i] is the vector that casts out the entries in the
        i-th column in the MOG. A table for computing col[i] is prepared
        in member function make_aux_table.
        """
        colored = self.colored
        basis_list = []
        value_list = []
        colored_6 = []
        self.basis_augment = []
        done = False
        for pos in lmap( bits2list,iter_ascending_bitweight(1 << len(colored))):
            basispos = [ colored[i] for i in pos]
            value = reduce(__xor__, [self.basis[p] for p in basispos], 0)
            if bw24(value) == 12:
                 basis_list.append(basispos)
                 value_list.append(value)
                 for i, val in enumerate(value_list):
                     if bw24(value ^ val) == 12:
                         for j in [i, -1]:
                              b,v = basis_list[j], value_list[j]
                              if len(b) > 1:
                                  self.basis_augment.append(b)
                                  b = [len(self.basis)] 
                                  self.basis.append(v)
                                  self.colored.append(b)
                              colored_6.append(b[0])
                         done = True
                         break
            if done: break

        c0, c1 = colored_6[0], colored_6[1]
        self.col_aux_table = [
             [(self.odd_v,0)],  [(self.odd_v,1)],
             [(c0,0), (c1,0)], [(c0,0), (c1,1)], [(c0,1), (c1,0)]
        ]
        col_aux_values = lmap(self.and_basis, self.col_aux_table)
        self.out_table = []
        for i in range(24):
            mask = 0xf << (i & -4)
            for j, v in enumerate(col_aux_values):
                 if v & mask == 1 << i:
                     self.out_table.append(j)
                     break
        assert len(self.out_table) == 24
        if self.verbose:
            print("augmented basis for Golay code to perm rep of Mat24") 
            print(lmap(hex, self.basis))
            print("augment:", self.basis_augment, "col", colored_6)
            print("col_aux_table:", self.col_aux_table)
            print("values:", [hex(x & 0xffffff) for x in col_aux_values])
            print("t_out:", self.out_table)



    def make_aux_table(self):
        """Prepare table for casting out the i-th column in the MOG

        Let col[i], i=0,...,5 is the vector that casts out the entries 
        in the i-th column in the MOG. col[i] is computed from the 
        basis vectors in self.basis with function and_basis().
         
        self.bw_aux_table is a table with 6 entries, where entry i
        is used to compute  col[i] by member function make_aux_table.
        """
        blackwhite = self.blackwhite
        lbw = len(blackwhite) 
        bw_combinations = [None] * 6
        bw_dict = dict(zip([0xf << (4*i) for i in range(6)], range(6))) 
        done = False
        for k in range(1,lbw+1):
            for sgn in iter_ascending_bitweight((1 << k)-1):
                signs = [ (sgn >> i) & 1 for i in range(k)]
                for pos in lmap(bits2list,iter_bitweight(k,1 << lbw)):
                    basispos = [ blackwhite[i] for i in pos]
                    z = list(zip(basispos,signs))
                    res = self.and_basis(z)
                    try:
                        bw_combinations[bw_dict[res]] = z
                        del  bw_dict[res]
                        done = not None in bw_combinations
                        if done: break
                    except:
                        pass
                if done: break
            if done: break
        self.bw_aux_table = bw_combinations 
        if self.verbose:
            print("bw_comb",  bw_combinations ) 


    def compute(self, basis):
        """Calculate the permutation from a (transformed) basis.

        Here we calculate the singletons from the transformed basis.
        The same calculations with the standard basis self.basis yields 
        the singletons  0,...,23 in that order. So the singletons from 
        the transformed basis are the images of 0,...,23. 

        The tables self.basis_augment, self.col_aux_table, 
        self.bw_aux_table, self.out_table[i] are used.

        Singleton i is represented as 2**i, so at the end we have to take
        binary logarithms.
        """
        basis = basis[:12]
        for l in self.basis_augment:
            basis.append(reduce(__xor__, [basis[x] for x in l], 0))
        t = []
        for z in self.col_aux_table + self.bw_aux_table:
            t.append(reduce(__and__, [int(basis[i]) ^ -s for (i,s) in z]))
        out = []
        st = len(self.col_aux_table)
        for i in range(24):
            m =  t[st + (i>>2)] & t[self.out_table[i]]
            out.append( v2(m) )
        return out

    def test(self, verbose = 0):
        """A simple self test. 

        self.compute(self.basis) must return [0,...,23].
        """
        basis = self.basis[:12]
        result = self.compute(basis)
        assert result == list(range(24)), result 
        if verbose:
            print("Test of class Matrix_to_perm passed." )


    def generate_and(self, z, basis_name):
        """Return C code for the expession self.and_basis(z)

        basis_name is the name of the array for the basis         
        """
        v = [[], []]
        for i,s in z:
            v[s].append("%s[%d]"  % (basis_name, i))
        neg = " | ".join(v[1])
        if len(neg): v[0].append("~(%s)" % neg)
        return " & ".join(v[0])

    def generate(self, names, in_name, out_name, 
             basis_name = "ba",  temp_name = "t"):
        """Generate C code roughly corresponding to self.calculate()

        'in_name' is the name of the input array representing the transformed
        basis. Each entry of that input array represent a Golay code word.
        In contrast to function self.calculate(), the Golay code vectors in
        the input array are given in 'gcode' represention, see documentation
        of class GolayCode. They are converted to vector representation first.
        
        'out_name' is the name of the output array representing the 
        resulting element of the Mathieu group Mat24 as a permutation.

        The remaining parameters are optional names for some internal arrays.
        """
        self.names = names
        s =  """uint_fast32_t _i, {0}[14], {1}[{2}];
for (_i=0; _i < 12; ++_i) {0}[_i] =
    {3};
""".format(
           basis_name, temp_name, 
           len(self.col_aux_table + self.bw_aux_table),
           generate_golay24_decode(names, "%s[_i]" % in_name)
        )


        for i,l in enumerate(self.basis_augment):
            v = " ^ ".join( ["%s[%d]" % (basis_name, j)  for j in l] )
            s += "%s[%d] = %s;\n" % (basis_name, i+12, v)
        for i, z in enumerate(self.col_aux_table + self.bw_aux_table):
            v = self.generate_and(z, basis_name)
            s += "%s[%d] = %s;\n" % (temp_name, i, v)

        t = [sum( x << (4*j) for j, x in enumerate(self.out_table[i:24:4]) )
                  for i in range(4)  ]
        s1 = "uint_fast32_t _w, _k = %s[_i + %d];\n" % (
                               temp_name, len(self.col_aux_table))
        for j in range(4):
            s1 += "_w =  %s[(%s >> (_i << 2)) & 0xf] & _k;\n" % (
                  temp_name, hex(t[j]) ) 
            # Next put *<out_name>++ = lsbit24(_w)
            lsb =  Lsbit24Function.generate(names, "_w", arg_is_pwr2=True)
            s1 +=  "*%s++ = %s;\n" % (out_name, lsb);

        s += "for (_i = 0; _i < 6; ++_i) {\n"
        s += prepend_blanks(s1,3)
        s += "}\n"
        s += "%s -= 24;\n" % out_name
        return s  
       
    def tables(self):
        """This class generates no tables"""
        return {}

    def directives(self):
        """Code generation function has name "MAT24_MATRIX_TO_PERM". """
        return   {
            "MAT24_MATRIX_TO_PERM" : UserDirective(self.generate,"ssss", 1),
        }
 

