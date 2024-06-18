r"""Module ``mat24tables`` contains the class ``Mat24Tables``

Class ``Mat24Tables`` is a base for class ``Mat24`` in module
``mmgroup.dev.mat24.mat24_ref``. It contains the basis of the Golay code 
and of (a set of representatives of) the Golay cocode. It also contains 
tables for fast conversion of a vector in ``V = GF(2)**24`` from the 
standard to the internal representation and vice versa. That class also 
contains a table of the syndromes of all ``2048`` odd elements of the 
cocode. The ``759`` octads, i.e. Golay code word of length ``8`` are 
numbered from ``0`` to ``758``. Class ``Mat24Tables`` provides tables
for computing the number of a octad from a Golay code word representing
an octad and vice versa. 

Class ``Mat24Tables`` also contains python versions of some C functions
in file ``mat24_functions.c`` using these tables.

"""
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

# python 3 string formatting see:
# https://pyformat.info/


import types
import sys
import re
import os
from random import randint, shuffle
from operator import __or__, __xor__, __and__
from functools import reduce

import numpy
from numpy import array, zeros, uint8, uint16, uint32

from mmgroup.bitfunctions import pivot_binary_high, bit_mat_transpose
from mmgroup.bitfunctions import bit_mat_inverse, bit_mat_basis_spanned
from mmgroup.bitfunctions import reverse24, lin_table, bw24
from mmgroup.bitfunctions import iter_bitweight, lmap
from mmgroup.bitfunctions import bitweight, bitparity, bits2list

from mmgroup.dev.mat24.mat24aux import Lsbit24Function


############################################################################
############################################################################
# Lexicographical basis of the Golay code
############################################################################
############################################################################

lex_golay_code_basis = None

def make_lex_golay_code_basis():
    r"""Return lexicographical basis of the Golay code.

    This basis is computed from the definition of the Golay code as a 
    lexicographic code, see :ref:`CCN+85`, Chapter. 11,  section 17.  

    This computation takes a long time. So we better save the result
    in the global variable 'lex_golay_code_basis'.

    This function is mainly of historical interest. It has been run
    once for a proof that our definition of the Golay code agrees with 
    that in :ref:`CCN+85'. 
    """ 
    global lex_golay_code_basis 
    if not lex_golay_code_basis is None:
        return lex_golay_code_basis
 
    print("Computing the lexicographical basis of the Golay code")
    bitweights = numpy.zeros(1 << 24, dtype = numpy.uint8)
    for i in range(24):
        bitweights[1<<i:2<<i] = 1 + bitweights[0:1<<i]
    # basis of already generated part of the Golay code 
    basis = []
    # codes[:len_code] is the list of all code worde generated from basis
    codes = numpy.zeros(4096, dtype = numpy.uint32)
    len_codes = 1
    for d in range(7,24):
        # For 8<=d1<=24 loop through all d1-bit integers i of bit weight 8
        for i in range(1<<d, 2<<d):
            if bitweights[i] == 8:
                bad = False
                # set 'bad' if an element of codes[] has distance < 8 from i
                for c in codes[:len_codes]:
                    bad = bitweights[c ^ i] < 8
                    if bad: break
                # If not bad: append i to basis and update codes
                if not bad:
                    basis.append(i)
                    if d < 24:
                        for j in range(len_codes):  
                            codes[len_codes+j] = i ^ codes[j]
                        len_codes <<= 1
                    break
    lex_golay_code_basis = basis
    print("Lexicographical basis of the Golay code has been computed")
    del bitweights
    del codes
    return basis

############################################################################
############################################################################
# The hexacode (used for making up a basis of the Golay code)
############################################################################
############################################################################


class HexacodeToGolay(object):
    r"""This class deals with the mapping from the Hexacode to the MOG
    
    This mapping is described in :ref:`CCN+85`, ch. 11, section 5. Rows of 
    the MOG are  numbered from 0 to 5 (left to right). Columns of the MOG 
    are numbered from 0 to 3 (top to bottom). MOG entries are represented 
    as integers with bit 4*i+j (of valence 2**(4*i+j)) corresponding to the 
    MOG entry in row i, column j::


                         +------------------------+
                         |  0   4   8  12  16  20 |                 
                         |  1   5   9  13  17  21 |                 
                         |  2   6  10  14  18  22 |                 
                         |  3   7  11  15  19  23 |                 
                         +------------------------+

    A hexacode word (or more precisely, an element of GF(4)^6) is 
    represented as a list of six entries of value 0,...,3, so that 0  
    and 1 are in the subfield GF(2) of GF(4) and an entry in column j
    has 'score' j, see  :ref:`CCN+85`, ch. 11, section 5.
    """
    HEXA_ENCODE = [0x0,0xc,0xa,0x6] 
         # Encoding of an even interpretation with score j:
         # HEXA_ENCODE[j] << (4*i) is the unique MOG vector with an 
         # even number of nonzero entries in column i (none in any other
         # column or in row 0) and with 'score' j in column i.

    HBASIS =  [
        [1,0,0,1,3,2], [0,1,0,1,2,3], [0,0,1,1,1,1],
        [2,0,0,2,1,3], [0,2,0,2,3,1], [0,0,2,2,2,2],
    ] 
         # This the GF(2) basis of the hexacode in [SPLG], ch. 11.

    HISTORICAL_REVERSED_BASIS = [
        0xff, 0xf0f, 0x3333, 0x5555, 0x9669, 0x30356, 
        0x50563, 0x9063a, 0x111178, 0x21121d, 0x41144e, 0x811724
    ]
         # Golay code basis taken from an older implementation for 
         # compatibility checks.


    basis_ = None
         # store basis here after the 1st call to cls.basis()

    @classmethod
    def show_basis(cls,basis):
        """display a basis of the Golay code or cocode"""
        for x in basis:
            print( "%06x" % x, end = " ")
        print( "" )

    @classmethod
    def hexacode_vector_to_mog(cls, h):
        """Maps a hexacode word h to a Golay code word.

        Entry h[i] is mapped to an even interpretation of MOG column i, 
        with score h[i], and the top bit of column i being zero, 
        see [SPLG], ch.11, section 5.         
        """  
        l = [ cls.HEXA_ENCODE[x]<<(4*i) for i,x in enumerate(h) ]
        return sum(l)

    @classmethod
    def mog_std_octad(cls, i):
        """Returns the MOG octad with ones in row 0 and column i

        The entry at the intersection of row 0 and column i is zero.
        """  
        return 0x111111 ^ (0xf << (4*i))

    @classmethod
    def historical_check(cls, GolayBasis, generate = False):
        """Checks the Golay code basis against basis a stored older basis

        The basis is checked against basis cls.HISTORICAL_REVERSED_BASIS.
        If 'generate' is set, the basis cls.HISTORICAL_REVERSED_BASIS
        is recomputed by function  make_lex_golay_code_basis().
        """
        ref = lmap(reverse24, cls.HISTORICAL_REVERSED_BASIS)  
        ref = pivot_binary_high(ref)[0] 
        assert  GolayBasis == ref
        if generate:
            assert HISTORICAL_REVERSED_BASIS ==  make_lex_golay_code_basis()


    @classmethod
    def any_golay_basis(cls, verbose = False):
        """Returns a basis of the Golay code.

        This basis is directly calculated from the hexacode.
        """
        basis = list(map(cls.hexacode_vector_to_mog, cls.HBASIS))
        basis += [cls.mog_std_octad(i) for i in range(6)]
        basis = pivot_binary_high(basis)[0]
        cls.historical_check(basis)
        return basis


    @classmethod
    def check_cocode_basis(cls, cocode_basis):
        """Checks correctness of a basis of a transversal of the Golay cocode"""
        assert len(cocode_basis) == 12, cocode_basis
        bit_mat_inverse( cls.any_golay_basis() + cocode_basis )
            # this inverse must exist!!

    @classmethod
    def better_cocode_basis(cls, verbose=False):
        """Returns a good basis of a traversal of the Golay cocode.

        The asserted properties of the cocode basis are as in the
        main comment of module mat24.py.
        """
        def e(i,a):
            v = [0]*6
            v[i] = a
            return v
        hexa = lmap(e, [2,2,1,1,0,0], [2,1,2,1,2,1])
        #basis = [0x111111,0x011111,0x110000,0x101000,0x100100,0x100010,]
        basis = [0x111111,0x1,0x110,0x1010,0x10010,0x100010,]
        basis +=  lmap(cls.hexacode_vector_to_mog, hexa)
        basis = basis[2:] + basis[:2]        
        cls.check_cocode_basis(basis)
        return basis



    @classmethod
    def code_basis_from_cocode_basis(cls, cocode_basis):
        """return orthogonal Golay code basis of the given cocode basis"""
        a = bit_mat_inverse(cls.any_golay_basis() + cocode_basis)
            # this inverse must exist!!
        res =  bit_mat_transpose(a)[12:]
        assert cls.any_golay_basis() == pivot_binary_high(res)[0]
        return res

    @classmethod
    def better_code_basis(cls):
        """returns a good basis of the Golay code

        The asserted properties of the Golay code basis are as in the
        main comment of module mat24.py.
        """
        if getattr(cls, "_BETTER_CODE_BASIS", None) is None:
            cls._BETTER_CODE_BASIS = cls.code_basis_from_cocode_basis(
                                                  cls.better_cocode_basis())
        return cls._BETTER_CODE_BASIS
        
    @classmethod
    def basis(cls):
        """Return basis for the GF(2)**24 to be used in the Golay code

        basis [0:12] is the basis for a transversal of the Golay code
        given by method better_cocode_basis().
        basis [12:24] is the basis of the Golay code in GF(2)**24 
        given by method better_code_basis().
        """
        if cls.basis_ is None:
            basis =  cls.better_code_basis()
            cobasis =  cls.better_cocode_basis()
            cls.basis_ = numpy.array(cobasis + basis, dtype = uint32)
        return cls.basis_

   
############################################################################
############################################################################
# Conversion tables for the Golay code
############################################################################
############################################################################


def encoding_tables(basis):
    table0 = bit_mat_basis_spanned(basis[0:8], uint32)
    table1 = bit_mat_basis_spanned(basis[8:16], uint32)
    table2 = bit_mat_basis_spanned(basis[16:], uint32)
    return table0, table1, table2


def make_syndrome_table(recip_basis):
     """Generates a table for fast calculation of syndrome of a bit vector

     Table input is a 11-bit number representing an odd cocode word
     ``c``. Since basis vector 11 is the only vector with odd parity, 
     exactly  the cocode words with number ``c  >= 0x800`` have odd 
     parity. In table entry ``c & 0x7ff``  we list the unique cocode 
     word of weight 1 or 3 equivalent to the cocode word with number
     ``c``in the following format:

              i   +  (j << 5)  +  (k << 10) ,

     representing the word with ones at position i, j and k. If that
     cocode word has weight 1, we put j = k = 24.

     We set bit 15 of entry ``e`` in the table to one iff the (even)
     cocode word with number ``e`` has bit weight 2.

     Input recip_basis is the reciprocal basis of the chosen basis
     for the Golay code. (Lower 12 basis vectors must span a 
     transversal of the Golay code. All but the least of these 
     lower basis vector must have even parity.)
     """ 
     syndrome_table = zeros(0x800, dtype=uint16)
     C1 = (24 << 5) | (24 << 10) 
     for i in range(24):
         assert recip_basis[i] & 0x800
     rb = recip_basis & 0x7ff
     for i in range(24):
         b_i = rb[i]
         syndrome_table[b_i] ^=  i | C1 
         for j in range(i+1,24):
             b_j = b_i ^ rb[j]
             syndrome_table[b_j] ^= 0x8000
             for k in range(j+1,24):
                 b_k = b_j ^ rb[k]
                 syndrome_table[b_k] ^= i | (j<<5) | (k<<10)  
     assert min(syndrome_table & 0x7fff) > 0
     return syndrome_table
                 

ODD_OCTADS_SPECIAL = 1


ODD_OCTADS_DICT = {
    7: [0, 1, 2], 11: [1 ,0, 2], 13: [1, 2, 0], 14: [2, 1, 0]
}

def octad_to_bitlist(vector):
    """Convert an octad to a bit list

    Here vector must be a vector of bitwight 8 representing an octad.
    The function returns the corresponding list of bits.
    """
    vector = int(vector)
    assert bitweight(vector) == 8
    global _shuffled_octad_list
    if ODD_OCTADS_SPECIAL > 1:
        # Then we shuffle the entries of an octad for test purposes.
        s = 1 + (vector % 19) % 7
        lst = bits2list(vector)
        return lst[s:] + lst[:s]
    if not ODD_OCTADS_SPECIAL or bitweight(vector & 15) & 1 == 0:
        return bits2list(vector)
    else:
        for i in range(0, 24, 4):
            v3 = vector & (15 << i)
            if v3 >> i in ODD_OCTADS_DICT:
                first, last = bits2list(v3), bits2list(vector & ~v3)
                seq = ODD_OCTADS_DICT[v3 >> i]
                return [first[seq[i]] for i in range(3)] + last 
        ERR = "Error in converting odd octad to bit list"
        raise ValueError(ERR)



def make_octad_tables(basis):
    """Return tables for encoding an decoding octads.

    Octads are numbered in the lexicographical order as given by the 
    numbers of the corresponding Golay codewords in 'gcode' 
    representation. 

    The function returns a triple  

       (oct_enc_table, oct_dec_table, suboctad_table)

    Given octad o, the corresponding Golay code word in gcode
    representation is oct_dec_table(o), for 0 <= 0 < 759.

    Given an octad (or a complement of an octad) v in 'gcode' 
    representation, the number of the corresponding octad is:

       (oct_enc_table[v & 0x7ff] >> 1) ,


    The vector v is a (possibly complemented) octad if the following 
    condition holds:

       oct_enc_table[v & 0x7ff] & 0x8000 == 0  .

    Vector v is not a complemented octad if in addition we have:

       ((v >> 11) ^ oct_enc_table[v & 0x7ff]) & 1 == 0

    suboctad_table[8*o+j] is the j-th entry of the octad with
    number o, for 0 <= o < 759. Here the entries of octad are the 
    the postions of the bits being set in that octad o.
    """
    codewords = lin_table(basis[:11])
    oct_dec_table = numpy.zeros(759, dtype = uint16)
    octad_table = numpy.zeros(759 * 8, dtype = uint8)
    octad = 0
    oct_enc_table = numpy.full(2048, 0xffff, dtype = uint16)
    for gcode, vector in enumerate(codewords[:2048]):
        weight = bw24(vector)
        if weight in [8, 16]:
            oct_dec_table[octad] = gcode + ((weight & 16) << 7)
            oct_enc_table[gcode] = ((weight - 8) >> 3) + 2 * octad
            oct_vector = vector if weight == 8 else vector ^ 0xffffff
            blist = octad_to_bitlist(oct_vector)
            octad_table[8 * octad : 8 * octad + 8] = blist
            octad += 1
    assert octad == 759
    return oct_enc_table, oct_dec_table, octad_table




def make_octad_index_table():
    """Return numpy array A for indexing suboctads octads

    Let 0 < i < 64 be the number of a suboctad. Let b[0],...,b[7]
    be the 8 basis vectors of the Golay cocode the make up the
    octad to which the suboctad refers, in the order given by
    the table created by function make_suboctad_table().

    Then b[4*i] + b[4*i+1] + b[4*i+2] + b[4*i+3] is cocode 
    vector corresponding to the suboctad with number i.
    """
    bl = []
    for i in range(64):
         j = (i << 1) + bitparity(i)
         if (bitweight(j) > 4):
             j ^= 0xff
         blist = bits2list(j)
         while len(blist) < 4:
              blist += [0, 0]
         assert len(blist) == 4
         bl += blist
    return numpy.array(bl, dtype = numpy.uint8)
         

class Mat24Tables(object):
    basis = HexacodeToGolay.basis()
    recip_basis = numpy.array(bit_mat_inverse(basis), dtype=uint32) 
    enc_table0, enc_table1, enc_table2 = encoding_tables(recip_basis)
    dec_table0, dec_table1, dec_table2 = encoding_tables(basis)
    syndrome_table = make_syndrome_table(recip_basis)
    B = basis[12:]
    oct_enc_table, oct_dec_table, octad_table = make_octad_tables(B)
    octad_index_table = make_octad_index_table()
    suboctad_weights = sum([1 << i for i in range(64)
                         if bitweight(i) in (1,2,5,6)])

    ###########################################################################
    # Some general bit operations
    ###########################################################################

    @staticmethod
    def lsbit24(v1):
        """Return position of least significant bit of the integer v1

        In case v1 & 0xffffff == 0 we return 24.
        """
        return Lsbit24Function.compute(int(v1))
        
    @staticmethod
    def bw24(v1):
        """Return the bit weight of the lowest 24 bits of v1"""
        v1 = int(v1)
        v1 = (v1 & 0x555555) + ((v1 & 0xaaaaaa) >> 1)
        v1 = (v1 & 0x333333) + ((v1 & 0xcccccc) >> 2)
        v1 = (v1 + (v1 >> 4)) & 0xf0f0f
        return (v1 + (v1 >> 8) + (v1 >> 16)) & 0x1f 

    ###########################################################################
    # Basis of Golay cocode and code
    ###########################################################################

    @classmethod
    def golay_basis(cls):
        """Return basis of Golay cocode followed by basis of Golay code

        The basis is returned as a list of 12 + 12 integers
        """
        return [int(cls.basis[i]) for i in range(24)]


    @classmethod
    def basis_weights_8(cls):
        """Return weights of basis of Golay code

        Let ``w[i]`` be the bit weight of the ``i``-th basis
        vector of the Golay code.
 
        The function returns a integer ``v`` such that bit ``i``  
        of ``v`` is equal to ``(w[i] >> 3) & 1``.
        """
        w = [cls.bw24(cls.basis[i + 12]) for i in range(12)]
        assert reduce(__or__, w) & 7 == 0
        bits = [(b >> 3) & 1 for b in w]
        return sum((b << i for i, b in enumerate(bits)))



    ###########################################################################
    # Conversion between bit vectors of GF(2)**24
    ###########################################################################
 
    @staticmethod
    def vect_to_bit_list(v1):
        """Return lists of positions of 1-bits and 0-bits of v1 

        Let w be the bit weight of v1 & 0xffffff, i.e. number n of bits 
        of v1 at positions < 24 equal  to one. Then the ordered bit
        positions where the corresponding bit of v1 is 1 are stored in
        a list  l1.
       
        Then v1 & 0xffffff has  24 - w  zero bits. The ordered list of 
        the positions of these zero bits is stored in l0. 

        The function returns the pair len(l1), l1 + l0.
        This return value is for compatibility with the corresponding
        C function.
        """
        v1 = int(v1)
        l1 = [ i for i in range(24) if v1 & (1 << i) ]
        l0 = [ i for i in range(24) if not v1 & (1 << i) ]
        return len(l1), l1 + l0


    @staticmethod
    def extract_b24(v1, u_mask):
        """Extract the bits of 24-bit vector v1 given by the mask u_mask
   
        If u_mask has bits equal to one at positions i_0, i_1, ..., i_k
        (in ascending order) then the bit of v1 at position i_j is copied 
        to the bit at position j of the return value for j = 0,...,k.
        """
        v1 = int(v1)
        u_mask = int(u_mask)
        l = [(v1 >> i) & 1 for i in range(24) if (u_mask >> i) & 1] 
        return sum(x << i for i, x in enumerate(l))

    @staticmethod
    def spread_b24(v1, u_mask):
        """Spread the bits of 24-bit vector v1 according to the mask u_mask

        If u_mask has bits equal to one at positions i_0, i_1, ..., i_k,
        (in ascending order) then the bit of v1 at position i is copied 
        to the bit at position i_j of the return value for j = 0,...,k.
        """
        v1 = int(v1)
        u_mask = int(u_mask)
        l = [i for i in range(24) if (u_mask >> i) & 1] 
        return sum( ((v1 >> i) & 1) << x for i, x in enumerate(l) )




    ############################################################################
    # Conversion between representations of GF(2)**24, Golay code, etc.
    ############################################################################


    @classmethod
    def vect_to_vintern(self, v1):        
        return  int(    self.enc_table0[v1 & 0xff]
                   ^ self.enc_table1[(v1 >> 8) & 0xff]
                   ^ self.enc_table2[(v1 >> 16) & 0xff]  )


    @classmethod
    def vintern_to_vect(self, v1 ):
        return  int(    self.dec_table0[v1 & 0xff]
                   ^ self.dec_table1[(v1 >> 8) & 0xff]
                   ^ self.dec_table2[(v1 >> 16) & 0xff]  )

    @classmethod
    def vect_to_cocode(self, v1 ):
        return self.vect_to_vintern(int(v1)) & 0xfff

    @classmethod
    def gcode_to_vect(self, v1):
        v1 = int(v1)
        return  int(   self.dec_table1[(v1 << 4) & 0xf0]
                   ^ self.dec_table2[(v1 >> 4) & 0xff]  )

    @classmethod
    def cocode_to_vect(self, c1):
        return self.vintern_to_vect(c1  & 0xfff)

    @classmethod
    def vect_to_gcode(self, v1):
        cn = self.vect_to_vintern(v1)
        assert cn & 0xfff == 0, "not a Golay code word"
        return cn >> 12


    @classmethod
    def vect_to_octad(self, v1):
        return self.gcode_to_octad(self.vect_to_gcode(v1))


    @classmethod
    def gcode_to_octad(self, v1, u_strict = 1):
       v1 = int(v1)
       y = int(self.oct_enc_table[v1 & 0x7ff])
       if y >> 15 or ((v1 >> 11) ^ y) & u_strict & 1:
           raise ValueError("Golay code vector is not an octad")
       return y >> 1 


    @classmethod
    def octad_to_gcode(self, u_octad):
        u_octad = int(u_octad)
        assert 0 <= u_octad < 759, "not a Golay code octad"
        return int(self.oct_dec_table[u_octad]) & 0xfff

    @classmethod
    def octad_to_vect(self,u_octad):
        u_octad = int(u_octad)
        return self.gcode_to_vect(self.octad_to_gcode(u_octad))

    ###########################################################################
    # Golay code syndromes and weights
    ###########################################################################


    @classmethod
    def cocode_syndrome(cls, c1, u_tetrad=24):
        """Return equivalent Golay cocode to 'c1' of length at most four.

        Same as mat24_syndrome(mat24_cocode_to_vect(c1, u_tetrad)).  
        """
        c1 = int(c1)
        u_tetrad = int(u_tetrad)
        if not 0 <= u_tetrad <= 24:
            raise ValueError("Bad argument in function cocode_syndrome()")
        bad = (u_tetrad >= 24) & ((c1 >> 11) + 1) & 1
        u_tetrad -= (u_tetrad + 8) >> 5          # change u_tetrad == 24 to 23
        y = - (((c1 >> 11) + 1) & 1)             # y = 0 if c1 is odd else -1
        c1 ^= int(cls.recip_basis[u_tetrad]) & y # make c1 odd, if c1 is odd:
        y  &=   1 << u_tetrad                    #     set bit u_tetrad in y
        syn = int(cls.syndrome_table[ c1 & 0x7ff ])
        syn = (1 << (syn & 31)) | (1 << ((syn>>5) & 31)) | (1 << ((syn>>10) & 31))
        # bit 24 in syn is also set if syn has weight 1
        bad &= not bool(syn & (y | 0x1000000)) 
        if bad:
            raise ValueError("Golay code syndrome is not unique")
        syn ^= y
        return syn & 0xffffff


    @classmethod
    def syndrome(cls, v1, u_tetrad=24):
        """Return equivalent Golay cocode word to 'v1' of length at most four.

        If that cocode word has length four, 0 <= 'u_tetrad' < 24 is a bit 
        position where the (not uniquely defined) equivalent cocode word will 
        have a 'one' bit. Otherwise, that cocode word is uniquely defined, and 
        'tetrad' is ignored. in any case 0 <= 'u_tetrad' <= 24 is required. 
        Here 24 means don't care.  Thus u_tetrad = Mat24.lsbit24(x) is legal 
        for any x. 
        """
        v1 = int(v1)
        u_tetrad = int(u_tetrad)
        return cls.cocode_syndrome(cls.vect_to_vintern(v1), u_tetrad)


    @classmethod
    def vect_type(cls, v1):
        """Return type of bit vector v1 in GF(2)**24. 

        This is either (w,d,y), 0 <= d <= 3,  or (w,4,t)

        Here w is the bit weight of vector v1. In the first case, there
        is a unique syndrome s(v1) of v1 of weight d, 0 <= d <= 3, and
        the weight of v & s(v1) is equal to y in case w <= 12 and to
        d - y in case w >= 12.

        In the second case the syndrome may be any of a tetrad of 
        six vectors s_i of length 4 with disjoint support. In case 
        w = 12, the value  t is the number of vectors s_i completely 
        contained in v1. For other values of w the value t is zero.
  
        It is well known that all vectors of the same type are
        in one orbit of the automorphism group of the Golay code.
        """
        v1 = int(v1)
        x = v1
        w = cls.bw24(x)
        s = cls.syndrome(x)
        d = cls.bw24(s)
        if d < 4: 
           y =  bw24(x & s) 
           if w >= 12: y = d - y
        elif  w != 12: 
           y = 0  
        else:
           y, rest = 0 + (s == x & s),  x & ~s & 0xffffff 
           while rest:
               s =  cls.syndrome(x, cls.lsbit24(rest))
               y +=  s == x & s  
               rest &= ~s
        return ( w, d, y )


   
    @classmethod
    def gcode_weight(cls, v1):
        """Return bit weight of Golay code word v1 divided by 4

        Here 0 <= v1 < 4096 is the number of a Golay code word.
        """
        v1 = int(v1)
        t = - ((v1 >> 11) & 1)
        w = (int(cls.theta_table[v1 & 0x7ff]) >> 12) & 7
        return ((w & 7) ^ t) + (t & 7)  

    @classmethod
    def gcode_to_bit_list(cls, v1):
        """Return Golay code word 'v1' as a list of bit positions.

        Here 0 <= v1 < 4096 is the number of a Golay code word.
        """
        v = cls.gcode_to_vect(v1)
        length, bit_list = cls.vect_to_bit_list(v)
        return bit_list[:length]

    @classmethod
    def cocode_weight(cls, c1):
        """Return minimum possible weight of cocode vector c1.
  
        Here c1 is a cocode vector in cocode representation. 
        """
        syn = cls.cocode_syndrome(c1, u_tetrad=0)
        return cls.bw24(syn)

    @classmethod
    def cocode_to_bit_list(cls, c1, u_tetrad=24):
        """Return cocode syndrome of 'c1' as list of bit positions.

        Here c1 is a cocode vector in cocode representation. The 
        function returns the sorted list of the bit positions of 
        the syndrome of c1.  That syndrome is unique if it has 
        length < 4. If the syndrome has length four then we return 
        the following syndrome:
          - A syndrome containing u_tetrad if 0 <= u_tetrad < 24, 
          - We raise ValueError if  'u_tetrad' == 24.  
        """
        v = cls.cocode_syndrome(c1, u_tetrad)
        length, bit_list = cls.vect_to_bit_list(v)
        return bit_list[:length]


    @classmethod
    def cocode_to_sextet(cls, c1):
        """Return a cocode word 'c1' in as a sextet.
 
        Here 'c1' is an cocode word in cocode representation. The 
        function stores concatenation of the 6 tetrads that make up 
        the sextet 'c1' is a list of length 24. The (ordered) 
        tetrads are stored in lexical order. 
        It raises ValueError if 'c1' is not a tetrad.
        """
        c1 = int(c1)
        sextet = cls.cocode_to_bit_list(c1, 0)
        if len(sextet) != 4:
            raise ValueError("Golay cocode word is not a sextet") 
        v = set(range(24)) - set(sextet)
        while len(v):
            tetrad = cls.cocode_to_bit_list(c1, min(v)) 
            v = v - set(tetrad)
            sextet += tetrad
        return sextet
        


    ############################################################################
    # Scalar product of Golay code and cocode
    ############################################################################


    @classmethod
    def scalar_prod(cls, v1, c1):
        """Return scalar product (v1,c1) of Golay code and cocode vector

        v1 is a Golay code vector in 'gcode' representation, c1 is a
        cocode vector in cocode representation.
        """
        v1 = int(v1)
        c1 = int(c1)
        r = v1 & c1
        r ^= r >> 6
        r ^= r >> 3
        return (0x96 >> (r & 7)) & 1


    ############################################################################
    # Conversion from and to suboctads
    ############################################################################

    @classmethod
    def suboctad_to_cocode(cls, u_sub, octad):
        r"""Convert even suboctad u_sub of octad v to cocode representation

        Each octad ``octad`` has 64 even subsets, when a subset x of
        ``octad`` is identified with its complement ``octad \ x``.
        These subsets are called suboctads.Let b_0,...,b_7 be the
        elements of the octad ``octad`` in order specified by function
        ``mat24_cocode_to_suboctad``.
        Then the even subset (b_0 , b_i)  has suboctad number 2**(i-1) for 
        i = 1,...,6. Combining suboctads by disjoint union corresponds to 
        combining their numbers by xor. The empty subocatad has number 0.
        This yields a 1:1 correspondence between the integers 0,...,63 and
        the suboctads of a fixed octad ``octad``, when identifying a
        suboctad with its complement.

        The function returns the suboctad of ``octad`` with number u_sub
        in 'cocode' representation. Octad ``octad`` must be given in
        'gcode' representation. An exception is raised if ``octad`` is
        not the nunber of an octad.
        """
        u_sub = int(u_sub)
        octad = int(octad)
        parity =  (0x96 >> ((u_sub ^ (u_sub >> 3)) & 7)) & 1
        sub = parity + ((u_sub & 0x3f) << 1)
        octad_entries = [int(x) for x in cls.octad_table[8*octad:8*octad+8]]
        vector =  sum(1 << o for i,o in enumerate(octad_entries) if
                      (1 << i) & sub)
        return cls.vect_to_cocode(vector) 
        
    @classmethod
    def cocode_to_suboctad(cls, c1, v1, u_strict = 0):
        """Convert cocode element c1 to suboctad of octad v1.

        This is an inverse of member function suboctad_to_cocode().
        The function returns the suboctad number corresponding to the
        cocode element c1, if v1 is an octad and c1 is an even subset 
        of v1. c1 and v1 must be given in 'cocode' and 'gcode'
        representation, respectively. The function raises an exception
        if v1 is not an octad or c1 cannot be represented as an even
        subset of v1.

        If bit 0 of ``u_strict`` is set then the suboctad must also
        correspond to a short Leech lattice vector; otherwise
        ValueError is raised. 
        """ 
        v1 = int(v1)
        c1 = int(c1)
        octad = cls.gcode_to_octad(v1, u_strict = 0)
        octad_entries = [int(x) for x in cls.octad_table[8*octad:8*octad+8]]
        syn = cls.cocode_syndrome(c1, octad_entries[0])
        v = sum(1 << i for i in octad_entries)
        if c1 & 0x800 or syn & v != syn:
            raise ValueError("Cocode word is not a suboctad")
        suboctad = sum(1 << i for i, o in enumerate(octad_entries) if
                      (1 << int(o)) & syn)
        if u_strict & 1:
            w = cls.bw24(cls.gcode_to_vect(v1))
            if not(((w >> 3) ^ cls.suboctad_weight(suboctad)) &  1):
                print("Octad = %s, suboctad = %s" % (hex(v), hex(suboctad)))
                E = "Leech vector corresponding to suboctad is not short"
                raise ValueError(E) 
        suboctad ^= 0xff if suboctad & 0x80 else 0
        return (octad << 6) + (suboctad >> 1)

    @classmethod
    def octad_entries(self, u_octad):
        assert 0 <= u_octad < 759, "not a Golay code octad"
        return [int(x) for x in octad_table[8*u_octad : 8*u_octad+8]]

    @staticmethod
    def suboctad_weight(u_sub):
        """Return parity of halved bit weight of suboctad u_sub
    
        See member function suboctad_to_cocode for documentation of 
        suboctads. The function returns 0 is the bit weight of a
        suboctad is divisible by four and 1 othewise.
        """
        u_sub = int(u_sub)
        w = bw24(u_sub & 0x3f)
        return ((w + 1) >> 1) & 1
  
    @staticmethod
    def suboctad_scalar_prod(u_sub1, u_sub2):
        """Return scalar product of suboctads u_sub1, u_sub2
    
        Here the scalar product is the parity of the  vector  
        u_sub1 & u_sub2 when u_sub1 and u_sub2 are given as subsets
        of an octad in vector notation.

        But in this functions parameters u_sub1, u_sub2 are suboctad 
        numbers as described in function suboctad_to_cocode.
        """
        u_sub = int(u_sub1)
        u_sub = int(u_sub1)
        wp = (0x96 >> ((u_sub1 ^ (u_sub1 >> 3)) & 7)) 
        wp &=  (0x96 >> ((u_sub2 ^ (u_sub2 >> 3)) & 7))
        u_sub1 &= u_sub2
        wp ^= (0x96 >> ((u_sub1 ^ (u_sub1 >> 3)) & 7)) 
        return wp & 1


