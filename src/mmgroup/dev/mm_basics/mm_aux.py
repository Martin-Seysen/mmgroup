"""Accessing vectors of the representation of the monster

We deal with the rational representation [Seys19]  of the monster 
group MM, which is  based on the Conway's construction [Conw85] of
the monster, modulo various small integers p = 2**n-1, 2 <= n <= 8.
Here the integer p is called the modulus. 

All representations R_p of MM are 196884-dimensional matrix 
representations over the integers modulo p. Here the matrices
operate on vectors which are represented as arrays of integers.
The matrix operations on these vectors are performed by C functions.

For vectors there is an external representation as an array of 
integers of type uint8_t and there is also an internal representation. 
The internal representation is fairly close to the description of the 
basis vectors of the representation in [Seys19] , so that matrix 
operations can take advantage of the internal structure of that 
repesentation.

The purpose of this module is access to the internal representation 
of such a vector and conversion between the external and the 
internal representation.  

There are versions optimized for 32-bit and for 64-bit computers. 
This version is optimized for {INT_BITS} bits. It supports 
moduli p = {P_LIST}.


References
----------
see file references.txt
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
import collections
import warnings
from numbers import Integral
import numpy as np


from mmgroup.dev.mm_basics import mm_basics   


from mmgroup.dev.mm_basics.mm_basics import MM_Basics, MM_Const
from mmgroup.dev.generate_c.make_c_tables import c_snippet
from mmgroup.dev.generate_c.make_c_tables import TableGenerator
from mmgroup.dev.generate_c.make_c_tables import format_item

from mmgroup.bitfunctions import bitlen, hibit


INT_BITS =  MM_Basics.INT_BITS 


internal_rep = """Internal representation of a vector

Accoding to [Seys19], the 196884 basis vectors of the representation 
R_p are labeled as tuples (tag, index1, index) follows:

   Name      Tag   Entries                               Remarks
   (ii)       A    i, i;  0 <= i < 24                    (1,4)
   (ij)       A    i, j;  0 <= j < i < 24                (1)
   X^+(ij)    B    i, j;  0 <= j < i < 24                (1)
   X^-(ij)    C    i, j;  0 <= j < i < 24                (1)
   X_o,s      T    o, s;  0 <= o < 759, 0 <= s < 64      (2)
   X_d,j      X    d, i;  0 <= d < 2**11, 0 <= j < 24    (1,3)
   d^+(*)j    Z    d, i;  0 <= d < 2**11, 0 <= j < 24    (1,3)
   d^-(*)j    Y    d, i;  0 <= d < 2**11, 0 <= j < 24    (1,3)


The tag is used in some python modules (e.g. autogroup.mvector)
to distinguish between different types of basis vectors.

Remarks
(1)  i and j, 0 <= i,j < 24  refer to basis vectors of the
     Boolean vector space in which the Golay code is defined.
(2)  o is one of 759 octads, s is one of 64 even subsets of
     octad d (modulo its complement in d). In the documentation
     mat24_functions_c_doc.txt we refer to the 'octad' 
     representation for o and to function 
     mat24_suboctad_to_cocode() for s.
(3)  d, 0 <= d < 2**11 refers to a Golay code vector modulo the 
     word Omega of weight 24. More specifically, if d0 is a 
     Golay code vector in 'gcode' representation (see
     documentation mat24_functions_c_doc.txt) then we 
     put d = d0 >> 1.
(4)  (ii) is the diagonal entry of a 24 times 24 matrix A

Basis vectors (ij) correspond to the entries of a symmetric
24 x 24 matrix. So it is natural to store the coefficients for (ij) 
in A[i,j] and also in A[j,i], with A an array of 24 small arrays of 
24 entries. Basis vectors (ii) correspond to the diagonal entries of
matrix A, see [Conw85] or [Seys19]  for background. Similarly, the 
coefficients for basis vectors X^+(ij) and X^-(ij) are stored as two 
symmetric 24 x 24 matrices. This module contains functions for 
dealing with that kind of redundancy. 

So it is natural to arrange a vector in R_p as an array of 759 small 
arrays with 64 entries, followed by 3*24 + 3*2048 small arrays with 
24 entries. 

In summary a vector of R_p in internal representation contains 
(in that order):

    24 small arrays of 24 entries 'A' for  (ij)
       with entries (ii) stored in the diagonal positions
    24 small arrays of 24 entries 'B' for  X^+(ij)
    24 small arrays of 24 entries 'C' for  X^-(ij)
   759 small arrays of 64 entries 'T' for  X_o,s
  2048 small arrays of 24 entries 'X' for  X_d,j   
  2048 small arrays of 24 entries 'Z' for  d^+(*)j  
  2048 small arrays of 24 entries 'Y' for  d^-(*)j   

Diagonal positions in the arrays for X^+(ij) and X^-(ij) are unused. 

So an entry of a small array corresponds to an entry of a vector in 
R_p. Such an entry is always coded in a bit field inside an integer 
of type  uint_mmv_t with  uint_mmv_t =  uint{INT_BITS}_t. The number 
of bits in a bit field is always a power of two. So e.g. for p = 3 
we use 2 bits; for p = 7 we use 4 bits with the highest bit unused. 
In case p = 2**k - 1, legal values for an entry are  0,...,2**k - 1, 
with 2**k - 1 equal to 0. Thus negation of a value can be done by 
complementing all k bits of that value. Apart from negation, the 
matrices in R_p perform the operations addition, subtraction and 
halfing on the entries of a vector, which can easily be done on 
several entries simultaneously by manipulating a just single integer 
of type uint_mmv_t.


Actually, we reserve 32 entries for small array with 24 etries. Thus 
a vector in R_p has 759*24 + (3*42 + 3*2048)*32 = 247488 entries. 
Compared to the 196884 entries actually needed for a vector in R_p 
we require about 25.7% more memory than necessary. In some cases we 
need this extra memory anyway. E.g. for p=3 an uint32_t may store 16 
and an uint64_t may store 32 entries, so that there will always be 
a slack of 8 entries for a small array with 24 entries. 

Function mm_aux_mmv_size(p) returns the number of integers of type 
uint_mmv_t required for storing a vector of R_p. 

When reading or calculating entries corresponding to (ij), X^+(ij) 
or X^+(ij), the entry with index i > j should be read or calculated 
first. Procedures for copying all such entries with index i > j to 
their transposed positions are provided.

Direct access to an entry of a vector in R_p is dififcult. Function 
mm_aux_get_mmv(p, mv, i) returns entry i of the vector mv of R_p.
mm_aux_put_mmv(p, value, mv, i) sets that entry with index i to the
given value. Here valid indices are 0 <= i < 196884, with the order
compatible to the order in section 'External representation'.


"""





doc_odd = """The order of the entries d^+(*)i, d^-(*)i and X_d,i  

Let Pl be the Parker loop as described in [Seys19] . Each of the 
arrays for  d^+(*)i, d^-(*)i and X_d,i  is organized as a 2048 x 24 
matrix, with the 2048 rows corresponding to a certain transversal 
of the center Z(Pl) of Pl in Pl.  The 24 columns of these matrices 
correspond to the 24 basis vectors of the space GF(2)^24  in which 
the binary Golay code is defined. 

We define a certain basis of the Golay code C in file 
mat24.functions.c. In that file we also fix a cocycle 
theta :  C  x  C  ->  GF(2), and we represent the elements of Pl 
as pairs (d, l) with d in C, l in GF(2), so that 

      (d1, l1) * (d2, l2) = (d1 + d2, l1 + l2 + theta(d1, d2)). 
 
Elements (d,0) of Pl are called positive. (d, 1) is the negative of 
(d, 0). Pl has a unit element 1 = (0, 0), its negative is -1 = (0,1),
and the quotient of Pl by {{1,-1}} is the binary Golay code C. The 
center Z(Pl) consists of 1, -1 and two more elements Omega, -Omega
represented as ('Omega, 0), ('Omega, 1). Here 'Omega is the all-one 
element of the Golay code C. More generally, for any vd in Pl we
write 'v for the Golay code word corresponding to v.  

A basis of the Golay code is given in file mat24_functions_c_doc.txt.
Basis vectors b[k] of C are numbered from 0 to 11, and basis vector 0 
in the chosen basis is 'Omega. So we obtain a transversal of Pl/Z(Pl) 
by mapping the integer d = sum d[k] * 2**k, k = 0,...,10, d[k] = 0,1,
to the element (d', 0) of the Parker loop Pl, where d' is the Golay 
code word given by  d' = sum d[k] * b[k+1].

This mapping from the interval [0,2048) to Pl yields the order of 
the rows in the arrays for d^+(*)i and d^-(*)j. Note that
e^+(*)j = d^+(*)j  and  e^-(*)j = - d^-(*)j  for  e = Omega*d.

We remark that the standard numbering of Golay code vectors in file
mat24.functions.c maps the integer (sum d[k] * 2**k)  to the Golay
code word (sum d[k] * b[k]). Thus for  X_d,i,  d^+(*)i,  and  d^-(*)i 
the index d corresponds to the Golay code word with number 2*d. 

Basis vector X_d,j is associated with the pair (('d, 0), c_j) of the 
Cartesic product Pl times C*, where c_j is the Golay cocode vector 
represented by the singleton {{j}}. Here the pair ('d, c_j) must 
correspond to a short Leech lattice vector, see [Seys19], section 6.2.
('d, c_j) corresponds to a short vector if and only if 
<'d, j> = |'d|/4 (mod 2), with  |'d| is the Hamming weight of 'd
and <.,.> the scalar product. If this is not the case, we decree
that  X_d,j is associated with the pair (('d + 'Omega, 0), c_j).
In any case, exactly one of the pairs  ('d, c_j) and
('d + 'Omega, c_j) corresponds to a short Leech lattice vector.
"""

#234567890123456789012345678901234567890123456789012345678901234567890



doc_octads = """The order of the entries X_o,s

The array for X_o,s is organized as a 759 x 64 matrix, with the 759 
rows corrsponding to the octads in the Golay code, and the 64 
entries of a row columns corresponding to the 64 suboctads of that
octad. Here an octad is a code word of the Golay code of Hamming 
weight 8 and a suboctad of an octad is an element of the Golay cocode
represented by a subset of the octad with even Hamming weight. A 
subset of an octad represents the same suboctad as its complement 
with respect to that octad.

We refer to the description of the Golay code and the Parker loop 
in the previous section. Row o, 0 <= 0 < 759, of the matrix 
corresponds to the element ('o, 0) of the Parker loop, where o is 
the number of the octad and 'o is the octad itself considered as a 
Golay code word. The numbering of the octads is given by function 
mat24_octad_to_gcode(o) in file mat24_functions.c. 

It remains to specify the order of the suboctads within the row 
corresponding to octad 'o. Notice that Golay code word 'o can be in
terpreted as a subset {{i0,...,i7}} of GF(2)^24. We assume that the
indices i0,...,i7 are in natural order. Any vector in GF(2)^24 may 
be interpreted as subsets of GF(2)^24 with vector addition correspon-
ding to the symmetric difference of subsets. The even subsets of 'o 
(modulo  the whole set 'o) have a basis

     b_0 = {{i0,i1}},  b_1 = {{i0,i2}},  b_2 = {{i0,i3}},   
     b_3 = {{i0,i4}},  b_4 = {{i0,i5}},  b_5 = {{i0,i6}}.   

Then column j with  j = sum j[k] * 2**k, k = 0,...,5, j[k] = 0,1, 
of that row corresponds to the subsoctad (sum j[k] * b_k) of 'o.
 
For the sake of accuracy we remark that X_o,s must be interpreted
as X_e,s, with 'e = 'Omega * 'o, in all cases where |s|/2 is odd. 
"""
#234567890123456789012345678901234567890123456789012345678901234567890

external_rep = """External representation of a vector 

There is also a so-called external representation of a vector in R_p. 
This is used to facilitate the access to vectors by external modules. 
Here the vector is represented as an array of 196884 integers of type
uint8_t. Basis vectors are ordered similar to the ordering for the
internal representation, but here the entries are in 1:1 correspondence 
with the basis vectors. In the external representation there are no 
unused or duplicated entries.

More precisely, the order of the entries is:

         24 entries with tag 'A' for  (ii)
        276 entries with tag 'A' for  (ij), i > j
        276 entries with tag 'B' for  X^+(ij), i > j
        276 entries with tag 'C' for  X^-(ij), i > j
   759 * 64 entries with tag 'T' for  X_d,delta
  2048 * 24 entries with tag 'X' for  X_d,j   
  2048 * 24 entries with tag 'Z'  for  d^+(*)j  
  2048 * 24 entries with tag 'Y'  for  d^-(*)j  

Indices (i,j) for entries (ij), X^+(ij), X^-(ij) are ordered as follows:

     (1,0),
     (2,0), (2,1),
     (3,0), (3,1), (3,2),
     ...
     (i,0), (i,1), ..., (i,i-1), 
     ...
     (24,0), (24,1), ..., (24,23).

Function mm_aux_bytes_to_mmv() converts a vector from external to
internal representation, Function mm_aux_mmv_to_bytes() does the 
inverse conversion.
"""


doc_sparse = """Sparse format of a vector of the 196884-dimensional rep R_p

A vector of the rep R_p can also be stored in sparse format. Here
a vector is stored as an array of 32-bit integers, where each entry 
stands for a multiple of a baiss vector. Basically, component of
a vector is stored in the bit fields of an integer as a pair
(label, coordinate). Here the 'label' is the label of a basis vector 
as given by the internal representation of a vector. Thus the label
is a tuple of the form (tag, index1, index2). 

An entry of teh sparse representation has the following bit fields:

   Bits     Meaning
   27..25   Tag: A = 1, B = 2, C = 3, T = 4, X = 5, Z = 6, Y = 7
   24..15   Index1: i or d with respect to the tag
   13.. 8   Index2: j or delta with respect to the tag
    7.. 0   Coordinate of basis vector; if the moulus p of the 
            representation is k bits long, then only the lowest
            k bits are evaluated.
            
The length of a sparse representation of a  vector must be given as 
a extra parameter to the functions in this module. A  sparse 
representation may contain several entries with he same label; in 
this case the corresponding coordinates are added. The order of
the entries is irrelevant in a sparse representation. Function
in this module generate sparse representations with at most one
entry for each label, but they accept sparse representations 
with several entries for each label.

In an entry with tag 'A', 'B' or 'C' generated by this module we
always have i >= j. An sparse rep generated by this module has
entries with coefficient at least  and at most p - 1. 

When reading an entry, coefficients c with 0 <= c <= p are accepted.
Entries with tag 'A', 'B' or 'C' and i < j are also accepted.
Illegal tags or indices are usually ignored on input.
"""


#2345678901234567890123456789012345678901234567890123456789012345678901

class MM_IO24(MM_Const):
    """Generate code to access vectors of R_p in internal representation

    This class provides code generation functions for fast data
    transfer between a small array in a vector in R_p and an
    array of bytes. 

    We use the code generation mechanism in class  
    make_c_tabels.TableGenerator. The generated functions take
    an input parameter 'log_field_bits' coresponding to
    LOG_FIELD_BITS. This indicates the a bit field containing an
    entry contains 2**log_field_bits.
    """
    #inv_triangle = []
    #for i in range(24): 
    #    inv_triangle += [i] * i
    #assert len(inv_triangle) == 276 

    @staticmethod
    def iter_ext_ABC():
        """Yield triples (tag, i, j) for external rep for tags A, B, C

        Triples are given in the order of the external representation
        for tags A, B and C.
        tag = 0, 1, 2 means 'A', 'B', 'C', respectively
        """
        for i in range(24):
            yield 0, i, i
        for tag in range(3):
            for i in range(24):
                for j in range(i):
                     yield tag, i, j

    @classmethod
    def abc_table(cls):
        """Table for expanding entries for tags 'A', 'B', 'C'.

        Entry k0 of the external representation of the monster is
        mapped to location k1 in the internal representation with
        k1 = (Table[k0] & 0x7ff) + k0 - 24. Entry k0 is also copied 
        to location k1 - 31 * (Table[k0] >> 11) of the internal rep.

        See comments 'Internal representation of a vector' and
        'External representation of a vector' for internal and
        external representation of a vector.
        """
        a = np.zeros(24 + 3 * 276, dtype = np.int16)
        for k0, (tag, i, j) in enumerate(cls.iter_ext_ABC()):
            k1 = (24 * tag + i) * 32 + j
            k1_entry = k1 - k0 + 24
            d_t  = i - j
            assert 0 <= k1_entry <  2048
            a[k0] = k1_entry  + (d_t << 11)
        return a
               


    @classmethod
    def reduce_table(cls):
        """Masks for reducing the fields of an uint_mmv_t modulo p

        Table entries 2*i-2 and 2*i-1 refer to modulus p = 2**i-1 for 
        i = 2,...,8. In the standard case, when i is not a power of 
        two, entries 2*i-2 and  2*i-1 have the following values:

        Index  Table entry
        2*i-2: A mask containg the value 1 in each field. 
        2*i-1: A mask containg the value 2**i-1 in each field. 

        If i is a power of two then table entry 2*i-1 is mask containing 
        the  value 2**(i/2)-1 in each field.

        The reason for that special case is that there is no space for 
        a carry bit between two adjacent fields if i is a power of two.
        In that case we need a special trick for detecting the value 
        2**i-i in a field.        
        """ 
        tbl = []
        for p_bits in range(2, 9):
            field_bits = p_bits
            while (field_bits & (field_bits - 1)): 
                field_bits += 1 
            mask_1 = sum(1 << i for i in range(0, INT_BITS, field_bits))
            q = p_bits >> (field_bits == p_bits)
            mask_q = ((1 << q) - 1) * mask_1
            tbl += [mask_1, mask_q]
        return tbl      
                    
    def __init__(self):
        """For code generation with class make_c_tables.TableGenerator

        """
        super(MM_IO24, self).__init__()
        my_new_tables = {
            "MM_AUX": sys.modules[__name__],
            "MM_AUX_IO24": self,
            "MM_AUX_TBL_ABC" : self.abc_table(), 
            "MM_AUX_TBL_REDUCE": self.reduce_table(),
            #"MM_TBL_INV_TRIANGLE24" : self.inv_triangle,
        }
        self.make_tables(my_new_tables, {})





