/////////////////////////////////////////////////////////////////////////////
// This C file has been created automatically. Do not edit!!!
/////////////////////////////////////////////////////////////////////////////

// %%COMMENT
// %%PY_DOCSTR MM_AUX, 1
// Accessing vectors of the representation of the monster
// 
// We deal with the rational representation [Seys19]  of the monster
// group MM, which is  based on the Conway's construction [Conw85] of
// the monster, modulo various small integers p = 2**n-1, 2 <= n <= 8.
// Here the integer p is called the modulus.
// 
// All representations R_p of MM are 196884-dimensional matrix
// representations over the integers modulo p. Here the matrices
// operate on vectors which are represented as arrays of integers.
// The matrix operations on these vectors are performed by C functions.
// 
// For vectors there is an external representation as an array of
// integers of type uint8_t and there is also an internal representation.
// The internal representation is fairly close to the description of the
// basis vectors of the representation in [Seys19] , so that matrix
// operations can take advantage of the internal structure of that
// repesentation.
// 
// The purpose of this module is access to the internal representation
// of such a vector and conversion between the external and the
// internal representation.
// 
// There are versions optimized for 32-bit and for 64-bit computers.
// This version is optimized for 64 bits. It supports
// moduli p = [3, 7, 15, 127].
// 
// 
// References
// ----------
// see file references.txt

// %%COMMENT
// %%PY_DOCSTR MM_AUX.internal_rep, 1
// Internal representation of a vector
// 
// Accoding to [Seys19], the 196884 basis vectors of the representation
// R_p are labeled as tuples (tag, index1, index) follows:
// 
//    Name      Tag   Entries                               Remarks
//    (ii)       A    i, i;  0 <= i < 24                    (1,4)
//    (ij)       A    i, j;  0 <= j < i < 24                (1)
//    X^+(ij)    B    i, j;  0 <= j < i < 24                (1)
//    X^-(ij)    C    i, j;  0 <= j < i < 24                (1)
//    X_o,s      T    o, s;  0 <= o < 759, 0 <= s < 64      (2)
//    X_d,j      X    d, i;  0 <= d < 2**11, 0 <= j < 24    (1,3)
//    d^+(*)j    Z    d, i;  0 <= d < 2**11, 0 <= j < 24    (1,3)
//    d^-(*)j    Y    d, i;  0 <= d < 2**11, 0 <= j < 24    (1,3)
// 
// 
// The tag is used in some python modules (e.g. autogroup.mvector)
// to distinguish between different types of basis vectors.
// 
// Remarks
// (1)  i and j, 0 <= i,j < 24  refer to basis vectors of the
//      Boolean vector space in which the Golay code is defined.
// (2)  o is one of 759 octads, s is one of 64 even subsets of
//      octad d (modulo its complement in d). In the documentation
//      mat24_functions_c_doc.txt we refer to the 'octad'
//      representation for o and to function
//      mat24_suboctad_to_cocode() for s.
// (3)  d, 0 <= d < 2**11 refers to a Golay code vector modulo the
//      word Omega of weight 24. More specifically, if d0 is a
//      Golay code vector in 'gcode' representation (see
//      documentation mat24_functions_c_doc.txt) then we
//      put d = d0 >> 1.
// (4)  (ii) is the diagonal entry of a 24 times 24 matrix A
// 
// Basis vectors (ij) correspond to the entries of a symmetric
// 24 x 24 matrix. So it is natural to store the coefficients for (ij)
// in A[i,j] and also in A[j,i], with A an array of 24 small arrays of
// 24 entries. Basis vectors (ii) correspond to the diagonal entries of
// matrix A, see [Conw85] or [Seys19]  for background. Similarly, the
// coefficients for basis vectors X^+(ij) and X^-(ij) are stored as two
// symmetric 24 x 24 matrices. This module contains functions for
// dealing with that kind of redundancy.
// 
// So it is natural to arrange a vector in R_p as an array of 759 small
// arrays with 64 entries, followed by 3*24 + 3*2048 small arrays with
// 24 entries.
// 
// In summary a vector of R_p in internal representation contains
// (in that order):
// 
//     24 small arrays of 24 entries 'A' for  (ij)
//        with entries (ii) stored in the diagonal positions
//     24 small arrays of 24 entries 'B' for  X^+(ij)
//     24 small arrays of 24 entries 'C' for  X^-(ij)
//    759 small arrays of 64 entries 'T' for  X_o,s
//   2048 small arrays of 24 entries 'X' for  X_d,j
//   2048 small arrays of 24 entries 'Z' for  d^+(*)j
//   2048 small arrays of 24 entries 'Y' for  d^-(*)j
// 
// Diagonal positions in the arrays for X^+(ij) and X^-(ij) are unused.
// 
// So an entry of a small array corresponds to an entry of a vector in
// R_p. Such an entry is always coded in a bit field inside an integer
// of type  uint_mmv_t with  uint_mmv_t =  uint64_t. The number
// of bits in a bit field is always a power of two. So e.g. for p = 3
// we use 2 bits; for p = 7 we use 4 bits with the highest bit unused.
// In case p = 2**k - 1, legal values for an entry are  0,...,2**k - 1,
// with 2**k - 1 equal to 0. Thus negation of a value can be done by
// complementing all k bits of that value. Apart from negation, the
// matrices in R_p perform the operations addition, subtraction and
// halfing on the entries of a vector, which can easily be done on
// several entries simultaneously by manipulating a just single integer
// of type uint_mmv_t.
// 
// 
// Actually, we reserve 32 entries for small array with 24 etries. Thus
// a vector in R_p has 759*24 + (3*42 + 3*2048)*32 = 247488 entries.
// Compared to the 196884 entries actually needed for a vector in R_p
// we require about 25.7% more memory than necessary. In some cases we
// need this extra memory anyway. E.g. for p=3 an uint32_t may store 16
// and an uint64_t may store 32 entries, so that there will always be
// a slack of 8 entries for a small array with 24 entries.
// 
// Function mm_aux_mmv_size(p) returns the number of integers of type
// uint_mmv_t required for storing a vector of R_p.
// 
// When reading or calculating entries corresponding to (ij), X^+(ij)
// or X^+(ij), the entry with index i > j should be read or calculated
// first. Procedures for copying all such entries with index i > j to
// their transposed positions are provided.
// 
// Direct access to an entry of a vector in R_p is dififcult. Function
// mm_aux_get_mmv(p, mv, i) returns entry i of the vector mv of R_p.
// mm_aux_put_mmv(p, value, mv, i) sets that entry with index i to the
// given value. Here valid indices are 0 <= i < 196884, with the order
// compatible to the order in section 'External representation'.


// %%COMMENT
// %%PY_DOCSTR MM_AUX.doc_odd, 1
// The order of the entries d^+(*)i, d^-(*)i and X_d,i
// 
// Let Pl be the Parker loop as described in [Seys19] . Each of the
// arrays for  d^+(*)i, d^-(*)i and X_d,i  is organized as a 2048 x 24
// matrix, with the 2048 rows corresponding to a certain transversal
// of the center Z(Pl) of Pl in Pl.  The 24 columns of these matrices
// correspond to the 24 basis vectors of the space GF(2)^24  in which
// the binary Golay code is defined.
// 
// We define a certain basis of the Golay code C in file
// mat24.functions.c. In that file we also fix a cocycle
// theta :  C  x  C  ->  GF(2), and we represent the elements of Pl
// as pairs (d, l) with d in C, l in GF(2), so that
// 
//       (d1, l1) * (d2, l2) = (d1 + d2, l1 + l2 + theta(d1, d2)).
// 
// Elements (d,0) of Pl are called positive. (d, 1) is the negative of
// (d, 0). Pl has a unit element 1 = (0, 0), its negative is -1 = (0,1),
// and the quotient of Pl by {1,-1} is the binary Golay code C. The
// center Z(Pl) consists of 1, -1 and two more elements Omega, -Omega
// represented as ('Omega, 0), ('Omega, 1). Here 'Omega is the all-one
// element of the Golay code C. More generally, for any vd in Pl we
// write 'v for the Golay code word corresponding to v.
// 
// A basis of the Golay code is given in file mat24_functions_c_doc.txt.
// Basis vectors b[k] of C are numbered from 0 to 11, and basis vector 0
// in the chosen basis is 'Omega. So we obtain a transversal of Pl/Z(Pl)
// by mapping the integer d = sum d[k] * 2**k, k = 0,...,10, d[k] = 0,1,
// to the element (d', 0) of the Parker loop Pl, where d' is the Golay
// code word given by  d' = sum d[k] * b[k+1].
// 
// This mapping from the interval [0,2048) to Pl yields the order of
// the rows in the arrays for d^+(*)i and d^-(*)j. Note that
// e^+(*)j = d^+(*)j  and  e^-(*)j = - d^-(*)j  for  e = Omega*d.
// 
// We remark that the standard numbering of Golay code vectors in file
// mat24.functions.c maps the integer (sum d[k] * 2**k)  to the Golay
// code word (sum d[k] * b[k]). Thus for  X_d,i,  d^+(*)i,  and  d^-(*)i
// the index d corresponds to the Golay code word with number 2*d.
// 
// Basis vector X_d,j is associated with the pair (('d, 0), c_j) of the
// Cartesic product Pl times C*, where c_j is the Golay cocode vector
// represented by the singleton {j}. Here the pair ('d, c_j) must
// correspond to a short Leech lattice vector, see [Seys19], section 6.2.
// ('d, c_j) corresponds to a short vector if and only if
// <'d, j> = |'d|/4 (mod 2), with  |'d| is the Hamming weight of 'd
// and <.,.> the scalar product. If this is not the case, we decree
// that  X_d,j is associated with the pair (('d + 'Omega, 0), c_j).
// In any case, exactly one of the pairs  ('d, c_j) and
// ('d + 'Omega, c_j) corresponds to a short Leech lattice vector.


// %%COMMENT
// %%PY_DOCSTR MM_AUX.doc_octads, 1
// The order of the entries X_o,s
// 
// The array for X_o,s is organized as a 759 x 64 matrix, with the 759
// rows corrsponding to the octads in the Golay code, and the 64
// entries of a row columns corresponding to the 64 suboctads of that
// octad. Here an octad is a code word of the Golay code of Hamming
// weight 8 and a suboctad of an octad is an element of the Golay cocode
// represented by a subset of the octad with even Hamming weight. A
// subset of an octad represents the same suboctad as its complement
// with respect to that octad.
// 
// We refer to the description of the Golay code and the Parker loop
// in the previous section. Row o, 0 <= 0 < 759, of the matrix
// corresponds to the element ('o, 0) of the Parker loop, where o is
// the number of the octad and 'o is the octad itself considered as a
// Golay code word. The numbering of the octads is given by function
// mat24_octad_to_gcode(o) in file mat24_functions.c.
// 
// It remains to specify the order of the suboctads within the row
// corresponding to octad 'o. Notice that Golay code word 'o can be in
// terpreted as a subset {i0,...,i7} of GF(2)^24. We assume that the
// indices i0,...,i7 are in natural order. Any vector in GF(2)^24 may
// be interpreted as subsets of GF(2)^24 with vector addition correspon-
// ding to the symmetric difference of subsets. The even subsets of 'o
// (modulo  the whole set 'o) have a basis
// 
//      b_0 = {i0,i1},  b_1 = {i0,i2},  b_2 = {i0,i3},
//      b_3 = {i0,i4},  b_4 = {i0,i5},  b_5 = {i0,i6}.
// 
// Then column j with  j = sum j[k] * 2**k, k = 0,...,5, j[k] = 0,1,
// of that row corresponds to the subsoctad (sum j[k] * b_k) of 'o.
// 
// For the sake of accuracy we remark that X_o,s must be interpreted
// as X_e,s, with 'e = 'Omega * 'o, in all cases where |s|/2 is odd.


// %%COMMENT
// %%PY_DOCSTR MM_AUX.external_rep, 1
// External representation of a vector
// 
// There is also a so-called external representation of a vector in R_p.
// This is used to facilitate the access to vectors by external modules.
// Here the vector is represented as an array of 196884 integers of type
// uint8_t. Basis vectors are ordered similar to the ordering for the
// internal representation, but here the entries are in 1:1 correspondence
// with the basis vectors. In the external representation there are no
// unused or duplicated entries.
// 
// More precisely, the order of the entries is:
// 
//          24 entries with tag 'A' for  (ii)
//         276 entries with tag 'A' for  (ij), i > j
//         276 entries with tag 'B' for  X^+(ij), i > j
//         276 entries with tag 'C' for  X^-(ij), i > j
//    759 * 64 entries with tag 'T' for  X_d,delta
//   2048 * 24 entries with tag 'X' for  X_d,j
//   2048 * 24 entries with tag 'Z'  for  d^+(*)j
//   2048 * 24 entries with tag 'Y'  for  d^-(*)j
// 
// Indices (i,j) for entries (ij), X^+(ij), X^-(ij) are ordered as follows:
// 
//      (1,0),
//      (2,0), (2,1),
//      (3,0), (3,1), (3,2),
//      ...
//      (i,0), (i,1), ..., (i,i-1),
//      ...
//      (24,0), (24,1), ..., (24,23).
// 
// Function mm_aux_bytes_to_mmv() converts a vector from external to
// internal representation, Function mm_aux_mmv_to_bytes() does the
// inverse conversion.

// %%COMMENT
// %%PY_DOCSTR MM_AUX.doc_sparse, 1
// Sparse format of a vector of the 196884-dimensional rep R_p
// 
// A vector of the rep R_p can also be stored in sparse format. Here
// a vector is stored as an array of 32-bit integers, where each entry
// stands for a multiple of a baiss vector. Basically, component of
// a vector is stored in the bit fields of an integer as a pair
// (label, coordinate). Here the 'label' is the label of a basis vector
// as given by the internal representation of a vector. Thus the label
// is a tuple of the form (tag, index1, index2).
// 
// An entry of teh sparse representation has the following bit fields:
// 
//    Bits     Meaning
//    27..25   Tag: A = 1, B = 2, C = 3, T = 4, X = 5, Z = 6, Y = 7
//    24..15   Index1: i or d with respect to the tag
//    13.. 8   Index2: j or delta with respect to the tag
//     7.. 0   Coordinate of basis vector; if the moulus p of the
//             representation is k bits long, then only the lowest
//             k bits are evaluated.
// 
// The length of a sparse representation of a  vector must be given as
// a extra parameter to the functions in this module. A  sparse
// representation may contain several entries with he same label; in
// this case the corresponding coordinates are added. The order of
// the entries is irrelevant in a sparse representation. Function
// in this module generate sparse representations with at most one
// entry for each label, but they accept sparse representations
// with several entries for each label.
// 
// In an entry with tag 'A', 'B' or 'C' generated by this module we
// always have i >= j. An sparse rep generated by this module has
// entries with coefficient at least  and at most p - 1.
// 
// When reading an entry, coefficients c with 0 <= c <= p are accepted.
// Entries with tag 'A', 'B' or 'C' and i < j are also accepted.
// Illegal tags or indices are usually ignored on input.

#include <stdlib.h>
#include "mm_basics.h"



// %%USE_TABLE
static const uint32_t MMV_CONST_TABLE[] = {
// %%TABLE MMV_CONST_TAB, uint32
0x00044643UL,0x00000000UL,0x00034643UL,0x00011305UL,
0x0003c643UL,0x0002c643UL,0x00022484UL,0x0001a484UL
};


// %%COMMENT
// %%PY_DOCSTR MM_AUX_IO24.abc_table, 1
// Table for expanding entries for tags 'A', 'B', 'C'.
// 
// Entry k0 of the external representation of the monster is
// mapped to location k1 in the internal representation with
// k1 = (Table[k0] & 0x7ff) + k0 - 24. Entry k0 is also copied
// to location k1 - 31 * (Table[k0] >> 11) of the internal rep.
// 
// See comments 'Internal representation of a vector' and
// 'External representation of a vector' for internal and
// external representation of a vector.
static const uint16_t MM_AUX_TBL_ABC[] = {
// %%TABLE MM_AUX_TBL_ABC, uint16
0x0018,0x0038,0x0058,0x0078,0x0098,0x00b8,0x00d8,0x00f8,
0x0118,0x0138,0x0158,0x0178,0x0198,0x01b8,0x01d8,0x01f8,
0x0218,0x0238,0x0258,0x0278,0x0298,0x02b8,0x02d8,0x02f8,
0x0820,0x103f,0x083f,0x185d,0x105d,0x085d,0x207a,0x187a,
0x107a,0x087a,0x2896,0x2096,0x1896,0x1096,0x0896,0x30b1,
0x28b1,0x20b1,0x18b1,0x10b1,0x08b1,0x38cb,0x30cb,0x28cb,
0x20cb,0x18cb,0x10cb,0x08cb,0x40e4,0x38e4,0x30e4,0x28e4,
0x20e4,0x18e4,0x10e4,0x08e4,0x48fc,0x40fc,0x38fc,0x30fc,
0x28fc,0x20fc,0x18fc,0x10fc,0x08fc,0x5113,0x4913,0x4113,
0x3913,0x3113,0x2913,0x2113,0x1913,0x1113,0x0913,0x5929,
0x5129,0x4929,0x4129,0x3929,0x3129,0x2929,0x2129,0x1929,
0x1129,0x0929,0x613e,0x593e,0x513e,0x493e,0x413e,0x393e,
0x313e,0x293e,0x213e,0x193e,0x113e,0x093e,0x6952,0x6152,
0x5952,0x5152,0x4952,0x4152,0x3952,0x3152,0x2952,0x2152,
0x1952,0x1152,0x0952,0x7165,0x6965,0x6165,0x5965,0x5165,
0x4965,0x4165,0x3965,0x3165,0x2965,0x2165,0x1965,0x1165,
0x0965,0x7977,0x7177,0x6977,0x6177,0x5977,0x5177,0x4977,
0x4177,0x3977,0x3177,0x2977,0x2177,0x1977,0x1177,0x0977,
0x8188,0x7988,0x7188,0x6988,0x6188,0x5988,0x5188,0x4988,
0x4188,0x3988,0x3188,0x2988,0x2188,0x1988,0x1188,0x0988,
0x8998,0x8198,0x7998,0x7198,0x6998,0x6198,0x5998,0x5198,
0x4998,0x4198,0x3998,0x3198,0x2998,0x2198,0x1998,0x1198,
0x0998,0x91a7,0x89a7,0x81a7,0x79a7,0x71a7,0x69a7,0x61a7,
0x59a7,0x51a7,0x49a7,0x41a7,0x39a7,0x31a7,0x29a7,0x21a7,
0x19a7,0x11a7,0x09a7,0x99b5,0x91b5,0x89b5,0x81b5,0x79b5,
0x71b5,0x69b5,0x61b5,0x59b5,0x51b5,0x49b5,0x41b5,0x39b5,
0x31b5,0x29b5,0x21b5,0x19b5,0x11b5,0x09b5,0xa1c2,0x99c2,
0x91c2,0x89c2,0x81c2,0x79c2,0x71c2,0x69c2,0x61c2,0x59c2,
0x51c2,0x49c2,0x41c2,0x39c2,0x31c2,0x29c2,0x21c2,0x19c2,
0x11c2,0x09c2,0xa9ce,0xa1ce,0x99ce,0x91ce,0x89ce,0x81ce,
0x79ce,0x71ce,0x69ce,0x61ce,0x59ce,0x51ce,0x49ce,0x41ce,
0x39ce,0x31ce,0x29ce,0x21ce,0x19ce,0x11ce,0x09ce,0xb1d9,
0xa9d9,0xa1d9,0x99d9,0x91d9,0x89d9,0x81d9,0x79d9,0x71d9,
0x69d9,0x61d9,0x59d9,0x51d9,0x49d9,0x41d9,0x39d9,0x31d9,
0x29d9,0x21d9,0x19d9,0x11d9,0x09d9,0xb9e3,0xb1e3,0xa9e3,
0xa1e3,0x99e3,0x91e3,0x89e3,0x81e3,0x79e3,0x71e3,0x69e3,
0x61e3,0x59e3,0x51e3,0x49e3,0x41e3,0x39e3,0x31e3,0x29e3,
0x21e3,0x19e3,0x11e3,0x09e3,0x0a0c,0x122b,0x0a2b,0x1a49,
0x1249,0x0a49,0x2266,0x1a66,0x1266,0x0a66,0x2a82,0x2282,
0x1a82,0x1282,0x0a82,0x329d,0x2a9d,0x229d,0x1a9d,0x129d,
0x0a9d,0x3ab7,0x32b7,0x2ab7,0x22b7,0x1ab7,0x12b7,0x0ab7,
0x42d0,0x3ad0,0x32d0,0x2ad0,0x22d0,0x1ad0,0x12d0,0x0ad0,
0x4ae8,0x42e8,0x3ae8,0x32e8,0x2ae8,0x22e8,0x1ae8,0x12e8,
0x0ae8,0x52ff,0x4aff,0x42ff,0x3aff,0x32ff,0x2aff,0x22ff,
0x1aff,0x12ff,0x0aff,0x5b15,0x5315,0x4b15,0x4315,0x3b15,
0x3315,0x2b15,0x2315,0x1b15,0x1315,0x0b15,0x632a,0x5b2a,
0x532a,0x4b2a,0x432a,0x3b2a,0x332a,0x2b2a,0x232a,0x1b2a,
0x132a,0x0b2a,0x6b3e,0x633e,0x5b3e,0x533e,0x4b3e,0x433e,
0x3b3e,0x333e,0x2b3e,0x233e,0x1b3e,0x133e,0x0b3e,0x7351,
0x6b51,0x6351,0x5b51,0x5351,0x4b51,0x4351,0x3b51,0x3351,
0x2b51,0x2351,0x1b51,0x1351,0x0b51,0x7b63,0x7363,0x6b63,
0x6363,0x5b63,0x5363,0x4b63,0x4363,0x3b63,0x3363,0x2b63,
0x2363,0x1b63,0x1363,0x0b63,0x8374,0x7b74,0x7374,0x6b74,
0x6374,0x5b74,0x5374,0x4b74,0x4374,0x3b74,0x3374,0x2b74,
0x2374,0x1b74,0x1374,0x0b74,0x8b84,0x8384,0x7b84,0x7384,
0x6b84,0x6384,0x5b84,0x5384,0x4b84,0x4384,0x3b84,0x3384,
0x2b84,0x2384,0x1b84,0x1384,0x0b84,0x9393,0x8b93,0x8393,
0x7b93,0x7393,0x6b93,0x6393,0x5b93,0x5393,0x4b93,0x4393,
0x3b93,0x3393,0x2b93,0x2393,0x1b93,0x1393,0x0b93,0x9ba1,
0x93a1,0x8ba1,0x83a1,0x7ba1,0x73a1,0x6ba1,0x63a1,0x5ba1,
0x53a1,0x4ba1,0x43a1,0x3ba1,0x33a1,0x2ba1,0x23a1,0x1ba1,
0x13a1,0x0ba1,0xa3ae,0x9bae,0x93ae,0x8bae,0x83ae,0x7bae,
0x73ae,0x6bae,0x63ae,0x5bae,0x53ae,0x4bae,0x43ae,0x3bae,
0x33ae,0x2bae,0x23ae,0x1bae,0x13ae,0x0bae,0xabba,0xa3ba,
0x9bba,0x93ba,0x8bba,0x83ba,0x7bba,0x73ba,0x6bba,0x63ba,
0x5bba,0x53ba,0x4bba,0x43ba,0x3bba,0x33ba,0x2bba,0x23ba,
0x1bba,0x13ba,0x0bba,0xb3c5,0xabc5,0xa3c5,0x9bc5,0x93c5,
0x8bc5,0x83c5,0x7bc5,0x73c5,0x6bc5,0x63c5,0x5bc5,0x53c5,
0x4bc5,0x43c5,0x3bc5,0x33c5,0x2bc5,0x23c5,0x1bc5,0x13c5,
0x0bc5,0xbbcf,0xb3cf,0xabcf,0xa3cf,0x9bcf,0x93cf,0x8bcf,
0x83cf,0x7bcf,0x73cf,0x6bcf,0x63cf,0x5bcf,0x53cf,0x4bcf,
0x43cf,0x3bcf,0x33cf,0x2bcf,0x23cf,0x1bcf,0x13cf,0x0bcf,
0x0bf8,0x1417,0x0c17,0x1c35,0x1435,0x0c35,0x2452,0x1c52,
0x1452,0x0c52,0x2c6e,0x246e,0x1c6e,0x146e,0x0c6e,0x3489,
0x2c89,0x2489,0x1c89,0x1489,0x0c89,0x3ca3,0x34a3,0x2ca3,
0x24a3,0x1ca3,0x14a3,0x0ca3,0x44bc,0x3cbc,0x34bc,0x2cbc,
0x24bc,0x1cbc,0x14bc,0x0cbc,0x4cd4,0x44d4,0x3cd4,0x34d4,
0x2cd4,0x24d4,0x1cd4,0x14d4,0x0cd4,0x54eb,0x4ceb,0x44eb,
0x3ceb,0x34eb,0x2ceb,0x24eb,0x1ceb,0x14eb,0x0ceb,0x5d01,
0x5501,0x4d01,0x4501,0x3d01,0x3501,0x2d01,0x2501,0x1d01,
0x1501,0x0d01,0x6516,0x5d16,0x5516,0x4d16,0x4516,0x3d16,
0x3516,0x2d16,0x2516,0x1d16,0x1516,0x0d16,0x6d2a,0x652a,
0x5d2a,0x552a,0x4d2a,0x452a,0x3d2a,0x352a,0x2d2a,0x252a,
0x1d2a,0x152a,0x0d2a,0x753d,0x6d3d,0x653d,0x5d3d,0x553d,
0x4d3d,0x453d,0x3d3d,0x353d,0x2d3d,0x253d,0x1d3d,0x153d,
0x0d3d,0x7d4f,0x754f,0x6d4f,0x654f,0x5d4f,0x554f,0x4d4f,
0x454f,0x3d4f,0x354f,0x2d4f,0x254f,0x1d4f,0x154f,0x0d4f,
0x8560,0x7d60,0x7560,0x6d60,0x6560,0x5d60,0x5560,0x4d60,
0x4560,0x3d60,0x3560,0x2d60,0x2560,0x1d60,0x1560,0x0d60,
0x8d70,0x8570,0x7d70,0x7570,0x6d70,0x6570,0x5d70,0x5570,
0x4d70,0x4570,0x3d70,0x3570,0x2d70,0x2570,0x1d70,0x1570,
0x0d70,0x957f,0x8d7f,0x857f,0x7d7f,0x757f,0x6d7f,0x657f,
0x5d7f,0x557f,0x4d7f,0x457f,0x3d7f,0x357f,0x2d7f,0x257f,
0x1d7f,0x157f,0x0d7f,0x9d8d,0x958d,0x8d8d,0x858d,0x7d8d,
0x758d,0x6d8d,0x658d,0x5d8d,0x558d,0x4d8d,0x458d,0x3d8d,
0x358d,0x2d8d,0x258d,0x1d8d,0x158d,0x0d8d,0xa59a,0x9d9a,
0x959a,0x8d9a,0x859a,0x7d9a,0x759a,0x6d9a,0x659a,0x5d9a,
0x559a,0x4d9a,0x459a,0x3d9a,0x359a,0x2d9a,0x259a,0x1d9a,
0x159a,0x0d9a,0xada6,0xa5a6,0x9da6,0x95a6,0x8da6,0x85a6,
0x7da6,0x75a6,0x6da6,0x65a6,0x5da6,0x55a6,0x4da6,0x45a6,
0x3da6,0x35a6,0x2da6,0x25a6,0x1da6,0x15a6,0x0da6,0xb5b1,
0xadb1,0xa5b1,0x9db1,0x95b1,0x8db1,0x85b1,0x7db1,0x75b1,
0x6db1,0x65b1,0x5db1,0x55b1,0x4db1,0x45b1,0x3db1,0x35b1,
0x2db1,0x25b1,0x1db1,0x15b1,0x0db1,0xbdbb,0xb5bb,0xadbb,
0xa5bb,0x9dbb,0x95bb,0x8dbb,0x85bb,0x7dbb,0x75bb,0x6dbb,
0x65bb,0x5dbb,0x55bb,0x4dbb,0x45bb,0x3dbb,0x35bb,0x2dbb,
0x25bb,0x1dbb,0x15bb,0x0dbb
};



// %%COMMENT
// %%PY_DOCSTR MM_AUX_IO24.reduce_table, 1
// Masks for reducing the fields of an uint_mmv_t modulo p
// 
// Table entries 2*i-2 and 2*i-1 refer to modulus p = 2**i-1 for
// i = 2,...,8. In the standard case, when i is not a power of
// two, entries 2*i-2 and  2*i-1 have the following values:
// 
// Index  Table entry
// 2*i-2: A mask containg the value 1 in each field.
// 2*i-1: A mask containg the value 2**i-1 in each field.
// 
// If i is a power of two then table entry 2*i-1 is mask containing
// the  value 2**(i/2)-1 in each field.
// 
// The reason for that special case is that there is no space for
// a carry bit between two adjacent fields if i is a power of two.
// In that case we need a special trick for detecting the value
// 2**i-i in a field.
static const uint_mmv_t MM_AUX_TBL_REDUCE[] = {
// %%TABLE MM_AUX_TBL_REDUCE, uint{INT_BITS}
0x5555555555555555ULL,0x5555555555555555ULL,
0x1111111111111111ULL,0x7777777777777777ULL,
0x1111111111111111ULL,0x3333333333333333ULL,
0x0101010101010101ULL,0x1f1f1f1f1f1f1f1fULL,
0x0101010101010101ULL,0x3f3f3f3f3f3f3f3fULL,
0x0101010101010101ULL,0x7f7f7f7f7f7f7f7fULL,
0x0101010101010101ULL,0x0f0f0f0f0f0f0f0fULL
};




// %%COMMENT
// The order of the parameters of functions in this file is:
//
// 1. Modulus p, if present
// 2. The input data array
// 3. Any parameters that do not affect the positions in the output array
// 4. The output data array
// 5. Parameters (e.g. lengths, indices) that affect the positions of the 
//    data being modified in the output array


/**********************************************************************
*** Some macros
**********************************************************************/

// %%GEN h
// // Return nonzero value if p is a bad modulus,  
// // i.e. not p = 2**k - 1 for some 2 <= k <= 8
// #define mm_aux_bad_p(p) (((p) & ((p)+1)) | (((p)-3) & ((0UL-256UL))))
// 
// // Offsets for tags A,B,C,T,X,Z,Y in the internal representation
// #define MM_AUX_OFS_A       0UL
// #define MM_AUX_OFS_B     768UL    //    24*32
// #define MM_AUX_OFS_C    1536UL    //  2*24*32
// #define MM_AUX_OFS_T    2304UL    //  3*24*32
// #define MM_AUX_OFS_X   50880UL    //  MM_AUX_OFS_T +    759*64
// #define MM_AUX_OFS_Z  116416UL    //  MM_AUX_OFS_X +   2048*32
// #define MM_AUX_OFS_Y  181952UL    //  MM_AUX_OFS_X + 2*2048*32
// #define MM_AUX_OFS_E  247488UL    //  MM_AUX_OFS_X + 3*2048*32. i.e
//                                   //  total length of internal rep
// 
// // Offsets for tags A,B,C,T,X,Z,Y in the external representation
// #define MM_AUX_XOFS_A      24UL
// #define MM_AUX_XOFS_B     300UL    //  24 + 1*276
// #define MM_AUX_XOFS_C     576UL    //  24 + 2*276
// #define MM_AUX_XOFS_T     852UL    //  24 + 3*276
// #define MM_AUX_XOFS_X   49428UL    //  MM_AUX_XOFS_T +    759*64
// #define MM_AUX_XOFS_Z   98580UL    //  MM_AUX_XOFS_X +   2048*24
// #define MM_AUX_XOFS_Y  147732UL    //  MM_AUX_XOFS_X + 2*2048*24
// #define MM_AUX_XOFS_E  196884UL    //  MM_AUX_XOFS_X + 3*2048*24. i.e
//                                    //  total length of external rep
// 
// 
// // Tags for labels and values of vectors in the representation space
// // A multiple of a unit vector with coordinate 'coord' is encoded
// // in the bit fields of a 32-bit integers in the form. 
// //   coord (tag, par1, par2) 
// #define MM_SPACE_TAG_A      0x2000000
// #define MM_SPACE_TAG_B      0x4000000
// #define MM_SPACE_TAG_C      0x6000000
// #define MM_SPACE_TAG_T      0x8000000
// #define MM_SPACE_TAG_X      0xA000000
// #define MM_SPACE_TAG_Z      0xC000000
// #define MM_SPACE_TAG_Y      0xE000000 
// // Mask for all tags:
// // Use y = (x & MM_SPACE_MASK_PAR1) << MM_SPACE_SHIFT_PAR1
// // to set parameter par1 in y to the value x.
// #define MM_SPACE_MASK_TAG     0xE000000 
// // Mask and shift factor for parameter par1  
// // Use y = (x << MM_SPACE_SHIFT_PAR1) & MM_SPACE_MASK_PAR1
// // to set parameter par1 in y to the value x.
// #define MM_SPACE_MASK_PAR1    0x1FFC000   
// #define MM_SPACE_SHIFT_PAR1          14   
// // Mask and shift factor for parameter par12  
// // Use y = (x << MM_SPACE_SHIFT_PAR2) & MM_SPACE_MASK_PAR2
// // to set parameter par2 in y to the value x.
// #define MM_SPACE_MASK_PAR2       0x3F00   
// #define MM_SPACE_SHIFT_PAR2           8 
// // Mask for coordinate:  
// // Use y = x  & MM_SPACE_MASK_COORD
// // to set the coordiante in y to the value x.
// // Caution: some special routines for modulus p = 2**k - 1
// // use only th lowest k bits of the coordinate.
// #define MM_SPACE_COORD_PAR1    0x1FFC000   
// 
// %%GEN c


// %%COMMENT
/**********************************************************************
*** Low-level functions supporting vectors of type uint_mmv_t[]
**********************************************************************/




// %%EXPORT p
void mm_aux_read_mmv1(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len)
// Read entries of vector mv with modulus p and store these entries
// in the array b. len is the number of entries to be read.
// len must be a multiple of the number of entries in an integer
// of type uint_mmv_t. It is ok if len is a multiple of 32. 
// Vector b is reduced modulo p.
{
    uint_fast32_t i, sh, tmp;  
    uint_mmv_t source;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    sh = ((i >> 15) & 15);        // This is P_BITS
    i = ((i >> 9) & 3); // This is LOG_FIELD_BITS
    len >>= 6 - i;

    switch (i) {
        // %%FOR LOG_F in [1, 2, 3]
        case 1:
            while (len--) {
                source = *mv++;
                // %%FOR jj in range(0, INT_BITS, 1 << LOG_F)
                tmp = (source >> 0) & p;
                b[0] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 2) & p;
                b[1] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 4) & p;
                b[2] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 6) & p;
                b[3] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[4] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 10) & p;
                b[5] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 12) & p;
                b[6] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 14) & p;
                b[7] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[8] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 18) & p;
                b[9] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 20) & p;
                b[10] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 22) & p;
                b[11] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[12] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 26) & p;
                b[13] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 28) & p;
                b[14] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 30) & p;
                b[15] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 32) & p;
                b[16] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 34) & p;
                b[17] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 36) & p;
                b[18] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 38) & p;
                b[19] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 40) & p;
                b[20] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 42) & p;
                b[21] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 44) & p;
                b[22] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 46) & p;
                b[23] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 48) & p;
                b[24] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 50) & p;
                b[25] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 52) & p;
                b[26] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 54) & p;
                b[27] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 56) & p;
                b[28] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 58) & p;
                b[29] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 60) & p;
                b[30] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 62) & p;
                b[31] = (tmp + ((tmp + 1) >> sh)) & p;
                // %%END FOR
                b += 32;
            }
            break;
        case 2:
            while (len--) {
                source = *mv++;
                // %%FOR jj in range(0, INT_BITS, 1 << LOG_F)
                tmp = (source >> 0) & p;
                b[0] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 4) & p;
                b[1] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[2] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 12) & p;
                b[3] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[4] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 20) & p;
                b[5] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[6] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 28) & p;
                b[7] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 32) & p;
                b[8] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 36) & p;
                b[9] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 40) & p;
                b[10] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 44) & p;
                b[11] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 48) & p;
                b[12] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 52) & p;
                b[13] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 56) & p;
                b[14] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 60) & p;
                b[15] = (tmp + ((tmp + 1) >> sh)) & p;
                // %%END FOR
                b += 16;
            }
            break;
        case 3:
            while (len--) {
                source = *mv++;
                // %%FOR jj in range(0, INT_BITS, 1 << LOG_F)
                tmp = (source >> 0) & p;
                b[0] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[1] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[2] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[3] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 32) & p;
                b[4] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 40) & p;
                b[5] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 48) & p;
                b[6] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 56) & p;
                b[7] = (tmp + ((tmp + 1) >> sh)) & p;
                // %%END FOR
                b += 8;
            }
            break;
        // %%END FOR
    }
}

// %%EXPORT p
void mm_aux_read_direct_mmv1(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len)
// Same operation as mm_aux_read_mmv1(), but vector b is not reduced.
// This is for debugging and less optimized than mm_aux_read_mmv1(). 
{
    uint_fast32_t i, j;  
    uint_mmv_t source;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    i = ((i >> 9) & 3); // This is LOG_FIELD_BITS
    len >>= 6 - i;
    i = 1 << i;                       // This is FIELD_BITS

    while (len--) {
        source = *mv++;
        for (j = 0; j < 64; j += i) *b++ = (source >> j) & p;
    }
}



// %%EXPORT p
void mm_aux_write_mmv1(uint32_t p, uint8_t *b, uint_mmv_t *mv, uint32_t len)
// Write data from the array b to the vector mv with modulus p. 
// len is the number of entries to be written.
// len must be a multiple of the number of entries in an integer
// of type uint_mmv_t. It is ok if len is a multiple of 32. 
{
    uint_fast32_t i;  
    uint_mmv_t dest;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    i = ((i >> 9) & 3); // This is LOG_FIELD_BITS
    len >>= 6 - i;
    
    switch(i) {
        // %%FOR LOG_F in [1, 2, 3]
        case 1:
            while (len--) {
                dest =  0;
                // %%FOR jj in range(0, INT_BITS, 1 << LOG_F)
                dest +=  (uint_mmv_t)(b[0]) << 0;
                dest +=  (uint_mmv_t)(b[1]) << 2;
                dest +=  (uint_mmv_t)(b[2]) << 4;
                dest +=  (uint_mmv_t)(b[3]) << 6;
                dest +=  (uint_mmv_t)(b[4]) << 8;
                dest +=  (uint_mmv_t)(b[5]) << 10;
                dest +=  (uint_mmv_t)(b[6]) << 12;
                dest +=  (uint_mmv_t)(b[7]) << 14;
                dest +=  (uint_mmv_t)(b[8]) << 16;
                dest +=  (uint_mmv_t)(b[9]) << 18;
                dest +=  (uint_mmv_t)(b[10]) << 20;
                dest +=  (uint_mmv_t)(b[11]) << 22;
                dest +=  (uint_mmv_t)(b[12]) << 24;
                dest +=  (uint_mmv_t)(b[13]) << 26;
                dest +=  (uint_mmv_t)(b[14]) << 28;
                dest +=  (uint_mmv_t)(b[15]) << 30;
                dest +=  (uint_mmv_t)(b[16]) << 32;
                dest +=  (uint_mmv_t)(b[17]) << 34;
                dest +=  (uint_mmv_t)(b[18]) << 36;
                dest +=  (uint_mmv_t)(b[19]) << 38;
                dest +=  (uint_mmv_t)(b[20]) << 40;
                dest +=  (uint_mmv_t)(b[21]) << 42;
                dest +=  (uint_mmv_t)(b[22]) << 44;
                dest +=  (uint_mmv_t)(b[23]) << 46;
                dest +=  (uint_mmv_t)(b[24]) << 48;
                dest +=  (uint_mmv_t)(b[25]) << 50;
                dest +=  (uint_mmv_t)(b[26]) << 52;
                dest +=  (uint_mmv_t)(b[27]) << 54;
                dest +=  (uint_mmv_t)(b[28]) << 56;
                dest +=  (uint_mmv_t)(b[29]) << 58;
                dest +=  (uint_mmv_t)(b[30]) << 60;
                dest +=  (uint_mmv_t)(b[31]) << 62;
                // %%END FOR
                *mv++ = dest;
                b += 32;
            }
            break;
        case 2:
            while (len--) {
                dest =  0;
                // %%FOR jj in range(0, INT_BITS, 1 << LOG_F)
                dest +=  (uint_mmv_t)(b[0]) << 0;
                dest +=  (uint_mmv_t)(b[1]) << 4;
                dest +=  (uint_mmv_t)(b[2]) << 8;
                dest +=  (uint_mmv_t)(b[3]) << 12;
                dest +=  (uint_mmv_t)(b[4]) << 16;
                dest +=  (uint_mmv_t)(b[5]) << 20;
                dest +=  (uint_mmv_t)(b[6]) << 24;
                dest +=  (uint_mmv_t)(b[7]) << 28;
                dest +=  (uint_mmv_t)(b[8]) << 32;
                dest +=  (uint_mmv_t)(b[9]) << 36;
                dest +=  (uint_mmv_t)(b[10]) << 40;
                dest +=  (uint_mmv_t)(b[11]) << 44;
                dest +=  (uint_mmv_t)(b[12]) << 48;
                dest +=  (uint_mmv_t)(b[13]) << 52;
                dest +=  (uint_mmv_t)(b[14]) << 56;
                dest +=  (uint_mmv_t)(b[15]) << 60;
                // %%END FOR
                *mv++ = dest;
                b += 16;
            }
            break;
        case 3:
            while (len--) {
                dest =  0;
                // %%FOR jj in range(0, INT_BITS, 1 << LOG_F)
                dest +=  (uint_mmv_t)(b[0]) << 0;
                dest +=  (uint_mmv_t)(b[1]) << 8;
                dest +=  (uint_mmv_t)(b[2]) << 16;
                dest +=  (uint_mmv_t)(b[3]) << 24;
                dest +=  (uint_mmv_t)(b[4]) << 32;
                dest +=  (uint_mmv_t)(b[5]) << 40;
                dest +=  (uint_mmv_t)(b[6]) << 48;
                dest +=  (uint_mmv_t)(b[7]) << 56;
                // %%END FOR
                *mv++ = dest;
                b += 8;
            }
            break;
        // %%END FOR
    }
}


// %%EXPORT p
void mm_aux_read_mmv24(uint32_t p, uint_mmv_t *mv, uint8_t *b, uint32_t len)
// Read entries of vector mv with modulus p and store these entries in
// the array b. mv is a vector of contigous small arrays with 24 
// entries. len is the number of small arrays to be read. So
// altogether 24 * len entries are read from mv and written to array b.
// Vector b is reduced modulo p.
{
    uint_fast32_t i, sh, tmp;  
    uint_mmv_t source;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    sh = ((i >> 15) & 15);        // This is P_BITS
    i = ((i >> 9) & 3); // This is LOG_FIELD_BITS

    switch(i) {
        // %%FOR LOG_F in [1, 2, 3]
        case 1:
            while (len--) {
                // %%FOR j in range(0, 24)         
                source = *mv++;                 
                tmp = (source >> 0) & p;
                b[0] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 2) & p;
                b[1] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 4) & p;
                b[2] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 6) & p;
                b[3] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[4] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 10) & p;
                b[5] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 12) & p;
                b[6] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 14) & p;
                b[7] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[8] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 18) & p;
                b[9] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 20) & p;
                b[10] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 22) & p;
                b[11] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[12] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 26) & p;
                b[13] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 28) & p;
                b[14] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 30) & p;
                b[15] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 32) & p;
                b[16] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 34) & p;
                b[17] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 36) & p;
                b[18] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 38) & p;
                b[19] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 40) & p;
                b[20] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 42) & p;
                b[21] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 44) & p;
                b[22] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 46) & p;
                b[23] = (tmp + ((tmp + 1) >> sh)) & p;
                // %%END FOR                           
                b += 24;
            }
            break;
        case 2:
            while (len--) {
                // %%FOR j in range(0, 24)         
                source = *mv++;                 
                tmp = (source >> 0) & p;
                b[0] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 4) & p;
                b[1] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[2] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 12) & p;
                b[3] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[4] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 20) & p;
                b[5] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[6] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 28) & p;
                b[7] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 32) & p;
                b[8] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 36) & p;
                b[9] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 40) & p;
                b[10] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 44) & p;
                b[11] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 48) & p;
                b[12] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 52) & p;
                b[13] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 56) & p;
                b[14] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 60) & p;
                b[15] = (tmp + ((tmp + 1) >> sh)) & p;
                source = *mv++;                 
                tmp = (source >> 0) & p;
                b[16] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 4) & p;
                b[17] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[18] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 12) & p;
                b[19] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[20] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 20) & p;
                b[21] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[22] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 28) & p;
                b[23] = (tmp + ((tmp + 1) >> sh)) & p;
                // %%END FOR                           
                b += 24;
            }
            break;
        case 3:
            while (len--) {
                // %%FOR j in range(0, 24)         
                source = *mv++;                 
                tmp = (source >> 0) & p;
                b[0] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[1] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[2] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[3] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 32) & p;
                b[4] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 40) & p;
                b[5] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 48) & p;
                b[6] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 56) & p;
                b[7] = (tmp + ((tmp + 1) >> sh)) & p;
                source = *mv++;                 
                tmp = (source >> 0) & p;
                b[8] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[9] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[10] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[11] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 32) & p;
                b[12] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 40) & p;
                b[13] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 48) & p;
                b[14] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 56) & p;
                b[15] = (tmp + ((tmp + 1) >> sh)) & p;
                source = *mv++;                 
                tmp = (source >> 0) & p;
                b[16] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 8) & p;
                b[17] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 16) & p;
                b[18] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 24) & p;
                b[19] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 32) & p;
                b[20] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 40) & p;
                b[21] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 48) & p;
                b[22] = (tmp + ((tmp + 1) >> sh)) & p;
                tmp = (source >> 56) & p;
                b[23] = (tmp + ((tmp + 1) >> sh)) & p;
                // %%END FOR                           
                b += 24;
                mv += 1;
            }
            break;
        // %%END FOR
    }
}    



// %%EXPORT p
void mm_aux_write_mmv24(uint32_t p, uint8_t *b, uint_mmv_t *mv, uint32_t len)
// Write data from the array b to the vector mv with modulus p. 
// mv is a vector of contiguous small arrays with 24 entries. 
// len is the number of small arrays to be written. So altogether 
// 24 * len entries are read from b and written to the vector mv.
{
    uint_fast32_t i;  
    uint_mmv_t dest;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    i = ((i >> 9) & 3); // This is LOG_FIELD_BITS


    switch(i) {
        // %%FOR LOG_F in [1, 2, 3]
        case 1:
            while (len--) {
                // %%FOR j in range(0, 24)             
                dest = 0;
                dest += (uint_mmv_t)(b[0] & p) << 0;
                dest += (uint_mmv_t)(b[1] & p) << 2;
                dest += (uint_mmv_t)(b[2] & p) << 4;
                dest += (uint_mmv_t)(b[3] & p) << 6;
                dest += (uint_mmv_t)(b[4] & p) << 8;
                dest += (uint_mmv_t)(b[5] & p) << 10;
                dest += (uint_mmv_t)(b[6] & p) << 12;
                dest += (uint_mmv_t)(b[7] & p) << 14;
                dest += (uint_mmv_t)(b[8] & p) << 16;
                dest += (uint_mmv_t)(b[9] & p) << 18;
                dest += (uint_mmv_t)(b[10] & p) << 20;
                dest += (uint_mmv_t)(b[11] & p) << 22;
                dest += (uint_mmv_t)(b[12] & p) << 24;
                dest += (uint_mmv_t)(b[13] & p) << 26;
                dest += (uint_mmv_t)(b[14] & p) << 28;
                dest += (uint_mmv_t)(b[15] & p) << 30;
                dest += (uint_mmv_t)(b[16] & p) << 32;
                dest += (uint_mmv_t)(b[17] & p) << 34;
                dest += (uint_mmv_t)(b[18] & p) << 36;
                dest += (uint_mmv_t)(b[19] & p) << 38;
                dest += (uint_mmv_t)(b[20] & p) << 40;
                dest += (uint_mmv_t)(b[21] & p) << 42;
                dest += (uint_mmv_t)(b[22] & p) << 44;
                dest += (uint_mmv_t)(b[23] & p) << 46;
                // %%END FOR                           
                *mv++ = dest;                
                b += 24;
            }
            break;
        case 2:
            while (len--) {
                // %%FOR j in range(0, 24)             
                dest = 0;
                dest += (uint_mmv_t)(b[0] & p) << 0;
                dest += (uint_mmv_t)(b[1] & p) << 4;
                dest += (uint_mmv_t)(b[2] & p) << 8;
                dest += (uint_mmv_t)(b[3] & p) << 12;
                dest += (uint_mmv_t)(b[4] & p) << 16;
                dest += (uint_mmv_t)(b[5] & p) << 20;
                dest += (uint_mmv_t)(b[6] & p) << 24;
                dest += (uint_mmv_t)(b[7] & p) << 28;
                dest += (uint_mmv_t)(b[8] & p) << 32;
                dest += (uint_mmv_t)(b[9] & p) << 36;
                dest += (uint_mmv_t)(b[10] & p) << 40;
                dest += (uint_mmv_t)(b[11] & p) << 44;
                dest += (uint_mmv_t)(b[12] & p) << 48;
                dest += (uint_mmv_t)(b[13] & p) << 52;
                dest += (uint_mmv_t)(b[14] & p) << 56;
                dest += (uint_mmv_t)(b[15] & p) << 60;
                *mv++ = dest;                
                dest = 0;
                dest += (uint_mmv_t)(b[16] & p) << 0;
                dest += (uint_mmv_t)(b[17] & p) << 4;
                dest += (uint_mmv_t)(b[18] & p) << 8;
                dest += (uint_mmv_t)(b[19] & p) << 12;
                dest += (uint_mmv_t)(b[20] & p) << 16;
                dest += (uint_mmv_t)(b[21] & p) << 20;
                dest += (uint_mmv_t)(b[22] & p) << 24;
                dest += (uint_mmv_t)(b[23] & p) << 28;
                // %%END FOR                           
                *mv++ = dest;                
                b += 24;
            }
            break;
        case 3:
            while (len--) {
                // %%FOR j in range(0, 24)             
                dest = 0;
                dest += (uint_mmv_t)(b[0] & p) << 0;
                dest += (uint_mmv_t)(b[1] & p) << 8;
                dest += (uint_mmv_t)(b[2] & p) << 16;
                dest += (uint_mmv_t)(b[3] & p) << 24;
                dest += (uint_mmv_t)(b[4] & p) << 32;
                dest += (uint_mmv_t)(b[5] & p) << 40;
                dest += (uint_mmv_t)(b[6] & p) << 48;
                dest += (uint_mmv_t)(b[7] & p) << 56;
                *mv++ = dest;                
                dest = 0;
                dest += (uint_mmv_t)(b[8] & p) << 0;
                dest += (uint_mmv_t)(b[9] & p) << 8;
                dest += (uint_mmv_t)(b[10] & p) << 16;
                dest += (uint_mmv_t)(b[11] & p) << 24;
                dest += (uint_mmv_t)(b[12] & p) << 32;
                dest += (uint_mmv_t)(b[13] & p) << 40;
                dest += (uint_mmv_t)(b[14] & p) << 48;
                dest += (uint_mmv_t)(b[15] & p) << 56;
                *mv++ = dest;                
                dest = 0;
                dest += (uint_mmv_t)(b[16] & p) << 0;
                dest += (uint_mmv_t)(b[17] & p) << 8;
                dest += (uint_mmv_t)(b[18] & p) << 16;
                dest += (uint_mmv_t)(b[19] & p) << 24;
                dest += (uint_mmv_t)(b[20] & p) << 32;
                dest += (uint_mmv_t)(b[21] & p) << 40;
                dest += (uint_mmv_t)(b[22] & p) << 48;
                dest += (uint_mmv_t)(b[23] & p) << 56;
                // %%END FOR                           
                *mv++ = dest;                
                *mv++ = 0;
                b += 24;
            }
            break;
        // %%END FOR
    }

} 



// %%EXPORT p
uint8_t mm_aux_get_mmv1(uint32_t p, uint_mmv_t *mv, uint32_t i)
// Return the entry with given index of vector mv with modulus p.
// The result is reduced modulo p.
{
    uint_fast32_t  c, j, res;

    // %%MMV_LOAD_CONST p, c
    // Store constant table for p to c
    c = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    j = ((c) & 7); // This is LOG_INT_FIELDS
    mv += i >> j;
    i &= (1 << j) - 1;
    j = 6 - j;           // This is LOG_FIELD_BITS
    res = (mv[0] >> (i << j)) & p;    // This is the result
    j = ((c >> 15) & 15);         // This is P_BITS
    // return result reduced modulo p
    return (uint8_t) ( (res + ((res + 1) >> j)) & p );
}


// %%EXPORT p
void mm_aux_put_mmv1(uint32_t p, uint8_t value, uint_mmv_t *mv, uint32_t i)
// Set the entry of the vector mv with modulus p at the  index i  
// to the given value. 0 <= value <= p must hold.
{
    uint_fast32_t  j;

    // %%MMV_LOAD_CONST  p, j
    // Store constant table for p to j
    j = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    j = ((j) & 7); // This is LOG_INT_FIELDS
    mv += i >> j;
    i &= (1 << j) - 1;
    j = 6 - j;           // This is LOG_FIELD_BITS
    i <<= j;
    mv[0] &= ~(((uint_mmv_t)p) << i);
    mv[0] |= ((uint_mmv_t)(value & p)) << i;
}


// %%COMMENT
/**********************************************************************
*** Functions for data transfer from and to vectors in R_p.
*** Here such a vector is given in internal representation and of type 
*** uint_mmv_t[]. For modulus p, p + 1 must be a power of two.
**********************************************************************/


// %%EXPORT p
uint32_t mm_aux_mmv_size(uint32_t p)
// Return number of integers of type uint_mmv_t required to store
// a vector of the representation R_p for a given p.
{
    uint_fast32_t tbl;
    // %%MMV_LOAD_CONST  p, tbl
    // Store constant table for p to tbl
    tbl = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    // return the value MMV_INTS for the specified p
    return(247488 >> ((tbl) & 7));
}


// %%EXPORT p
void mm_aux_zero_mmv(uint32_t p, uint_mmv_t *mv)
// Zero the vector of the representation R_p referred by mv,
{
    uint_fast32_t j;
    if (mm_aux_bad_p(p)) return;
    // %%MMV_LOAD_CONST  p, j
    // Store constant table for p to j
    j = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    j = (247488 >> ((j) & 7));
    do {
        *mv++ = 0;
    } while(--j);
}



// %%EXPORT p
uint8_t mm_aux_get_mmv(uint32_t p, uint_mmv_t *mv, uint32_t i)
// Return the entry of the vector mv in R_p at  index  i.
{
    uint_fast32_t  j, c, res;
    if (mm_aux_bad_p(p)) return 0;
    if (i <  MM_AUX_XOFS_X) {
        if (i <  MM_AUX_XOFS_T) {
            // Tags A, B, C
            i = (MM_AUX_TBL_ABC[i] & 0x7ff) + i - 24;
        } else {
            // Tag T
            i += MM_AUX_OFS_T - MM_AUX_XOFS_T;
        } 
    } else {
        if (i >=  MM_AUX_XOFS_E) return 0;
        // Tags X, Z, Y
        i -=  MM_AUX_XOFS_X;
        // Put i += 8 * floor(i/24), for i <  3 * 2048 * 24
        i += (((i >> 3) * 0xaaab) >> 17) << 3; 
        i += MM_AUX_OFS_X;
    }

    // %%MMV_LOAD_CONST  p, c
    // Store constant table for p to c
    c = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    j = ((c) & 7); // This is LOG_INT_FIELDS
    mv += i >> j;
    i = (i & ((1 << j) - 1)) << (6 - j);
    res = (mv[0] >> i) & p;             // This is the result
    //Reduce reult modulo p
    c = ((c >> 15) & 15); // This is P_BITS
    return (uint8_t)((res + ((res + 1) >> c)) & p);
}











// %%EXPORT p
void mm_aux_put_mmv(uint32_t p, uint8_t value, uint_mmv_t *mv, uint32_t i)
// Set the entry of the vector mv in R_p at index i to the 
// given value. 0 <= value <= p must hold.
{

    uint_fast32_t  j, sh, diff;
    if (mm_aux_bad_p(p)) return;
    value &= p;

    // %%MMV_LOAD_CONST  p, j
    // Store constant table for p to j
    j = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    j = ((j) & 7); // This is LOG_INT_FIELDS

    if (i <  MM_AUX_XOFS_X) {
        if (i <  MM_AUX_XOFS_T) {
            // Tags A, B, C
            diff = 31 * (MM_AUX_TBL_ABC[i] >> 11);
            i = (MM_AUX_TBL_ABC[i] & 0x7ff) + i - 24;
            sh = (i & ((1 << j) - 1)) << (6 - j);
            mv[i >> j] = (mv[i >> j] &  ~(((uint_mmv_t)p) << sh))
                  |   ((uint_mmv_t)(value)) << sh;
            i -= diff;
        } else {
            // Tag T
            i += MM_AUX_OFS_T - MM_AUX_XOFS_T;
        } 
    } else {
        if (i >=  MM_AUX_XOFS_E) return;
        // Tags X, Z, Y
        i -=  MM_AUX_XOFS_X;
        // Put i += 8 * floor(i/24), for i <  3 * 2048 * 24
        i += (((i >> 3) * 0xaaab) >> 17) << 3; 
        i += MM_AUX_OFS_X;
    }

    mv += i >> j;
    sh = (i & ((1 << j) - 1)) << (6 - j);
    mv[0] = (mv[0] &  ~(((uint_mmv_t)p) << sh))
              |   ((uint_mmv_t)(value)) << sh;
}



// %%EXPORT p
void mm_aux_random_mmv(uint32_t p, uint8_t *seed, uint_mmv_t *mv)
// Randomize the vector mv in R_p with the internal random generator.
// The random generator must be seeded by function mm_rng_seed().
{
    uint8_t b1[3072];
    uint_fast32_t i, c;

    if (mm_aux_bad_p(p)) return;
    // %%MMV_LOAD_CONST p, c
    // Store constant table for p to c
    c = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    c = ((c) & 7); // This is LOG_INT_FIELDS

    // Do the small part
    mm_rng_gen_modp(p, seed, (uint8_t *)(mv), 24 + 3 * 276);
    mm_aux_small24_expand((uint8_t *)(mv), b1);
    mm_aux_write_mmv24(p, b1, mv, 72);
    mv += MM_AUX_OFS_T >> c;

    // Do the 759 * 64 vector; note that 759 = 11 * 69
    for (i = 0; i < 22; ++i) {
        mm_rng_gen_modp(p, seed, b1, 69 * 32);
        mm_aux_write_mmv1(p, b1, mv, 69 * 32);
        mv += (69 * 32) >> c;
    } 

    // Do the 6144 * 24 vector
    for (i = 0; i < 48; ++i) {  
        mm_rng_gen_modp(p, seed, b1, 3072);
        mm_aux_write_mmv24(p, b1, mv, 128);
        mv += 4096 >> c;
    } 
}


// %%COMMENT
/**********************************************************************
*** Reducing and checking a vector in R_p
**********************************************************************/


// %%EXPORT p
int32_t mm_aux_reduce_mmv(uint32_t p, uint_mmv_t *mv)
// Reduce the vector mv in R_p modulo p. Note that a zero entry
// is represented either by 0...0 or by 1...1. This functions sets 
// all zero entries of the vector mv to 0...0.
// The function returns 0 if it detects no error.
// It may return the following error codes:
//  -1: Bad input value p
//  -2: A one bit outside a field has been found
{
    uint_fast32_t i, sh;
    uint_mmv_t data, cy, mask_1, mask_p, acc;

    if (mm_aux_bad_p(p)) return -1;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    sh = ((i >> 15) & 15);          // This is P_BITS
    i = ((i) & 7);   // This is LOG_INT_FIELDS
    i = MM_AUX_OFS_E >> i;               // No of unit_mmv_t to process
    mask_1 =  MM_AUX_TBL_REDUCE[2*sh-4];
    mask_p =  MM_AUX_TBL_REDUCE[2*sh-3];
    if (sh & (sh - 1)) {
        // case P_BITS is not a power of two
        acc = 0;      // use acc for accumulating error bits
        do {
            data = *mv;
            acc |= data;
            data &= mask_p;
            cy = (data + mask_1) & ~mask_p;
            data += (cy >> sh) - cy;
            *mv++ = data;
        } while (--i);
        if (acc & ~mask_p) return -2;
    } else {
        // case P_BITS is a power of two
        // use acc for  (<high half> & <low half>) of fields
        sh >>= 1;   // halved P_BITS
        do {
            data = *mv;
            acc = data & (data >> sh) & mask_p;
            cy = (acc + mask_1) & ~mask_p;
            data += (cy >> sh) - (cy << sh);
            *mv++ = data;
        } while (--i);
    }
    return 0;
}




static int32_t check24(uint32_t p, uint_mmv_t *mv, uint32_t length)
{
    uint_fast32_t d;
    uint_mmv_t  acc = 0, mask;

    // %%MMV_LOAD_CONST  p, d
    // Store constant table for p to d
    d = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    // Put d =  5 - LOG_INT_FIELDS
    d = 5 - ((d) & 7); 
    switch (d) {
        // %%IF {INT_BITS} == 64
        case 0:
            mask = 0xffff000000000000ULL;
            while (length--) acc |= *mv++ & mask;
            break;
        // %%END IF
        case 1:
            mask = 0xffffffff00000000ULL;
            while (length--) {
                acc |= mv[1] & mask; mv += 2;
            }
            break;
        case 2:
            while (length--) {
                acc |= mv[3]; mv += 4;
            }
            break;
        // %%IF {INT_BITS} == 32
        // %%END IF
        break;
    }
    return acc ? -3 : 0;
}


static int32_t check_sym(uint32_t p, uint_mmv_t *mv, uint8_t * buffer)
// buffer is a work buffer of size at least 72 * 32
// Side effect: Tmv entries with tags A, B, C are read to the buffer
{
    uint_fast32_t i, acc = 0;
    uint8_t  *p_row, *p_col;
    mm_aux_read_direct_mmv1(p, mv, buffer, 72*32);
    for(i = 768; i < 1536; i += 33)  
        acc |= buffer[i] | buffer[i + 768];
    if (acc) return -4;
    p_row = buffer;
    acc = 0;
    for (p_col = buffer; p_col < buffer + 24; ++p_col) {
         for (i = 0; i < 24; ++i)
             acc |= (p_row[i] ^ p_col[i << 5]) 
                  | (p_row[i + 768] ^ p_col[(i << 5) + 768])
                  | (p_row[i + 1536] ^ p_col[(i << 5) + 1536]);
         p_row+= 32;
    }
    return acc ? -5 : 0;
}



static int32_t check_mmv_buffer(uint32_t p, uint_mmv_t *mv, uint8_t * buffer)
// Workhorse for function mm_aux_check_mmv. buffer must have size 72*32.
// Side effect: mv entries with tags A, B, C are read to the buffer
{
    uint_fast32_t i;
    i = mm_aux_reduce_mmv(p, mv);
    if (i) return i;                 // Errors -1, -2 may occur here
    i = check24(p, mv, 72);          // check tags A,B,C
    if (i) return i;                 // Error -3 may occur here
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    i = ((i) & 7); //  LOG_INT_FIELDS
    i = check24(p, mv + (MM_AUX_OFS_X >> i), 6144); // check tags X,Y,Z
    if (i) return i - 100;                 // Error -3 may occur here
    return check_sym(p, mv, buffer); // Errors -4, -5 may occur here
}


// %%EXPORT p
int32_t mm_aux_check_mmv(uint32_t p, uint_mmv_t *mv)
// Check the vector mv in R_p modulo p. 
// The function returns 0 if it detects no error.
// It may return the following error codes:
//  -1: Bad input value p
//  -2: A one bit outside a field has been found
//  -3: A subfield has an illegal nonzero entry at index >= 24
//  -4: Illegal nonzero diagonal entry 
//  -5: Symmetric part of vector is not symmetric 
// As a side effect, mv is reduced with mm_aux_reduce_mmv(p, mv).
{
    uint8_t buffer[72*32];
    return check_mmv_buffer(p, mv, buffer);
}

// %%COMMENT
/**********************************************************************
*** Low-level functions supperting external rep of vectors in R_p
**********************************************************************/

// %%EXPORT p
void mm_aux_small24_expand(uint8_t *b_src, uint8_t *b_dest)
//
{
    uint_fast16_t j0, j1t, j1e;
    uint8_t *b_transpose = b_dest;
    for (j0 = 0; j0 < 24 * 25; j0 += 25) {
        b_dest[j0] = *b_src++;
        b_dest[j0 + 1152] =  b_dest[j0 + 576] = 0;
    }

    for (j0 = 0; j0 < 24; ++j0)  {
        j1e = 24 * j0;
        for (j1t = 0; j1t < j1e; j1t += 24) {
            b_transpose[j1t] = b_dest[0] = b_src[0];
            b_transpose[j1t + 576] = b_dest[576] = b_src[276];
            b_transpose[j1t + 1152] = b_dest[1152] = b_src[552];
            ++b_dest; ++b_src;
        }
        b_dest += 24 - j0;
        ++b_transpose;
    }
}



// %%EXPORT p
void mm_aux_small24_compress(uint8_t *b_src, uint8_t *b_dest)
//
{
    uint_fast16_t  j0, j1;
    for (j0 = 0; j0 < 24 * 25; j0 += 25) 
        *b_dest++ = b_src[j0];
    for (j0 = 0; j0 < 24; ++j0)  {
        for (j1 = j0; j1; --j1) {
            b_dest[0] = b_src[0];
            b_dest[276] = b_src[576];
            b_dest[552] = b_src[1152];
            ++b_dest; ++b_src;
        } 
        b_src += 24 - j0;
    }
}






// %%COMMENT
/**********************************************************************
*** Conversion between internal and external rep of vectors in R_p
**********************************************************************/








// %%EXPORT p
void mm_aux_mmv_to_bytes(uint32_t p, uint_mmv_t *mv, uint8_t *b)
// Convert the vector mv in R_p to an array b of bytes. p is the
// modulus for mv. b must have length 196884.  
{
    uint8_t b1[3*576]; 
    uint_fast32_t c;

    // %%MMV_LOAD_CONST  p, c
    // Store constant table for p to c
    c = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    c = ((c) & 7); // This is LOG_INT_FIELDS
    
    mm_aux_read_mmv24(p, mv, b1, 72);
    mm_aux_small24_compress(b1, b);
    mv += MM_AUX_OFS_T >> c;
    b +=  MM_AUX_XOFS_T;
    mm_aux_read_mmv1(p, mv, b,  759*64);
    mv += (MM_AUX_OFS_X - MM_AUX_OFS_T) >> c;
    b +=  (MM_AUX_XOFS_X - MM_AUX_XOFS_T);
    mm_aux_read_mmv24(p, mv, b, 6144);
}



// %%EXPORT p
void mm_aux_bytes_to_mmv(uint32_t p, uint8_t *b, uint_mmv_t *mv)
// Convert the array b of 196884 bytes to the vector mv in R_p to.
// p is the modulus for mv.  
{
    uint8_t b1[3*576];
    uint_fast32_t  c;

    // %%MMV_LOAD_CONST  p, c
    // Store constant table for p to c
    c = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    c = ((c) & 7); // This is LOG_INT_FIELDS

    mm_aux_small24_expand(b, b1);   
    mm_aux_write_mmv24(p, b1, mv, 72);
    mv += MM_AUX_OFS_T >> c;
    b +=  MM_AUX_XOFS_T;
    mm_aux_write_mmv1(p, b, mv, 759*64);
    mv += (MM_AUX_OFS_X - MM_AUX_OFS_T) >> c;
    b +=  (MM_AUX_XOFS_X - MM_AUX_XOFS_T);
    mm_aux_write_mmv24(p, b, mv, 6144);
}




// %%COMMENT
/**********************************************************************
*** Conversion between internal and sparse rep of vectors in R_p
**********************************************************************/

// %%EXPORT p
int32_t mm_aux_mmv_to_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp)
// Convert the vector 'mv' in the representation R_p of the monster, 
// with 'mv' given in internal representation, to sparse format. The
// function writes the entries of the vector in sparse format to the
// array referred by 'sp' and returns the number of entries written
// to 'sp'. The array 'sp' may have up to 196884 entries.
{
    int32_t status;
    uint_fast32_t row, row_end,i, j, isp = 0, value;
    uint_fast32_t field_bits, lg_int_fields, ofs, sh;
    uint_mmv_t source;
    uint8_t b[72*32], *p_row;

    if ((status = check_mmv_buffer(p, mv, b)) != 0) return status;
    
    // %%MMV_LOAD_CONST  p, j
    // Store constant table for p to j
    j = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    field_bits = ((j >> 11) & 15); // This is FIELD_BITS
    lg_int_fields = ((j) & 7); // This is LOG_INT_FIELDS
    sh = 8 - 6 + lg_int_fields; // This is 8 - LOG_FIELD_BITS; 
 
    // Do tags A, B, C
    p_row = b;
    for (row = 0; row < 3; ++row) for (i = 0; i < 24; ++i) {
         for (j = 0; j <= i; ++j) {
            if ((value = p_row[j]) != 0)  sp[isp++] =
                0x2000000 + (row << 25) + (i << 14) + (j << 8) + value; 
        } 
        p_row += 32;
    }
    
    // Do tag T
    mv += MM_AUX_OFS_T >> lg_int_fields;
    row_end = (MM_AUX_OFS_X - MM_AUX_OFS_T) >> lg_int_fields;
    for (row = 0; row < row_end; ++row) if ((source = *mv++) != 0) {
        ofs = 0x8000000 + (row << (8 + lg_int_fields));
        for (j = 0; j < 64; j += field_bits) {
            if ((value = (source >> j) & p) != 0)  {
                sp[isp++] = ofs + (j << sh) + value;
            } 
        }           
    }

    row_end = (MM_AUX_OFS_E - MM_AUX_OFS_X) >> lg_int_fields;
    for (row = 0; row < row_end; ++row) if ((source = *mv++) != 0) {
        ofs = 0x5000000 + (row << (8 + lg_int_fields));
        ofs += ofs & 0xfffe000;
        for (j = 0; j < 64; j += field_bits) {
            if ((value = (source >> j) & p) != 0)  {
                 sp[isp++] = ofs + (j << sh) + value;
            } 
        }           
    }

    return (int32_t)isp; 
}

// %%EXPORT p
void mm_aux_mmv_extract_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp, uint32_t length)
// Extract certain entries of the vector  'mv' depending on the 
// vector 'sp'. The entries of vector 'sp' are updated with the
// corresponding entries of 'mv'. Here 'mv'  and 'sp' are vectors 
// in the representation R_p of the monster, with 'mv' given in 
// internal representation and 'sp' given in sparse format. 'sp' 
// has  length 'length'. If 'sp' has an entry with a certain label
// then the coordinate of that entry is set to the corresponding 
// coordinate of vector 'mv'. If several entries of 'sp' have the 
// same label then the same coordinate is taken from 'mv' 
// several times.
{
    uint_fast32_t i, lg_int_fields, lg_field_bits, p_bits, index_mask;

    if (mm_aux_bad_p(p)) return;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    p_bits = ((i >> 15) & 15); // This is P_BITS
    lg_field_bits = ((i >> 9) & 3); // This is LOG_FIELD_BITS
    lg_int_fields = ((i) & 7); // This is LOG_INT_FIELDS
    index_mask = (1 << lg_int_fields) - 1;

    for ( ;length--; ++sp) {
        uint_fast32_t v = *sp & ~0xff, index;
        uint_fast32_t tag = v >> 25,i = (v >> 14) & 0x7ff, j = (v >> 8) & 0x3f;
        switch (tag) {
            case 2:  // tag B
            case 3:  // tag C
                if (i == j) continue;
                // Fall trough to case tag A
            case 1:  // tag A
                if (i >= 24 || j >= 24) continue;
                index = (tag - 1) * 768  + (i << 5) + j;
                break;
            case 4:  // tag T
                if (i >= 759) continue;
                index = MM_AUX_OFS_T + (i << 6) + j;
                break;
            case 5:  // tag X
            case 6:  // tag Z
            case 7:  // tag Y
                if (j >= 24) continue;
                index = ((v >> 14) << 5) + j - 0x50000 + MM_AUX_OFS_X;
                break;
            default:
                continue;
        }
        i = ((index) & index_mask) << lg_field_bits; 
        i = (mv[index >> lg_int_fields] >> i) & p; 
        // i is the entry of vector mb. reduce i modulo p
        i = (i + ((i + 1) >> p_bits)) & p;
        *sp = v + i; 
    }
}




#define add_sparse_value(index, value) \
    { \
        uint_fast32_t idx = (index) >> lg_int_fields; \
        uint_fast32_t sh = ((index) & index_mask) << lg_field_bits; \
        uint_mmv_t old_value = (mv[idx] >> sh)  & p; \
        uint_mmv_t new_value = old_value + (value & p); \
        new_value = (new_value + (new_value >> p_bits)) & p; \
        mv[idx]  ^=  (old_value ^ new_value) << sh; \
    }
        

// %%EXPORT p
void mm_aux_mmv_add_sparse(uint32_t p, uint32_t *sp, uint32_t length, uint_mmv_t *mv)
// Add vector 'sp' to vector 'mv'. Here 'mv' in and 'sp' are vectors 
// in the representation R_p of the monster, with 'mv' given in internal
// representation and 'sp' given in sparse format. 'sp' has length 
// 'length'. Entries in 'sp' with the same index are added up.
{
    uint_fast32_t i, lg_int_fields, lg_field_bits, p_bits, index_mask;

    if (mm_aux_bad_p(p)) return;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    lg_field_bits = ((i >> 9) & 3); // This is LOG_FIELD_BITS
    lg_int_fields = ((i) & 7); // This is LOG_INT_FIELDS
    p_bits = ((i >> 15) & 15);               // This is P_BITS
    index_mask = (1 << lg_int_fields) - 1;

    for ( ;length--; ++sp) {
        uint_fast32_t v = *sp, index;
        uint_fast32_t tag = v >> 25,i = (v >> 14) & 0x7ff, j = (v >> 8) & 0x3f;
        switch (tag) {
            case 2:  // tag B
            case 3:  // tag C
                if (i == j) continue;
                // Fall trough to case tag A
            case 1:  // tag A
                if (i >= 24 || j >= 24) continue;
                index = (tag - 1) * 768 + (i << 5) + j;
                if (i != j) add_sparse_value(index, v);
                index += 31 * (j - i);
                break;
            case 4:  // tag T
                if (i >= 759) continue;
                index = MM_AUX_OFS_T + (i << 6) + j;
                break;
            case 5:  // tag X
            case 6:  // tag Z
            case 7:  // tag Y
                if (j >= 24) continue;
                index = ((v >> 14) << 5) + j - 0x50000 + MM_AUX_OFS_X;
                break;
            default:
                continue;
        }
        add_sparse_value(index, v);
    }
}



    

#define set_sparse_value(index, v) \
    { \
        uint_fast32_t idx = (index) >> lg_int_fields; \
        uint_fast32_t sh = ((index) & index_mask) << lg_field_bits; \
        uint_mmv_t value = ((mv[idx] >> sh)  ^ v) & p; \
        mv[idx]  ^= value << sh; \
    }




// %%EXPORT p
void mm_aux_mmv_set_sparse(uint32_t p, uint_mmv_t *mv, uint32_t *sp, uint32_t length)
// Set certain entries of the vector  'mv' depending on the vector
// 'sp'. Here 'mv' in and 'sp' are vectors in the representation R_p 
// of the monster, with 'mv' given in internal representation and 
// 'sp' given in sparse format. 'sp' has  length 'length'.
// If 'sp' has an entry with label 'l' then the corresponding entry 
// of 'mv' is set to to the value coded in the entry of 'sp'.
{
    uint_fast32_t i, lg_int_fields, lg_field_bits, index_mask;

    if (mm_aux_bad_p(p)) return;
    // %%MMV_LOAD_CONST  p, i
    // Store constant table for p to i
    i = MMV_CONST_TABLE[((((p) + 1) * 232) >> 8) & 7];
    lg_field_bits = ((i >> 9) & 3); // This is LOG_FIELD_BITS
    lg_int_fields = ((i) & 7); // This is LOG_INT_FIELDS
    index_mask = (1 << lg_int_fields) - 1;

    for ( ;length--; ++sp) {
        uint_fast32_t v = *sp, index;
        uint_fast32_t tag = v >> 25, i = (v >> 14) & 0x7ff, j = (v >> 8) & 0x3f;
        switch (tag) {
            case 2:  // tag B
            case 3:  // tag C
                if (i == j) continue;
                // Fall trough to case tag A
            case 1:  // tag A
                if (i >= 24 || j >= 24) continue;
                index = (tag - 1) * 768 + (i << 5) + j;
                set_sparse_value(index, v);
                index += 31 * (j - i);
                break;
            case 4:  // tag T
                if (i >= 759) continue;
                index = MM_AUX_OFS_T + (i << 6) + j;
                break;
            case 5:  // tag X
            case 6:  // tag Z
            case 7:  // tag Y
                if (j >= 24) continue;
                index = ((v >> 14) << 5) + j - 0x50000 + MM_AUX_OFS_X;
                break;
            default:
                continue;
        }
        set_sparse_value(index, v);
    }
}




// %%COMMENT
/**********************************************************************
*** Index conversion between external and sparse rep of vectors in R_p
**********************************************************************/



// %%EXPORT p
uint32_t mm_aux_index_extern_to_sparse(uint32_t i)
// Convert external index i to sparse index.
// Return 0 if index i is bad
{
    if (i <  MM_AUX_XOFS_X) {
        if (i <  MM_AUX_XOFS_T) {
            // Tags A, B, C
            i = (MM_AUX_TBL_ABC[i] & 0x7ff) + i - 24;
            // put i += (i / 0x300) * 0x100; assuming 0 <= i < 0x900 
            i += (0x2A54000 >> ((i >> 8) << 1)) & 0x300;
            // now 0 <= i < 0xc00. output bits of old i as 
            // (tag - 1) = bits 11..10, i = bits 9..5, j = bits 4..0
            return 0x2000000 + ((i & 0xc00) << 15) +
                   ((i & 0x3e0) << 9) + ((i & 0x1f) << 8);
        } else {
            // Tag T
            i += 0x80000 - MM_AUX_XOFS_T;
            return i << 8;
        } 
    } else if (i <  MM_AUX_XOFS_E) {
        // Tags X, Z, Y
        i -=  MM_AUX_XOFS_X;
        // Put i += 8 * floor(i/24), for i <  3 * 2048 * 24
        i += (((i >> 3) * 0xaaab) >> 17) << 3; 
        // shift bits 17..5 of i to bit positions 18...6
        i += i & 0x3ffe0;
        i += 0xA0000;
        return i << 8;
    } else return 0;
}



// %%EXPORT p
void mm_aux_array_extern_to_sparse(uint32_t *a, uint32_t len)
// Convert array a of external indices to sparse indices in place.
// Each index is converted using mm_aux_index_extern_to_sparse() 
{
    for(; len--; ++a) *a = mm_aux_index_extern_to_sparse(*a); 
}



// %%EXPORT p
int32_t mm_aux_index_sparse_to_extern(uint32_t i)
// Convert sparse index i to external index.
// Return -1 if index i is bad
{
    uint_fast32_t tag = i >> 25, j = (i >> 8) & 0x3f;
    i = (i >> 14) & 0x7ff;
    switch (tag) {
        case 2:  // tag B
        case 3:  // tag C
            if (i == j) return -1;
            // Fall trough to case tag A
        case 1:  // tag A
            if (i >= 24 || j >= 24) return -1;
            if (i == j) return i;
            return  MM_AUX_XOFS_A - 276 + tag * 276 
                  + ((i * i - i) >> 1) + j;
        case 4:  // tag T
            if (i >= 759) return -1;
            return MM_AUX_XOFS_T + (i << 6) + j;
        case 5:  // tag X
        case 6:  // tag Z
        case 7:  // tag Y
            if (j >= 24) return -1;
            return  MM_AUX_XOFS_X - 0x3c000
                + 24 * ((tag << 11) + i) + j; 
        default:
            return -1;
    }
}
