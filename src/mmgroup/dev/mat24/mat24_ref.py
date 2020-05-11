"""Implement binary (24,12,8) Golay code and Mathieu group Mat24 



"""
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals
from functools import reduce


from random import randint
import types
import sys
import re
import os
import subprocess
from operator import __or__, __xor__

import numpy


from mmgroup.bitfunctions import bw24, bit_mat_mul, bitparity

from mmgroup.dev.mat24.mat24tables import Mat24Tables

from mmgroup.dev.mat24.mat24aux import Lsbit24Function, MatrixToPerm

from mmgroup.dev.mat24.mat24heptad import HeptadCompleter

from mmgroup.dev.mat24.mat24theta import make_augmented_theta_table
from mmgroup.dev.mat24.mat24theta import make_autpl_qf_table

from mmgroup.dev.mat24.make_addition_table import generate_c_bitmatmul

from mmgroup.dev.mat24.make_mul_transp import generate_c_mul_transp

from mmgroup.dev.generate_c.make_c_tables import TableGenerator, make_doc
from mmgroup.dev.generate_c.generate_functions import UserDirective
from mmgroup.dev.generate_c.generate_functions import UserFormat



TABLE_NAMES = None

###########################################################################
###########################################################################
# The Golay code (a pure python implementation)
###########################################################################
###########################################################################


class Mat24(Mat24Tables): 
    """Provide functions for the Mathieu group Mat24 and the Parker loop.

    The Golay code C and its codode C*
    ---------------------------------- 

    The Mathieu group Mat24 operates as a permutation group on a set of 
    24 elements which we label with numbers 0,...,23 for use in Python 
    and C. So it also operates on a vector space V = GF(2)**24, with 
    GF(2) = {0,1}. 

    A vector v in a vector space over GF(2) is called a bit vector. We
    represent a bit vector as an integer, so that the i-th bit of v 
    (with valence 2**i) is the i-th component of v. 

    Mat24 is the automorphism group of the binary [24,12,8] Golay code C
    So C is a 12-dimensional subspace of V. Code words in C have weight
    0, 8, 12, 16 or 24. Up to isomorphism there is only one such Golay
    code. Our Golay code in V is compatible to the 'MOG' in [CoSl99], 
    Ch. 11 with numbering of the unit vectors in V as follows:

                         +------------------------+
                         |  0   4   8  12  16  20 |                 
                         |  1   5   9  13  17  21 |                 
                         |  2   6  10  14  18  22 |                 
                         |  3   7  11  15  19  23 |                 
                         +------------------------+

   
    There are member functions for checking and completing codewords and 
    for getting the syndrome or the type (i.e. the orbit under Mat24) of 
    a 24-bit vector.

    All relevant functions in this python class are implemented as 
    class methods. There are also corresponding funtions implemented
    in C. 

    We internally use a basis of GF(2**24) such that the first 12 basis
    vectors are a transversal of the Golay cocode and the last 12 basis 
    vectors span the Golay code. Here '**' means exponentiation. Our
    cocode basis is the reciprocal basis of the Golay code basis.

    The assured properties of the (transversal of the) cocode basis are:

        - All basis vectors except vector 0 have an even number of ones.
        - Basis vector 0 has a one in row 0, column 0 only.
        - Basis vectors 0,...,4 and 11 have ones in row 0 only.
        - The other basis have two nonzero entries in the last 3 rows 
          of exactly one column and zero entries everywhere else.
        - Basis vector 11 (i.e. the last one) has ones precisely in all 
          six columns of row 1

    The assured properties of the basis of the Golay code are:

        - Basis vector 0 has 24 ones, all other basis vectors have 
          weight 8 or 16.
        - All but basis vector 11 are even interpretations of the MOG, 
          i.e. the first MOG row and all MOG columns have even bit weight. 
        - In the first five basis vectors the entries in each column  
          are equal.
        - Basis vectors 5,...,10  have zero entries in row 0.
        - Basis vector 11 has ones precisely in row and column 0, except 
          for intersection of row and column 0, where it is zero.
     
    The user may rely on these properties for operations on the vectors.
      
    We represent vector in V, C and V* and in sthe subspace of 
    octads of C as follows:

    The 759 octads are numbered from 0 to 758. This is not a vector space.
    The 2**12 code words are represented as binary numbers 0 to 4095.
    The 2**12 cocode words are represented as binary numbers 0 to 4095.

    As usual, binary numbers representing bit vectors are added with the 
    XOR operation '^'. Unused high bits are ignored. The 2**12 code words
    and the 2**12 cocode words mentioned above can be added with '^'.

    Member functions changing from one representation xxx to another
    representation yyy are named xxx_to_yyy, where xxx, yyy is as follows:

    vect :     standard representation of a bit vector in V = GF(2)**24
               coded as a 24-bit integer. 
    vintern:   internal representation a bit vector in V as a vector in
               the basis given above, coded as 24-bit integer. 
    gcode:     representation of a Golay code word in the basis given
               given above, coded as 12-bit integer. Here we simply take 
               the upper 12 bits of the internal representation.
    octad:     representation as an octad numbered from 0 to 758
               (in lexical order given by representation  'gcode')
    cocode:    representation as a cocode word in the basis given above, 
               coded as 12-bit integer. Here we simply take the lower 
               12 bits of the internal representation.
    

    All these representations are given as integers.
     
    We implement the following conversion functions: 

        vect_to_vintern, vintern_to_vect, vect_to_cocode, 
        vintern_to_vect, gcode_to_vect, cocode_to_vect.

    Here irrelevant bits of the input are ignored. cocode_to_vect
    returns one of many possible solutions.

    In the following functions the input is checked and an exception
    is raised in case of an error:

        vect_to_gcode,   vect_to_octad,   gcode_to_octad,    
        octad_to_vect, octad_to_gcode
    
    Here we raise an exception if the given input is not a code word
    or an octad, as required by the member function. The correspinding
    C functions returns -1 (casted to the appropriate unsigned integer    
    type) in case of any error.  
       
    Function syndrome() takes a vector v and calculates its syndrome,
    which is a vector of minimum weight equivalent to v modulo the 
    Golay code. Function cocode_syndrome() takes a 'cocode'
    representation as a cocode word instead.

    Function scalar_prod() returns the scalar product of a Golay code
    vector in 'gcode' and a cocode vector in 'cocode' representation.


    The Mathieu group Mat24
    ----------------------- 

    This class also contains support for the Mathieu group Mat24.
    An element of Mat24 can be represented in one of the following ways:

    perm:     Representation as a array of length 24 encoding a 
              permutation of the integers 0,...,23 as a mapping.

    m24num:   Representation as an integer 0 <= i < 244823040. The
              identity permutation is coded as 0. Other codes are
              more or less arbitrary, see function m24num_to_perm().

    matrix:   Representation as a 12 x 12 bit matrix acting on the Golay
              code by right multiplication. This matrix acts on a Golay 
              code vectors (given in the 'gcode' representation) by
              right multiplication. 
              Such a matrix is implemented as an array of integers with
              each integer corresponding to a row vector of the matrix. 
              The raison d'etre of this representation is to support 
              the Parker loop and its automorphism group. Therefore a
              row vector is implemented as a 32-bit integer.

    We implement the following conversion functions

        m24num_to_perm, perm_to_m24num, perm_to_matrix, matrix_to_perm.

    There is a function perm_check() for checking if an array of 
    length 24 really represents an element of the Mathieu group Mat24.
    All other function operating on Mat24 in any way do not check if
    their inputs are really in Mat24. They will output garbage on bad
    input, but they are not suppused to crash.

    The easiest way to create a random element of Mat24 is to create 
    a random integer 0 <= n < 244823040, and to call function
    m24num_to_perm(x). You can use functions perm_from_heptads() or
    perm_complete_heptad() to create specific elements, see 
    documentation of these function for details. The group Mat24 is
    discussed in detail in [CoSl99] Ch. 11.


    Operation of the group Mat24 on vectors
    --------------------------------------- 

    Elements of Mat24 operate from the right on vectors in V = (2)**24
    or on Golay code or cocode vectors, as usual in finite group theory. 
    A function performing such an operation has the name

            op_<vector>_<group>

    where <vector> indicates the representation of the vector space and
    <group> indicates the representation of the group. We implement the
    functions

        op_vect_perm, op_gcode_matrix, op_gcode_perm, op_cocode_perm.

    E.g. op_gcode_matrix operates on a Golay code word (in 'gcode'
    representation) by right multiplying an element m of Mat24 with it.
    Here element m is a 12 x 12 matrix (in 'matrix' representation').      

    
    Group operation in the group Mat24
    ----------------------------------

    Multiplication and inversion in the group mat2 is supported for 
    the permutation representation 'perm'. Therefore we have functions
 
         mul_perm, inv_perm



    The Parker loop Pl
    ------------------

    We support the Parker loop Pl and also its automorphism group.

    The Parker loop Pl is nonassociative loop with a central element 
    -1 such that Pl/{1,-1} = C, with C the Golay code. We choose a
    transversal of C in Pl and write the elements of Pl as pairs
    (v, s) with v in C and s in {1,-1} and multiplication rule
 
      (v1, s1) * (v2, s2) = (v1 + v2, s1 * s2 * (-1)**theta(v1, v2)),  

    where the cocycle theta is quadratic in the first and linear in 
    the second element. So theta can also be considered as a 
    quadratic form on the Golay code C with values in the cocode C*.

    An element (v, (-1)**i) of the Parker loop Pl is represented as 
    a 13-bit integer, with bits 0,...,11 the Golay code word in
    'gcode' representation, and bit 12 the sign bit s. We call this
    representation of the Parker loo the 'ploop' representation.
    So we can convert and element of C in 'gcode' representation
    to an element pf Pl in 'ploop' representation by adjusting the
    sign in bit 12.
    
    Function ploop_theta(v) returns the element theta(v) of C* in
    'cocode' representation. Parameter v must be an element of C 
    or Pl in 'gcode' or 'ploop' representation. Similarly, function
    ploop_cocode(v1, v2) returns the value of the coycle 
    theta(v1, v2), which is 0 or 1. 

    Function mul_ploop() returns the product of two elements of 
    the Parker Loop. Function inv_ploop() returns the inverse of
    ab element of the Parker loop.


    The group AutPl of standard automorphisms of the Parker loop Pl
    ---------------------------------------------------------------

    An automorphism of the Parker loop is implemented as an array a
    of twelve 32-bit integers. The lowest 13 bits of a[i] encode the 
    image of the i-th basis vector of the Parker loop. Here the basis
    of Parker Loop corresponds to the basis of the Golay code, and 
    each basis vector has positive sign.

    The bits 13..24 of the vectors a[i] ancode a quadratic form which
    faciliates computations in AutPl. A description of the quadratic
    form is out of the scope of this document. 

    This representation of AutPl is called the 'autpl' representation.
    We only use the 'autpl' representaion for elements of AutPl.

    Function perm_to_autpl(c, p) returns an automorphism of the
    Parker loop created from a element p of Mat24 in 'perm' 
    representation and a cocode element c in 'cocode' representation.

    For m = perm_to_autpl(c, p) we can get back p and c be computing
    p = autpl_to_perm(m) and c = autpl_to_cocode(m). 

    m = cocode_to_autpl(c) returns the same result as 
    m = perm_to_autpl(c, p) with the identity permutation p.
    Note that

     perm_to_autpl(c, p) = cocode_to_autpl(c) * perm_to_autpl(0, p).   

    The automorphism perm_to_autpl(0, p) maps the basis vectors 
    (b_i, 0) of the Parker loop to some elements (b'_i, 0) of the
    Parker loop, with the Golay code automorphism b_i -> b'_i given
    by the permutation p.
    
    Function op_ploop_autpl(v, m) applies Parker loop automorphism m
    to cocode vector v and returns the result. Here all cocode 
    vectors are given in 'cocode' representation. We assume that
    Parker loop automorphisms operate by right multiplication on Pl.

    Function mul_autpl(m1, m2) returns the product m1 * m2 of the
    Parker loop automorphisms m1 and m2. Function inv_autpl(m1) 
    returns the inverse of the Parker loop automorphism m1. 


    Auxiliary functions
    ------------------
    Here is an overview of some auxiliary functions in this class. 
    They are described in the corresponding function documentation.

    bw24          bit weight of the lowest 24 bits of an integer
    lsbit24       min(24, least significant bit pos.) for an integer 
    gcode_weight  weight of a Golay code word in 'gtype' representation
    vect_type     orbit of a vector in V under the group Mat24
    vect_to_bit_list
                  given a bit vector in V, it returns the lists of
                  the positions of the 0 bits and of the 1 bits of v.
    extract_b24   extract bits from bit vector using a 24-bit mask
    spread_b24    spread bit vector according to a 24-bit mask


    Internal operation
    ------------------

    For switching from the standard representation to the internal
    representation we use 3 tables with 2**8 entries of 24 bit length.
    For switching back from internal to standard representation we use
    3 other tables of the same format. There are also tables for 
    computing the syndrome of a vector in V with respect to the Golay
    code.


    Abbreviations for functions and parameters in this class
    --------------------------------------------------------

    The following list of abbreviations used in names of functions
    allows to guess the action of most functions in this module:

    abbreviation  meaning                                   data type

    assoc         associator (in Golay code or Pl)
    autpl         autmorphism of the Parker loop Pl         uint32_t[12]
    bw24          bit weight of the lowest 24 bits of an int
    cap           intersection (of Golay code elements)
    cocode        element of Golay cocode C*                uint32_t
    cocycle       cocycle:  Pl times Pl  ->  {0,1}
    comm          commutator (in Golay code or Pl)
    gcode         element of Golay code C                   uint32_t
    inv           inversion (in Mat24, Pl, or AutPl)
    lsbit24       least significant bit of an integer, 
                  counting bits 0,...,23 only 
    m24num        number of an element of Mat24             uint32_t         
    matrix        element of Mat24 as binary matrix 
                  acting on the Golay code C                uint32_t[12] 
    mul           multiplication (in Mat24, Pl, or AutPl)
    net           Benes network for an element of Mat24     uint32_t[9]
    octad         number of an octad, i.e. a Golay code
                  element of weight 8                       uint32_t
    op            op_<vector>_<operation> means: 
                  apply <operation> to <vector>  
    op_all        apply operation to all vectors
    perm          element of Mat24 as a permutation         uint8_t[24]
    ploop         element of the Parker loop Pl             uint32_t    
    pow           power operator (in Pl)
    scalar_prod   scalar product (of Golay code and cocode)
    suboctad      suboctad, see function suboctad_to_cocode
    syndrome      syndrome (after decoding Golay code)      uint32_t 
    theta         cocycle theta: Pl -> C^* in Parker loop
    to            <x>_to_<y> means: return representation <y>
                  of an object given in representation <x>
    vect          vector in V = GF(2)**24                   uint32_t
    vintern       vector in V, in internal representation   uint32_t
    

    Parameters of functions are either integers or arrays of integers.
    Here all integer types are unsigned and of fixed length, such as
    uint8_t, uint16_t or uint32_t. 

    The type of a parameter is given by a single letter in the name
    of the parameter:

    name  meaning                                           type
    a     array specified in documentation of function      unspecified 
    c     Golay cocode element, represented as 'cocode'     uint32_t
    m     permutation in Mat24 or automorphism of Pl
          represented as a bit matrix                       uint32_t[12]
    p     permutation in Mat24 represented as 'perm'        uint8_t[24] 
    u_<x> unsigned integer, e.g.                            unspecified
           u_exp:   integer denoting an exponent
           u_m24:   number of a permutation in Mat24
           u_octad: number of octad, 0 < u_octad < 259
           u_width: integer denoting a bit width
    v     vector in V, Golay code C or Parker loop Pl 
          represented as vect, vintern, gcode or ploop      uint32_t           

    Integer input parameters have shape u_<x>, e.g. u_m24, u_exp.
    An integer computed by a function is returned as return value.
    Input array parameters have a digit as a suffix, e.g.: v1, v2, m1.
    Output array parameters have the suffix _out, e.g.: p_out.
    Input/output array parameters  have the suffix _io, e.g.: m_io.

    References
    ----------

    See file references.txt
    """
    ## Create tables for faster computetion
    MAT24_ORDER =  244823040 
    matrix_to_perm_ = MatrixToPerm(Mat24Tables.basis[12:])
    heptad_completer = HeptadCompleter(Mat24Tables)
    theta_table = make_augmented_theta_table()
    autpl_qf_table = make_autpl_qf_table(theta_table)
    verbose = False

    # recip_basis(i & 31) shall not fail in C 
    recip_basis_c = numpy.append(Mat24Tables.recip_basis.copy(), [0]*8)

    ## Collect tables and coding functions for generating C code
    tables = {
           "Mat24_enc_table0"    : Mat24Tables.enc_table0,
           "Mat24_enc_table1"    : Mat24Tables.enc_table1,
           "Mat24_enc_table2"    : Mat24Tables.enc_table2,
           "Mat24_dec_table0"    : Mat24Tables.dec_table0,
           "Mat24_dec_table1"    : Mat24Tables.dec_table1,
           "Mat24_dec_table2"    : Mat24Tables.dec_table2,
           "Mat24_basis"         : Mat24Tables.basis,
           "Mat24_recip_basis"   : recip_basis_c,
           "Mat24_syndrome_table"    : Mat24Tables.syndrome_table,
           "Mat24_oct_enc_table"     : Mat24Tables.oct_enc_table,
           "Mat24_oct_enc_offset"    : Mat24Tables.oct_enc_offset,
           "Mat24_oct_dec_table"     : Mat24Tables.oct_dec_table, 
           "Mat24_theta_table"       : theta_table,
           "Mat24_autpl_qf_table"    : autpl_qf_table,
    }

    directives = {
           "BITMATMUL" :  UserDirective(generate_c_bitmatmul, "psss"),
           "BITVMULTRANSP" :  UserDirective(generate_c_mul_transp, "ss"),
    }

    _subgenerators = [
            Lsbit24Function, 
            matrix_to_perm_,
            heptad_completer
    ]

    for gen in _subgenerators:
        tables.update(gen.tables())
        directives.update(gen.directives())


    @classmethod
    def str_basis(cls, with_reciprocal_basis = False):
        def show(text, basis):
            s = text + "\n"
            for i,x in enumerate(basis): 
                if i % 6 == 0: s += "    "
                s +=  "%06x " % x
                if i % 6 == 5: s += "\n"
            return s
        s = "We list the basis vectors of the Golay code and of its cocode."
        s += """

Basis vectors have indices 0,...,11. Each basis vector is displayed 
as a hexadecimal number with bit i (of valence 2**i) corresponding  
to component i of the basis vector in GF(2)^24 for i = 0,...,23.
Golay cocode vectors are to be understood modulo the Golay code. 

"""
        s +=  show( "Golay cocode basis", cls.basis[:12] )
        s +=  show( "Golay code basis" , cls.basis[12:] )
        if with_reciprocal_basis:
             s += show( "Reciprocal basis", cls.recip_basis )
        return s
      
    @classmethod
    def show_basis(cls, with_reciprocal_basis = True):
        print(cls.str_basis(with_reciprocal_basis))

    ###########################################################################
    # Parker Loop
    ###########################################################################

    @classmethod
    def ploop_theta(cls, v1):
        """Return the theta function for the Parker loop.

        theta is a quadratic from from the Golay code C to the cocode C*.
        Here parameter v1 of function theta is represented as a Golay code 
        word. The result of the function is represented as a Golay cocode 
        word. The cocycle of the Parker loop is given by:

             cocycle(v1,v2) =   scalar(theta(v1), v2)
        
        with  scalar(.,.) the scalar product.
        """
        return cls.theta_table[(v1 >> 1) & 0x7ff] & 0xfff

    @classmethod
    def ploop_cocycle(cls, v1, v2):
        """Return the cocycle of the Parker loop.

        Then the Parker Loop product is given by

             v1 (*) v2  =  v1 ^ v2 * (-1)**cocycle(v1, v2) . 
        """
        s = cls.ploop_theta(v1) & v2 & 0xfff
        s ^= s >> 6
        s ^= s >> 3
        s = 0x96 >> (s & 7)
        return s & 1 


    @classmethod
    def mul_ploop(cls, v1, v2):
        """Return the Parker loop product v1 (*) v2

        Here v1 and v2 are integers coded as follows:
        bit 0,...,11:   representation as Golay code word
        bit 12:         Parker loop sign
        otther bits:    ignored
        """
        return v1 ^ v2 ^ (cls.ploop_cocycle(v1, v2) << 12)


    @classmethod
    def pow_ploop(cls, v1, u_exp):
        """Return power v1 ** u_exp of the Parker loop element v1."""
        return (v1 & -(u_exp & 1)) ^ (
             cls.theta_table[(v1 >> 1) & 0x7ff] & ((u_exp & 2) << 11) )


    @classmethod
    def ploop_comm(cls, v1, v2):
        """Return commutator of Golay code word v1 and v2.

        This is 0 if the intersection of the vectors v1 and v2 has
        bit weight 0 mod 4 and 1 is that intersection has bit weight 
        2 mod 4. v1 and v2 are in 'gvect' or 'ploop' representation.
        """
        cap = cls.gcode_to_vect(v1) & cls.gcode_to_vect(v2)
        return (cls.bw24(cap) >> 1) & 1

    @classmethod
    def ploop_assoc(cls, v1, v2, v3):
        """Return associator of Golay code words v1, v2 and v3
 
        This the parity of the intersection of the vectors v1, v2 and 
        v3.  v1, v2 and v3 are in 'gvect' or 'ploop' representation.
        """
        assoc = (cls.gcode_to_vect(v1) & cls.gcode_to_vect(v2)
                & cls.gcode_to_vect(v3))
        return cls.bw24(assoc) & 1


    @classmethod
    def ploop_cap(cls, v1, v2):
        """Return intersection of two Golay code words as cocode word.

        v1 and v2 are in 'gvect' or 'ploop' representation, the result
        is returned in 'cocode' representation.
        """
        return ( cls.theta_table[((v1 ^ v2) >> 1) & 0x7ff] 
                  ^  cls.theta_table[(v1 >> 1) & 0x7ff] 
                  ^  cls.theta_table[(v2 >> 1) & 0x7ff]  )


    ###########################################################################
    # Mathieu group M24: conversion between representations
    ###########################################################################
    


    @classmethod
    def perm_complete_heptad(cls, p_io):
        """Complete a permutation p given by p_io to an element of Mat24.

        p must be a list of length 24. Entries p[i], i = 0,1,2,3,4,5,8.
        must make up a valid umbral heptad, i.e. a heptad not contained 
        in an octad. p[0],...,p[5] must be contained in an octad, p[8]
        must not be contained in that octad. The other entries of 
        input p are ignored.

        It can be shown that such a permutation p can be completed to 
        a unique element of Mat24.

        The function returns 0 in case of success and a nonzero value
        otherwise. In case of success, p_io is completed to an element of
        the Mathieu group Mat24. 
        """
        return cls.heptad_completer.compute(p_io)


    @classmethod
    def perm_check(cls, p1):
        """Check if permutation p1 is in in the Mathieu group Mat24.

        The function returns zero iff this is the case.
        """
        p2 = list(p1[:])
        if  cls.perm_complete_heptad(p2):
             return -1
        return  p2 != list(p1[:])

    @classmethod
    def perm_from_heptads(cls, h1, h2):
        """Try to find a permutation p that maps heptad h1 to h2

        h1 and h2 must be lists of length 7 defining two umbral heptads,
        i.e. heptads not contained in an octad. If a permutation p in
        the Mathieu group Mat24 that maps h1 to h2 exists, it is unique. 

        The function returns p if such a p exists an is unique and
        it returns None otherwise.
        """
        return cls.heptad_completer.perm_from_heptads(h1, h2)

    @classmethod
    def m24num_to_perm(cls, u_m24):
        """Return permutation with number u_m24 in the Mathieu group Mat24.

        The inverse of this function is member function Mat24_perm_to_int()
        This is just a short and convenient way to number elements of Mat24.
        Input u_m24 = 0 gives the identity permutation.

        0 <= u_m24 < 244823040 = order(Mat24) must hold.

        For internal operation see
        mat24Heptad.HeptadCompleter.int_to_perm
        """
        p = cls.heptad_completer.int_to_perm(u_m24)
        return p

    @classmethod
    def perm_to_m24num(cls, p1):
        """Convert permutation p1 in the Mathieu group Mat24 to an integer.

        This reverses member function int_to_perm(). The input permutation
        is not checked.
        """
        return cls.heptad_completer.perm_to_int(p1[:])
            

    @classmethod
    def perm_to_matrix(cls, p1):
        """Convert a permutation p1 in  Mat24 to a matrix.

        The matrix is a 12 x 12 bit matrix acting on the Golay code
        vectors by right multiplication.

        Permutation p is not checked to be a member of the Mathieu group.
        """
        mat = [cls.recip_basis[p1[i]] >> 12 for i in range(24)]
        return bit_mat_mul(cls.basis[12:], mat)


    @classmethod
    def matrix_to_perm(cls, m1):
        """Convert element of Mat24 from matrix to permutation.

        The matrix m1 is a 12 x 12 bit matrix acting on the Golay code
        vectors by right multiplication. The matrix is not checked.
    
        The permutation is represented as a list.
        """
        basis = [cls.gcode_to_vect(v & 0xfff) for v in m1]
        return cls.matrix_to_perm_.compute(basis)


    #####################################################################
    # Mathieu group M24: operation of group elements
    #####################################################################

    @classmethod
    def op_vect_perm(cls, v1, p1):
        """Apply a permutation p1 to a vector v1 in GF(2)**24

        Here p1 is the permutation that maps i to p1[i]  for i=0,...,23.
        """      
        return sum([ ((v1 >> i) & 1) << p1[i] for i in range(24) ])


    @classmethod
    def op_gcode_matrix(cls, v1, m1):
        """Apply the 12 x 12 bit matrix m1 to a Golay code vector

        Here application means right multiplication v1 * m1. The code
        vector v1 is given in 'gcode' representation.
        """    
        return reduce(
                __xor__, [x for i,x in enumerate(m1) if v1 & (1<<i)], 0 )


    @classmethod
    def op_gcode_perm(cls, v1, p1):
        """Apply a permutation p1 to a Golay code vector v1
  
        Here p1 is the permutation that maps i to p1[i] for i=0,...,23,
        representing an element of the Mathieu group M24.
        
        Golay code vector v1 is given in gcode representation.
        """    
        v1 = cls.gcode_to_vect(v1)
        v1 = cls.op_vect_perm(v1, p1)
        return  cls.vect_to_vintern(v1) >> 12


    @classmethod
    def op_cocode_perm(cls, c1, p1):
        """Apply a permutation p to a Golay cocode vector c1

        Here p1 is the permutation that maps i to p1[i]  for i=0,...,23,
        representing an element of the Mathieu group M24.

        Golay cocode vector c1 is given in cocode representation.
        """
        v = c1
        y = - ((v + 1) & 1)                  # y = 0 if v is odd else -11
        v ^= cls.recip_basis[0] & y          # make v odd
        res = cls.recip_basis[p1[0]] & y     # .. and adjust result
        syn = cls.syndrome_table[ (v >> 1) & 0x7ff ]  # get syndrome
        res ^=  cls.recip_basis[p1[syn & 31]] # do 1st syndrome vector
        syn = (syn >> 5) & 0x3ff             # mask out 2nd, 3rd synd. vector
        syn &=  ((syn + 0x100) >> 10) - 1    # kill syn if v >= 24*32
        res ^=  cls.recip_basis[p1[syn & 31]] ^ cls.recip_basis[p1[syn >> 5]]
                                             #  do 1st and 2nd syndrome vector
        return res & 0xfff


    @classmethod
    def mul_perm(cls, p1, p2):
        """Return p1 * p2, with p1, p2 in the Mathieu group M24

        p1, p2 are represented as permutations
        """
        return [ p2[p1[i]] for i in range(24) ]
 
    @classmethod
    def inv_perm(cls, p1):
        """Return inverse of p1, with p1 in the Mathieu group M24

        p1, is represented as a permutations
        """
        l = [None] * 24
        for i, x in enumerate(p1): l[x] =  i
        return l

             
    ###########################################################################
    # Automorphisms of the Parker Loop
    ###########################################################################



    @classmethod
    def autpl_set_qform(cls, m_io):
        """Recompute quadratic form on Parker loop automorphism m_io

        This functions augments the Parker loop automorphism m_io by
        a quadratic form qf. The form qf simplifies the application
        of m_io to Parker loop elements and also the multiplication
        of  Parker loop automorphisms. The form qf is stored in bits
        13,...,24 of the entries of m_io. 

        At present the documatation of the mathematical background
        behind that form is out of the scope of the module.
        """
        m_io[0] &= 0x1fff 
        m_io[1] &= 0x1fff
        for i in range(2,12):
            qq = cls.ploop_theta(1<<i) & ((1 << i) - 1)
            s = cls.ploop_theta(m_io[i])
            for j in range(1,i):
                t =  s & m_io[j] 
                t ^= t >> 6
                t ^= t >> 3
                qq ^= ((0x96 >> (t & 7)) & 1) << j
            m_io[i] = (m_io[i] & 0x1ffff) ^ (qq << 13)
        return m_io

    @classmethod
    def perm_to_autpl(cls, c1, p1):
        """Combine Mat24 and cocode element to Parker loop automorphism

        Given an element p1 of the Mathieu group Mat24 (in permutation 
        representation) and a Golay cocode element c1 (in cocode 
        representation), the function returns a Parker loop automorphism 
        m as a 12 x (12+13) matrix, i.e. in 'autpl' representation.
        m contains the 12 images of the basis vectors of the Parker loop
        and a quadratic form for simplfying its operation on Pl.
        """
        m = cls.perm_to_matrix(p1)
        for i in range(12):
            m[i] ^= ((c1 >> i) & 1) << 12
        return  cls.autpl_set_qform(m)

    @classmethod
    def cocode_to_autpl(cls,  c1):
        """Convert cocode element c1 to Parker loop automorphism

        Same as perm_to_autpl(c1, p), with p the identity permutation
        """
        return [(1 << i) + (((c1 >> i) & 1) << 12) for i in range(12)]

    autpl_to_perm = matrix_to_perm 

    @classmethod
    def autpl_to_cocode(cls, m1):
        """Extract cocode vector c from Parker loop automorphism m1

        Then m1 = perm_to_autpl(c, p), where p is the permutation
        obtained by calling autpl_to_perm(m1).

        Note that m1 = cocode_to_autpl(c) *  perm_to_autpl(0, p).
        """
        return sum( ((m1[i] >> 12) & 1) << i for i in range(12) )
       

    @classmethod
    def op_ploop_autpl(cls, v1, m1):
        """Apply Parker loop automorphism m1 to Parker Loop element v1

        Here m1 is a Parker loop autmorphism (in autpl representation)
        and v1 is a cocde vector (in cocode representation).
        The function returns the resluting cocode vector  v1 * m1.
        """
        v = v1
        t = v & 0x1000
        for i in range(12):
            t ^= m1[i] & -((v >> i) & 1) 
        v = (t >> 13) & v 
        v ^= v >> 6
        v ^= v >> 3 
        v = (0x96 >> (v & 7)) & 1
        return (t & 0x1fff) ^ (v << 12)

    @classmethod
    def mul_autpl(cls, m1, m2): 
        """Return product m1 * m2 of Parker loop automorphisms m1, m2

        Here  m1, m2  and the result is in 'autpl' representation.
        """
        m3 = [cls.op_ploop_autpl(m1[i], m2) for i in range(12)]
        return cls.autpl_set_qform(m3)

    @classmethod
    def inv_autpl(cls, m1): 
        """Return inverse of Parker loop automorphism m1

        Here m1 and the inverse of it is in 'autpl' representation.
        """
        p = cls.matrix_to_perm(m1)
        p = cls.inv_perm(p)
        mi = cls.perm_to_matrix(p)
        m0 = cls.mul_autpl(mi, m1)
        for i in range(12):
            assert mi[i] & -0x1000 == 0, map(hex, mi)
            assert m0[i] & 0xfff == 1 << i
            mi[i] ^= m0[i] & 0x1000
        return cls.autpl_set_qform(mi)

    @classmethod
    def perm_to_iautpl(cls, c1, p1):
        """Return (inverse of p1, inverse of element of perm_to_autpl(c1, p1))
     
        This will be optimized in the fast version.  
        """
        return cls.inv_perm(p1), cls.inv_autpl(cls.perm_to_autpl(c1, p1)) 


