r"""Module ``mat24_ref`` contains the class ``Mat24.`` 

Class ``Mat24`` in module ``mmgroup.dev.mat24.mat24_ref`` is the 
table-providing class used by the code generator to generate the C 
file ``mat24_functions.c``. That C file contains the functionality
of the python extension ``mmgroup.mat24``.
Class ``Mat24`` may also be used as a table-providing class for 
generating other C files. 

Class ``Mat24`` also exports the same functions as the ``mmgroup.mat24``
extension as class methods. So it can be used for testing that extension 
and also as a pure-python substitute for that extension. 

Class ``Mat24`` is based on class ``mat24tables.Mat24Tables``. It also
uses classes, functions and tables from the following modules 
in ``mmgroup.dev.mat24``:

 *  ``make_addition_table``, ``make_mul_transp``,
    ``mat24aux``, ``mat24heptad``, ``mat24theta``

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
from mmgroup.bitfunctions import pivot_binary_low

from mmgroup.dev.mat24.mat24tables import Mat24Tables

from mmgroup.dev.mat24.mat24aux import Lsbit24Function, MatrixToPerm

from mmgroup.dev.mat24.mat24heptad import HeptadCompleter

from mmgroup.dev.mat24.mat24theta import make_augmented_theta_table
from mmgroup.dev.mat24.mat24theta import make_autpl_qf_table

from mmgroup.dev.mat24.make_addition_table import BitMatrixMulFix

from mmgroup.dev.mat24.make_mul_transp import BitMatrixMulTransp

from mmgroup.dev.mat24.mat24_doc import __doc__ as Mat24__doc__


from mmgroup.generate_c import UserDirective, UserFormat



TABLE_NAMES = None

###########################################################################
###########################################################################
# The Golay code (a pure python implementation)
###########################################################################
###########################################################################


class Mat24(Mat24Tables): 
    r"""Provide functions for the Mathieu group Mat24 and the Parker loop.

    These functions are exported as class methods. They are considered as
    pure python substitutes for the functions in the ``mmgroup.mat24``
    extension.

    So they can be used in an early stage of the build process where that
    extension is not yet available. The are also used for testing the
    ``mmgroup.mat24`` package.

    See module ``mmgroup.dev.mat24.mat24_doc`` for an overview of the
    functionality of the ``mmgroup.mat24`` extension. 

    This class is also a table-providing class used by the code
    generator for generating C code containing the functionality
    of the  ``mmgroup.mat24`` extension. For details, see section 
    *How the code generator is used* in 
    *The mmgroup guide for developpers*. 
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
           "Mat24_doc"               : Mat24__doc__ 
    }

    directives = {}

    _subgenerators = [
            Lsbit24Function, 
            matrix_to_perm_,
            heptad_completer,
            BitMatrixMulFix(),
            BitMatrixMulTransp(),
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
        return cls.theta_table[v1 & 0x7ff] & 0xfff

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
             cls.theta_table[v1 & 0x7ff] & ((u_exp & 2) << 11) )


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
        return ( cls.theta_table[(v1 ^ v2) & 0x7ff] 
                  ^  cls.theta_table[v1 & 0x7ff] 
                  ^  cls.theta_table[v2 & 0x7ff]  ) & 0xfff

    @classmethod
    def ploop_solve(cls, a):
        """Return cocode element that kills signs of Parker loop elements
 
        Here 'a' is an array of Parker loop elements. The function tries 
        to find a cocode element that makes all these Parker loop 
        elements positive, when operating on them as a diagonal 
        automorphism. The function returns the least cocode element in 
        lexical order satisfying that condition. For that order we 
        assume that lower bits have higher valence.
        If no such cocode element exists, ValueError is raised.
        """
        a1 = [x & 0x1fff for x in a]
        basis, columns = pivot_binary_low(a1)
        res = 0
        for b, col in zip(basis, columns):
            res |= ((b >> 12) & 1) << col 

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
        y = - (((v >> 11) + 1) & 1)          # y = 0 if v is odd else -1
        v ^= cls.recip_basis[0] & y         # make v odd
        res = cls.recip_basis[p1[0]] & y     # .. and adjust result
        syn = cls.syndrome_table[v & 0x7ff ]  # get syndrome
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
        for i in range(12):
            qq = cls.ploop_theta(1 << i) & ((1 << i) - 1)
            s = cls.ploop_theta(m_io[i])
            for j in range(i):
                t =  s & m_io[j] 
                t ^= t >> 6
                t ^= t >> 3
                qq ^= ((0x96 >> (t & 7)) & 1) << j
            m_io[i] = (m_io[i] & 0x1fff) ^ (qq << 13)
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


