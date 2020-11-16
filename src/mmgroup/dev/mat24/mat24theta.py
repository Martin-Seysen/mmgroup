"""The module computes a table for the Parker Loop cocycle ``theta``.

"""


from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


from random import randint
import types
import sys
import re
import os
from operator import __or__, __xor__

import numpy
from numpy import array, zeros, uint8, uint16, uint32, uint64


from mmgroup.bitfunctions import bw24, v2

from mmgroup.dev.mat24.mat24tables import Mat24Tables


from mmgroup.dev.mat24.mat24aux import split_golay_codevector





def theta_to_basis_vector(v):
    """Return suitable theta(v) value for Golay code basis vector v

    Here theta(v) is a cocode word, which is (loosely speeking)
    the cocycle of the Parker loop. Both, v and theta(v) are 
    vectors given in 'vect' representation.

    We assume that vector v is either 'grey' or 'colored' in the
    sense of [Seys19]. If v is a 'grey' vector, we return the
    cocycle theta(v) as given by Lemma 3.9. For each 'colored'
    vasis vector v, we assume that theta(v) is 'grey'; then the
    'grey' part of theta(v) is uniquely defined by Lemma 3.9.  
    """
    bw, col = split_golay_codevector(v, check=Mat24Tables)
    if col and not bw:
        assert bw24(col) == 8
        return ((col >> 1) | (col >> 2) | (col >> 3)) & 0x111111
    assert col == 0 
    bw &= 0x111111
    weight = bw24(bw)
    if weight < 2 or weight > 4:   return 0
    if weight == 2: return bw ^ 0x111111
    if weight == 4: return bw
    if weight == 3: return 0x111111  
    
    
def theta_basis_vectors(verbose = 0):
    """Return cocycle of the standard basis of the Golay code.

    Here the standard basis of the Golay code is given by
    Mat24Tables.basis[12:]. We return list containing theta(b_i) 
    for the basis vectors b_0,...,b_11. Here all values 
    theta(b_i) are computed by function theta_to_basis_vector(v).

    We also check that these values theta(b_i) are consisten.
    """
    basis = Mat24Tables.basis[12:]
    thetas = [theta_to_basis_vector(v) for v in basis] 
    for i in range(12):
        for j in range(i):
            cocycle = bw24(thetas[i] & basis[j]) 
            cocycle = (cocycle  ^ (bw24(thetas[j] & basis[i]))) & 1
            expected = (bw24(basis[i] & basis[j]) >> 1) & 1
            assert cocycle == expected & 1, (i, j, cocycle, expected)
        cocycle = bw24(thetas[i] & basis[i]) & 1
        assert cocycle == (bw24(basis[i]) >> 2) & 1
    if verbose:
        print( "Cocycle test for basis of Parker loop passed" )
    return thetas



def make_theta_table():
    r"""Creates the theta function for the Parker Loop
    theta() is a function from the Golay code C to the Golay cocode C*.

    Multiplication '(*)' in the Parker loop is defined by:

           a (*) b = (a ^ b) * (-1) ** scalar(theta(a),b) ,

    Where '^' is the vector addition in GF(2)**12, and  scalar(.,.)
    is the scalar product defined on C x C*.     

    More specifically, theta is a quadratic form from C onto C*. The 
    basis cc_i, i=0,...,11 of the cocode defined here is reciprocal to 
    the basis c_i, i=0,...,11 of the cocode defined here. 
    Define   theta(i)[j] = theta(i) * c_j.   Then 

            theta(i)  = sum[1 << theta(i)[j] for j in range(12) ]. 
 
    To define theta as a quadratic form, it suffices to define theta 
    on the basis vectors. We define theta on the basis vectors as
    given by function theta_to_basis_vector(v).

    There is a symmetric bilinear form B : from C x C onto C*
    associated with theta on all pairs of basis vectors satisfying:

             B(x,y) = theta(x+y) + theta(x) + theta(y)   .        
             B(x,y) = x :math:`\cap` Y

    Here :math:`\cap means intersection (i.e. bitwise 'and' in 
    GF(2)**24) of two Golay code words. This allows us to compute 
    theta(x), for all Golay code words, if theta is given on a 
    basis of the Golay code.
    """
    assert Mat24Tables.gcode_to_vect(0x800) == 0xffffff
    theta_table = numpy.zeros(0x800, dtype = uint16)

    thetas = theta_basis_vectors()
    for i in range(11):
         theta_table[1 << i] = Mat24Tables.vect_to_cocode(thetas[i])
 
    for i in range(0x800):
        if (i & (i-1)) != 0:    # i.e. if i is not 0 or a power of 2
            i0 = 1 << v2(i) 
            i1 = i ^ i0
            cap = Mat24Tables.gcode_to_vect(i0) 
            cap &= Mat24Tables.gcode_to_vect(i1)
            cc = Mat24Tables.vect_to_cocode(cap)
            theta_table[i] = theta_table[i0] ^ theta_table[i1] ^ cc
    return theta_table




def augment_theta_table(theta_table):
    """store bitweight of Golay code word i in table[i/2], bits 14..12"""
    for i in range(0x800):
        w = bw24(Mat24Tables.gcode_to_vect(i))
        assert w in [0,8,12,16]
        theta_table[i] |=  (w >> 2) << 12
    return  theta_table


def make_augmented_theta_table():
    theta_table = make_theta_table()
    return augment_theta_table(theta_table)
            



    

def make_autpl_qf_table(theta_table, bitwidth=32):
    qf_table = [theta_table[(1 << i)] & 0x7ff for i in range(11)]
    mask_table =  [(1 << i) - 1 for i in range(11)]
    if  bitwidth==32:
        tab = numpy.zeros(10, dtype = uint32)
        for i in range(0, 10, 2):
            tab[i] = qf_table[i+1] + (qf_table[i+2] << 16) 
            tab[i+1] = mask_table[i+1] + (mask_table[i+2] << 16) 
        return tab
    if bitwidth == 64:
        t = [sum((int(qf_table[5*i+j+1]) << (12*j) for j in range(5)))
            for i in range(2)]
        return ["0x%xULL" % x for x in t]
    raise ValueError("Illegal  bitwidth for table 'qf_table'")
        



