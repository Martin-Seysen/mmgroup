"""Create a table for the Parker Loop cocycle theta


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
from numpy import array, zeros, uint8, uint16, uint32


from mmgroup.bitfunctions import bw24, v2

from mmgroup.dev.mat24.mat24tables import Mat24Tables


from mmgroup.dev.mat24.mat24aux import split_golay_codevector





def theta_to_basis_vector(v):
    """returns suitable theta(v) value for Golay code basis vector v

    Here theta(v) is a cocode word, which is (loosely speeking)
    the cocycle of the Parker loop.

    The reader may take this as black magic and study the checking
    function theta_basis_vectors() instead.
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
    """yet to be commented!!!!!!!!!!!!!!!!!"""
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
    """Creates the theta function for the Parker Loop
    theta() is a function from the Golay code C to the Golay cocode C*.

    Multiplication '(*)' in the Parker loop is defined by:

           a (*) b = (a ^ b) * (-1) ** scalar(theta(a),b) ,

    Where '^' is the vector addition in GF(2)**12, and  scalar(.,.)
    is the scalar product defined on C x C*.     

    More specifically, theta is a quadratic form from C onto C*. The 
    basis cc_i, i=0,...,11 of the cocode defined here is reciprocal to 
    the basis c_i, i=0,...,11 of the cocode defined here. 
    Define   theta(i)[j] = theta(i) * c_j.   Then 

            theta(i)  = sum[ theta(i)[j] for j in range(12) ]. 
 
    To define Theta as a quadratic form, it suffices to define theta 
    on the basis vectors, and to define the symmetric bilinear form B 
    (from C x C onto C*) associated with theta on all pairs of basis 
    vectors. We have:
 
             B(x,y) = theta(x+y) + theta(x) + theta(y)   .        

    We define: 
            theta(c_i)[j] = bw(c_i & c_j) >> 1  (mod 2)    ; j < i
            theta(c_i)[j] = bw(c_i) >> 2        (mod 2)    ; j == i
            theta(c_i)[j] = 0                              ; j > i
            B(c_i, c_j)   = 0                              ; i = j,
            B(c_i, c_j)   = c_i & c_j           (mod 2)    ; i != j 
  
    Here c_i & c_j is the  intersection (i. e. bitwise product in GF(2)) 
    of c_i and c_j, interpreted as a cocode word in C*.
    """
    assert Mat24Tables.gcode_to_vect(1) == 0xffffff
    tab = numpy.zeros(0x1000, dtype = uint16)
    for i in range(12):
        t =  0
        for j in range(i):
            v = Mat24Tables.gcode_to_vect(1<<i) 
            v &= Mat24Tables.gcode_to_vect(1<<j)
            w = (bw24(v) >> 1) & 1
            t += w << j
        v = Mat24Tables.gcode_to_vect(1<<i)
        w = (bw24(v) >> 2) & 1
        t +=  w << i
        tab[1 << i] = t 

    # Exprimental new cocycle !!!!!!!!!!!
    thetas = theta_basis_vectors()
    for i in range(12):
         tab[1 << i] = Mat24Tables.vect_to_cocode(thetas[i])
    # End of exprimental new cocycle !!!!!!!!!!!

    for i in range(0x1000):
        if (i & (i-1)) != 0:    # i.e. if i is not 0 or a power of 2
            i0 = 1 << v2(i) 
            i1 = i ^ i0
            cap = Mat24Tables.gcode_to_vect(i0) 
            cap &= Mat24Tables.gcode_to_vect(i1)
            cc = Mat24Tables.vect_to_cocode(cap)
            tab[i] = tab[i0] ^ tab[i1] ^ cc
    theta_table =  numpy.zeros(0x800, dtype = uint16)

    for i in range(0x800):
        theta_table[i] = tab[2*i]
        assert tab[2*i] == tab[2*i+1]
    del tab
    return theta_table




def augment_theta_table(theta_table):
    """store bitweight of Golay code word i in table[i/2], bits 14..12"""
    for i in range(0x800):
        w = bw24(Mat24Tables.gcode_to_vect(2*i))
        assert w in [0,8,12,16]
        theta_table[i] |=  (w >> 2) << 12
    return  theta_table


def make_augmented_theta_table():
    theta_table = make_theta_table()
    return augment_theta_table(theta_table)
            


def make_autpl_qf_table(theta_table, bitwidth=32):
    mask_table = []
    qf_table =  []
    for i in range(2,12):
        qf_table.append((theta_table[(1 << (i-1))] & 0xfff) >> 1)
        mask_table.append ((1 << (i-1)) - 1)
    assert  bitwidth==32
    tab = numpy.zeros(12, dtype = uint32)
    for i in range(0, 10, 2):
        tab[i] = qf_table[i] + (qf_table[i+1] << 16) 
        tab[i+1] = mask_table[i] + (mask_table[i+1] << 16) 
    return tab
    




