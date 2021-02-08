import sys
import os
import warnings

from random import randint
from collections import OrderedDict

import numpy as np
import pytest


from mmgroup.clifford12 import leech3matrix_echelon
from mmgroup.clifford12 import leech3matrix_sub_diag
from mmgroup.clifford12 import leech3matrix_load_3
from mmgroup.clifford12 import leech3matrix_load_15

def rand_matrix_mod_n(n, dtype = np.int64):
    """Create random 24 x 24 matrix of intgers modulo n"""
    m =  np.zeros((24, 24), dtype = dtype)
    for i in range(24):
        for j in range(24):
            m[i,j] = randint(0,n-1)
    return m

def print_mod3(m):
    """Print 24 x 24 matrix of intgers modulo 3"""
    m1 = m % 3
    for i in range(24):
        print("".join(map(str, m1[i] % 3)))
    print("")



def matrix_to_rep_modp(m, p):
    """Convert 24 x 24 matrix to monster rep, part A, modulo p

    The function returns the first block of a vector of the monster
    rep modulo p. Ths block corresponds to the entries of the
    vector with tag 'A'. 
    """
    assert p in (3, 15)
    if p == 3:
        m0 = (m - 1) % 3 + 1
        b = np.zeros(24, dtype = np.uint64)
        for i in range(24):
            b[i] = sum(m0[i,j] << (2*j) for j in range(24))
    elif p == 15:
        m0 = m % 3 + 3 * rand_matrix_mod_n(5)
        b = np.zeros(2*24, dtype = np.uint64)
        for i in range(24):
            r = sum((int(m[i,j]) % 3 ) << (4*j) for j in range(24))
            b[2*i] = r & 0xffffffffffffffff
            b[2*i + 1] = r >> 64
    else:
        err = "Module Leech3matrix does not support p=%s"
        raise ValueError(err, str(p))
    return b


load_functions = {
     3: leech3matrix_load_3,
    15: leech3matrix_load_15,
}

def from_array(m, load_mode = None):
    """Load 24 x 24 matrix 'm' to array of type 'leech3matrix'

    If load_mode is 3 or 5 then matrix is converted to a vector of
    of the monter rep using function ``matrix_to_rep_modp(m, p)``
    with ``p = load_mode``.

    By default we convert the matrix directly.

    The function returnes an numpy array that can be processed
    with the function in file ``leech3matrix.c``.
    """
    a = np.zeros(72, dtype = np.uint64)
    if load_mode in (3, 15):
        b = matrix_to_rep_modp(m, load_mode)
        f = load_functions[load_mode]
        f(b, a)
    else:
        for i in range(24):
            s = int(sum(int(m[i,j]) % 3 << 4*j for j in range(24)))
            a[3*i] = s & 0xffffffffffffffff
            a[3*i+1] = (s >> 64) & 0xffffffff
    return a


def as_array(a, high = 0, dtype = np.int64):
    """Convert  array of type 'leech3matrix' to 24 x 24 matrix"""
    m =  np.zeros((24, 24), dtype = dtype)
    jj = [96 * bool(high) + 4 * j for j in range(24)]
    jj = [(j >> 6, j & 0x3f) for j in jj]
    for i in range(24):
        for j, (ofs, shift) in enumerate(jj):
             m[i,j] = ((int(a[3*i + ofs]) >> shift) & 3) % 3
    return m
 
#######################################################################
# Test function monomial_to_word leech3matrix_echelon
#######################################################################


@pytest.mark.qstate
def test_leech3matrix_echelon(verbose = 0):
    unit = np.eye(24)
    load_modes = [None, 3, 15]
    for i in range(20):
        load_mode = load_modes[i % len(load_modes)]
        m = rand_matrix_mod_n(3, dtype = np.int64)
        a =  from_array(m, load_mode)
        leech3matrix_sub_diag(a, 2, 24)
        if verbose:
            print("load_mode =", load_mode)
            print_mod3(as_array(a,0))
            print_mod3(as_array(a,1))
        leech3matrix_echelon(a)
        if verbose:
            print_mod3(as_array(a,0))
            print_mod3(as_array(a,1))
        mi = as_array(a, high=1)
        prod = as_array(a, high=0)
        if verbose:
            print_mod3(prod)
            print_mod3(mi @ m)
            #print((mi @ m - prod) % 3)
        assert ((mi @ m - prod) % 3 == 0).all()







        
         
         