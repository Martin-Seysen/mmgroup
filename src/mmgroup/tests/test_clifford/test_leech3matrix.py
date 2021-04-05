import sys
import os
import warnings

from random import randint
from collections import OrderedDict

import numpy as np
import pytest

from mmgroup import MM, MMSpace
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.clifford12 import leech3matrix_echelon
from mmgroup.clifford12 import leech3matrix_sub_diag
from mmgroup.clifford12 import leech3matrix_load
from mmgroup.clifford12 import leech3matrix_kernel_vector
from mmgroup.clifford12 import leech3matrix_compress
from mmgroup.clifford12 import xsp2co1_from_vect_mod3
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import leech2matrix_solve_eqn
from mmgroup.clifford12 import leech3matrix_watermark
from mmgroup.clifford12 import leech3matrix_watermark_perm_num
from mmgroup.mat24 import MAT24_ORDER
from mmgroup.mm15 import op_word_tag_A

def rand_matrix_mod_n(n, dtype = np.int64):
    """Create random 24 x 24 matrix of intgers modulo n"""
    m =  np.zeros((24, 24), dtype = dtype)
    for i in range(24):
        for j in range(24):
            m[i,j] = randint(0,n-1)
    return m

def print_mod3(m):
    """Print 24 x 24 matrix of integers modulo 3"""
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



def from_array(m, load_mode = None):
    """Load 24 x 24 matrix 'm' to array of type 'leech3matrix'

    If load_mode is 3 or 15 then matrix is converted to a vector of
    of the monter rep using function ``matrix_to_rep_modp(m, p)``
    with ``p = load_mode``.

    By default we convert the matrix directly.

    The function returnes an numpy array that can be processed
    with the function in file ``leech3matrix.c``.
    """
    a = np.zeros(72, dtype = np.uint64)
    if load_mode in (3, 15):
        b = matrix_to_rep_modp(m, load_mode)
        leech3matrix_load(load_mode, b, a)
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
# Test function leech3matrix_echelon
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




#######################################################################
# Test function leech3matrix_kernel_vector
#######################################################################

INT_BITS = 64

def kernel_testmatrix(p, tags = "dxypt"*3, diag = 0):
    space = MMSpace(p)
    a = space()
    diag %= 3
    for i in range(24):
        d = diag + randint(1, 2) + 3 * randint(0,4)
        a['A', i, i] = d % p if i else diag
    # Now ``a`` encodes the tag A part of the rep of the monster mod p
    # It encodes a symmetric matrix with kernel :math:`\Omega`.
    v3 = 1  # \Omega in Leech lattice mod 3 encoding
    g = MM.rand_word(tags, len(tags))
    space.mm.op_word_tag_A(a.data, g.data, len(g.data), 1)
    v3_kernel = gen_leech3_op_vector_word(v3, g.data, len(g.data))
    return a, v3_kernel

def neg_leech_mod3(x):
    t = (x ^ (x >> 24)) & 0xffffff
    return x ^ t ^ (t << 24)

def str_v3(v3):
    d = {0:0, 1:1, 0x1000000:2, 0x1000001:0}
    l = [str(d[(v3 >> i) & 0x1000001]) for i in range(24)]
    return "".join(l)

@pytest.mark.qstate
def test_leech3matrix_kernel_vector(verbose = 0):
    print("Testing function leech3matrix_kernel_vector()")
    for i in range(10):
        for p in (3, 15):
            diag = randint(0,2)
            v_a, v3_kernel =  kernel_testmatrix(p, "pxydl" * 4, diag)
            a = v_a.data
            a1 = np.zeros(72, dtype = np.uint64)
            a2 = np.zeros(48, dtype = np.uint64)
            if verbose:
                print("p = %d, diagonal = %d" % (p, diag))
                print("Part of tag A matrix")
                print(v_a["A", :16, :16])
                print("Kernel vector mod 3 expected:", str_v3(v3_kernel))
            v3_obt = leech3matrix_kernel_vector(p, a, diag)
            if verbose:
                print("Kernel vector mod 3 obtained:", str_v3(v3_obt))
                print("")
            leech3matrix_load(p, a, a1)
            leech3matrix_sub_diag(a1, diag, 0);
            leech3matrix_sub_diag(a1, 2, 24);
            leech3matrix_echelon(a1);
            leech3matrix_compress(a1, a2)
            #for i in range(24): print(i, hex(a2[2*i]))
            #for i in range(24): print(i, hex(a2[2*i+1]))
            v3_obt1 = xsp2co1_from_vect_mod3(a2[47])
            hh = xsp2co1_from_vect_mod3
            #print(hex(v3_obt), hex(v3_obt1) ) 
            data = [v3_obt1, v3_kernel, neg_leech_mod3(v3_kernel)]
            assert v3_obt in data[1:], list(map(hex, data))
            assert v3_obt1 in data[1:], list(map(hex, data))
            


#######################################################################
# Test functions leech2matrix_add_eqn and leech2matrix_solve_eqn
#######################################################################
       

def one_test_leech2matrix_eqn(n, verbose = 0):
    ml = [1 << i for i in range(n)]
    for k in range(2 * n * n):
        i = randint(0, n-1)
        j = randint(0, n-1)
        if i != j:
            ml[i] ^= ml[j]
    if verbose:
        print("Test function leech2matrix_add_eqn")
        print("m =", [hex(x) for x in ml])
    m, rows = np.zeros(n, dtype = np.uint64), 0
    for i, l in enumerate(ml):
        found = leech2matrix_add_eqn(m, rows, n, l)
        assert found == 1
        rows += found
        if i == 2:
            l_bad = ml[0] ^ ml[1]
            assert leech2matrix_add_eqn(m, rows, n, l_bad) == 0
        if i == 10:
            l_bad = ml[0] ^ ml[3] ^ ml[6] ^ ml[9]
            assert leech2matrix_add_eqn(m, rows, n, l_bad) == 0
    v = randint(0, (1 << n) - 1) 
    w = leech2matrix_solve_eqn(m, n, v)
    if verbose:
        print("inv =", [hex(int(x)) for x in m])
        print("v =", hex(v), ";", hex(w), " * m = v" )
    check = 0
    for i, l in enumerate(ml):
        if (w >> i) & 1:
            check ^=  l
    assert v == check


@pytest.mark.qstate
def test_leech2matrix_eqn():
    for n in range(2,33):
        one_test_leech2matrix_eqn(n)


       
#######################################################################
# Test leech3matrix_watermark and leech3matrix_watermark_perm_num
#######################################################################
  

MMV15 = MMSpace(15)
MM = MMV15.group    
         
def one_test_watermark(verbose = 0):
    v0 = MMV15.rand_uniform()
    v = v0.data[:48]
    w0 = np.zeros(24, dtype = np.uint32)
    result = leech3matrix_watermark(15, v, w0)
    if result < 0:
        return 0
    pi_num = randint(0, MAT24_ORDER-1)
    TAG_y, TAG_p = 0x40000000, 0x20000000
    y1, y2 = randint(0, 0xfff), randint(0, 0xfff)
    op = np.array([TAG_y + y1, TAG_p + pi_num, TAG_y + y2], 
         dtype = np.uint32)
    op_word_tag_A(v, op, len(op), 1)
    pi_num_obt = leech3matrix_watermark_perm_num(15, w0, v)   
    assert pi_num_obt == pi_num 
    return 1


WATERMARK_TESTS = 10
WATERMARK_MIN_SUCCESS = 4

@pytest.mark.qstate
def test_watermark():
    success = 0
    for n in range(WATERMARK_TESTS):
        success += bool(one_test_watermark())
    if success < WATERMARK_MIN_SUCCESS:
        err = "%d of %s permutation watermark tests failed"
        raise ValueError(err, WATERMARK_TESTS-success, WATERMARK_TESTS)