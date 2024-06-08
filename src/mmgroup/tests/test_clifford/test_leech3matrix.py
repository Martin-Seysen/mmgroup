import sys
import os
import warnings

from random import randint, sample
from collections import OrderedDict
import time

import numpy as np
import pytest

from mmgroup import MM0, MMSpace, MMV
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.clifford12 import leech3matrix_echelon
from mmgroup.clifford12 import leech3matrix_reduced_echelon
from mmgroup.clifford12 import leech3matrix_kernel_image
from mmgroup.clifford12 import leech3matrix_compress
from mmgroup.clifford12 import xsp2co1_from_vect_mod3
from mmgroup.clifford12 import leech3matrix_sub_diag
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import leech2matrix_solve_eqn
from mmgroup.clifford12 import leech3matrix_rank
from mmgroup.clifford12 import leech3matrix_vmul
from mmgroup.mat24 import MAT24_ORDER
from mmgroup.mm_op import mm_op_word_tag_A, mm_op_load_leech3matrix


def rand_matrix_mod_n(n, dtype = np.int64):
    """Create random 24 x 24 matrix of intgers modulo n"""
    return np.random.randint(0, n-1, size = (24, 24), dtype = dtype)

def print_mod3(m1, m2 = None):
    """Print 24 x 24 matrix of integers modulo 3"""
    m1a = m1 % 3
    if m2 is not None:
        m2a = m2 % 3
    for i in range(24):
        s = "".join(map(str, m1a[i] % 3))
        if m2 is not None:
            s += " " + "".join(map(str, m2a[i] % 3))
        print(s)
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
        err = "C functions do not support modulus p=%s"
        raise ValueError(err, str(p))
    return b



def from_array(m, load_mode = None):
    """Load 24 x 24 matrix 'm' to array of type 'leech3matrix'

    If load_mode is 3 or 15 then matrix is converted to a vector of
    of the monster rep using function ``matrix_to_rep_modp(m, p)``
    with ``p = load_mode``.

    By default we convert the matrix directly.

    The function returnes an numpy array that can be processed
    with the function in file ``leech3matrix.c``.
    """
    a = np.zeros(72, dtype = np.uint64)
    if load_mode == 15:
        b = matrix_to_rep_modp(m, load_mode)
        mm_op_load_leech3matrix(15, b, a)
    elif load_mode == 3:
        b = matrix_to_rep_modp(m, load_mode)
        mm_op_load_leech3matrix(3, b, a)
    elif not load_mode:
        for i in range(24):
            s = int(sum(int(m[i,j]) % 3 << 4*j for j in range(24)))
            a[3*i] = s & 0xffffffffffffffff
            a[3*i+1] = (s >> 64) & 0xffffffff
    else:
        err = "Loading from representation modulo % not suppoerted"
        raise ValueError(err % load_mode)
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


def subunit_matrix(rank, isect, dtype = np.int64):
    a = np.zeros((24, 24), dtype = dtype)
    assert isect <= rank
    assert rank + isect <= 24
    dif = rank - isect
    for i in range(dif):
        a[i, i] = randint(1,2)
        for j in range(i + 1, dif):
            a[i, j] =  randint(0, 2)
    for i in range(dif, rank):
        a[i, i + rank - dif] = randint(1,2)
    return a


def randomize_matrix(a):
    rows = list(range(a.shape[0]))
    cols = list(range(a.shape[1]))
    rows = list(range(4))
    for i in range(60):
        j1, j2 = sample(rows, 2)
        d = randint(1, 2)
        a[j1] += d * a[j2]
        a[:, j2] += (3 - d) * a[:, j1]
        a %= 3
    return a


TEST_MATRIX_0 = [ [2], [1, 1, 2, 2], [1, 1, 2, 2], [2] ]

TEST_MATRICES = [
    (TEST_MATRIX_0, (2, 1)),
]


def leech3matrix_echelon_testdata():
    for mat, data in TEST_MATRICES:
        m = np.zeros((24, 24), dtype = np.int64)
        for i, row in enumerate(mat):
            for j, entry in enumerate(row):
                 m[i, j] = entry
        yield m, data
    for rank in range(0, 24):
      for _ in range(2):
        for isect in range(0, min(rank, 24 - rank)):
            m = subunit_matrix(rank, isect, dtype = np.int64)
            m = randomize_matrix(m)
            yield m, (rank, isect)
    for i in range(20):
        yield rand_matrix_mod_n(3, dtype = np.int64), None


@pytest.mark.qstate
def test_leech3matrix_echelon(verbose = 0):
    if verbose:
        print("\nTest echelonization of 24 x 24 matrix mod 3\n")
    unit = np.eye(24)
    load_modes = [None, 15]
    for i, (m, rk) in enumerate(leech3matrix_echelon_testdata()):
        load_mode = load_modes[i % len(load_modes)]
        a =  from_array(m, load_mode)
        leech3matrix_sub_diag(a, 2, 24)
        if verbose:
            print("Test %d, load_mode = %s, rk = %s" %
                 (i+1, load_mode, rk))
            print_mod3(as_array(a,0), as_array(a,1))
        a1 = np.copy(a)
        leech3matrix_echelon(a1)
        if verbose:
     #       print("echelon:")
     #       print_mod3(as_array(a1,0), as_array(a1,1))
             pass
        mi = as_array(a1, high=1)
        prod = as_array(a1, high=0)
        assert ((mi @ m - prod) % 3 == 0).all()

        a2 = np.copy(a)
        cols = leech3matrix_reduced_echelon(a2, 0 if i % 3 else 24)
        mi2 = as_array(a2, high=1)
        prod2 = as_array(a2, high=0)
        if verbose:
            print("reduced echelon, columms: ", hex(cols))
            print_mod3(prod2, mi2)
        assert ((mi2 @ m - prod2) % 3 == 0).all()

        a3 = np.copy(a)
        l0 = leech3matrix_kernel_image(a3)
        assert l0 >= 0, l0
        l_isect, l_ker, l_image = l0 >> 16, (l0 >> 8) & 255, l0 & 255
        mi3 = as_array(a3, high=1)
        prod3 = as_array(a3, high=0)
        if verbose:
            S = "Intersection kernel/image, dim image = %d, dim rad = %d"
            print(S % (l_image, l_isect))
            print_mod3(prod3, mi3)
        assert ((mi3 @ m - prod3) % 3 == 0).all()

        for i in range(l_isect):
            assert (prod3[i] == mi3[i + l_image]).all()
        for i in range(l_image, 24):
            assert (prod3[i] == 0).all()
        if rk:
            ref_l_image, ref_l_isect = rk
            assert l_image == ref_l_image, (l_image == ref_l_image)
            assert l_isect == ref_l_isect, (l_isect == ref_l_isect)
            assert l_ker + l_image == 24, (l_ker, l_image)


@pytest.mark.qstate
@pytest.mark.bench
def test_bench_leech3matrix_kernel_image():
    matrices = []
    N_TESTS = 50000
    N_SAMPLES = 64
    MASK = N_SAMPLES - 1
    a = np.zeros(72, dtype = np.uint64)
    for i in range(N_SAMPLES):
        m = subunit_matrix(24 - i % 8 - i % 3, i % 3, dtype = np.int64)
        m = randomize_matrix(m)
        matrices.append(from_array(m))
    t_start = time.process_time()
    for i in range(N_TESTS):
        np.copyto(a, matrices[i & MASK])
        leech3matrix_kernel_image(a)
    t = time.process_time() - t_start
    S = "Runtime of function leech3matrix_kernel_image: %.3f us"
    print("\n" + S % (1.0e6 * t / N_TESTS))



#######################################################################
# Test function leech3matrix_kernel_vector
#######################################################################

INT_BITS = 64

def kernel_testmatrix(p, tags = "dxypt"*3, diag = 0):
    space = MMV(p)
    a = space()
    diag %= 3
    for i in range(24):
        d = diag + randint(1, 2) + 3 * randint(0,4)
        a['A', i, i] = d % p if i else diag
    # Now ``a`` encodes the tag A part of the rep of the monster mod p
    # It encodes a symmetric matrix with kernel :math:`\Omega`.
    v3 = 1  # \Omega in Leech lattice mod 3 encoding
    g = MM0([(t,'r') for t in tags])
    mm_op_word_tag_A(p, a.data, g.mmdata, len(g.mmdata), 1)
    v3_kernel = gen_leech3_op_vector_word(v3, g.mmdata, len(g.mmdata))
    return a, v3_kernel

def neg_leech_mod3(x):
    t = (x ^ (x >> 24)) & 0xffffff
    return x ^ t ^ (t << 24)

def str_v3(v3):
    d = {0:0, 1:1, 0x1000000:2, 0x1000001:0}
    l = [str(d[(v3 >> i) & 0x1000001]) for i in range(24)]
    return "".join(l)


def load_leech3matrix(p, v):
    a = np.zeros(3*24, dtype = np.uint64)
    if p in [3, 15]:
        mm_op_load_leech3matrix(p, v, a)
    else:
        err = "Computation of kernel vector not supporte for modulus %s"
        raise ValueError(err % p)
    return a


def kernel_vector(p, v, diag):
    a = load_leech3matrix(p, v)
    return leech3matrix_rank(a, diag) & 0xffffffffffff

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
            v3_obt = kernel_vector(p, a, diag)
            if verbose:
                print("Kernel vector mod 3 obtained:", str_v3(v3_obt))
                print("")
            a1 = load_leech3matrix(p, a)
            leech3matrix_sub_diag(a1, diag, 0);
            leech3matrix_sub_diag(a1, 2, 24);
            leech3matrix_echelon(a1);
            leech3matrix_compress(a1, a2)
            #for i in range(48): print(i, hex(a2[i]))
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
    m32 = np.array(m, dtype = np.uint32)
    w = leech2matrix_solve_eqn(m32, n, v)
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
# Test function leech3matrix_vmul
#######################################################################


TBL = [0, 1, 0x1000000,  0x1000001]
I_TBL_DICT = {0:0, 1:1,  0x1000000:2,  0x1000001:0}

def vector24_to_mod3(v):
    r = 3 * randint(0,1)
    w = 0
    for i in range(min(24, len(v))):
        x = v[i] % 3
        if x == 0:
            x = r
            r ^= 3
        w ^= TBL[x] << i
    return w

def vector_mod_3_to_24(v):
    l = []
    for i in range(24):
        l.append(I_TBL_DICT[(v >> i) & 0x1000001])
    return np.array(l, dtype = np.uint8)


def rand_v24():
    return np.random.randint(3, size=24, dtype=np.uint8)

def leech3matrix_vmul_testdata():
    """Yield test data for function leech3matrix_vmul()

    The function yields pairs (v, m), where v is a numpy array of
    24 integers between 0 and 23. m is a random vector in the rep
    :math:`\rho_3` of the Monster; so that m['A'] is a 24 times 24
    numpy matrix. One can use function vector24_to_mod3() to convert
    vector v to  **Leech lattice mod 3** encoding, as expected for
    the input of function leech3matrix_vmul(). function
    vector_mod_3_to_24() reverses the effect of function
    vector24_to_mod3(), to that the output of function
    leech3matrix_vmul() can be converted to a numpy vector again.
    """
    V3 = MMV(3)
    v3 = V3('R')
    for i in range(20):
        yield rand_v24(), v3
        v3 *= MM0(
            [('p', 'r'), ('t', 'n'), ('l', 'n'), ('p', 'r'), ('l', 'n')]
        )



@pytest.mark.qstate
def test_leech3matrix_vmul():
    for v, m in leech3matrix_vmul_testdata():
        v3 = vector24_to_mod3(v)
        v3prod = leech3matrix_vmul(v3, m.data)
        v3prod_vector = vector_mod_3_to_24(v3prod)
        v3_ref = (v @ m['A']) % 3
        assert (v3prod_vector == v3_ref).all()


