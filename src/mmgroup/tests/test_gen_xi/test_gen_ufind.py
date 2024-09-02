"""Test C functions for the union-find algorithm.

In this module we test the python function leech2_orbits_raw() in
module mmgroup.structures.xleech2, and the low-level C functions
in file gen_union_find.c used by that python function.

Function leech2_orbits_raw() deals with a huge data set; so we'd
better test the C functions used by that function on a smaller data
set first.

For testing the C functions we deal with the affine group H of
structure 2^8.L_3(2) acting naturally as a permutation group on the
affine space A8 correspondng to the linear space GF(2)^3. A8 has
size 8. Let 2**A8 be the power set of A8. For testing the C functions
we want to compute the orbits on the power set 2**A8 of A8 under the
action of H. Therefore we use the union-find algorithm implemented in
file gen_union_find.c.

These orbits are well known, so that that we can check the correctness
of the union-find algorithm. We generate H as a group of permutation
matrices acting the set GF(2)^8, which is isomorphic to 2**A8. Using
the union-find algorithm over the set 2**A8 of size 256, we obtain
the orbits on 2**A8 under the action of H.
"""

from collections import defaultdict
import pytest
import numpy as np
from copy import deepcopy
from random import sample

from mmgroup import Cocode, XLeech2, PLoop, Xsp2_Co1
from mmgroup import leech2_orbits_raw

from mmgroup.generators import gen_ufind_init
from mmgroup.generators import gen_ufind_find
from mmgroup.generators import gen_ufind_union
from mmgroup.generators import gen_ufind_find_all_min
from mmgroup.generators import gen_ufind_partition

from mmgroup.generators import gen_ufind_lin2_init
from mmgroup.generators import gen_ufind_lin2_size
from mmgroup.generators import gen_ufind_lin2_dim
from mmgroup.generators import gen_ufind_lin2_n_gen
from mmgroup.generators import gen_ufind_lin2_gen
from mmgroup.generators import gen_ufind_lin2_n_orbits
from mmgroup.generators import gen_ufind_lin2_orbits
from mmgroup.generators import gen_ufind_lin2_get_map
from mmgroup.generators import gen_ufind_make_map
from mmgroup.generators import gen_ufind_lin2_check
from mmgroup.generators import gen_ufind_lin2_len_orbit_v
from mmgroup.generators import gen_ufind_lin2_orbit_v
from mmgroup.generators import gen_ufind_lin2_rep_v
from mmgroup.generators import gen_ufind_lin2_map_v_gen
from mmgroup.generators import gen_ufind_lin2_map_v
from mmgroup.generators import gen_ufind_lin2_finalize
from mmgroup.generators import gen_ufind_lin2_check_finalized
from mmgroup.generators import gen_ufind_lin2_representatives

#####################################################################
# Bit matrix operations
#####################################################################


def chk(res, expected = None):
    """Check result 'res' of a C function

    If ``expected`` is None, the result ``res`` must be a nonnegative
    number. Otherwise ``res`` must be equal to ``expected``.
    """
    ok = (expected is None and res >= 0) or  res == expected
    if ok:
        return res
    if not ok:
        if expected is None:
            err = "Result of C function is %d" % res
        else:
            err = "Result obtained: %d, expected: %d" % (res, expected)
        raise ValueError(err)

#####################################################################
# Bit matrix operations
#####################################################################

def vmatmul(v, m, n):
    """Multiply bit vector v with n times n bit matrix m"""
    w = 0
    for i in range(n):
        if (v >> i) & 1:
            w ^= m[i]
    return w & ((1 << n) - 1)

def is_inv(m, mi, n):
    """Return True if mi is the inverse of the n times n bit matrix m"""
    acc = 0
    for i in range(n):
        acc |= vmatmul(m[i], mi, n) ^ (1 << i)
    return acc == 0

def v_mul_g(a, v, g):
    """Multiply vector v with group word g stored in a"""
    n = chk(gen_ufind_lin2_dim(a))
    m = np.zeros(n, dtype = np.uint32)
    for i in g:
        chk(gen_ufind_lin2_gen(a, i, m, n))
        vmatmul(v, m, n)
    return v

#####################################################################
# Define 8 times 8 bit batrices M7, M2, MT generating the group H
#####################################################################

def make_L3_2_matrix(t):
    data = [0, t[0], t[1], t[0] ^ t[1]]
    data += [t[2] ^ x  for x in data]
    return np.array([1 << x for x in data], dtype = np.uint32)

# Elements M2 and M7 generate the linear subgroup of H
M7 = make_L3_2_matrix([2,4,3])   # M7 is of order 7
M2 = make_L3_2_matrix([1,2,5])   # M7 is of order 7

# define a translation MT in the affine group h
MT = np.array([1 << (x ^ 1) for x in range(8)],  dtype = np.uint32)

# 8 times 8 unit bit matrix
M_UNIT8 = [1 << x for x in range(8)]

#print(M2, M7, MT)


#####################################################################
# Functions for converting internal data structures to a partition
#####################################################################



def partition_as_llist(ind, data):
    """Convert output of function leech2_orbits_raw() to partition

    Function leech2_orbits_raw() returns a triple
    (n_sets, ind, data) that describes a partition of a set of
    integers. This function takes the pair (ind, data) from that
    output and converts it to a list of list of integers
    decribing a partition of a set of integers. Here each list
    of integers is sorted; and the list of lists is sorted
    by the smallest elements of the lists of integers.
    """
    lst = []
    for i in range(len(ind) - 1):
        lst.append([int(x) for x in data[ind[i]:ind[i+1]]])
    return lst


def union_linear_low_level(generators):
    def vmatmul(v, m):
        v1 = 0;
        for j, x in enumerate(m):
            if (v >> j) & 1: v1 ^= x
        return v1

    gen = np.array(generators, dtype = np.uint32)
    n_gen, dim = gen.shape
    t_len = 1 << dim
    t =  np.zeros(t_len, dtype = np.uint32)
    map =  np.zeros(t_len, dtype = np.uint32)
    data =  np.zeros(t_len, dtype = np.uint32)
    chk(gen_ufind_init(t, t_len))
    for m in gen:
        for i in range(t_len):
           gen_ufind_union(t, t_len, i, vmatmul(i, m))
    n_orbits = chk(gen_ufind_find_all_min(t, t_len))
    chk(n_orbits)
    assert n_orbits > 0   
    ind = np.zeros(n_orbits + 1, dtype = np.uint32)
    chk(gen_ufind_partition(t, t_len, data, ind, n_orbits + 1), n_orbits)
    chk(gen_ufind_make_map(t, t_len, map), 0)
    del t
    return n_orbits, ind, data, map


def union_linear_high_level(generators):
    gen = np.array(generators, dtype = np.uint32)
    n_gen, dim = gen.shape
    len_a = chk(gen_ufind_lin2_size(dim, n_gen))
    assert len_a > 0, len_a
    a = np.zeros(len_a, dtype = np.uint32)
    chk(gen_ufind_lin2_init(a, len_a, dim, gen.ravel(), n_gen))
    t_len = 1 << chk(gen_ufind_lin2_dim(a))
    n_orbits = chk(gen_ufind_lin2_n_orbits(a))
    data = np.zeros(t_len, dtype = np.uint32)
    ind = np.zeros(n_orbits + 1, dtype = np.uint32)
    status = chk(gen_ufind_lin2_orbits(a, data, t_len, ind, n_orbits + 1))
    assert status >= 0, (2, status, t_len, n_orbits)
    map = np.zeros(t_len, dtype = np.uint32) 
    assert gen_ufind_lin2_get_map(a, map, t_len) == t_len
    assert gen_ufind_lin2_n_gen(a) == len(generators)
    assert gen_ufind_lin2_check(a, len(a)) >= 0
    return (n_orbits, ind, data, map), a


def equ_union_linear(u1, u2):
    """assert that two functions return equal data

    Here we check that the return values ``u1`` and ``u2`` of
    two calls to functions ``union_linear_high_level`` (or
    ``union_linear_low_level`) are equal.`  
    """
    assert u1[0] == u2[0]
    for a1, a2 in zip(u1[1:], u2[1:]):
        assert a1.shape == a2.shape
        assert (a1 == a2).all()


#####################################################################
# Check the orbits of GF(2)^8 under the action of H
#####################################################################

def check_properties_a(a, generators):
    """Check array ``a`` generated by ``gen_ufind_lin2_init``

    The function check some elementary properties of the array ``a``
    obtained by applying function ``gen_ufind_lin2_init`` to the
    set ``generators`` of geneators of a group.
    """
    gen = np.array(generators, dtype = np.uint32)
    n_gen, dim = gen.shape
    assert gen_ufind_lin2_n_gen(a) == n_gen
    assert gen_ufind_lin2_dim(a) == dim
    for i in range(n_gen):
        m = np.zeros(dim, dtype = np.uint32)
        mi = np.zeros(dim, dtype = np.uint32)
        chk(gen_ufind_lin2_gen(a, 2*i, m, dim))
        chk(gen_ufind_lin2_gen(a, 2*i + 1, mi, dim))
        assert (m == gen[i]).all()
        assert is_inv(m, mi, dim)
    chk(gen_ufind_lin2_check(a, len(a)))


def check_properties_a_llist(a, llist):
    for orbit in llist:
        v_set = sample(orbit, min(2, len(orbit)))
        for v in v_set:
            o = np.zeros(len(orbit), dtype = np.uint32)
            l_o = gen_ufind_lin2_len_orbit_v(a, v)
            assert l_o == len(orbit)
            assert gen_ufind_lin2_orbit_v(a, v, o, len(o)) == l_o
            assert list(o) == list(orbit), (l_o, len(orbit), len(o))
            assert gen_ufind_lin2_rep_v(a, v) == orbit[0]
    n_orbits = len(llist)
    r = np.zeros(n_orbits, dtype = np.uint32)
    assert n_orbits == gen_ufind_lin2_representatives(a, r, len(r))
    assert [x[0] for x in llist] == list(r)

def check_orbits_H(llist):
    """Check orbits on GF(2)^8 under the action of H

    Parameter ``llist`` must encode the partition defined by
    these orbits as described in function partition_as_llist().

    We check some well-known properties of this partition.

    H is transitive on the sets of GF(2)^8 of a given bit
    weight, except for weight 4, where H has 2 orbits.
    """
    from mmgroup.bitfunctions import bitweight
    SIZES =   [1, 8, 28, 56, 14, 56, 56, 28, 8, 1]
    WEIGHTS = [0, 1,  2,  3,  4,  4,  5,  6, 7, 8]
    # SIZES[i] is equal to the binomial coefficient
    # (8 choose WEIGHTS[i]), for WEIGHTS[i] != 4.
    assert len(llist) == len(SIZES)
    last = -1
    all_ = set()
    for i, lst in enumerate(llist):
        assert len(lst) == SIZES[i]
        assert lst[0] > last
        weight = WEIGHTS[i]
        assert bitweight(lst[0]) == weight
        for j in range(1, len(lst)):
            assert lst[j] > lst[j-1]
            assert bitweight(lst[j]) == weight
        last = lst[0]
        current = set(lst)
        assert len(current & all_) == 0
        all_ |= current
    assert all_ == set(range(256))

def check_properties_a_group(a, llist):
    g0 = np.zeros(30, dtype = np.uint8)
    #print([hex(x) for x in a[:20]])
    chk(gen_ufind_lin2_finalize(a))
    assert a[0] == 2, a[0]
    print("aaa", a[len(a) - 2])

    n_orbits = gen_ufind_lin2_n_orbits(a)
    assert n_orbits > 0, n_orbits
    r = np.zeros(n_orbits, dtype = np.uint32)
    assert n_orbits == gen_ufind_lin2_representatives(a, r, len(r))
    reps = set(r)
    table = a[4:]  # main table in internal structure
    n = chk(gen_ufind_lin2_dim(a))
    n_gen = chk(gen_ufind_lin2_n_gen(a))
    #print(a[0], n, n_gen, n_orbits)

    map_v_tmp = np.zeros(n, dtype = np.uint32)
    def map_v(a, v):
        g = a[v + 4] >> 24
        if g == 0xfe:
            return v
        if g < 2*n:
            assert gen_ufind_lin2_gen(a, g, map_v_tmp, n) == n
            return vmatmul(v, map_v_tmp, n)
        raise ValueError("Error in main table")

    for v in range(1 << n):
        g = table[v] >> 24
        if g == 0xfe:
            assert v in reps
        else:
            assert g < 2 * n_gen, (v, g)
            group_elem = np.zeros(n, dtype = np.uint32)
            chk(gen_ufind_lin2_gen(a, g, group_elem, n))
            w = vmatmul(v, group_elem, n)
            ref = chk(gen_ufind_lin2_rep_v(a, v))
            assert gen_ufind_lin2_rep_v(a, w) == ref

    a_g = np.zeros(20, dtype = np.uint8)
    for orbit in llist:
        m = np.zeros(n, dtype = np.uint32)
        v_set = sample(orbit, min(2, len(orbit)))
        for v in v_set:
            w = v
            g0 = chk(gen_ufind_lin2_map_v_gen(a, v))
            #print(v,g0)
            while g0 < 0xfe:
                chk(gen_ufind_lin2_gen(a, g0, m, n))
                w = vmatmul(w, m, n)
                g0 = chk(gen_ufind_lin2_map_v_gen(a, w))
            #print(v, w)
            assert w == chk(gen_ufind_lin2_rep_v(a, v))

            l_g = chk(gen_ufind_lin2_map_v(a, v, a_g, len(a_g)))
            assert l_g >= 0, (v, l_g, a_g[:10])
            a_g0 = a_g[:l_g]
            w = v_mul_g(a, w, a_g0)
            #print(v, w, orbit[0], a_g0)
            assert w == orbit[0], ("error", v, w, orbit[0], a_g0)

    assert gen_ufind_lin2_check_finalized(a, len(a)) == len(a)



@pytest.mark.gen_xi
def test_ufind_L3_2(verbose = 0):
    r"""Test the union-find algorithm on the goup H

    Here the group H is acting as a permutation group of the basis
    vectors of GF(2)^8 as described in the documentation string
    of this file.

    This function computes the orbits of GF(2)^8 under the action
    of H using the C functions in file gen_union_find.c.
    Then it checks some well-known properties of these orbits.
    """
    print("Testing C functions for union-find algorithm")

    generators = (M2, M7, MT, M_UNIT8)
    gen_hi, a = union_linear_high_level(generators)
    gen_lo = union_linear_low_level(generators)
    equ_union_linear(gen_hi, gen_lo)
    n_sets, ind, data, map = gen_hi
    assert n_sets == 10
    llist = partition_as_llist(ind, data)
    check_properties_a(a, generators)
    check_orbits_H(llist)
    check_properties_a_llist(a, llist)
    check_properties_a_group(a, llist)

#####################################################################
#  testing function eech2_orbits_raw()
#####################################################################


# This test is yet under construction!

@pytest.mark.gen_xi
def test_ufind_leech2(verbose = 0):
     g_list = [Xsp2_Co1('r', 'G_x0') for i in range(4)]
     n_sets, indices, data, mapping  = leech2_orbits_raw(g_list, map=True)
     lengths = [indices[i+1] - indices[i] for i in range(n_sets)]
     
