"""Test C functions for the union-find algorithm.

In this module we test the C functions in module ``gen_ufind_lin2.c``
using their python wrappers in the ``mmgroup.generators`` extension.
In that file we define an *orbit array* containing orbits of the
action of a group on a linear space GF(2)^n. We also test the python
class ``Orbit_Lin2`` for wrapping such an orbit array. On the way we
also test the C functions in module ``gen_ufind_lin2.c`` implementing
a union-find algorithm used by the functions in ``gen_ufind_lin2.c``.

For testing the C functions we deal with the affine group H of
structure 2^8.L_3(2) acting naturally as a permutation group on the
affine space A8 corresponding to the linear space GF(2)^3. A8 has
size 8. Let P := 2**A8 be the power set of A8. P has a natural
structure an a vector space GF(2)^8. For our tests we compute the
orbits on the space P under the action of H using the functions in
module ``gen_ufind_lin2.c`` and also the python class ``Orbit_Lin2``.
These orbits are well known, so that that we can check the
correctness of the functions to be tested. We generate H as a
group of permutation matrices acting the vector space P.
"""

from collections import defaultdict
from numbers import Integral
import pytest
import numpy as np
from copy import deepcopy
from random import sample, randint

from mmgroup import Cocode, XLeech2, PLoop, Xsp2_Co1
from mmgroup import leech2_orbits_raw
from mmgroup.bitfunctions import bitweight


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
from mmgroup.generators import gen_ufind_lin2_orbit_lengths
from mmgroup.generators import gen_ufind_lin2_get_table
from mmgroup.generators import gen_ufind_lin2_add
from mmgroup.generators import gen_ufind_lin2_compressed_size
from mmgroup.generators import gen_ufind_lin2_compress



from mmgroup.general import Orbit_Lin2
from mmgroup.tests.test_general.bitmatrix import chk, BitMatrix
from mmgroup.tests.test_general.bitmatrix import vmatmul, is_inv
from mmgroup.tests.test_general.bitmatrix import  v_mul_g



a = None




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
    global a
    a = np.zeros(len_a, dtype = np.uint32)
    chk(gen_ufind_lin2_init(a, len_a, dim, n_gen))
    for g in gen:
        chk(gen_ufind_lin2_add(a, g, dim))
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
# Generate the main orbit array 'a' for a certain group H
#####################################################################

"""Bit matrices generating a group H

For testing the C functions we deal with the affine group H of
structure 2^8.L_3(2) acting naturally as a permutation group on the
affine space A8 correspondng to the linear space GF(2)^3. A8 has
size 8. Let P := 2**A8 be the power set of A8. P has a natural
structure an a vector space GF(2)^8. H has a natural action on
as a linear group on the vector space P. In this section we define
8 times 8 bit matrices M2, M7, MT (acting on P by right
multiplication) that generate the group H.
"""

def make_L3_2_matrix(t):
    data = [0, t[0], t[1], t[0] ^ t[1]]
    data += [t[2] ^ x  for x in data]
    return np.array([1 << x for x in data], dtype = np.uint32)

# Elements M2 and M7 generate the linear subgroup of H
M7 = make_L3_2_matrix([2,4,3])   # M7 is of order 7
M2 = make_L3_2_matrix([1,2,5])   # M2 is of order 7

# define a translation MT in the affine group h
MT = np.array([1 << (x ^ 1) for x in range(8)],  dtype = np.uint32)

# 8 times 8 unit bit matrix
M_UNIT8 = [1 << x for x in range(8)]

#print(M2, M7, MT)



def generate_main_orbit_array(generators):
    r"""Return main orbit array and orbits of group H

    Given a list of generators of the group ``H`` described above, the
    function returns a pair ``(a, llist)``. Here ``a`` is an (opaque)
    **orbit array** for the group H and its natural action on
    :math:`\mbox{GF}_2^3` as desribed in file ``gen_ufind_lin2.c``.
    Object ``llist`` is a list of lists describing these orbits.

    The function makes some checks that are specific for the
    group ``H``.
    """
    global a
    gen_hi, a = union_linear_high_level(generators)
    gen_lo = union_linear_low_level(generators)
    equ_union_linear(gen_hi, gen_lo)
    n_sets, ind, data, map = gen_hi
    assert n_sets == 10
    llist = partition_as_llist(ind, data)
    return a, llist


#####################################################################
# Check the orbits of GF(2)^8 under the action of H
#####################################################################


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

#####################################################################
# Check the main orbit array ``a`` for the group H
#####################################################################


def check_properties_a(a, generators):
    """Check orbit array ``a`` generated by ``generate_main_orbit_array``

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
    """Check orbits in the orbit array ``a``

    This function checks the orbits in the orbit array ``a`` generated
    by  function ``generate_main_orbit_array``. That function also
    generates a list ``llist`` of these orbits. Function
    ``check_properties_a_llist`` check the orbits stored in the
    array ``a`` against the orbits o ``llist``.
    """
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


#####################################################################
# Check the Schreier vector in the main orbit array ``a``
#####################################################################



def check_Schreier_vector(a, llist):
    """Check Schreier vectors in the orbit array ``a``

    This function checks the Schreier vectors in the orbit array
    ``a`` generated  by  function ``generate_main_orbit_array``.
    There are functions in module ``gen_ufind_lin2_c`` that use a
    Schreier vector in the orbit array ``a`` for the group ``H``
    for finding an element of ``H`` that maps an arbitrary vector
    in an orbit the representative of that orbit stored in ``a``.

    This function tests these C functions in
    module ``gen_ufind_lin2_c``.
    """
    g0 = np.zeros(30, dtype = np.uint8)
    #print([hex(x) for x in a[:20]])
    chk(gen_ufind_lin2_finalize(a))
    assert a[0] == 2, a[0]

    n_orbits = gen_ufind_lin2_n_orbits(a)
    assert n_orbits > 0, n_orbits
    r = np.zeros(n_orbits, dtype = np.uint32)
    assert n_orbits == gen_ufind_lin2_representatives(a, r, len(r))
    reps = set(r)
    n = chk(gen_ufind_lin2_dim(a))
    table =np.zeros(1 << n, dtype = np.uint32)
    # Copy main table from internal structure ``a`` to ``table``
    chk(gen_ufind_lin2_get_table(a, table, len(table)))
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
            assert w == orbit[0], ("error", v, w, orbit[0], a_g0)

            l_g = chk(gen_ufind_lin2_map_v(a, v, a_g, len(a_g)))
            assert l_g >= 0, (v, l_g, a_g[:10])
            a_g0 = a_g[:l_g]
            w = v_mul_g(a, v, a_g0)
            #print(v, w, orbit[0], a_g0)
            assert w == orbit[0], ("error", v, w, orbit[0], a_g0)

    assert gen_ufind_lin2_check_finalized(a, len(a)) == len(a)



#####################################################################
# Check the Python class ``Orbit_Lin2`` for wrapping an orbit array
#####################################################################


def mul_v_word_G(orbits, v, word):
    gen = orbits.generators()
    for i, e in word:
        v *= gen[i] ** e
    return v

def check_orbit_object(orbits, generators, a, llist):
    reps, lengths = orbits.representatives()
    assert orbits.dim == gen_ufind_lin2_dim(a)
    assert len(orbits.generators()) == len(generators)
    #print(orbits.generators())
    py_generators = [BitMatrix(g) for g in generators]
    #print(py_generators)
    assert orbits.generators() == py_generators
    assert orbits.n_orbits() == len(llist)

    for i, o in enumerate(llist):
        orbit = orbits.orbit(o[0])
        assert list(orbit) == list(o)
    for n in range(10):
        v = randint(0, (1 << orbits.dim) - 1)
        w = orbits.orbit_rep(v)
        assert w == gen_ufind_lin2_rep_v(a, v)
        g = orbits.map_v_G(v)
        assert w == v * g, (v, w)
        g_list = [2 * g + int(s < 0) for g, s in orbits.map_v_word_G(v)]
        b = np.zeros(50, dtype = np.uint8)
        lb = chk(gen_ufind_lin2_map_v(a, v, b, len(b)))
        assert g_list == list(b[:lb])
        w = sample(list(orbits.orbit(v)), 1)[0]
        g = orbits.map_v_G(v, w)
        assert w == v * g, (v, w)
        assert w == orbits.mul_v_g(v, g)
        word = orbits.map_v_word_G(v, w)
        assert w == mul_v_word_G(orbits, v, word)
    for lst in llist:
        for v in sample(lst, min(2, len(lst))):
            assert orbits.orbit_size(v) == len(lst)
            assert list(orbits.orbit(v)) == lst
            assert orbits.orbit_rep(v) == lst[0]
    rep, lengths = orbits.representatives()
    assert len(rep) == len(lengths) == len(llist)
    for i, lst in enumerate(llist):
        assert rep[i] == lst[0]
        assert lengths[i] == len(lst)


def check_orbits_py_class(generators, a, llist):
    """Return Orbit_Lin2 object corresponding to main orbit array ``a``

    Given a list ``generators`` of the group ``H``, the function
    returns ob object ``orbits`` of class ``Orbit_Lin2`` containing
    the orbit array for the action of the group ``H`` on the vector
    space ``P``.

    Arguments ``a`` and ``llist`` must be provided as returned by
    function ``generate_main_orbit_array(generators)``. This test
    function also tests the returned python object ``orbits``
    against the orbit array ``a`` and the list ``llist`` of orbits.
    """
    orbits = Orbit_Lin2(lambda g : g.matrix,
                            [BitMatrix(g) for g in generators])
    check_orbit_object(orbits, generators, a, llist)
    data, functions = orbits.pickle()
    new_orbits = Orbit_Lin2(data, functions)
    check_orbit_object(new_orbits, generators, a, llist)
    return orbits

#####################################################################
# Check a compressed version of  orbit array ``a``
#####################################################################


def gen_weight4():
    """Genereate integers 0 <= i < 256 of bit weight 4"""
    for i in range(256):
        if bitweight(i) == 4:
            yield i



def check_properties_a_compressed(a):
     """Check a compressed version of the main orbit array ``a``

     Function ``gen_ufind_lin2_compress`` returns a compressed
     version ``c`` of the main orbit array ``a`` that contains fewer
     orbits the the original array. This function computes such a
     comprressed orbit array ``c`` for the group ``H`` containing
     orbits of vectors of bit weight 4 only. This function tests
     the functions in module ``gen_ufind_lin2.c`` acting on the
     compressed array ``c``.
     """
     o = np.array([0x0f, 0x17], dtype = np.uint32)
     sizes = {0x0f : 14, 0x17 : 56}
     l_o = len(o)
     l_c = chk(gen_ufind_lin2_compressed_size(a, o, l_o))
     c = np.zeros(l_c, dtype = np.uint32)
     status = gen_ufind_lin2_compress(a, o, l_o, c, l_c)
     chk(status)
     chk(gen_ufind_lin2_check_finalized(c, l_c))
     dim = gen_ufind_lin2_dim(a)
     assert gen_ufind_lin2_dim(c) == dim
     n_gen = gen_ufind_lin2_n_gen(c)
     assert n_gen == gen_ufind_lin2_n_gen(a)
     assert gen_ufind_lin2_n_orbits(c) == len(o)
     q = np.zeros(len(o), dtype = np.uint32)
     gen_ufind_lin2_representatives(c, q, len(q))
     assert list(q) == list(o)
     gen_ufind_lin2_orbit_lengths(c, q, len(q))
     assert list(q) == [sizes[x] for x in o]
     for i in range(n_gen):
          g = np.zeros((dim,2), dtype = np.uint32)
          chk(gen_ufind_lin2_gen(a, i, g[0], dim))
          chk(gen_ufind_lin2_gen(c, i, g[1], dim))
          assert (g[0] == g[1]).all()
     #print(c[:10])
     #print("cccc", c[10:26])
     for i, v in enumerate(gen_weight4()):
         w = chk(gen_ufind_lin2_rep_v(c, v))
         assert w == gen_ufind_lin2_rep_v(a, v)
         size = chk(gen_ufind_lin2_len_orbit_v(c, v))
         #print(hex(v), hex(w), size, sizes[w])
         assert size == sizes[w], (hex(v), hex(w), size, sizes[w])
         if i & 3 == 0:
            b = np.zeros(50, dtype = np.uint8)
            lb = chk(gen_ufind_lin2_map_v(c, v, b, len(b)))
            w1 = v_mul_g(c, v, b[:lb])
            #print("bbbb", b[:lb], hex(w), hex(w1))
            assert w1 == w
     return c

#####################################################################
# Check a compressed version of Python class ``Orbit_Lin2``
#####################################################################

def check_compressed(c_py, a_py, c, orbits):
    assert c_py.dim == a_py.dim
    assert c_py.n_orbits() == len(orbits)
    assert (c_py.generators() == a_py.generators())
    reps, sizes = [list(x) for x in c_py.representatives()]
    assert reps == orbits
    ref_sizes_a = [a_py.orbit_size(v) for v in orbits]
    ref_sizes_c = [c_py.orbit_size(v) for v in orbits]
    assert sizes == ref_sizes_a == ref_sizes_c
    v1 = 0
    # Actual test is yet to be implemented!
    for i, v in enumerate(gen_weight4()):
        if i % 3 != 2:
            continue
        w = a_py.orbit_rep(v)
        assert  w == c_py.orbit_rep(v)
        assert a_py.orbit_size(v) == c_py.orbit_size(v)
        # print(v, w, orbits[1], v1)
        if w == orbits[1] and v1:
            assert a_py.map_v_word_G(v, v1) == c_py.map_v_word_G(v, v1)
            g = a_py.map_v_G(v, v1)
            assert g == c_py.map_v_G(v, v1)
            assert c_py.mul_v_g(v, g) == v1 
        if w == orbits[1]:
            v1 = v

def check_properties_a_py_compressed(a_py, c):
    """Check Orbit_Lin2 object corresponding to compressed orbit array

    Let ``a_py`` be an object of class ``Orbit_Lin2`` corresponding
    to the main orbit array ``a`` for the group ``H``; and let ``c``
    be the compresses orbit array returned by function
    ``check_properties_a_group_compressed``.

    This function computes a python object ``c_py`` of class
    ``Orbit_Lin2`` corresponding to the compressed orbit array ``c``.
    Therefore it uses method ``compress`` of class ``Orbit_Lin2``. It
    tests the methods of the python  object ``c_py`` against the
    corresponding functions in module ``gen_ufind_lin2.c`` acting
    on the compressed array ``c``.
    """
    orbits = [0x0f, 0x17]
    c_py = a_py.compress(orbits)
    check_compressed(c_py, a_py, c, orbits)
    data, functions = c_py.pickle()
    new_c_py = Orbit_Lin2(data, functions)
    check_compressed(new_c_py, a_py, c, orbits)
    return c_py

#####################################################################
# Main test function
#####################################################################


@pytest.mark.general
def test_ufind_L3_2(verbose = 0):
    r"""Test the union-find algorithm on the goup H

    Here the group H is acting as a permutation group of the basis
    vectors of GF(2)^8 as described in the documentation string
    of this file.

    This function tests the functions in module ``gen_ufind_lin2.c``
    and also the python class ``Orbit_Lin2`` with the action of
    the group H on the vector space GF(2)^8.
    """
    print("Testing C functions in module gen_ufind_lin2.c")
    # Generators for the group H
    generators = (M2, M7, MT, M_UNIT8)
    a, llist = generate_main_orbit_array(generators)
    check_orbits_H(llist)
    check_properties_a(a, generators)
    check_properties_a_llist(a, llist)
    check_Schreier_vector(a, llist)
    a_py = check_orbits_py_class(generators, a, llist)
    c = check_properties_a_compressed(a)
    check_properties_a_py_compressed(a_py, c)

