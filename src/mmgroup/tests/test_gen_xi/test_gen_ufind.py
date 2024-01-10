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
from mmgroup.generators import gen_ufind_union_affine
from mmgroup.generators import gen_ufind_union_leech2

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




def table_as_llist(w):
    """Convert internal union-find data structure to partition

    The union-find algorithm works with an internal data structure w
    created by function gen_ufind_init(). This function returns
    the partition encoded in data structure w in the same format
    as function partition_as_llist().
    """
    w = deepcopy(w)
    d = defaultdict(list)
    lw = len(w)
    for i in range(lw):
        f = gen_ufind_find(w, lw, i)
        if (f >= 0):
            d[f].append(i)
    lst = []
    for l1 in d.values():
        lst.append(l1)
    lst.sort()
    return lst




def union_affine(w, m, v):
    """Test function for debugging, which is now deprecated!

    Equivalent to gen_ufind_union_affine(w, len(w), m, len(m), v)
    """
    def vmatmul(v, m):
        v1 = 0;
        for j, x in enumerate(m):
            if (v >> j) & 1: v1 ^= x
        return v1
            
    for i in range(len(w)):
        im = v ^ vmatmul(i, m)
        gen_ufind_union(w, len(w), i, im)


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
    w = np.zeros(256, dtype = np.uint32)
    gen_ufind_init(w, len(w))
    w_copy = deepcopy(w)
    for m in (M2, M7, MT):
        gen_ufind_union_affine(w, len(w), m, len(m), 0)
        union_affine(w_copy, m, 0)
    assert (w == w_copy).all()
    ref_llist = table_as_llist(w)
    if verbose:
        for x in ref_llist:
            print(len(x), [hex(y) for y in x])
    gen_ufind_find_all_min(w, len(w))
    llist1 = table_as_llist(w)
    assert llist1 == ref_llist
    indices = np.zeros(11, dtype = np.uint32)
    n_sets = gen_ufind_partition(w, len(w), indices, len(indices))
    assert n_sets == 10
    llist2 = partition_as_llist(indices, w)
    assert llist2 == ref_llist
    check_orbits_H(ref_llist)

#####################################################################
#  testing function eech2_orbits_raw()
#####################################################################


# This test is yet under construction!

@pytest.mark.gen_xi
def test_ufind_leech2(verbose = 0):
     g_list = [Xsp2_Co1('r', 'G_x0') for i in range(4)]
     n_sets, indices, data = leech2_orbits_raw(g_list)
     lengths = [indices[i+1] - indices[i] for i in range(n_sets)]
     print(lengths)
     
