"""Test C functions for union-find algorithm

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


def make_L3_2_matrix(t):
    data = [0, t[0], t[1], t[0] ^ t[1]]
    data += [t[2] ^ x  for x in data]
    return np.array([1 << x for x in data], dtype = np.uint32)


M7 = make_L3_2_matrix([2,4,3])
M2 = make_L3_2_matrix([1,2,5])
MT = np.array([1 << (x ^ 1) for x in range(8)],  dtype = np.uint32)

#print(M2, M7, MT)


def table_as_llist(w):
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


def partition_as_llist(ind, data):
    lst = []
    for i in range(len(ind) - 1):
        lst.append([int(x) for x in data[ind[i]:ind[i+1]]])
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

RR = list(range(256))

@pytest.mark.gen_xi
def test_ufind_L3_2(verbose = 0):
    r"""Test generation union-find algorithm

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


@pytest.mark.gen_xi
def test_ufind_leech2(verbose = 0):
     g_list = [Xsp2_Co1('r', 'G_x0') for i in range(4)]
     n_sets, indices, data = leech2_orbits_raw(g_list)
     lengths = [indices[i+1] - indices[i] for i in range(n_sets)]
     print(lengths)
     
