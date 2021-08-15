from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample

import numpy as np
import pytest

from mmgroup import MM, MMSpace
from mmgroup.mm import mm_aux_index_extern_to_sparse
from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup.mm import mm_aux_index_leech2_to_sparse
from mmgroup.mm15 import op_find_short as mm_op15_find_short

LEECH_SHORT = np.zeros(98280, dtype = np.int32)
for i in range(98280):
    sparse = mm_aux_index_extern_to_sparse(300 + i)
    l1 = LEECH_SHORT[i] = mm_aux_index_sparse_to_leech2(sparse)
    sparse1 = mm_aux_index_leech2_to_sparse(l1)
    assert sparse == sparse1, [hex(x) for x in (i, sparse, l1, sparse1)]

V = MMSpace(15)

v = V.rand_uniform()




def from_v(v_abs, value, value1 = 0):
    if not value1:
        a = [LEECH_SHORT[i] for i in range(98280) if v_abs[i] == value]
        return np.array(a, dtype = np.uint32)
    else:
        a0, a1 = from_v(v_abs, value), from_v(v_abs, value1)
        return np.concatenate([a0, a1  + 0x1000000])


def do_test_find_short15(v, value, value1= 0, verbose = 0):
    v_mon = v['E', 300:300+98280] 
    v_abs = np.where(v_mon > 7, 15 - v_mon, v_mon)
    a_ref = from_v(v_abs, value, value1)
    a = np.zeros(98280, dtype = np.uint32)
    n = mm_op15_find_short(v.data, a, 98280, value, value1)
    a = a[:n]
    if verbose:
        print("\nlen(a)=", len(a), ", len(a_ref)=", len(a_ref))
        print(a, a_ref)
    assert len(a) == len(a_ref)
    assert (a == a_ref).all()




testdata = [
    [V((1, "B", 21, 17), (3, "B", 2, 1)), 1, 3],
    [V((4, "C", 3, 0), (11, "T", 753, 53), (4, "X", 1753, 23)), 4, 0],
    [V.rand_uniform(), 1, 0], 
    [V.rand_uniform(), 7, 2], 
#   [V.rand_uniform(), 5, 6], 
]




@pytest.mark.involution
def test_find_short15(verbose = 0):  
    for v, y0, y1 in testdata:
        v_mon = v['E', 300:300+98280] 
        v_abs = np.where(v_mon > 7, 15 - v_mon, v_mon)
        do_test_find_short15(v, y0, y1, verbose)

