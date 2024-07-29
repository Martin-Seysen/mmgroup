
from random import randint #, shuffle, sample

import numpy as np
import pytest

from mmgroup import MM0, MMSpace, MMV
from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
from mmgroup.mm_op import mm_aux_index_leech2_to_sparse
from mmgroup.mm_op import mm_op_eval_X_find_abs




LEECH_SHORT =  None

def load_leech_short():
    global LEECH_SHORT
    LEECH_SHORT = np.zeros(98280, dtype = np.int32)
    for i in range(98280):
        sparse = mm_aux_index_extern_to_sparse(300 + i)
        l1 = LEECH_SHORT[i] = mm_aux_index_sparse_to_leech2(sparse)
        sparse1 = mm_aux_index_leech2_to_sparse(l1)
        assert sparse == sparse1, [hex(x) for x in (i, sparse, l1, sparse1)]

def test_extern_sparse():
    load_leech_short()


V = MMV(15)

v = V('R')




def from_v(v_abs, value, value1 = 0):
    load_leech_short()
    if not value1:
        a = [LEECH_SHORT[i] for i in range(98280) if v_abs[i] == value]
        return np.array(a, dtype = np.uint32)
    else:
        a0, a1 = from_v(v_abs, value), from_v(v_abs, value1)
        return np.concatenate([a0, a1  + 0x1000000])


def do_test_eval_X_find_abs(v, value, value1= 0, verbose = 0):
    v_mon = v['E', 300:300+98280] 
    v_abs = np.where(v_mon > 7, 15 - v_mon, v_mon)
    a_ref = from_v(v_abs, value, value1)
    a = np.zeros(98280, dtype = np.uint32)
    n = mm_op_eval_X_find_abs(15, v.data, a, 98280, value, value1)
    a = a[:n]
    if verbose:
        print("\nlen(a)=", len(a), ", len(a_ref)=", len(a_ref))
        print(a, a_ref)
    assert len(a) == len(a_ref)
    assert (a == a_ref).all()




testdata = [
    [V([(1, "B", 21, 17), (3, "B", 2, 1)]), 1, 3],
    [V([(4, "C", 3, 0), (11, "T", 753, 53), (4, "X", 1753, 23)]), 4, 0],
    [V(), 1, 0], 
    [V('R'), 7, 2], 
#   [V('R'), 5, 6], 
]




@pytest.mark.axes
def test_eval_X_find_abs(verbose = 0):  
    for v, y0, y1 in testdata:
        v_mon = v['E', 300:300+98280] 
        v_abs = np.where(v_mon > 7, 15 - v_mon, v_mon)
        do_test_eval_X_find_abs(v, y0, y1, verbose)

