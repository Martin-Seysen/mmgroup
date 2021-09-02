from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample


import numpy as np
import pytest

from mmgroup import MM, MMSpace, Cocode
from mmgroup.clifford12 import leech2matrix_eval_A_odd_mod15_aux
from mmgroup.clifford12 import leech2matrix_eval_A
from mmgroup.mm import mm_aux_index_extern_to_sparse
from mmgroup.mm import mm_aux_index_leech2_to_sparse
from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup.mm import mm_aux_index_sparse_to_leech
V = MMSpace(15)



##################################################################

def eval_a_odd(v, signs):
    a = np.zeros(2, dtype = np.uint64)
    x = leech2matrix_eval_A_odd_mod15_aux(v.data, signs, a)
    v1 = [((int(a[i >> 4]) >> (4 * (i & 15)))) & 15 for i in range(24)]
    w = np.array(v1, dtype = np.int32) % 15
    return x, w


def eval_a_odd_ref(v, signs):
    A = np.array(v["A"], dtype = np.int32)
    v_list = [(-1)**((signs >> i) & 1) for i in range(24)]
    v1 = np.array(v_list, dtype =  np.int32)
    w = (v1 @ A) % 15
    return (w @ v1) % 15, w
  


def eval_a_odd_testdata():
    data = [
        (V(("A",2,2)), 0x3),
    ]    
    for d in data:
        yield d
    for i in range(5):
        v = V.rand_uniform()
        for j in range(10):
            yield v, randint(0, 0xffffff)




@pytest.mark.involution
def test_eval_a_odd(verbose = 0):   
    for n, (v, b) in enumerate(eval_a_odd_testdata()):
        if verbose:
            print("Test %d" % (n+1))
            print("b =", hex(b))
            print(v["A"])
        x, w = eval_a_odd(v, b)
        x_ref, w_ref =  eval_a_odd_ref(v, b)  
        assert (w == w_ref).all(), (w, w_ref)
        assert x == x_ref, (x, x_ref)


##################################################################


def eval_a_ref(v, v2, e = 0):
    if (e % 3):
        v = v * MM(('t', e))
    A = np.array(v["A"], dtype = np.int32)
    v2_sp = mm_aux_index_leech2_to_sparse(v2) 
    v_2 = np.zeros(24, dtype = np.int32)
    res = mm_aux_index_sparse_to_leech(v2_sp, v_2)
    return (v_2 @ A @ v_2) % 15



def rand_leech2():
    ext = 300 + randint(0, 98279)
    sp = mm_aux_index_extern_to_sparse(ext)
    return mm_aux_index_sparse_to_leech2(sp)


def eval_a_testdata():
    data = [
        (V(("A",2,2)), Cocode([2,3]).ord, 0),
    ]
    for d in data:
        yield d
    for i0 in range(24):
        for i1 in range(i0):
            yield V.rand_uniform(),  Cocode([i0,i1]).ord, 0
    for k in range(100):
        yield V.rand_uniform(), rand_leech2(), 0 
    for e in (1,2):
        for k in range(100):
            yield V.rand_uniform(), rand_leech2(), e 




@pytest.mark.involution
def test_eval_a(verbose = 0):   
    for n, (v, v2, e) in enumerate(eval_a_testdata()):
        if verbose:
            print("Test %d" % (n+1))
            print("v['A'] =")
            print(v['A'])
            print("v2 =", hex(v2), ", e =", e)
        m = v.eval_A(v2, e)
        m_ref =  eval_a_ref(v, v2, e)  
        assert m == m_ref, (hex(v2), e, m, m_ref)
        if e == 0:
            m_orig = leech2matrix_eval_A(15, v.data, v2)
            assert m == m_orig


