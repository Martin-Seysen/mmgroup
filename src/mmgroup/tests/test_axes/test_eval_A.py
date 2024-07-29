

from random import randint #, shuffle, sample


import numpy as np
import pytest

from mmgroup import MM0, MMV, MMSpace, Cocode
from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_index_leech2_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
from mmgroup.mm_op import mm_aux_index_sparse_to_leech
from mmgroup.mm_op import mm_op_eval_A_aux
from mmgroup.mm_op import mm_op_eval_A

V = MMV(15)



##################################################################

def eval_a_aux(v, masks, signs, row = 24):
    x = mm_op_eval_A_aux(15, v.data, masks, signs, row)
    res, row = x & 0xffff, x >> 16
    return res % 15,  row % 15


def eval_a_aux_ref(v, masks, signs, row = 24):
    A = np.array(v["A"], dtype = np.int32)
    v_list = [(-1)**((signs >> i) & 1) * ((masks >> i) & 1)   
        for i in range(24)] 
    v1 = np.array(v_list, dtype =  np.int32)
    w = (v1 @ A) % 15
    row_out = w[row] * v1[row] if 0 <= row < 24 else 0
    return (w @ v1) % 15, row_out % 15
    
 

  


def eval_a_aux_testdata():
    data = [
        (V("A",2,2), 0xffffff, 0x3, 24),
    ]    
    for d in data:
        yield d
    for i in range(5):
        v = V('R')
        for j in range(10):
            yield v, randint(0x0, 0xffffff), randint(0, 0xffffff), randint(0,24)




@pytest.mark.axes
def test_eval_a_aux(verbose = 0):   
    for n, (v, masks, signs, row) in enumerate(eval_a_aux_testdata()):
        if verbose:
            print("Test %d" % (n+1))
            print("masks = %s, signs = %s, row = %d" % (
                hex(masks), hex(signs), row))
            print(v["A"])
        x, w = eval_a_aux(v, masks, signs, row)
        x_ref, w_ref =  eval_a_aux_ref(v, masks, signs, row)  
        assert x == x_ref, (x, x_ref)
        assert (w == w_ref), (w, w_ref)


##################################################################


def eval_a_ref(v, v2, e = 0):
    if (e % 3):
        v = v * MM0('t', e)
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
        (V("A",2,2), Cocode([2,3]).ord, 0),
    ]
    for d in data:
        yield d
    for i0 in range(24):
        for i1 in range(i0):
            yield V('R'),  Cocode([i0,i1]).ord, 0
    for k in range(100):
        yield V('R'), rand_leech2(), 0 
    for e in (1,2):
        for k in range(100):
            yield V('R'), rand_leech2(), e 




@pytest.mark.axes
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
            m_orig = mm_op_eval_A(15, v.data, v2)
            assert m == m_orig


