"""Test computing the scalar prdict in class MMVector"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint
import numpy as np

import pytest

from mmgroup.mm_space import characteristics
from mmgroup import MMV, MM0, MMVector, mmv_scalprod
from mmgroup.mm_op import mm_aux_mmv_size, mm_aux_index_mmv
from mmgroup.mm_op import mm_op_scalprod_ind


PRIMES = characteristics()


def py_scalprod(v1, v2):
    """Scalar product of vectors v1 and v2 in class MMVector"""
    assert isinstance(v1, MMVector) and isinstance(v2, MMVector)
    p = v1.p
    assert p == v2.p
    b1 = np.array(v1.as_bytes(), dtype = np.int64)
    b2 = np.array(v2.as_bytes(), dtype = np.int64)
    assert len(b1) == 196884
    res = np.dot(b1, b2) + np.dot(b1[24:300], b2[24:300])
    return int(res % p)




def dense_v(p):
    """Generate vector mod p containing many one bits"""
    V = MMV(p)
    v1 = V('R')
    v2 = V('R')
    if p <= 3:
        offset = 1
    elif p <= 15:
        offset = 2
    else:
        offset = 4
    offset *= 72
    v1.data[offset:] |= v2.data[offset:]
    v1.check()
    #print(v1['T', 3, 4:10])
    return v1




@pytest.mark.mm_op
def test_scalprod(ntests=10):
    """Test function mmgroup.mmv_scalprod()"""
    RAND_G = [MM0('r', 3) for i in range(4)]
    def rand_g():
        return RAND_G[randint(0, len(RAND_G)-1)]
    for p in PRIMES:
        ind_buf = None
        if p in [3, 15]:
            ind_buf = np.zeros(mm_aux_mmv_size(p) + 1, dtype = np.uint16)
        V = MMV(p)
        for i in range(ntests):
            v1, v2 = V('R'), V('R') 
            prod = mmv_scalprod(v1, v2)
            ref_prod = py_scalprod(v1, v2)
            assert prod == ref_prod, p
            if i >= 2:
                continue
            g = rand_g()
            v1 *= g
            v2 *= g
            assert mmv_scalprod(v1, v2) == ref_prod
            v1, v2 = dense_v(p), dense_v(p)  
            prod = mmv_scalprod(v1, v2)
            ref_prod = py_scalprod(v1, v2)
            assert prod == ref_prod, p
            if ind_buf is not None:
                st = mm_aux_index_mmv(p, v1.data, ind_buf, len(ind_buf))
                assert st > 0
                prod = mm_op_scalprod_ind(p, v1.data, v2.data, ind_buf)
                assert prod == ref_prod, p
