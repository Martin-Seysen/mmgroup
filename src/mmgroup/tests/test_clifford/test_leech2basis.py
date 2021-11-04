import sys
import os
import warnings

from random import randint
from collections import OrderedDict

import numpy as np
import pytest
from operator import __xor__
from functools import reduce

from mmgroup import MM0, MMSpace, MMV
from mmgroup.generators import gen_leech3_op_vector_word
from mmgroup.generators import gen_leech2_scalprod
from mmgroup.clifford12 import bitmatrix64_echelon_h
from mmgroup.clifford12 import leech2_matrix_basis
from mmgroup.clifford12 import leech2_matrix_orthogonal
from mmgroup.clifford12 import bitmatrix64_xch_bits
from mmgroup.clifford12 import bitmatrix64_rot_bits
from mmgroup.clifford12 import bitmatrix64_cap_h
from mmgroup.clifford12 import leech2_matrix_radical
from mmgroup.clifford12 import leech2_matrix_expand



def leech2_matrix_basis_testcases():
    for i in range(3, 500, 17):
        l = [randint(0, 0xffffff) for j in range(i)]
        yield np.array(l, dtype = np.uint32)
   

def leech_echelon(a):
    a = np.copy(a)
    bitmatrix64_xch_bits(a, 12, 0x800);
    bitmatrix64_rot_bits(a, 1, 12,0);
    k=bitmatrix64_echelon_h(a, len(a), 24, 24);
    bitmatrix64_rot_bits(a, 11, 12,0);
    bitmatrix64_xch_bits(a, 12, 0x800);
    return a[:k]


@pytest.mark.qstate
def test_leech2_matrix_basis():
    for a in leech2_matrix_basis_testcases():
        a1 = np.array(a, dtype = np.uint64)
        k1 = bitmatrix64_echelon_h(a1, len(a1), 24, 24)
        a1 = leech_echelon(a1)
        a2 = np.zeros(24, dtype = np.uint64)
        k2 = leech2_matrix_basis(a, len(a), a2, 24)
        assert k1 == k2 # == k2a
        assert (a1[:k1] == a2[:k1]).all()



def leech2_matrix_orthogonal_testcases():
    for i in [1, 3, 8, 11, 12, 15, 19, 23, 23, 24]:
        l = [randint(0, 0xffffff) for j in range(i)]
        yield np.array(l, dtype = np.uint64)
   

@pytest.mark.qstate
def test_leech2_matrix_orthogonal():
    for a in leech2_matrix_orthogonal_testcases():
        a1 = np.copy(a)
        o1 = np.zeros(24, dtype = np.uint64)
        ko1 = leech2_matrix_orthogonal(a, o1, len(a))
        o1 = o1[ko1:]
        ka1 = bitmatrix64_echelon_h(a1, len(a1), 24, 24)
        assert ka1 ==  ko1
        for v in a1:
            for w in o1:
                assert gen_leech2_scalprod(v, w) == 0
                 



@pytest.mark.qstate
def test_leech2_matrix_radical():
    for a in leech2_matrix_orthogonal_testcases():
        a1 = np.copy(a)
        k1 = bitmatrix64_echelon_h(a1, len(a1), 24, 24)
        a1 = a1[:k1]
        o1 = np.zeros(24, dtype = np.uint64)
        ko1 = leech2_matrix_orthogonal(a1, o1, len(a))
        o1 = o1[ko1:]
        ko1 = 24 - ko1
        ko1 = bitmatrix64_echelon_h(o1, len(o1), 24, 24)
        kcap =  bitmatrix64_cap_h(a1, o1, 24, 24);
        #print(k1, ko1, kcap)
        kcap1 = kcap[0]
        radical = a1[kcap1:k1]
        radical = leech_echelon(radical)

        a_copy = np.array(a, dtype = np.uint32)
        a2 = np.zeros(24, dtype = np.uint64)
        k2 = leech2_matrix_radical(a_copy, len(a_copy), a2, 24)
        a2 = a2[:k2]

        assert len(radical) == k2
        assert (a2 == radical).all()
        

@pytest.mark.qstate
def test_leech2_matrix_radical():
    for a in leech2_matrix_orthogonal_testcases():
        k = len(a)
        if (k > 14):
            continue
        v = np.zeros(1 << k, dtype = np.uint32)
        leech2_matrix_expand(a, k, v)
        for j in range(5):
            index = randint(0, (1 << k) - 1)
            terms = [int(x) for i, x in enumerate(a) 
                       if (1 << (k - 1 - i)) & index]
            ref = reduce(__xor__, terms, 0)
            assert v[index] == ref
        
 