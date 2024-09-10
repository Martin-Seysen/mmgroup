import sys
import os
import warnings

from random import randint
from collections import OrderedDict

import numpy as np
import pytest
from operator import __xor__
from functools import reduce

from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import leech2matrix_prep_eqn
from mmgroup.clifford12 import leech2matrix_solve_eqn
from mmgroup.clifford12 import leech2matrix_subspace_eqn
from mmgroup.clifford12 import uint64_parity



def rand_bit_matrix(rows, cols):
    imax = (1 << cols) - 1
    return np.array([randint(0, imax) for i in range(rows)], 
         dtype = np.uint64)


def leech2_matrix_basis_testcases():
    for rows, cols in [(3, 3), (3,8), (8,3), (29, 29)]:
        yield rand_bit_matrix(rows, cols), rows, cols
  
def mul_matrix_vector(a, v):
    y = 0
    for i, row in enumerate(a):
        y += uint64_parity(int(row) & int(v)) << i
    return y

 
@pytest.mark.qstate
def test_leech2_matrix_solve():
    for aa, rows, cols in leech2_matrix_basis_testcases():
        m = np.zeros(min(rows, cols), dtype = np.uint64)
        a_list = []
        nrows = 0
        for a in aa:
            a_in_space = leech2matrix_subspace_eqn(m, nrows, cols, a)
            status = leech2matrix_add_eqn(m, nrows, cols, a)
            assert 0 <= status <= 1
            if status:
                a_list.append(a)
                assert a_in_space < 0
            else:
                assert a_in_space >= 0 
            nrows += status
        assert 0 <= nrows <= cols
        a_valid = np.array(a_list, dtype = np.uint64)
        b = np.zeros(nrows, dtype = np.uint32)
        status = leech2matrix_prep_eqn(m, nrows, cols, b)
        assert status >= 0
        for k in range(3):
            v = randint(0, (1 << nrows) - 1)
            w = leech2matrix_solve_eqn(b, nrows, v)
            #print(a_valid, nrows, cols,  v, w)
            assert mul_matrix_vector(a_valid, w) == v
            leech2matrix_subspace_eqn(m, nrows, cols, v) == w
 