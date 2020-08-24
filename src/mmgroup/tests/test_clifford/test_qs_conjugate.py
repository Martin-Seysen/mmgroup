from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix, rand_qs_matrix
from mmgroup.structures.qs_matrix import qstate_column_monomial_matrix
from mmgroup.tests.test_clifford.test_qs_product import create_m 
from mmgroup.tests.test_clifford.test_qs_product import rand_mul_scalar 

#####################################################################
# Create test matrices
#####################################################################




qs_matrix_data = [
    [ 0, 0, (1,0), []],
    [ 0, 4, (0,0), [0b00_0101, 0b00_0011]  ],
]



def create_testmatrices():
    """yield instances of class ``QStateMatrix`` for tests """
    for m_data in qs_matrix_data:
        m = create_m(*m_data)
        yield m
    yield QStateMatrix(4, 4,  1)  
    for i in range(20):
        yield rand_qs_matrix(2, 2, 4)
    for cols in range(0,7):
        for rows in range(0, 7):
            if rows + cols > 9:
                 break
            yield QStateMatrix(rows, cols)
            for data_rows in range(1, rows + cols + 3):
                for n in range(2):
                    m = rand_qs_matrix(rows, cols, data_rows)
                    rand_mul_scalar(m)  
                    yield m   
 

def bitrange_mask(imin, imax):
    """Return bit mask with bits set in ``range(imin, imax)``"""
    return (1 << imax) - (1 << imin)

def is_leading_low(x, j):
    """Return True if bit ``j`` is the lowest bit set in ``x``."""
    return x & ((1 << (j+1)) - 1) == 1 << j

def zero_column(d, j, i):
    """Return True if ``d[k,j] = 0`` for ``0 <= k < i``

    Here ``d`` is a bit matrix implemented as an array of integers.
    """
    s = 0
    for k in range(i):
        s |= d[k]
    return (s & (1 << j)) == 0

EPS = 1.0e-8

@pytest.mark.qstate
def test_reduce_matrix(verbose = 0):
    for ntest, m in enumerate(create_testmatrices()):
        # Check method m.reduce_matrix() for matrix m
        m1 = m.copy()
        if verbose:
            print("Test %d: reduce matrix, m =" % (ntest + 1), m)
            print("Standard reduced m =", m.copy().reduce())
        row_tab = m1.reduce_matrix()
        m1.check()
        if verbose:
            print("Reduced matrix, m1 =", m1)
            print("table", row_tab)
        if  m != m1:
            print("m reduced", m.copy().reduce())
            print("m1", m1.copy())
            print("m1 reduced", m1.copy().reduce())
            err = "m and m1 are not equal"
            raise ValueError(err)
        if m1.nrows == 0: 
            assert  m.lb_rank() == -1
            continue
        assert len(row_tab) == m1.nrows + m1.ncols
        d = m1.data
        imin = 0
        fst_row = row_tab[m1.rows + m1.cols]
        for j in range(m1.cols):
            i =  row_tab[j]
            if i < 255: 
                assert i == imin + 1
                imin = i
                if not is_leading_low(d[i], j):
                    err = "error in row %d of m1"% i 
                    print("m1", m1)
                    print("table", row_tab)
                    raise ValueError(err)
        row_mask = bitrange_mask(m1.cols, m1.cols + m1.rows)
        row_set = set(range(1, m1.nrows))
        for j in range(m1.cols, m1.cols + m1.rows):
            i =  row_tab[j]
            ok = True
            if i < fst_row:
                ok = is_leading_low(d[i] & row_mask, j)
                ok = ok and zero_column(d, j, i)
                assert row_tab[m1.ncols + i] == j
                row_set.remove(i)
            elif i < 255: 
                assert i == imin + 1
                imin = i
                ok = is_leading_low(d[i], j)
                ok = ok and zero_column(d, j, fst_row)
                row_set.remove(i)
            if not ok:
                err = "error in column %d, row %d of m1" % (j, i) 
                print("m1", m1)
                print("table", row_tab)
                raise ValueError(err)
        assert imin == m1.nrows - 1
        for i in row_set:
            assert d[i] & row_mask == 0

        q_start, q_end = m1.ncols + fst_row, m1.ncols + m1.nrows
        q_mask = bitrange_mask(q_start, q_end)
        for j in range(q_start, q_end):
            i = row_tab[j]
            if i < 255:
                q_row = d[i] & q_mask
                assert i in row_set
                assert q_row == 1 << j
                row_set.remove(i)
        for i in row_set:
            assert d[i] & q_mask == 0
            
        # Finally, check the computation of the rank of matrix m              
        if m.ncols > 8:
            continue
        rk = np.linalg.matrix_rank(m[:], tol=EPS)
        lb_rk = m.lb_rank()
        assert rk == 1 << lb_rk, (m1[:], rk, lb_rk)   

        