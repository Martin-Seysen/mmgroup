from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from math import log
from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix, rand_qs_matrix
from mmgroup.structures.qs_matrix import qstate_column_monomial_matrix
from mmgroup.tests.test_clifford.test_qs_product import create_m 
from mmgroup.tests.test_clifford.test_qs_product import rand_mul_scalar 
from mmgroup.clifford12 import qstate12_bit_matrix_t

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


#####################################################################
# Test method reduce_matrix()
#####################################################################


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
            print("shape", m1.shape)
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
        fst_row = row_tab[m1.ncols]
        for j in range(m1.shape[1]):
            i =  row_tab[j]
            if i < 255: 
                assert i == imin + 1
                imin = i
                if not is_leading_low(d[i], j):
                    err = "error in row %d of m1"% i 
                    print("m1", m1)
                    print("table", row_tab)
                    raise ValueError(err)
        row_mask = bitrange_mask(m1.shape[1], m1.ncols)
        row_set = set(range(1, m1.nrows))
        for j in range(m1.shape[1], m1.ncols):
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
            
        # Check the computation of the rank of matrix m              
        if m.ncols > 8:
            continue
        rk = np.linalg.matrix_rank(m[:], tol=EPS)
        lb_rk = m.lb_rank()
        assert rk == 1 << lb_rk, (m1[:], rk, lb_rk)   

        # Finally, check the computation of the norm of matrix m              
        n2 = np.linalg.norm(m[:], ord = 2)
        log_n2 = -1 if n2 == 0 else 2.0 * log(n2) / log(2.0)
        ok = abs(log_n2 - m.lb_norm2()) < EPS
        if not ok:
            print("Matrix m =", m.reduce())
            print("log_2(squared norm2) = %f, expected = %.13f" %
                (m.lb_norm2(), log_n2 ))
            err = "Error in computing norm of matrix"
            raise ValueError(err)            


#####################################################################
# Test function qstate12_bit_matrix_t()
#####################################################################


def rand_bit_matrix(rows, cols):
    m = (1 << cols) - 1
    return [randint(0, m) for i in range(rows)]
    
def as_bit_array(m, cols):
    a = np.zeros( (len(m), cols), dtype = np.uint8)
    for i in range(len(m)):
        for j in range(cols):
            a[i,j] = (int(m[i]) >> j) & 1
    return a

            
def create_bitmatrices():
    """yield a bit matrix as a list of integers """
    for rows in range(8):
        for cols in range(8):
            m = rand_bit_matrix(rows, cols)
            yield m, cols
        
        
@pytest.mark.qstate
def test_bitmatrix_t(verbose = 0):
    """Test the transposition of a bit matrix"""
    for ntest, (m, cols) in enumerate(create_bitmatrices()):
        t = qstate12_bit_matrix_t(m, cols)
        m1 = as_bit_array(m, cols)
        #print(m1, "shp", m1.shape)
        t1 = as_bit_array(t, len(m))
        if verbose:
            print("Test %d: " % ntest)
            print("Transpose %s bit matrix, m =" % str(m1.shape))
            print(m1, "\nTransposed: %s\n%s" % (str(t1.shape), t1))
        assert m1.T.shape == t1.shape, (m1.shape, t1.shape)
        assert (m1.T == t1).all()



        