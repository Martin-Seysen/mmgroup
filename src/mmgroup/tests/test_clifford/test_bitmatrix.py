from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from math import log
from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.clifford12 import bitmatrix64_t
from mmgroup.clifford12 import bitmatrix64_echelon_h
from mmgroup.clifford12 import bitmatrix64_echelon_l
from mmgroup.clifford12 import bitmatrix64_cap_h
from mmgroup.clifford12 import bitmatrix64_error_pool

 


#####################################################################
# Test function bitmatrix64_t()
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
        t = bitmatrix64_t(m, cols)
        m1 = as_bit_array(m, cols)
        #print(m1, "shp", m1.shape)
        t1 = as_bit_array(t, len(m))
        if verbose:
            print("Test %d: " % ntest)
            print("Transpose %s bit matrix, m =" % str(m1.shape))
            print(m1, "\nTransposed: %s\n%s" % (str(t1.shape), t1))
        assert m1.T.shape == t1.shape, (m1.shape, t1.shape)
        assert (m1.T == t1).all()


#####################################################################
# Test functions bitmatrix64_echelon_h(), bitmatrix64_echelon_l()
#####################################################################



def create_echelon_matrices():
    """yield a bit matrix as a list of integers """
    for rows in range(8):
        for _ in range(8):
            m = rand_bit_matrix(rows, 64)
            i, j = randint(0,63), randint(0,63)
            yield m, min(i,j), max(i,j)

 
def ref_echelon(m, j0 = 64, n = 64, high = True):
    m = [int(x) for x in m]
    row = 0
    if high:
        assert n <= j0
        range_ = range(j0-1, j0-n-1, -1)
        mask =  ((1 << n) - 1) << (j0 - n)
    else:
        range_ = range(j0, j0+n)
        mask =  ((1 << n) - 1) << (j0)        
    for col in range_:
        col_mask = 1 << col
        for i in range(len(m)-1, row-1, -1):
            if m[i] & col_mask:
                for i1 in range(i-1, -1, -1):
                    if  m[i1] & col_mask: 
                        m[i1] ^= m[i]
                
                m[i], m[row] = m[row], m[i]
                row += 1
                break
    return row, m


@pytest.mark.qstate
def test_bitmatrix_echelon(verbose = 0):
    for ntest, (m, imin, imax) in enumerate(create_echelon_matrices()):
         m1h = np.array(m, dtype = np.uint64, copy = True)
         j0h, nh = imax + 1, imax - imin
         lenh = bitmatrix64_echelon_h(m1h, len(m1h), j0h, nh)
         m1h = list(m1h)
         lenh_ref, m1h_ref =  ref_echelon(m, j0h, nh)
         ok_h = m1h == m1h_ref and lenh == lenh_ref

         m1l = np.array(m, dtype = np.uint64, copy = True)
         j0l, nl = imin, imax - imin
         lenl = bitmatrix64_echelon_l(m1l, len(m1l), j0l, nl)
         m1l = list(m1l)
         lenl_ref, m1l_ref =  ref_echelon(m, j0l, nl, False)
         ok_l = m1l == m1l_ref and lenl == lenl_ref

         ok = ok_h and ok_l
         if verbose or not ok:
              print("Test ", ntest)
              print("m =", [hex(x) for x in m], j0h, nh)
              print("echelon h", [hex(x) for x in m1h], lenh)
              print("expected", [hex(x) for x in m1h_ref], lenh_ref)
              if not ok_h:
                  err = "Error in function bitmatrix64_echelon_h"
                  raise ValueError(err)
              print("m =", [hex(x) for x in m], j0l, nl)
              print("echelon l", [hex(x) for x in m1l], lenl)
              print("expected", [hex(x) for x in m1l_ref], lenl_ref)
              if not ok_l:
                  err = "Error in function bitmatrix64_echelon_l"
                  raise ValueError(err)


#####################################################################
# Test functions bitmatrix64_cap_h()
#####################################################################

def create_cap_matrices():
    """yield a bit matrix as a list of integers """
    test_matrices = [
        ([3,7,11], [3,1], 0,4)
    ]
    for t in test_matrices:
        yield t
    for rows1 in range(5):
        for rows2 in range(5):
            m1 = rand_bit_matrix(rows1, 5)
            m2 = rand_bit_matrix(rows2, 5)
            i, j = randint(0,6), randint(0,6)
            yield m1, m2, min(i,j), max(i,j)
    for rows1 in range(1,65,7):
        for rows2 in range(1,65,7):
            m1 = rand_bit_matrix(rows1, 64)
            m2 = rand_bit_matrix(rows2, 64)
            i, j = randint(0,63), randint(0,63)
            yield m1, m2, min(i,j), max(i,j)




def basis_to_set(b):
    if len(b) == 0:
        return set([0])
    s0 = basis_to_set(b[1:])
    b0 = b[0]
    return s0 | set((x ^ b0 for x in s0))

def as_set(m, mask):
    m = [int(x) & int(mask) for x in m]
    return basis_to_set(m)
        


@pytest.mark.qstate
def test_bitmatrix_cap(verbose = 0):
    for ntest, (m1, m2, imin, imax) in enumerate(create_cap_matrices()):
         if verbose:
             print("Test", ntest+1, "imin =", imin, "imax =", imax )
             print("m1 =", [hex(x) for x in m1])
             print("m2 =", [hex(x) for x in m2])
         m1h = np.array(m1, dtype = np.uint64, copy = True)
         m2h = np.array(m2, dtype = np.uint64, copy = True)
         j0h, nh = imax + 1, imax - imin + 1 
         l1h, l2h = bitmatrix64_cap_h(m1h, m2h, j0h, nh)
         m1h, m2h = list(m1h), list (m2h)
         if verbose:
             print("Non intersecting parts of m1 and m2")
             print("out1 =", [hex(x) for x in m1h[:l1h]])
             print("out2 =", [hex(x) for x in m2h[:l2h]])
             mask = (1 << (imax + 1)) - (1 << imin)
             print("Intersecting parts (mask = %s):" % hex(mask))
             print("out1 =", [hex(x) for x in m1h[l1h:]])
             print("out2 =", [hex(x) for x in m2h[l2h:]])
         if verbose > 1:    
             print("Intermediate results") 
             pool = np.zeros(20, dtype = np.uint64)
             bitmatrix64_error_pool(pool, 20)
             for i in range(20): print(i,  hex(pool[i]))

         assert len(m1) == len(m1h) and len(m2) == len(m2h)
         assert ref_echelon(m1) == ref_echelon(m1h)
         assert ref_echelon(m2) == ref_echelon(m2h)
         mask = (1 << (imax + 1)) - (1 << imin)
         l1, l2 = len(m1), len(m2)
         while l1 and int(m1h[l1-1]) & mask == 0: 
             l1 -= 1
         while l2 and int(m2h[l2-1]) & mask == 0: 
             l2 -= 1
         assert l1-l1h == l2-l2h, (l1, l1h, l2, l2h)
         if len(m1) < 5 and len(m2) < 5:
             set1 = as_set(m1h, mask)
             set2 = as_set(m2h, mask)
             cap  = set1 & set2
             set_cap1 = as_set(m1h[l1h:], mask)
             set_cap2 = as_set(m2h[l2h:], mask)
             assert cap == set_cap1 , (set1, set2, cap, set_cap1)
             assert cap == set_cap2 , (set1, set2, cap, set_cap2)
             if verbose:
                 print("Intersection testd successfully")
            

