from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from math import log
from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.bitfunctions import bit_mat_mul, bit_mat_inverse

from mmgroup.clifford12 import bitmatrix64_t
from mmgroup.clifford12 import bitmatrix64_echelon_h
from mmgroup.clifford12 import bitmatrix64_echelon_l
from mmgroup.clifford12 import bitmatrix64_cap_h
from mmgroup.clifford12 import bitmatrix64_mul
from mmgroup.clifford12 import bitmatrix64_inv
from mmgroup.clifford12 import bitmatrix64_error_pool
from mmgroup.clifford12 import bitmatrix64_find_low_bit


 


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
            

#####################################################################
# Test functions bitmatrix64_mul() and bitmatrix64_inv()
#####################################################################

# This tests also function bitmatrix64_mask_rows()
# and  bitmatrix64_add_diag().


def create_mul_inv_matrices():
    """yield a bit matrix as a list of integers """
    for i in [1, 3, 5, 7, 17, 24, 31, 32]:
        m2i = None
        for _ in range(1000):
            m2 = rand_bit_matrix(i, i)
            try:
                m2i = bit_mat_inverse(m2)
                break
            except ZeroDivisionError:       
                pass
        for j in [1,2,7,32, 63, 64]:
            m1 = rand_bit_matrix(j, i)
            yield m1, m2, m2i
            m2i = None
             
@pytest.mark.qstate
def test_bitmatrix_mul_inv(verbose = 0):
    for ntest, (m1, m2, m2i) in enumerate(create_mul_inv_matrices()):
         #print(m1, m2, m2i)
         m1a = np.array(m1, dtype = np.uint64, copy = True)
         m2a = np.array(m2, dtype = np.uint64, copy = True)
         m3 =  bit_mat_mul(m1, m2)
         m3a = bitmatrix64_mul(m1a, m2a)
         assert list(m3a) == m3
         if m2i is not None:
            m2ia = bitmatrix64_inv(m2)
            assert list(m2ia) == m2i 


   


#####################################################################
# Test function  bitmatrix64_find_low_bit()
#####################################################################


def int_to_bitmartix(x, excess = 0, bmax = 0):
    l = []
    assert x >= 0
    while x:
        l.append(x & 0xffffffffffffffff)
        x >>= 64
    l += [0] * excess
    nmax = (bmax + 63) >> 6
    if len(l) < nmax:
        l += [0] * (nmax - len(l))
    return np.array(l, dtype = np.uint64)


bvalues = [0,1,2, 62, 63, 64, 65, 87, 127, 128, 129]


def find_low_bit_testdata():
    lv = len(bvalues)
    for i in range(lv):
        for j in range(i, lv):
            for k in range(j, lv):
                t = tuple([bvalues[n] for n in (i,j,k)]) 
                yield t + (1,)
                r = randint(1, 1 << 65) | 1
                yield t + (r,)

@pytest.mark.qstate
def test_bitmatrix64_find_low_bit(verbose = 0):
    if verbose:
        print("\nTesting function bitmatrix64_find_low_bit")
    for n, (imin, b, imax, factor) in enumerate(find_low_bit_testdata()):
        b_int = factor << b
        b_int |= randint(0, (1 << imin) - 1)
        bm = int_to_bitmartix(b_int, bmax = imax)
        if verbose:
           print("Test %d:" % n, imin, b, imax, 
               ", factor =", hex(factor))
        b1 = bitmatrix64_find_low_bit(bm, imin, imax)
        assert b1 == b, (imin, b, imax, b1)
 
         
    
    

