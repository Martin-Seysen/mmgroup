from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix
from mmgroup.structures.qs_matrix import qs_rand_real_matrix
from mmgroup.structures.qs_matrix import qs_unit_matrix
from mmgroup.structures.qs_matrix import qs_column_monomial_matrix
from mmgroup.structures.qs_matrix import qs_row_monomial_matrix
from mmgroup.structures.qs_matrix import qs_from_signs
from mmgroup.structures.qs_matrix import FORMAT_REDUCED
from mmgroup.structures.qs_matrix import qstate12_mul_matrix_mod3


#####################################################################
# Display and check some matrices
#####################################################################




qs_matrix_data = [
     (0, 0, (1,0), []),
     (0, 4, (1,4), [0, 8, 4, 2, 1]),
     (0, 5, (1,4), [0, 16, 8, 4, 2, 1]),
     (0, 6, (1,4), [0, 32, 16, 8, 4, 2]),
     (0, 6, (2,0), [0, 32, 16, 8, 4, 2, 1]),
    ]


def create_testvectors():
    """yield instances of class ``QStateMatrix`` for tests """
    for rows, cols, factor, data in qs_matrix_data:
        m = QStateMatrix(rows, cols, data)
        m.mul_scalar(*factor) 
        yield m
    for i in (2,3,4,5, 6, 7):
        yield qs_unit_matrix(i) 
    for i in range(20):
        yield qs_rand_real_matrix(2, 0, 3)
    for i in range(20):
        yield qs_rand_real_matrix(1, 0, 2)
    for rows in range(6):
        for cols in range(5):
            for nr in [1,2] + list(range(rows+cols, rows+cols+3)):
                for _ in range(200):
                    m = qs_rand_real_matrix(rows, cols, nr)
                    m.mul_scalar(randint(-8, 8), randint(0,7) & 4)  
                    yield m                
                
#####################################################################
# Test conversion to array of signs
#####################################################################



def modify_sign_bitmap(bm, ncols):
    bm1 = np.copy(bm)
    r = randint(0,  (2 << ncols) - 1)
    index, bit = divmod(r, 64)
    bm1[index] = int(bm1[index]) ^ (1 << bit)
    return bm1



@pytest.mark.qstate
def test_qs_matrix(verbose = 0):
    """Basic test for class ``QStateMatrix`` 
    
    It tests the ``__str__ `` method, the conversion to a
    complex matrix and the echelonization and the reduction of 
    an instance of class ``QStateMatrix``.
    """
    FORMAT_REDUCED = False
    display_len =  len(qs_matrix_data)
    for ntest, m in enumerate(create_testvectors()):
        if verbose:
            print("TEST %s" % (ntest+1))
            print(m.copy().reshape((sum(m.shape), 0)))
        ac = m.complex()
        a = ac.real.ravel()
        fr, fphi = m.factor
        if 0 <= fr < 62 and fr & 1 == 0 and fphi & 3 == 0:
            assert (ac == m.int32()).all()
        b = [0 if x == 0 else 1 + ((x < 0) << 1) for x in a]
        bm = m.to_signs()
        blist = [((int(bm[i >> 5])) >> ((i & 31) << 1)) & 3
            for i in range(len(a)) ]
        if b != blist:
            print("TEST %s" % (ntest+1))
            print(m.copy().reshape((sum(m.shape), 0)))
            print("expected:", b)
            print("obtained:", blist)
            print("delta   :", [x ^ y for x, y in zip(b,blist)])
            err = "Error in computing array of signs"
            raise ValueError(err)
        assert m.compare_signs(bm)
        for i in range(2):
            bm1 = modify_sign_bitmap(bm, sum(m.shape))
            assert m.compare_signs(bm1) == False


        abs_m = m.maxabs()
        m1 = m / abs_m if abs_m > 0 else m * 1
        m1 = m1.reshape((-1, 0))
        #m1 = abs(m1)
        try:
            m_from_signs = qs_from_signs(bm, sum(m.shape))
        except:
            print(m_from_signs)
            #raise
            print("WTF")
        try:
            eq = m_from_signs == m1
        except:
            print("original      m", m1)
            print("reconstructed m", m_from_signs)
            raise
        if not eq:
            print("original      m", m1)
            print("reconstructed m", m_from_signs)
            if sum(m.shape) <= 8: print("signs", bm)
            err = "Reconstruction of state vector from signs failed"
            raise ValueError(err)


#####################################################################
# Test function qstate12_mul_matrix_mod3
#####################################################################

def create_rational_testmatrices():
    """yield instances of class ``QStateMatrix`` for tests """
    for rows, cols, factor, data in qs_matrix_data:
        if not 0 <= factor[0] < 62:
            continue
        m = QStateMatrix(rows, cols, data)
        factor = (factor[0] & -2, factor[1] & 4)
        m.mul_scalar(*factor) 
        yield m
    for i in (2,3,4,5, 6, 7):
        yield qs_unit_matrix(i) 
    for i in range(20):
        yield qs_rand_real_matrix(2, 0, 3)
    for i in range(20):
        yield qs_rand_real_matrix(1, 0, 2)
    for rows in range(6):
        for cols in range(5):
            for nr in [1,2] + list(range(rows+cols, rows+cols+3)):
                for _ in range(5):
                    m = qs_rand_real_matrix(rows, cols, nr)
                    m.mul_scalar(randint(-8, 8) & -2, randint(0,7) & 4)  
                    yield m  

MAX_UINT64 = 0xffffffffffffffff

def create_test_matrix_mod3(n):
    a = [randint(0, MAX_UINT64) for i in range(1 << n)]
    w = MAX_UINT64
    return a, w



def scalar_prod_mod3_uint64(a, b):
    s = sum(((a >> i) & 3) * ((b >> i) & 3) for i in range(0, 64, 2))
    return s % 3

def mul_a_w(a, w, max_check = 3):
    w_and = (w ^ (w >> 1)) & 0x5555555555555555
    w_and |= w_and << 1
    w_xor = w & 0xaaaaaaaaaaaaaaaa
    w_xor |= w_xor >> 1
    aw = [((x ^ w_xor) & w_and) % 3 for x in a]
    for i in range(min(max_check, len(a))):
        assert aw[i] == scalar_prod_mod3_uint64(a[i], w)
    return np.array(aw, dtype = np.int32)


def mul_m_a_w(m, a, w):
    f, phi = m.factor
    assert phi & 3 == 0 and f & 1 == 0
    m1 = m.copy().mul_scalar(-(f & -16) + 16)
    return m1.int32().ravel() @ mul_a_w(a, w)
          
@pytest.mark.qstate
def test_qstate12_mul_matrix_mod3(verbose = 0):
    for ntest, m in enumerate(create_rational_testmatrices()):
        a, w = create_test_matrix_mod3(m.ncols)
        if verbose:
            print("TEST %s" % (ntest+1))
            print(m.copy().reshape((sum(m.shape), 0)))
            if m.ncols < 4:
                print("a =", [hex(x) for x in a])
                print("w = ", hex(w))
        ref_prod = mul_m_a_w(m, a, w)
        if verbose:
            print("Product mod 3 =", ref_prod)
        prod = qstate12_mul_matrix_mod3(m, a, w)
        if prod != ref_prod:
            print("Product obtained =", prod)
            err = "Error in function qstate12_mul_matrix_mod3"
            raise ValueError(err)
