from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix, rand_qs_matrix
from mmgroup.structures.qs_matrix import prep_mul, flat_product

from mmgroup.tests.test_clifford.test_qs_matrix import compare_complex

from mmgroup.clifford12 import QState12, qstate12_product
from mmgroup.clifford12 import qstate12_prep_mul

#####################################################################
# 
#####################################################################




qs_matrix_data = [
   # [ (0,0, (1,0), []),  (0,0, (1,0), []), 0 ],
    [ (0,4, (0,0), [0b00_0101, 0b00_0011]),
        (0,4, (0,0), [0b00_0111, 0b_00_1001]), 2],
]


def create_m(rows, cols, factor, data):
    m = QStateMatrix(rows, cols, data)
    m.mul_scalar(*factor) 
    return m


def create_product_testvectors():
    """yield instances of class ``QStateMatrix`` for tests """
    for m1_data, m2_data, nqb in qs_matrix_data:
        m1 = create_m(*m1_data)
        m2 = create_m(*m2_data)
        yield m1, m2, nqb
    yield QStateMatrix(4, 4,  1) , QStateMatrix(4, 4,  1), 3 
    for i in range(20):
        nqb = randint(0,2)
        yield rand_qs_matrix(2, 0, 3), rand_qs_matrix(2, 0, 3), nqb
    for cols1 in range(2,6):
        for cols2 in range(2,6):
            for nqb in range(min(cols1, cols2)):
               for n in range(2):
                   m1 = rand_qs_matrix(0, cols1, cols1+3)
                   m1.mul_scalar(randint(-8, 8), randint(0,7))  
                   m2 = rand_qs_matrix(0, cols2, cols2+3)
                   m1.mul_scalar(randint(-8, 8), randint(0,7))  
                   yield m1, m2, nqb    
    # Sparse states           
    for cols1 in range(2,6):
        for cols2 in range(2,6):
            for nqb in range(min(cols1, cols2)):
               for r in range(3):
                   for n in range(2):
                       m1 = rand_qs_matrix(0, cols1,r)
                       m1.mul_scalar(randint(-8, 8), randint(0,7))  
                       m2 = rand_qs_matrix(0, cols2, r)
                       m1.mul_scalar(randint(-8, 8), randint(0,7))  
                       yield m1, m2, nqb    
                


def qs_complex_prod(a, b, nqb, nc):
    nb = nqb - nc
    assert nb >= 0
    a1 = a.reshape((-1, 1 << nb, 1 << nc))
    b1 = b.reshape((-1, 1 << nb, 1 << nc))
    prod = np.einsum("ikl,jkl->ijk", a1, b1)
    return np.ravel(prod)


def check_eq_cols(qs1, qs2, nqb):
    s = 0
    n = min(qs1.nrows, qs2.nrows)
    m1, m2 = qs1.data, qs2.data
    for i in range(n):
        s |= m1[i] ^ m2[i]
    for i in range(n, qs1.nrows):
        s |= m1[i]
    for i in range(n, qs2.nrows):
        s |= m2[i]
    s &= (1 << nqb) - 1
    assert s == 0
   

@pytest.mark.qstate
def test_qs_prep_mul(verbose = 0):
    """Basic test for class ``QStateMatrix`` 
    
    It tests the ``__str__ `` method, the conversion to a
    complex matrix and the echelonization and the reduction of 
    an instance of class ``QStateMatrix``.
    """
    for ntest, (m1, m2, nqb) in enumerate(create_product_testvectors()):
        m1a, m2a = prep_mul(m1, m2, nqb)
        if verbose:
            print("TEST %s" % (ntest+1))
            print(m1)
            print(m2)
            print("nqb =", nqb)
            print("Output states")
            print(m1a)
            print(m2a)
        c1 = qs_complex_prod(m1.complex(), m2.complex(), nqb, nqb)
        c1a = qs_complex_prod(m1a.complex(), m2a.complex(), nqb, nqb)
        try:
            err = "Function prep_mul() has failed"
            compare_complex(c1, c1a, err)
        except ValueError:
            print("\n" + err +"\nInput states:")
            print(m1); print(m2);  
            print("nqb =", nqb)      
            print("\nOutput states:")
            print(m1a); print(m2a);  
            raise      
        check_eq_cols(m1a, m2a, nqb)


@pytest.mark.qstate1
def test_qs_product(verbose = 0):
    for ntest, (m1, m2, nqb) in enumerate(create_product_testvectors()):
        nc = randint(0, nqb)
        if verbose:
            print("TEST %s" % (ntest+1))
            print(m1)
            print(m2)
            print("nqb =", nqb, " nc =", nc)
            print("Output states")
        try:
            err = "Execution of function flat_product() has failed"
            #m3 = flat_product(m1, m2, nqb, nc)
            qm1, qm2 = QState12(m1), QState12(m2)
            qstate12_product(qm1, qm2, nqb, nc)
            m3, m4 = QStateMatrix(qm1), QStateMatrix(qm2)
        except ValueError:
            print("\n" + err +"\nInput states:")
            print(m1); print(m2);  
            print("nqb =", nqb, ", nc =", nc)  
            pm1, pm2 = QState12(m1), QState12(m2)
            row_pos =qstate12_prep_mul(pm1, pm2, nqb)            
            print("prep_mul, row pos = ", row_pos) 
            print(QStateMatrix(pm1)); print(QStateMatrix(pm2))            
            print("\nOutput states:")
            print(QStateMatrix(qm1)); print(QStateMatrix(qm2)); 
            print(qm1.nrows, qm1.ncols, qm2.nrows, qm2.ncols)
            raise      
        c1, c2, c3 = m1.complex(), m2.complex(), m3.complex()
        c3_ref = qs_complex_prod(c1, c2, nqb, nc)
        try:
            err = "Function flat_product() has failed"
            compare_complex(c3_ref, c3, err)
        except ValueError:
            print("\n" + err +"\nInput states:")
            print(m1); print(m2);  
            print("nqb =", nqb, ", nc =", nc)      
            print("\nOutput states:")
            print(m3);  
            raise      
        if verbose:
            print("Output state")
            print(m3)

        