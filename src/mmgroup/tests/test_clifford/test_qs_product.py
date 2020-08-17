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
    [ (0,0, (1,0), []),  (0,0, (1,0), []), 0 ],
    [ (0,4, (0,0), [0b00_0101, 0b00_0011]),
        (0,4, (0,0), [0b00_0111, 0b_00_1001]), 2],
]


def create_m(rows, cols, factor, data):
    """Create ``QStateMatrix`` object from input data"""
    m = QStateMatrix(rows, cols, data)
    m.mul_scalar(*factor) 
    return m



def rand_mul_scalar(m):
    m.mul_scalar(randint(-16, 15), randint(-8,7))

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
                   rand_mul_scalar(m1)  
                   m2 = rand_qs_matrix(0, cols2, cols2+3)
                   rand_mul_scalar(m2)  
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
    """Complex equivalent to function ``qstate12_product``

    ``a`` and ``b`` must be complex numpy arrays of dimension 2.
    The function returns ``qstate12_product(a, b, nqb, nc)``
    as a complex umpy array of dimension 2.
    """
    nb = nqb - nc
    assert nb >= 0
    a1 = a.reshape((-1, 1 << nb, 1 << nc))
    b1 = b.reshape((-1, 1 << nb, 1 << nc))
    prod = np.einsum("ikl,jkl->ijk", a1, b1)
    return np.ravel(prod)


def check_eq_cols(qs1, qs2, nqb):
    """Check if the first ``nqb`` columns of states are equal
    
    Raise ValueError if this is not the case. 
    ``qs1`` and ``qs2`` must be instances of class ``QStateMatrix``.
    """
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
    """Test function ``qstate12_prep_mul``. """
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


@pytest.mark.qstate
def test_qs_product(verbose = 0):
    """Test function ``qstate12_product``. """
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
            m3 = flat_product(m1, m2, nqb, nc)
        except ValueError:
            print("\n" + err +"\nInput states:")
            print(m1); print(m2);  
            print("nqb =", nqb, ", nc =", nc)  
            pm1, pm2 = QState12(m1), QState12(m2)
            row_pos =qstate12_prep_mul(pm1, pm2, nqb)            
            print("prep_mul, row pos = ", row_pos) 
            print(QStateMatrix(pm1)); print(QStateMatrix(pm2))            
            print("\nOutput states:")
            qm1, qm2 = QState12(m1), QState12(m2)
            qstate12_product(qm1, qm2, nqb, nc)
            print(QStateMatrix(qm1)); print(QStateMatrix(qm2)); 
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
            
            
            
def create_matmul_testvectors():
    """Yield pairs of matrices for testing matrix multiplication

    These matrices are instances of class ``QStateMatrix``.
    """
    for i in range(2,6):
        for j in range(2,6):
            for k in range(2, 6):
                for r in range(3, 12, 2):
                    for n in range(1):
                        m1 = rand_qs_matrix(i, j, r)
                        rand_mul_scalar(m1)  
                        m2 = rand_qs_matrix(j, k, r)
                        rand_mul_scalar(m2)  
                        yield m1, m2   

@pytest.mark.qstate
def test_matmul(verbose = 0):
    """Test matrix multiplication. """
    for ntest,(m1, m2) in enumerate(create_matmul_testvectors()):
        if verbose:
            print("TEST %s QStateMatrix multiplication" % (ntest+1))
            print(m1)
            print(m2)
        m3 = m1 @ m2
        c1, c2 = m1.complex(), m2.complex()
        c3 = m3.complex()
        c3_ref = c1 @ c2        
        err = "Matrix multiplcation has failed"
        compare_complex(c3_ref, c3, err)
        if verbose:
            print("Product")
            print(m3)



def large_matmul_testvectors():
    """Yield test data for large matrix multiplication

    yield pairs ``m1``, ``m2`` of 2**12 times 2**12
    and also pairs random slices of about 40 elements
    in range(2**12)
    """
    nn = 12
    for r in range(10, 27):
        m1 = rand_qs_matrix(nn, nn, r)
        rand_mul_scalar(m1)  
        m2 = rand_qs_matrix(nn, nn, r)
        rand_mul_scalar(m2)
        r1 = slice(randint(0,100), 1<<nn, 2*randint(150,160)+1)       
        r2 = slice(randint(0,100), 1<<nn, 2*randint(150,160)+1)       
        yield m1, m2, r1, r2   


@pytest.mark.qstate
def test_large_matmul(verbose = 0):
    """Testlarge  matrix multiplication. 
    
    For the mmgroup project wee have to multiply marices
    of type QStateMatrix with shape (12,12). So we'd
    better test this operation here!
    
    We generate random matrices m1, m2. We test the 
    multiplication of a random  r times 4096 submatrix
    of m1 with a random 4096 times r submatrix of m2.
    Here r is about 14.
    """
    for ntest,(m1, m2, r1, r2) in enumerate(large_matmul_testvectors()):
        if verbose:
            print("TEST %s QStateMatrix multiplication" % (ntest+1))
            print(m1.echelon())
            print(m2.echelon())
        m3 = m1 @ m2
        if verbose:
            print("Product")
            print(m3)
        c1, c2 = m1[r1,:], m2[:,r2]
        c3 = m3[r1, r2]
        c3_ref = c1 @ c2        
        err = "Matrix multiplcation has failed"
        compare_complex(c3_ref, c3, err)




def create_mul_testvectors():
    """Yield pairs of matrices for testing matrix multiplication

    These matrices are instances of class ``QStateMatrix``.
    """
    for i in range(2,6):
        for j in range(2,6):
            for r in range(2, 12, 2):
                for n in range(2):
                    m1 = rand_qs_matrix(i, j, r)
                    rand_mul_scalar(m1)  
                    m2 = rand_qs_matrix(i, j, r)
                    rand_mul_scalar(m2)  
                    yield m1, m2   



@pytest.mark.qstate
def test_mul(verbose = 0):
    """Test elementwise and scalar multiplication for `QStateMatrix``"""
    for ntest, (m1, m2) in enumerate(create_mul_testvectors()):
        if verbose:
            print("TEST %s Matrix multiplication of entries" % (ntest+1))
            print(m1.echelon())
            print(m2.echelon())
        m3 = (m1 * m2)
        if verbose:
            print("Product")
            print(m3)
        c3 = m3[:,:]
        c3_ref = m1[:,:] * m2[:,:]
        err = "Matrix multiplcation of entries has failed"
        compare_complex(c3_ref, c3, err)
        err = "Matrix negation has failed"
        compare_complex(-(m1[:,:]), (-m1)[:,:], err)
        err = "Matrix scalar multiplication failed"
        f = (1+1j)**randint(-8,7) * 2**(0.5*randint(-10,10))
        compare_complex(f*(m1[:,:]), (f*m1)[:,:], err)
        compare_complex(f*(m1[:,:]), (m1*f)[:,:], err)
        compare_complex((m1[:,:]/f), (m1/f)[:,:], err)
        if ntest < 20:
            compare_complex(0*(m1[:,:]), (0*m1)[:,:], err)
   

def create_conj_testvectors():
    """Yield pairs of matrices for testing matrix multiplication

    These matrices are instances of class ``QStateMatrix``.
    """
    for i in range(2,6):
        for j in range(2,6):
            for r in range(2, 12, 2):
                for n in range(2):
                    m1 = rand_qs_matrix(i, j, r)
                    rand_mul_scalar(m1)  
                    yield m1  

   
@pytest.mark.qstate
def test_conj(verbose = 0):
    """Test conjugation and transposition for `QStateMatrix``"""
    for ntest, m1 in enumerate(create_conj_testvectors()):
        if verbose:
            print("TEST %s Matrix multiplication of entries" % (ntest+1))
            print(m1.echelon())
        if verbose:
            print("Conjugate")
            print(m1.conj())
        err = "Matrix conjugation failed" 
        compare_complex(m1.conj()[:,:], m1[:,:].conj(), err)
        if verbose:
            print("Transpose")
            print(m1.T)
        err = "Matrix transpostion failed" 
        compare_complex(m1.T[:,:], m1[:,:].T, err)


     