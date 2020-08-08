from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix, rand_qs_matrix



#####################################################################
# Convert matrix of type QStateMatrix to a complex matrix
#####################################################################

def eval_A_Q(ncols, data, v):
    """Auxiliary function for function ``slow_complex``.
    
    Evaluate parts ``A``and ``Q`` of a matrix for a vector ``v``.
    
    Given a state matrix ``qs`` type ``QStateMatrix``, parameter
    ``ncols`` should be equal to ``qs.rows + qs.cols``, and parameter
    ``data`` should be equal to ``qs.data``.  Then parts ``A``and ``Q``
    of ``qs`` are evaluated for a bot vector ``v`` of bit length
    ``qs.nrows``, where . ``qs.nrows`` is equal to ``len(data)``.
    
    ``A`` is the ``qs.nrows`` times ``qs.ncols`` bit matrix that
    consistts of the first `qs.ncols`` columns of bit matrix ``data``.
    ``Q`` is the ``qs.nrows`` times ``qs.nrows`` bit matrix that
    consistts of the following `qs.rows`` columns of matrix ``data``.
    
    The function returns a pair ``(a, q)``, where ``a`` is the matrix
    product :math:``v \cdot A``.  ``q`` is defined as follows:
    Let ``R`` be the bit matrix of the same shape as ``Q``
    with ``R[i,j] = v[i] * v[j] * Q[i,j]``. Then ``q`` is the
    sum of all bits in ``R`` modulo 4
    
    We pass the entries ``ncols,  data`` of the state matrix
    ``qs`` instead of ``qs`` itself, for reasons of efficiency.
    """
    # Clear irrelevant bits of v 
    v &= (1 << len(data)) - 1  
    # Sum up  :math:``v \cdot A`` in ``a``, the off-diagonal elements 
    # of ``R`` in ``q``, and the diagonal elements of ``R`` in ``diag``. 
    # Since ``R`` is symmetric, we XOR  the lower triagular off-diagonal 
    # elements of ``R`` in ``q`` and double the result.
    a = q = diag = 0
    for i, d in enumerate(data):
        # ``d`` = row ``i`` of ``(A,Q)`` if bit ``i`` of ``v``
        # is set an ``d = 0`` otherwise.        
        d &=  -((v >> i) & 1)
        # XOR  (the ``A`` part of)  ``d`` to ``A``.
        a ^= d
        # Next let ``d`` be row ``i`` of the ``R``.
        d = (d >> ncols) & v
        # Add diagonal bit of row ``i`` of ``R`` to ``diag``
        diag += (d >> i) & 1
        # XOR lower triangular part of that row to ``q``
        q ^= d &  ((1 << i) - 1)
    # Calculate the bit parity of ``q`` in bit 0 of ``q``
    for sh in [32, 16, 8, 4, 2, 1]:
        q ^= q >> sh 
    # Mask out the relevant bits of the result ``a``
    a &= (1 << ncols) - 1
    # Caclulate ``q`` as the doubled sum of the bits the lower 
    # triangular part of the bits of ``R``, plus the sum
    # of the diagonal bits of ``R``. 
    q = (2*q + diag) & 3
    return a, q   




def slow_complex(m):
    """Convert matrix ``m`` of type ``QStateMatrix`` to complex
    
    The function returns a complex ``numpy`` array of appropriate
    shape containing the matrix represented by the input ``m``.
    
    This is a pathon implementation for testing the fast
    method ``m.complex()``    
    """
    # Compute the complex scalar factor ``f``  from ``m.factor``.
    e, phi = m.factor 
    f = 2.0**(0.5*(e - phi)) * (1+1j)**phi
    # Put ``phases[i] = sqrt(-1)**i * f,  0 <= i < 4.``
    phases = [f * 1j**i for i in range(4)]
    # Create comples zero array ``a`` od appropriate length.
    a = np.zeros(1 << m.ncols, dtype = np.complex)
    # For all vectors of length ``m.nrows`` with bit 0
    # of ``v`` equal to 1 put ``index =  v * A``, and  let ``q``
    # be the sum of the matrix ``R`` derivend from ``v`` and
    # ``Q`` as in function ``eval_A_Q()``.
    # Then add ``phases[q mod 4]`` to entry ``a[index]``
    ncols, data = m.ncols, m.data
    for v in range(1, 1 << len(data), 2):
        index, q = eval_A_Q(ncols, data, v)
        a[index] += phases[q]
    # Change shape of array ``a`` as required
    return a.reshape((1 << m.rows, 1 << m.cols))


#####################################################################
# Check that states are approcimately equal
#####################################################################

EPS = 1.0e-8

def check_eq_qstate(m1, m2, text = None):
    """Raise ValueError unless matrices ``m1, m2``  are equal.
    
    Inputs ``m1, m2`` must be matrices of type ``QStateMatrix``.
    Both matrices are converted to complex matrices (with a
    method of class ``QStateMatrix`` that does not reduce the
    matrices) and testes for (approximate) equality.   

    ``text`` is an optional text to be displayed in case of
    an error.    
    """
    c1, c2 = m1.complex_unreduced(), m2.complex_unreduced()
    eq = np.amax(abs(c1 - c2)) < EPS
    if not eq:
        err = "Error in comparing instances of QStateMatrix"
        print("\n" + (text if text else err) + "\n")
        print("m1 = ", m1)
        print("m2 = ", m2)
        if max(m1.shape + m2.shape) < 4:
            print("complex m1 = \n", c1)
            print("complex m2 = \n", c2)
        raise ValueError(err)
         

def compare_complex(c1, c2, text, qstate = None):
    """Check that two complex matrices are equal
    
    ``c1`` is the expectd and ``c2`` is the obtained complex matrix.
    ``qstate`` is an optiona state representing the matrix.
    """
    diff = np.amax(abs(c1 - c2))
    if diff > EPS:
        err = text if text else "Error in comparing complex matrices"
        print("\n" + err + "\n")
        print("qstate  = ", qstate)
        if max(c1.shape + c2.shape) < 6:
            print("expected\n", c1)
            print("obtained\n\n", c2)
        else:
            print("Shape of expected matrix:", c1.shape)
            print("Shape of obtained matrix:", c2.shape)
        print("Max. absolute error:", diff)
        raise ValueError(err)


def check_complex(m):
    """Check conversion of ``QStateMatrix`` to complex matrix
    
    It test the methods for conversion of a matrix to a 
    complex matrix 8with and without previous reduction) against
    the python function ````QStateMatrix()``.    
    """
    c0 = slow_complex(m)
    c1 = m.complex_unreduced()
    err = "Complex matrix from unreduced state"
    compare_complex(c0, c1, err, m)
    c2 = m.complex()
    err = "Complex matrix from reduced state"
    compare_complex(c0, c2, err, m)


#####################################################################
# display some matrices
#####################################################################

qs_matrix_data = [
     (0,0, (1,0), []),
     (0,0, (-3,6), 0),
     (3,0, (-3,5), 5),
     (0,6, (9,2), 25),
     (3,1, (1,0), []),
     (2, 2, (0,0), [0b110_10_01, 0b101_01_11, 0b011_01_00]),
     (2, 0, (0,0), [0b000_10, 0b010_01, 0b100_01]),
    ]


def create_display_testvectors():
    """yield instances of class ``QStateMatrix`` for tests """
    for rows, cols, factor, data in qs_matrix_data:
        m = QStateMatrix(rows, cols, data)
        m.mul_scalar(*factor) 
        yield m
    yield QStateMatrix(4, 4,  1) 
    for i in range(20):
        yield rand_qs_matrix(2, 0, 3)
    for i in range(20):
        yield rand_qs_matrix(1, 0, 2)
    for rows in range(6):
        for cols in range(5):
            for nr in [1,2] + list(range(rows+cols, rows+cols+3)):
                m = rand_qs_matrix(rows, cols, nr)
                m.mul_scalar(randint(-8, 8), randint(0,7))  
                yield m                
                


@pytest.mark.qstate
def test_qs_matrix(verbose = 0):
    """Basic test for class ``QStateMatrix`` 
    
    It tests the ``__str__ `` method, the conversion to a
    complex matrix and the reduction of an instance of class
    ``QStateMatrix``.
    """
    for ntest, m in enumerate(create_display_testvectors()):
        if verbose or ntest < len(qs_matrix_data):
            print("TEST %s" % (ntest+1))
            print(m)
            m2 = QStateMatrix(m)
            print("Echelon")
            print(m2.echelon())
            print(m2)
        check_complex(m)
        m2 = QStateMatrix(m)
        if verbose: 
            print("Echelonize....")
            m2.echelon()
            print(m2)
            check_eq_qstate(m,m2)
            m2 = QStateMatrix(m)
        if verbose: print("Reducing....")
        m2.reduce()
        if verbose: print(m2)
        check_eq_qstate(m,m2)








        