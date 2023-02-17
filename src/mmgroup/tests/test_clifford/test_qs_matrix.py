from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix
from mmgroup.structures.qs_matrix import qs_unit_matrix, qs_rand_matrix
from mmgroup.structures.qs_matrix import qs_column_monomial_matrix
from mmgroup.structures.qs_matrix import qs_row_monomial_matrix
from mmgroup.structures.qs_matrix import FORMAT_REDUCED

#####################################################################
# Convert matrix of type QStateMatrix to a complex matrix
#####################################################################

def eval_A_Q(ncols, data, v):
    r"""Auxiliary function for function ``slow_complex``.
    
    Evaluate parts ``A``and ``Q`` of a matrix for a vector ``v``.
    
    Given a state matrix ``qs`` type ``QStateMatrix``, parameter
    ``ncols`` should be equal to ``qs.rows + qs.cols``, and parameter
    ``data`` should be equal to ``qs.data``.  Then parts ``A``and ``Q``
    of ``qs`` are evaluated for a bit vector ``v`` of bit length
    ``qs.nrows``, where ``qs.nrows`` is equal to ``len(data)``.
    
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
    a = np.zeros(1 << m.ncols, dtype = complex)
    # For all vectors of length ``m.nrows`` with bit 0
    # of ``v`` equal to 1 put ``index =  v * A``, and  let ``q``
    # be the sum of the matrix ``R`` derivend from ``v`` and
    # ``Q`` as in function ``eval_A_Q()``.
    # Then add ``phases[q mod 4]`` to entry ``a[index]``
    ncols, data = m.ncols, m.copy().check().raw_data[:m.nrows]
    for v in range(1, 1 << len(data), 2):
        index, q = eval_A_Q(ncols, data, v)
        a[index] += phases[q]
    # Change shape of array ``a`` as required
    rows, cols = m.shape
    return a.reshape((1 << rows, 1 << cols))



#####################################################################
# Check echelon form  matrix of type QStateMatrix
#####################################################################


def check_echelon(m, reduced = False):
    """Check if a matrix ``m`` is in echelon form

    ``m`` must be a matrix of type ``QStateMatrix``.
    If ``reduced`` is True them the matrix is checked to be
    in reduced echelon form.

    Raise ``ValueError`` if matrix ``m`` is not in correct
    echelon form.
    """
    err = "Bad vector in row % of part A of QSMatrix"
    data = m.copy().check().raw_data[:m.nrows]
    mask = (1 << m.ncols) - 1
    index = m.ncols + 1
    for i, v in enumerate(data[1:]):
        v &= mask
        new_index = v.bit_length()
        if new_index >= index:
            print("\nError: " + err % i)
            print(m)
            err1 = "QSMatrix is not in echelon form"
            raise ValueError(err1)
        index = new_index
        if reduced:
            v = reduce(__or__, data[:i], 0) 
            if (v & (1 << (index - 1))):
                print("\nError: " + err % i)
                print(m)
                err1 = "QSMatrix is not in reduced echelon form"
                raise ValueError(err1)
                
            
        
        


#####################################################################
# Check that states are approximately equal
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
    c1 = m1.copy().complex()
    c2 = m2.copy().complex()
    eq = np.amax(abs(c1 - c2), initial = 0.0) < EPS
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
    
    ``c1`` is the expected and ``c2`` is the obtained complex matrix.
    ``qstate`` is an optiona state representing the matrix.
    """
    try:
        diff = np.amax(abs(c1 - c2), initial = 0.0)
    except:
        err = "\nCannot compare arrays of shapes %s and %s\n"
        print(err % (c1.shape, c2.shape))
        raise
    if diff > EPS:
        err = text if text else "Error in comparing complex matrices"
        print("\n" + err + "\n")
        print("qstate  = ", qstate)
        if max(sum(c1.shape), sum(c2.shape)) <= 64:
            print("expected\n", c1)
            print("obtained\n\n", c2)
            print("expected - obtained\n\n", c1 - c2)
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
    c1 = m.copy().complex()
    err = "Complex matrix from unreduced state"
    compare_complex(c0, c1, err, m)


#####################################################################
# Create radnom indices
#####################################################################

def random_slice(length):
    """Return data for indexing a random slice of a vector

    Let ``v`` be a one-dimensional numpy array of a given ``length``.
    
    We return a random triple ``(start, stop, step)`` for indexing a 
    the array ``v`` in the form ``v[start:stop:step]``.
    """
    start = randint(-length, length-1)
    stop = randint(-length, length-1)
    down = stop % length < start % length
    step = max(1, randint(0, abs(start - stop)//2))
    step *= (-1)**down 
    return start, stop, step    

#####################################################################
# Display and check some matrices
#####################################################################




qs_matrix_data = [
     (0, 0, (1,0), []),
     (0, 0, (-3,6), 0),
     (3, 0, (-3,5), 5),
     (0, 6, (9,2), 25),
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
    yield qs_unit_matrix(4) 
    for i in range(20):
        yield qs_rand_matrix(2, 0, 3)
    for i in range(20):
        yield qs_rand_matrix(1, 0, 2)
    for rows in range(6):
        for cols in range(5):
            for nr in [1,2] + list(range(rows+cols, rows+cols+3)):
                m = qs_rand_matrix(rows, cols, nr)
                m.mul_scalar(randint(-8, 8), randint(0,7))  
                yield m                
                


@pytest.mark.qstate
def test_qs_matrix(verbose = 0):
    """Basic test for class ``QStateMatrix`` 
    
    It tests the ``__str__ `` method, the conversion to a
    complex matrix and the echelonization and the reduction of 
    an instance of class ``QStateMatrix``.
    """
    FORMAT_REDUCED = False
    display_len =  len(qs_matrix_data)
    for ntest, m in enumerate(create_display_testvectors()):
        if verbose or ntest < display_len:
            print("TEST %s" % (ntest+1))
            print(m.copy())
        check_complex(m)
        m2 = m.copy().echelon()
        check_echelon(m2)
        check_eq_qstate(m,m2, "Echelonized")
        if verbose:
            print("Echelonized")
            print(m2)
        m2 = m.copy().reduce()
        check_echelon(m2, reduced = True)
        check_eq_qstate(m,m2, "Reduced")
        if verbose: 
            print("Reduced")
            print(m2)
        assert m2._check_reduced(forced=True)
    FORMAT_REDUCED = True

 

@pytest.mark.qstate
def test_getitem(verbose = 0):
    FORMAT_REDUCED = False
    for ntest, m in enumerate(create_display_testvectors()):
        if verbose:
            print("TEST %s" % (ntest+1))
            print(m)
        c = m.complex()
        f0, t0, d0 = random_slice(1 << m.shape[0])
        f1, t1, d1 = random_slice(1 << m.shape[1])
        c1_ref = c[f0:t0:d0, f1:t1:d1] 
        c1 =  m[f0:t0:d0, f1:t1:d1] 
        s = "index m[%d:%d:%d, %d:%d:%d]" % (f0, t0, d0, f1, t1, d1)
        if verbose:
            print(s)
            print(c1)
        compare_complex(c1_ref, c1, s, m) 
        
        i0 = randint(-1 <<  m.shape[0], (1 <<  m.shape[0]) - 1 )     
        i1 = randint(-1 <<  m.shape[1], (1 <<  m.shape[1]) - 1 ) 
        c1_ref = c[f0:t0:d0, i1] 
        c1 =  m[f0:t0:d0, i1] 
        s = "index m[%d:%d:%d, %s]" % (f0, t0, d0, i1)        
        if verbose:
            print(s)
            print(c1)
        compare_complex(c1_ref, c1, s, m)  

        c1_ref = c[i0, f1:t1:d1] 
        c1 =  m[i0, f1:t1:d1] 
        s = "index m[%d, %d:%d:%d]" % (i0, f1, t1, d1)
        if verbose:
            print(s)
            print(c1)
        compare_complex(c1_ref, c1, s, m)  

        c1_ref = c[i0, i1] 
        c1 =  m[i0, i1] 
        s = "index m[%d, %d]" % (i0, i1)
        if verbose:
            print(s)
            print(c1)
        compare_complex(np.array(c1_ref), np.array(c1), s, m)  
    FORMAT_REDUCED = True

 
#####################################################################
# Test the generation of a monomial matrix
#####################################################################


def eval_nonomial(data):
    """Auxiliary function for function ``slow_complex``.
    
    Yet to be documented
    
    TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    This function sets map(0b000..01) to data[0] ^ data[1]
    The C function might set  map(0b000..01) to data[1]
    
    Find a reasonable definition of this function!!!!!!!!!!!!!!!
    """
    nqb = len(data) - 1
    res = np.zeros( (1 << nqb, 2), dtype = np.int32)
    data = data[:]
    mask = (2 << nqb) - 1
    data[0] &= mask
    for i in range(1, nqb + 1):
        data[i] &= mask
        mask += mask + 1
        
    mask =  (1 << nqb) - 1       
    for j in range(1 << nqb):
        v = j + j + 1    
        a = 0
        for i, d in enumerate(data): 
            a ^= d & -((v >> i) & 1)
        q = (a >> nqb) & v        
        for sh in [32, 16, 8, 4, 2, 1]:
            q ^= q >> sh 
        res[j,0] = a & mask
        res[j,1] = 1 - 2*(q & 1)
    return res


def a_binary(data, length):
    s = ""
    for n in data:
        n &= ((1 << length) - 1)
        b = format(n, "b")
        s +=  "0" * (length - len(b)) + b + "\n"
    return s


      


qs_monomial_matrix_data = [
     [1],
     [0b10, 0b00],
     [0b11, 0b11],
    ]
 
def rand_monomial_data(nqb):
    k = 1 << (2 * nqb)
    return [randint(0, k-1) for i in range(nqb+1)]
 
def create_monomial_testvectors():
    """yield instances of class ``QStateMatrix`` for tests """
    for data in qs_monomial_matrix_data:
        yield data
    for nqb in range(8):
        for k in range(3):
            yield(rand_monomial_data(nqb))
            
 
@pytest.mark.qstate
def test_monomial(verbose = 0):
    FORMAT_REDUCED = False
    for ntest, data in enumerate(create_monomial_testvectors()):
        nqb = len(data) - 1
        m = qs_column_monomial_matrix(data)
        if verbose:
            print("TEST %s, nqb = %d" % (ntest+1, nqb))
            s_data =  a_binary(data, 2*nqb)
            print("input:\n" + s_data)
            print(m)

        c = m[:,:].T
        s = np.sum(abs(c), axis = 1) - 1
        error = np.amax(abs(s)) > EPS
        if (verbose or error):
            if error:
                err = "Error in checking monomial matrix"
                print("\n" + err + "\n")
            s_data =  a_binary(data, 2*nqb)
            print("data  = \n" +  s_data)
            print("obtained:", m)
            if nqb <= 3:
                print("obtained matrix\n\n", c)
            print("Shape of obtained matrix:", m.shape, c.shape)
            if error:
                print("Max. absolute error:", np.amax(abs(s)))
                raise ValueError(err)
        res = eval_nonomial(data)
        for j, w in  enumerate(res):
            k = w[0]
            value = w[1]
            assert c[j,k] == value, (nqb, j, k, c[j,k], value, res)
        # Test qs_rown_monomial_matrix(data)
        m_r = qs_row_monomial_matrix(data)
        #print("m_r =", m_r, "m_c=", m ," m_T =", m.T)      
        assert m_r == m.T, (str(m_r), str(m.T))
    FORMAT_REDUCED = True


#####################################################################
# Test some error cases
#####################################################################
  
@pytest.mark.qstate
def test_qs_errors(verbose = 0):
    FORMAT_REDUCED = False
    with pytest.raises(ValueError):
        QStateMatrix(10000, 10000) 
    with pytest.raises(ValueError):
        QStateMatrix(1, 1, [0, 0, 0xf]) 
    m1, m2 = qs_unit_matrix(1),  qs_unit_matrix(2)
    with pytest.raises(ValueError):
        m1 @  m2
    with pytest.raises(ValueError):
        m1 *  m2
    m =  QStateMatrix(3, 3, [0, 0, 0xf], 2) 
    with pytest.raises(ValueError):
       m.power(-1) 
    FORMAT_REDUCED = True



