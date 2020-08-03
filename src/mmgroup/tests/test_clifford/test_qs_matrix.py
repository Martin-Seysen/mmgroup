from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix, rand_qs_matrix



#####################################################################
# display some matrices
#####################################################################

def eval_A_Q(ncols, data, v):
    v &= (1 << len(data)) - 1
    a = q = diag = 0
    for i, d in enumerate(data):
        d &=  -((v >> i) & 1)
        a ^= d
        d = (d >> ncols) & v
        diag += (d >> i) & 1
        q ^= d &  ((1 << i) - 1)
    for sh in [32, 16, 8, 4, 2, 1]:
        q ^= q >> sh 
    a &= (1 << ncols) - 1
    q = ((q << 1) + diag) & 3
    return a, q   




def slow_complex(m):
    e, phi = m.factor 
    f = 2.0**(0.5*(e - phi)) * (1+1j)**phi
    phases = [f * 1j**i for i in range(4)]
    a = np.zeros(1 << m.ncols, dtype = np.complex)
    ncols, data = m.ncols, m.data
    for v in range(1, 1 << len(data), 2):
        index, q = eval_A_Q(ncols, data, v)
        a[index] += phases[q]
    return a.reshape((1 << m.rows, 1 << m.cols))


#####################################################################
# Check that states are approcimately equal
#####################################################################

EPS = 1.0e-8

def check_eq_qstate(m1, m2, text = None):
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
         


def check_complex(m):
    c1 = m.complex_unreduced()
    c2 = slow_complex(m)
    eq = np.amax(abs(c1 -c2)) < EPS
    if not eq:
        err = "Error in computing complex matrix from QStateMatrix m"
        print("\n" + err + "\n")
        print("m = ", m)
        if max(m.shape) < 6:
            print("C result = \n", c1)
            print("Python result = \n", c2)
        raise ValueError(err)


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
    ]


def create_display_testvectors():
    """yield rows, cols, factor, data"""
    for rows, cols, factor, data in qs_matrix_data:
        m = QStateMatrix(rows, cols, data)
        m.mul_scalar(*factor) 
        yield m
    yield QStateMatrix(4, 4,  1) 
    for i in range(20):
        yield rand_qs_matrix(1, 0, 2)
    for rows in range(6):
        for cols in range(5):
            for nr in range(rows+cols, rows+cols+3):
                m = rand_qs_matrix(rows, cols, nr)
                m.mul_scalar(randint(-8, 8), randint(0,7))  
                yield m                
                


@pytest.mark.qstate
def test_display_qs_matrix(verbose = 0):
    for ntest, m in enumerate(create_display_testvectors()):
        if verbose or ntest < len(qs_matrix_data):
            print("TEST %s" % (ntest+1))
            print(m)
        check_complex(m)
        m2 = QStateMatrix(m)
        if verbose: print("Echelonize....")
        m2.echelon()
        if verbose: print(m2)
        check_eq_qstate(m,m2)
        if verbose: print("Reducing....")
        m2.reduce()
        if verbose: print(m2)
        check_eq_qstate(m,m2)








        