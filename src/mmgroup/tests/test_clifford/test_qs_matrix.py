from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix, rand_qs_matrix



#####################################################################
# display some matrices
#####################################################################

EPS = 1.0e-8

def check_eq_qstate(a, b):
    eq = np.amax(abs(a.complex_unreduced()-b.complex_unreduced())) < EPS
    if not eq:
        err = "Error in comparing instances of QStateMatrix"
        print("\n" + err + "\n")
        print("m1 = ", a)
        print("m2 = ", b)
        if max(a.shape + b.shape) < 4:
            print("complex m1 = \n", a.complex_unreduced())
            print("complex m2 = \n", b.complex_unreduced())
        raise ValueError(err)
         

#####################################################################
# display some matrices
#####################################################################


def create_display_testvectors():
    """yield rows, cols, factor, data"""
    qs_matrix_data = [
     (0,0, (1,0), []),
     (0,0, (-3,6), 0),
     (3,0, (-3,5), 5),
     (0,6, (9,2), 25),
     (3,1, (1,0), []),
     (2, 2, (0,0), [0b110_10_01, 0b101_01_11, 0b011_01_00]),
    ]
    for rows, cols, factor, data in qs_matrix_data:
        m = QStateMatrix(rows, cols, data)
        m.mul_scalar(*factor) 
        yield m
    for i in range(20):
        yield rand_qs_matrix(1, 0, 2)
    for rows in range(6):
        for cols in range(5):
            for nr in range(rows+cols, rows+cols+3):
                m = rand_qs_matrix(rows, cols, nr)
                m.mul_scalar(randint(-8, 8), randint(0,7))  
                yield m                
                


@pytest.mark.qstate
def test_display_qs_matrix():
    print( QStateMatrix(4, 4,  1) )
    for ntest, m in enumerate(create_display_testvectors()):
        print("TEST %s" % (ntest+1))
        print(m)
        #m1 = QStateMatrix(m)
        #m1.mul_scalar(8 + ((factor[0] + factor[1]) & 1))
        m2 = QStateMatrix(m)
        print("Echelonize....")
        m2.echelon()
        print(m2)
        check_eq_qstate(m,m2)
        print("Reducing....")
        m2.reduce()
        print(m2)
        check_eq_qstate(m,m2)








        