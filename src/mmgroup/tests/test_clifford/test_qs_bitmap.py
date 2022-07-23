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
        a = m.complex().real.ravel()
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



          
    
