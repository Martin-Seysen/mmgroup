from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint #, shuffle, sample

from mmgroup.structures.qs_matrix import QStateMatrix




#####################################################################
# display some matrices
#####################################################################


def create_display_testvectors():
    """yield rows, cols, factor, data"""
    yield  (0,0, (1,), [])
    yield  (0,0, (-3,6), 0)
    yield  (3,0, (-3,5), 5)
    yield  (0,6, (9,2), 25)
    yield  (3,1, (1,), [])
    yield  (2, 2, (0,), [0b110_10_01, 0b101_01_11, 0b011_01_00])




@pytest.mark.qstate
def test_display_qs_matrix():
    for rows, cols, factor, data in create_display_testvectors():
        m = QStateMatrix(rows, cols, data)
        m.mul_scalar(*factor)
        print(str(m))








        