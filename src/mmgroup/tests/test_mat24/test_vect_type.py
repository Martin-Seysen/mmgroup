from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from mmgroup import mat24, MAT24_ORDER
from mmgroup.mat24 import Mat24Sub_count_vtypes


@pytest.mark.mat24
def test_mat24_vect_type(verbose = 0):
    print("Testing function mat24_vect_type()")
    d, dc = Mat24Sub_count_vtypes(verbose)
    assert len(d) == len(dc) == 49
    for n in d.values():
        assert MAT24_ORDER % n == 0 
    #print("passed")

