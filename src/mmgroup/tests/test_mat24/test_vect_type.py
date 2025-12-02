from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import time
import pytest

from mmgroup import mat24, MAT24_ORDER
from mmgroup.mat24 import Mat24Sub_count_vtypes


@pytest.mark.mat24
def test_mat24_vect_type(verbose = 0):
    print("Testing function mat24_vect_type()")
    t_start = time.process_time()
    d = Mat24Sub_count_vtypes(verbose)
    t = time.process_time() - t_start
    assert len(d) == 49
    assert sum(d.values()) == 1 << 24
    assert max(d.keys()) < 160
    for n in d.values():
        assert MAT24_ORDER % n == 0 
    #print("Time = %.1f ms" % (1000*t)) 
    #print("passed")

