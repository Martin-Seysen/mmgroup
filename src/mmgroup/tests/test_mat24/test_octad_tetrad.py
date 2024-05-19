from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint, sample

from mmgroup import mat24

from mmgroup.bitfunctions import bit_mat_mul, v2, bitweight




#####################################################################
# Test function mat24.intersect_octad_tetrad()
#####################################################################



def mul_mat24(v, g):
    pi = mat24.m24num_to_perm(g % mat24.MAT24_ORDER)
    return mat24.op_vect_perm(v, pi)



def intersect_octad_tetrad_testcases(ntests = 10):
    # good cases
    for i in range(ntests):
        for size in range(1, 8):
            data = sample(range(8), size)
            o = 0xff0
            v2 = sum((1 << x for x in data))
            g = randint(0, 1 << 40)
            o, v2 = mul_mat24(o, g), mul_mat24(v2, g)
            yield o, v2, True
    BAD_SAMPLES =  [
        ([], []),
        (range(8), [8, 9, 10, 12]), 
    ]   
    for o, v2 in BAD_SAMPLES:
        o = sum((1 << x for x in o))
        v2 = sum((1 << x for x in v2))
        g = randint(0, 1 << 40)
        o, v2 = mul_mat24(o, g), mul_mat24(v2, g)
        yield o, v2, False





@pytest.mark.mat24
def test_intersect_octad_tetrad():
    FAIL = "Function mat24.intersect_octad_tetrad should fail but didn't"
    print("Testing function mat24.intersect_octad_tetrad_testcases")
    for o, v2, good in intersect_octad_tetrad_testcases():
        #print("octad, vector", hex(o), hex(v2), hex(o & v2))
        if good:
            o1 = mat24.intersect_octad_tetrad(o, v2)
            #print(hex(o1))
            assert mat24.bw24(o1) == 8
            assert mat24.bw24(o1 & o) == 4
            assert o1 & v2 == v2
        else:
            failed = False
            try:
               o1 = mat24.intersect_octad_tetrad(o, v2)
            except:
               failed = True
            assert failed