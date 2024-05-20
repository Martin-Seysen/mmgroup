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
            if randint(0, 1):
                o ^= 0xffffff
            yield o, v2, True
    BAD_SAMPLES =  [
        ([], [], None),
        (range(7), [1,2,3], None),
        (range(8), [8, 9, 10, 12], False), 
    ]   
    for o, v2, status in BAD_SAMPLES:
        o = sum((1 << x for x in o))
        v2 = sum((1 << x for x in v2))
        g = randint(0, 1 << 40)
        o, v2 = mul_mat24(o, g), mul_mat24(v2, g)
        yield o, v2, status

@pytest.mark.mat24
def test_intersect_octad_tetrad(verbose = 0):
    FAIL = "Function mat24.intersect_octad_tetrad should fail but didn't"
    print("Testing function mat24.intersect_octad_tetrad_testcases")
    for o, v2, good in intersect_octad_tetrad_testcases():
        # print("octad, vector", hex(o), hex(v2), hex(o & v2))
        if good:
            o1 = mat24.intersect_octad_tetrad(o, v2)
            ok = mat24.bw24(o1) == 8
            ok |= mat24.bw24(o1 & o) == 4
            ok |= o1 & v2 == v2
            if verbose or not ok:
                print("Testing o = %s, v2 = %s" % (hex(o), hex(v2)))
                print("Obtained: o1 =", hex(o1))
                print("o1 & v2 =", hex(o1 & v2))
                for var in ['o', 'v2', 'o1', 'o1 & o', '(o1 & v2) ^ v2']:
                    w = eval('mat24.bw24(%s)' % var)
                    print(var, "has weight", w)
                if not ok:
                    raise ValueError("Test failed")
        elif good is not None:
            o1 = mat24.intersect_octad_tetrad(o, v2)
            assert o1 == 0, [hex(x) for x in (o, v2, o1)]
        else:
            failed = False
            try:
                o1 = mat24.intersect_octad_tetrad(o, v2)
                # print([hex(x) for x in (o, v2, o1)])
            except:
               failed = True
            assert failed, [hex(x) for x in (o, v2, o1)] 
