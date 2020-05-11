"""Test auxiliary C functions in generated file  mat24_functions.c

These functions are wrapped in class mat24fast.Mat24 and used for the
representation of the monster group.
"""
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys 
from random import shuffle, randint
import pytest

from mmgroup import mat24
Mat24 = mat24

def op_all_autpl_testdata():
    def to_autpl(c, perm_num):
        p = Mat24.m24num_to_perm(perm_num)        
        m = Mat24.perm_to_autpl(c, p) 
        assert (m[11] >> 12) & 1 == (c >> 11) & 1
        return m
    testdata = [ 
        (1, 0, [1 << i for i in range(11)]),
        (3, 17, [1 << i for i in range(11)]),
    ]
    for d, n,  vlist in  testdata:
        yield to_autpl(d, n),  vlist
    for i in range(500):
        m = to_autpl(randint(0, Mat24.MAT24_ORDER - 1), randint(0,0xfff))
        yield (m ,  [randint(0, 0x7ff) for i in range(20)])





@pytest.mark.mat24
def test_op_all_autpl():
    print("Testing function op_all_autpl()...", end="")
    for m, vlist in op_all_autpl_testdata():
        v_all = Mat24.op_all_autpl(m)
        odd = (m[11] & 0x1000) >> 12
        for v in vlist:
            w0 = v_all[v]
            signs = [(w0 >> i) & 1 for i in (14,13,12)]
            w = (w0 & 0xfff)
            im_v = Mat24.op_ploop_autpl(v, m) 
            assert (im_v & 0xfff) ==  w, (v, im_v, hex(w0))
            assert im_v >> 12 == signs[1], (v, im_v, signs, hex(w0))
            assert ((im_v >> 12) ^ (im_v >> 11)) & 1 == signs[0], (
                  v, im_v, signs, hex(w0))
            pwr = signs[1] ^ signs[2]
            weight = Mat24.gcode_weight(im_v)  & 1
            assert pwr == weight & odd, (v, im_v, signs, weight, hex(w0))
    print(" passed")
  

@pytest.mark.mat24
def test_autpl_cocode():
    print("Testing function op_all_cocode()...", end="")
    for i in range(100):
        c = randint(0,0xfff)
        m = Mat24.perm_to_autpl(c, list(range(24))) 
        v_all = Mat24.op_all_autpl(m)
        v_c = Mat24.op_all_cocode(c)
        v_all = [(x >> 12) & 7 for x in  v_all]
        v_c =  [x & 7 for x in  v_c]
        if v_all != v_c:
            for j in range(len(v_all)):
                if v_all[j] != v_c[j]:
                    print("cocode = %s, v_all[%s] = %s, v_c[%s] = %s" % 
                        tuple(map(hex, (c, j, v_all[j], j, v_c[j]))))
                    print(v_all[:32])
                    print(v_c[:32])
                    err = "Testing of function op_all_cocode() failed"
                    raise ValueError(err) 
                    break

