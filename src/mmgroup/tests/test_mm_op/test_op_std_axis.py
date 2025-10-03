from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from random import randint
from importlib import import_module
import time
from functools import partial

from mmgroup import MM0, MMV, AutPL, XLeech2, MMVector, SubOctad
from mmgroup.mm_space import characteristics 
from mmgroup.mm_op import mm_op_mul_std_axis
from mmgroup.mm_op import mm_sub_table_octad_to_std_ax_op

import pytest




#######################################################################
# Test function mm_group_prepare_op_ABC
#######################################################################

STD_XL = XLeech2(0, [2,3])

def mul_std_axis(v):
    v1 = v.copy()
    mm_op_mul_std_axis(v1.p, v1.data)
    v1.check()
    return v1 



def eigenvector_from_xleech2(x0, V, verbose = 0):
    assert isinstance(x0, XLeech2)
    t0, t0sub = x0.subtype
    x1 = x0 * STD_XL
    t1, t1sub = x1.subtype
    assert t0 == 2
    v0 = V(x0)
    if (t0sub == 2) and verbose:
        o_, sub_ = x0.as_suboctad()
        octad_code = mm_sub_table_octad_to_std_ax_op(o_)
        print('cc', octad_code)
        print('oo', hex(x0.ord), o_, hex(sub_), (t0, t0sub, t1, t1sub))
    if t1 == 2:
        v1 = V(x1)
        yield v0 + v1, 0
        yield v0 - v1, 8
    elif t1 == 3:
        yield v0, 1
    elif t1 == 4:
        yield v0, 0
    else:
        assert t1 == 0
        yield V("A_2_2+A_3_3-A_2_3+2*B_2_3"), 0
        yield V("A_2_2+A_3_3-A_2_3-2*B_2_3"), 32


def x_leech2_cases(ntests = 5):
   cases = [ 
       (0, [2,3]),  (0, [2,1]),  (0, [0,1]), 
       ('T',561,0x16), ('T',561,0x3b), ('T',735,0x16),
       ('X',0x534,2), ('X',0x534,3), ('X',0x534,17),
       ('X',0x2f5,3),
   ]
   for data in cases:
       yield  XLeech2(*data)
   for i in range(ntests):
       yield XLeech2('r', 2)


def eigenvector_testcases(V, ntests = 5, verbose = 0):
    for x0 in x_leech2_cases(ntests):
        yield from eigenvector_from_xleech2(x0, V, verbose)



def linear_testcases(V, ntests=5, debug_cases = False):
    S_DATA = [
        ("A_18_2", "A_18_2"),
        ("A_17_13+A_18_2+A_19_0", "A_11_3+A_15_2+A_18_2+A_21_0+A_21_5"),
        ("X_18_2", "X_18_2"),
    ]
    if debug_cases:
       for s1, s2 in S_DATA:
           yield V(s1),  V(s2)
    for i in range(3):
        yield V("R"), V("R")



@pytest.mark.mm_op
def test_eigenvectors_linear(verbose = 1):
    for p in  reversed(characteristics()):
        V = MMV(p)
        for i, (v1, v2) in enumerate(linear_testcases(V, ntests=1000)):
            v3 = v1 + v2
            w1 = mul_std_axis(v1)
            w2 = mul_std_axis(v2)
            w3 = mul_std_axis(v3)
            ok = w3 == w1 + w2
            if not ok:
                print("\nv1=", v1)
                print("v2=", v2)
                print("v3=", v3)
                print("w1=", w1)
                print("w2=", w2)
                print("w3=", w3)
                print("w1+w2=", w1+w2)
                print("err=", w3 - w1 - w2)
                #for i in range(-4,10):
                #    print(i,  mul_std_axis(v1) * i, mul_std_axis(v1*i))
                ERR = "Linearity error in Griess algebra computation mod %d" 
                raise ValueError(ERR % p)
           




@pytest.mark.mm_op
def test_eigenvectors(verbose = 0):
    print("Testing function mm_op_mul_std_axis")
    n_all = 0
    for p in  reversed(characteristics()):
        n = 0
        V = MMV(p)
        v_sum, img_sum = V(),  V()
        for vector, value in eigenvector_testcases(V, verbose=verbose):
            vector *= MM0('r', 'B')
            vector *= randint(1, p-1)
            n += 1
            if verbose:
                print("p = %d, testcase %d" % (p, n))
            v1 = mul_std_axis(vector)
            ok = v1 == value * vector
            if verbose or not ok:
                print("v= ", vector)
                print("v1=", v1)
                print("eigenvalue =", value)
                if not ok:
                    ERR = "v1 should be equal to %d * v0 (mod %d)" 
                    raise ValueError(ERR % (value, p)) 
            v_sum +=  vector          
            img_sum +=  v1
        n_all += n 
        img_computed = v_sum.copy()
        mm_op_mul_std_axis(p, img_computed.data)
        img_computed.check()
        ok = img_sum == img_computed
        if not ok:
            diff = img_sum - img_computed
            print(diff)
            ERR = "Error in general Griess algebra computation mod %d" 
            raise ValueError(ERR % p)
    print(n_all, "test passed")


