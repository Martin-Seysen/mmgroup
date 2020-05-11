
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import collections
import re
import warnings
from random import randint
import pytest


from mmgroup import mat24 as m24
from mmgroup.tests.spaces.sparse_mm_space import SparseMmSpace





p = 241
rep_mm = SparseMmSpace(p)
grp = rep_mm.group


def vector_data(v):
    data = v.as_tuples()
    assert len(data) == 1
    return data[0]


def as_suboctad(v1, d):
    c = m24.ploop_cap(v1, d)
    return m24.cocode_to_suboctad(c, d)




def one_test_op_xy(v, eps, e, f, verbose = 0):
    pi_atom = grp.atom('d', eps)
    x_atom = grp.atom('x', e)**(-1)
    y_atom = grp.atom('y', f)**(-1)
    w = v * pi_atom * x_atom  * y_atom 
    if verbose:
        print("%s * %s * %s * %s = %s" % (v, pi_atom, x_atom, y_atom, w))
        
    v_value, v_tag, v_d, v_j = vector_data(v)
    w_value, w_tag, w_d, w_j = vector_data(w)
    neg = v_value != w_value
    if neg:
        assert (v_value + w_value) % v.space.p == 0
    err = sign_err = False
    if v_tag == 'X':
        err |= w_tag != 'X'
        err |= w_j != v_j 
        err |= w_d != v_d ^ ((f >> 1) & 0x7ff)
        # Next we check the sign
        d = (v_d << 1) & 0xfff
        c =  m24.vect_to_cocode(1 << v_j)
        sign =  m24.gcode_weight(e ^ f) 
        sign += m24.gcode_weight(f) 
        sign += m24.gcode_weight(d) * (eps + 1)
        sign += m24.gcode_weight(d ^ e ^ f)
        sign += m24.scalar_prod(e, c)
        cc = c if eps & 1 else 0 
        cc ^= eps ^ m24.ploop_theta(f) ^ m24.ploop_cap(e, f) 
        sign += m24.scalar_prod(d, cc)
        sign += f >> 12
    elif v_tag in 'YZ':
        err |= w_tag not in 'YZ'
        err |= w_j != v_j 
        tau = v_tag == 'Y'
        sigma = (tau + eps) & 1
        sigma_w = w_tag == 'Y'
        err |= sigma != sigma_w
        sigma_1 = 0 if sigma else 0x7ff
        err |= w_d != v_d ^ ((e >> 1) & 0x7ff) ^ ((f >> 1) & sigma_1)
        # Next we check the sign
        d = (v_d << 1) & 0xfff
        c =  m24.vect_to_cocode(1 << v_j)
        sign = m24.ploop_cocycle(f, e) * (sigma + 1)
        sign += m24.gcode_weight(f) * (sigma)
        sign += m24.gcode_weight(d ^ e)  
        sign += m24.gcode_weight(d ^ e ^ f)  
        sign += m24.scalar_prod(f, c)
        cc = eps ^ m24.ploop_theta(e) 
        cc ^=  m24.ploop_theta(f) * (sigma ^ 1) 
        sign += m24.scalar_prod(d, cc)
        sign += (e >> 12) + (f >> 12) * (sigma + 1)
        sign += (e & sigma)
    elif v_tag == 'T':
        d = m24.octad_to_gcode(v_d)
        err |= w_tag != 'T'
        err |= w_d != v_d 
        err |= w_j != v_j  ^ as_suboctad(f, d)
        sign = m24.gcode_weight(d ^ e) + m24.gcode_weight(e)
        sign +=  m24.scalar_prod(d, eps)
        sign += m24.suboctad_scalar_prod(as_suboctad(e ^ f, d), v_j)
        sign += m24.suboctad_weight(v_j) * eps
    elif v_tag in 'BC':
        err |= w_tag not in 'BC'
        m, n = v_tag == 'C', w_tag == 'C'
        v_i, v_j = max(v_d, v_j), min(v_d, v_j)
        w_i, w_j = max(w_d, w_j), min(w_d, w_j)
        err |= v_i != w_i or v_j != w_j
        c = m24.vect_to_cocode((1 << v_i) ^ (1 << v_j))
        err |= m ^ n ^ m24.scalar_prod(f, c)
        sign = m * eps + m24.scalar_prod(e ^ f, c)
    elif v_tag == 'A':
        err |= w_tag != 'A'
        v_i, v_j = max(v_d, v_j), min(v_d, v_j)
        w_i, w_j = max(w_d, w_j), min(w_d, w_j)
        err |= v_i != w_i or v_j != w_j
        c = m24.vect_to_cocode((1 << v_i) ^ (1 << v_j))
        sign = m24.scalar_prod(f, c)
    else:
        raise ValueError("Bad tag " + v_tag)
    sign_err = (sign & 1) != (neg & 1)
    err |= sign_err
    if verbose or err:
        err_str = "Sign error" if sign_err else "Error"
        if err: print(err_str, end = ": ")
        print("%s * d_%x * x_%x * y_%x = %s" % (v, eps, e, f, w))
        if err:
            raise ValueError("Testing operation xy failed")
            pass






def rand_v(tag):
    return rep_mm.unit(tag, "r")
   

def op_xy_testdata():
    data = [
       [ ("Y", 3, 6),  0, 0, 0x1171 ],
       [ ("Y", 3, 6),  12, 0, 0 ],
       [ ("Y", 3, 6),  12, 1111, 0 ],
       [ ("Y", 3, 6),  12, 0, 1111],
       [ ("Z", 3, 6),  0, 0, 0x1171 ],
       [ ("Z", 3, 6),  12, 0, 0 ],
       [ ("Z", 3, 6),  12, 1111, 0 ],
       [ ("Z", 3, 6),  12, 0, 1111],
       [ ("X", 3, 6),  0, 0, 0x1171 ],
       [ ("X", 3, 6),  12, 0, 0 ],
       [ ("X", 3, 6),  12, 1111, 0 ],
       [ ("X", 3, 6),  12, 0, 1111],
    ]
    data = []
    for v, f, e, eps in data:
        if isinstance(v, str): 
            v1 = rand_v(v)
        else :
            v1=  rep_mm.unit(*v)
        yield v1,  f, e, eps
    for v in "ZY":
        for i in range(500):
            v1 =  rand_v(v)
            eps = randint(0, 0xfff) 
            e = randint(0, 0x1fff)
            yield v1,  eps, e, 0 
    v_tags = "ABCTYZX"
    #v_tags = "ZYX"
    for v in v_tags:
        for i in range(1000):
            v1 =  rand_v(v)
            eps = randint(0, 0xfff)
            e = randint(0, 0x1fff)
            f = randint(0, 0x1fff)
            yield v1,  eps, e, f
             

@pytest.mark.space
def test_op_xy(verbose = 0):
    print("Testing group operations x, y ...")
    for v,  f, e, eps in op_xy_testdata():
        one_test_op_xy(v, f, e, eps, verbose = verbose)
    print("passed")






