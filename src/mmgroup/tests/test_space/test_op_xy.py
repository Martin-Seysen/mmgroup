
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
from mmgroup.tests.spaces.sparse_mm_space import SparseMmV
from mmgroup.tests.groups.mgroup_n import MGroupNWord







p = 241
rep_mm = SparseMmV(p)
grp = MGroupNWord


def vector_data(v):
    data = v.as_tuples()
    assert len(data) == 1
    return data[0]


def as_suboctad(v1, d):
    c = m24.ploop_cap(v1, d)
    return m24.cocode_to_suboctad(c, d) & 0x3f


1234567890123456789012345678901234567890123456789012345678901234567890
def op_xy(v, eps, e, f):
    """Multiply unit vector v with group element

    This function multplies a (multiple of a) unit vector v
    with the group element

        g  =  d_<eps> * (x_<e>)**(-1)  * (y_<f>)**(-1) .
   
    It returns the vector w = v * g. This function uses the same 
    formula for calculating  w = v * g,  which is used in the 
    implementation of the monster group for computing
    v = w * g ** (-1).

    Input vector v  must be given as a tuple as described in class
    mmgroup.structures.abstract_mm_rep_space.AbstractMmRepSpace.
    Output vector is returned as a tuple of the same shape.
    """
    if len(v) == 3:
        v = (1,) + v
    v_value, v_tag, v_d, v_j = v
    parity = (eps >> 11) & 1  # parity of eps
    if v_tag == 'X':
        w_d = v_d ^ (f & 0x7ff)
        sign = v_d >> 12
        d = v_d & 0x7ff
        c =  m24.vect_to_cocode(1 << v_j)
        sign +=  m24.gcode_weight(e ^ f) 
        sign += m24.gcode_weight(f) 
        sign += m24.gcode_weight(d) * (parity + 1)
        sign += m24.gcode_weight(d ^ e ^ f)
        sign += m24.scalar_prod(e, c)
        cc = c if parity else 0 
        cc ^= eps ^ m24.ploop_theta(f) ^ m24.ploop_cap(e, f) 
        sign += m24.scalar_prod(d, cc)
        sign += f >> 12
        return (-1)**sign * v_value, 'X', w_d, v_j
    elif v_tag in 'YZ':
        tau = v_tag == 'Y'
        sigma = (tau + parity) & 1
        w_tag = "ZY"[sigma]
        sign = (v_d >> 12) + ((v_d >> 11) & tau)
        w_d = (v_d ^ e ^ (f & ~-sigma)) & 0x7ff
        #print("YZ, w_d", hex(w_d), hex(v_d ^ (e & 0x7ff) ^ (f & sigma_1)), err)
        # Next we check the sign
        d = v_d & 0x7ff
        c =  m24.vect_to_cocode(1 << v_j)
        sign += m24.ploop_cocycle(f, e) * (sigma + 1)
        sign += m24.gcode_weight(f) * (sigma)
        sign += m24.gcode_weight(d ^ e)  
        sign += m24.gcode_weight(d ^ e ^ f)  
        sign += m24.scalar_prod(f, c)
        cc = eps ^ m24.ploop_theta(e) 
        cc ^=  m24.ploop_theta(f) * ((sigma ^ 1) & 1)
        sign += m24.scalar_prod(d, cc)
        sign += (e >> 12) + (f >> 12) * (sigma + 1)
        sign += ((e >> 11) & sigma)
        return (-1)**sign * v_value, w_tag, w_d, v_j
    elif v_tag == 'T':
        d = m24.octad_to_gcode(v_d)
        w_j = v_j  ^ as_suboctad(f, d)
        sign = m24.gcode_weight(d ^ e) + m24.gcode_weight(e)
        sign +=  m24.scalar_prod(d, eps)
        sign += m24.suboctad_scalar_prod(as_suboctad(e ^ f, d), v_j)
        sign += m24.suboctad_weight(v_j) * parity
        return (-1)**sign * v_value, 'T', v_d, w_j
    elif v_tag in 'BC':
        m = v_tag == 'C'
        c = m24.vect_to_cocode((1 << v_d) ^ (1 << v_j))
        n = m ^ m24.scalar_prod(f, c)
        w_tag = "BC"[n]
        w_i, w_j = max(v_d, v_j), min(v_d, v_j)
        sign = m * parity + m24.scalar_prod(e ^ f, c)
        return (-1)**sign * v_value, w_tag, w_i, v_j
    elif v_tag == 'A':
        w_i, w_j = max(v_d, v_j), min(v_d, v_j)
        c = m24.vect_to_cocode((1 << v_d) ^ (1 << v_j))
        sign = m24.scalar_prod(f, c)
        return (-1)**sign * v_value,'A', w_i, v_j
    else:
        raise ValueError("Bad tag " + v_tag)




def one_test_op_xy(v, eps, e, f, verbose = 0):
    pi_atom = grp('d', eps)
    x_atom = grp('x', e)**(-1)
    y_atom = grp('y', f)**(-1)
    w = v * pi_atom * x_atom  * y_atom 
    if verbose:
        print("%s * %s * %s * %s = %s" % (v, pi_atom, x_atom, y_atom, w))
    v_tuple = vector_data(v)   
    w_tuple = vector_data(w) 
    w1_tuple = op_xy(v_tuple, eps, e, f)  
    w1 = v.space(v.p, [w1_tuple])
    if w != w1:
        print("%s * %s * %s * %s = %s" % (v, pi_atom, x_atom, y_atom, w))
        print("with xy formula:", w1)
        raise ValueError("Calculation with xy formula failed")







def rand_v(tag):
    return rep_mm(tag, "r")
   

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
            v1=  rep_mm(*v)
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
    for v, eps, e, f in op_xy_testdata():
        one_test_op_xy(v, eps, e, f, verbose = verbose)
    print("passed")






