from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import numpy as np
from random import randint

import pytest

from mmgroup.mm import mm_sub_test_prep_xy
from mmgroup import mat24 as m24

from mmgroup.tests.spaces.spaces import MMTestSpace






def _as_suboctad(v1, o):
    d = m24.octad_to_gcode(o)
    c = m24.ploop_cap(v1, d)
    return m24.cocode_to_suboctad(c, d)


class prep_xy:
    def __init__(self, eps, e, f):
        self.f = f & 0x1fff
        self.e = e & 0x1fff
        self.eps = eps = eps & 0xfff
        self.odd = eps & 1
        lin = np.zeros(6, dtype = np.uint32) 
        mm_sub_test_prep_xy(f, e, eps, 1, lin)
        self.lin_i = lin[:3]
        self.lin_d = lin[3:6]
        self.sign_XYZ = np.zeros(2048, dtype = np.uint32) 
        mm_sub_test_prep_xy(f, e, eps, 2, self.sign_XYZ)
        self.s_T = np.zeros(759, dtype = np.uint32) 
        mm_sub_test_prep_xy(f, e, eps, 3, self.s_T)


    def inv_op_unit(self, tag, d, j):
        if tag == 'X':
            tag1 = 'X' 
            d1 = d ^ self.lin_d[0] 
            j1 = j
            sign =  (self.sign_XYZ[d] & 1)
            sign ^= (self.lin_i[0] >> j) & 1 
            if self.odd:
                cc = m24.vect_to_cocode(1 << j)
                sign ^=  m24.scalar_prod(d << 1, cc)
        elif tag in 'ZY':
            s = self.odd ^ (tag == 'Y')
            tag1 = 'ZY'[s]
            s += 1
            d1 = d ^ self.lin_d[s]
            j1 = j
            sign =  (self.sign_XYZ[d] >> s) & 1
            sign ^= (self.lin_i[s] >> j) & 1 
        elif tag == 'T':
            tag1 = 'T'
            d1 = d
            te = self.s_T[d] 
            so_exp =  _as_suboctad(self.f, d)
            assert  te & 0x3f == so_exp , (hex(te), hex(so_exp))
            j1 = j ^ (te  & 0x3f)
            sign = m24.suboctad_scalar_prod(j, (te >> 8) & 0x3f)
            sign ^= (te >> 14) & 1
            sign ^= m24.suboctad_weight(j) & self.eps & 1
            assert ((te >> 15) ^ self.eps) & 1 == 0
        else:
            raise ValueError("Illegal tag " + str(tag))
        return sign & 1, tag1, d1, j1

    def inv_op(self, v):
        w = v.space.zero()
        p = v.space.p
        #for sp, value in (v.data.items()): 
        #    tag, d, j = v.space.sparse_to_tag(sp)
        for value, tag, d, j in v.as_tuples():
            sign, tag, d, j = self.inv_op_unit(tag, d, j) 
            if sign & 1: 
                value = -value % p
            #w.add_monomial(value, (tag, d, j))
            w += value * w.space.unit(tag, d, j)
        return w

    def check_v(self, v, verbose = 0):
        grp = v.space.group
        delta_atom = grp.atom('d', self.eps)
        x_atom = grp.atom('x', self.e)**(-1)
        y_atom = grp.atom('y', self.f)**(-1)
        w_ref = v * delta_atom * x_atom  * y_atom 
        w = self.inv_op(v)
        error = w != w_ref
        if error or verbose:
            eps, e, f = self.eps, self.e, self.f 
            print("vector", v)
            print("operation", "d_%x * x_%x * y_%x" % (eps, e, f))
            print("obtained:", w)
            if error:
                print("expected:", w_ref)
                raise ValueError("x-y operation failed")
                print("Error: x-y operation failed!!!")


def as_vector(x, space):
    v = space.zero()
    if isinstance(x, str):
       for tag in x:
            v += space.unit(tag, 'r')
       return v
    if isinstance(x, tuple):
       x = [x]
    if isinstance(x, list):
       for y in x:
           value = 1 if len(y) <= 3 else y[3] 
           #v.add_monomial(value, tuple(y[:3]))  
           tag, d, j = tuple(y[:3])
           v += value * v.space.unit(tag, d, j)
       return v
    raise TypeError("Bad type for vector of rep")


p = 7
space = MMTestSpace(p)


def prep_xy_testcases():
    testcases = [
       [ [("X", 3, 6)],  0, 0, 0x1171 ],
       [ [("X", 3, 6)],  12, 0, 0 ],
       [ [("X", 3, 6)],  12, 1111, 0 ],
       [ [("X", 3, 6)],  12, 0, 1111],
    ] 
    for v, eps, e, f in testcases:
        yield as_vector(v, space), prep_xy(eps, e, f) 
    v_tags = "TZYX"
    for v in v_tags:
        for i in range(1000):
            v1 =  as_vector(v, space)
            eps = randint(0, 0xfff) 
            e = randint(0, 0x1fff)
            f = randint(0, 0x1fff)
            yield v1,  prep_xy(eps, e, f)


@pytest.mark.mm
def test_prep_xy(verbose = 0):
    print("Testing preparation of operation x-y...")
    for v, op in prep_xy_testcases():
        op.check_v(v, verbose = verbose)
        if verbose: print("")
    print("passed")
    

    

       



    