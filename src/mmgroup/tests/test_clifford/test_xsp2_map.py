from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral

import numpy as np
import pytest

from mmgroup import Xsp2_Co1, XLeech2
from mmgroup.clifford12 import xsp2co1_Co1_get_mapping






#####################################################################
# Creating test vectors
#####################################################################


STD_BASIS = [1 << i  for i in range(24)]


def create_test_bases_get_mapping():
    def ar(m):
        return np.array(m, dtype = np.uint64)
    def make_from(c, t):
        m =  c.xsp_conjugate(STD_BASIS)
        m1 = t.xsp_conjugate(STD_BASIS)
        m2 = c.xsp_conjugate(m1)
        #assert m2 == t.c.xsp_conjugate(m)
        return ar(m1), ar(m2), ar(m)


    testdata = [
        [[1 << i for i in range(24)]]*3
    ]
    for m1, m2, m in testdata:
        yield ar(m1), ar(m2), ar(m)
    for i in range(100):
        c = Xsp2_Co1('r', 'G_x0') 
        t = Xsp2_Co1('r', 'G_x0')
        yield make_from(c, t)
 



@pytest.mark.xsp2co1
def test_xsp2co1_Co1_get_mapping(verbose = 0):
    for n, (m1, m2, mref) in enumerate(create_test_bases_get_mapping()):
        m = np.zeros(24, dtype = np.uint64)
        res = xsp2co1_Co1_get_mapping(m1, m2, m)

        ok = res == 0 and (m == mref).all()
        if verbose or not ok:
            print("Test", n)
            print('m1='   , m1)
            print('m2='   , m2)
            print('m_obt=', m)
            print('m_exp=', mref)
            if not ok:
                if (res): 
                    s = "Function xsp2co1_Co1_get_mapping returns %d"
                    print(s % res)
                ERR = "Error in conjugating with element l"
                raise ValueError(ERR)


            
