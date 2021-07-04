from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral
import time
from math import floor
from random import randint, shuffle, sample
from collections import defaultdict

import numpy as np
import pytest

from mmgroup import MM, AutPL, PLoop, Cocode, Xsp2_Co1
from mmgroup.generators import gen_leech2_type

try:
    from mmgroup.tests.test_involutions.Co1_2A_involutions import Co1_class_2A_samples
except ImportError:
    from mmgroup.tests.test_involutions.make_Co1_2A_involutions import write_samples
    write_samples() 
    from mmgroup.tests.test_involutions.Co1_2A_involutions import Co1_class_2A_samples




def make_involution_samples():
    g0 = MM(("y", range(8)))
    for i in range(100):
        yield g0 ** MM.rand_G_x0() * MM(("x","r"), ("d","r")), 8
    for data in Co1_class_2A_samples:
        g = MM(data)
        for i in range(10): 
            t = MM.rand_G_x0()
            h = g**t
            h.in_G_x0()
            yield  h, 8 


def do_test_involution_invariants(g, n, verbose = 0):
    gg = Xsp2_Co1(g)
    invar, v1, v0 = gg._involution_invariants()
    if int(invar[0]) & 0x8000000:
        if (g**2).in_Q_x0():
            err = "Invariant error in involution"
            raise ValueError(err) 
    invar = [int(x) for x in invar if int(x)]
    for x in invar[2:]:
        assert x & 0xff000000ff000000 == 0
     
    data = [(x >> 32, x & 0x1ffffff) for x in invar if (x>>26) & 1 == 0]
    if n is not None:
        assert n == len(data)
    #print(data)
    preimages, images = zip(*data)
    g_images = [x for x in gg.xsp_conjugate(preimages)]
    g_images2 = [x for x in gg.xsp_conjugate(images)]
    for preim, im, g_im, g_im2 in zip(preimages, images, g_images, g_images2):
        t = preim ^ im ^ g_im
        if verbose:
            print([hex(x) for x in (preim, im, g_im, g_im2, t)])
        assert (preim ^ im ^ g_im) & 0xffffff == 0
        assert g_im2 == im & 0xffffff
    


@pytest.mark.mmm
def test_xx(verbose = 0):
    for i, (g, n) in enumerate(make_involution_samples()):
         if verbose:
             print(i)
         do_test_involution_invariants(g, n, verbose=verbose)
