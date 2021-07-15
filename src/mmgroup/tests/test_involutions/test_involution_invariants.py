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
from mmgroup.tests.test_involutions.make_involution_samples import invariant_count_type2

try:
    from mmgroup.tests.test_involutions import involution_samples
except ImportError:
    import os
    _DIR = os.path.split(__file__)[0]
    PY_FILENAME = os.path.join(_DIR, "involution_samples.py")
    from mmgroup.tests.test_involutions import make_involution_samples
    make_involution_samples.print_invariants(file = PY_FILENAME)
    from mmgroup.tests.test_involutions import involution_samples

INVOLUTION_SAMPLES = involution_samples.INVOLUTION_SAMPLES
#print(INVOLUTION_SAMPLES)





def make_involution_samples():
    for invar, _, g in INVOLUTION_SAMPLES:
        g = MM(g)
        yield g, invar
        for i in range(20): 
            t = MM.rand_G_x0()
            h = g**t
            h.in_G_x0()
            yield  h, invar


def do_test_involution_invariants(g, ref_invariants, verbose = 0):
    gg = Xsp2_Co1(g)
    ref_ord, ref_chi, ref_involution_invariants = ref_invariants
    invar, v1, v0  = gg._involution_invariants()
    assert ref_ord[0] == g.order()
    #print("square", (g**2).reduce())
    assert ref_ord[1] == (g**2).type_Q_x0()
    assert ref_chi == g.chi_G_x0()
    
    if int(invar[0]) & 0x8000000:
        if (gg**2).in_Q_x0():
            err = "Invariant error in involution"
            raise ValueError(err) 
    assert (int(invar[0]) >> 24) & 7 == ref_involution_invariants[1]
    assert (int(invar[1]) >> 24) & 7 == ref_involution_invariants[2]
    assert gen_leech2_type(v0) >> 4 == ref_involution_invariants[3]
    assert invariant_count_type2(invar) == ref_involution_invariants[4]
    invar = [int(x) for x in invar if int(x)]
    assert len(invar) == ref_involution_invariants[0]
    for x in invar[2:]:
        assert x & 0xff000000ff000000 == 0
    data =  [x for x in invar if (x >> 26) & 1 == 0]   
    preimages = [x >> 32 for x in data]
    images = [x & 0x1ffffff for x in data]
    #print("\n",  preimages, images )
    #print( list(images), list(preimages) )
    g_images = [x for x in gg.xsp_conjugate(list(preimages))]
    g_images2 = [x for x in gg.xsp_conjugate(list(images))]
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
