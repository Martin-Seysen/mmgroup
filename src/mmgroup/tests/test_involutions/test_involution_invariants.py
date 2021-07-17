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
from mmgroup.clifford12 import xsp2co1_elem_find_type4
from mmgroup.clifford12 import xsp2co1_involution_find_type4
from mmgroup.clifford12 import xsp2co1_elem_conj_G_x0_to_Q_x0
from mmgroup.clifford12 import chk_qstate12

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



def invariant_status(ref_invariants):
    ref_ord, ref_chi, ref_involution_invariants = ref_invariants
    length = ref_involution_invariants[0]
    if length <= 9:
        st = 1;
    elif length == 12:
        if ref_involution_invariants[1] & 2 == 0:
            st = 0
        else:
            st = 1
    else:
        st = 0
    if ref_chi[0] in (196883, 275, 4371):
        assert st == 1
        st = 2
    return st


def conj_G_x0_to_Q_x0(g):
    gg = Xsp2_Co1(g)
    a = np.zeros(7, dtype = np.uint32)
    lv = chk_qstate12(xsp2co1_elem_conj_G_x0_to_Q_x0(gg._data, a))
    length, q = lv >> 25, lv & 0x1ffffff
    h = MM.from_data(a[:length])
    gh = MM(g)**h
    assert gh.in_Q_x0()
    assert gh == MM(('q', q))
    return h


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

    # test conjugation
    istate = invariant_status(ref_invariants)
    if ref_chi[0] in (196883, 275, 4371):
        # The g is of type 2A, 2B or 2A in the monster
        assert istate == 2, (istate, ref_invariants)
    v = xsp2co1_elem_find_type4(gg._data)
    if istate == 0:
        ok = v <= 0
    else:
        ok = v > 0
        mv = MM(("c", v))**-1
        h = g**mv
        if istate > 1: 
            assert h.in_N_x0(), (g, h, ref_invariants)
            print(g)
            print("%-28s" % h, ref_invariants)
            h1 = conj_G_x0_to_Q_x0(g)
            """
            print(g, mv**-1)
            print("%-28s" % h, ref_invariants)
            h1 = conj_G_x0_to_Q_x0(g)
            h1.in_Q_x0()
            h2 = g **h1
            h2.in_Q_x0()
            print("%-28s" % h2, ref_invariants);print("")
            assert  h2.in_Q_x0(), (h2, h2.in_Q_x0())  
            """
            print("")
            
    if not ok:
        print("\nError in conjugating involution invariants")
        print("Invariants:", ref_invariants)
        print("Conjugtion status expected:", istate)
        print("Conjugation vector:", hex(v))
        print("Involution invariants")
        for x in invar: print(" 0x%014x" % x)
        length = ref_invariants[2][0]
        if 8 <= length <= 9: 
            coset = ref_invariants[2][4] > 0
            _invar = np.array(invar + [0]*12, dtype= np.uint64)
            v1 =  xsp2co1_involution_find_type4(_invar, coset)
            print("Low level conjugation vector:", hex(v1))
        err = "Search for conjugation vector failed"
        raise ValueError(err)


 



@pytest.mark.mmm
def test_xx(verbose = 0):
    for i, (g, n) in enumerate(make_involution_samples()):
         if verbose:
             print(i)
         do_test_involution_invariants(g, n, verbose=verbose)