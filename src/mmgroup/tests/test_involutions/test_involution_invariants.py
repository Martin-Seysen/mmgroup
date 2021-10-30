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

from mmgroup import MM0, AutPL, PLoop, Cocode, Xsp2_Co1, MMV, XLeech2
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
        g = MM0(g)
        g.in_G_x0()
        yield g, invar
        n_samples = 20 if invariant_status(invar) == 3 else 10
        for i in range(n_samples): 
            t = MM0('r', 'G_x0')
            h = g**t
            h.in_G_x0()
            yield  h, invar



def invariant_status(ref_invariants):
    """Return invariant status 

    Here parameter ``ref_invariants`` must be en entry of the list 
    ``INVOLUTION_SAMPLES`` defined in file ``involution_sample.py``.
    This list shows some invariants of all classes in the subgroup 
    :math:`G_x0` of the monster that square up to an element of 
    :math:`Q_x0`. Let :math:`g` be an element of the subgroup  
    :math:`G_{x0}` of the monster that has inariants as given by
    parameter  ``ref_invariants``. The function returns:

    0 if :math:`g` cannot be mapped to :math:`N_{x0}` by conjugation

    1 if :math:`g` can be mapped to :math:`N_{x0}`, but not to 
    :math:`Q_{x0}`.

    2 if math:`g` can be mapped to :math:`Q_{x0}`, but not to 
    the central involution of :math:`Q_{x0}`. Then  math:`g`
    is of type 1A, 2A, 2B or 4A in the monster.

    3 if math:`g` can be mapped to the central involution of 
    :math:`Q_{x0}`. Then  math:`g`  is of type 2B in the monster.
    """
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
        st = 3 if (ref_chi[0] == 275 and ref_ord[0] == 2) else 2
    return st


def conj_G_x0_to_Q_x0(g):
    gg = Xsp2_Co1(g)
    a = np.zeros(7, dtype = np.uint32)
    lv = chk_qstate12(xsp2co1_elem_conj_G_x0_to_Q_x0(gg._data, a))
    length, q = lv >> 25, lv & 0x1ffffff
    h = MM0('a', a[:length])
    gh = MM0(g) ** h
    assert gh.in_Q_x0()
    assert gh == MM0('q', q)
    return h

Z = MM0('x', 0x1000)

def do_test_involution_invariants(g, ref_invariants, verbose = 0):
    g.reduce()
    gg = Xsp2_Co1(g)
    ref_ord, ref_chi, ref_involution_invariants = ref_invariants
    invar, v1, v0  = gg._involution_invariants()
    assert ref_ord[0] == gg.order()
    #print("square", (g**2).reduce())
    assert ref_ord[1] == XLeech2(gg**2).type
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
        assert istate >= 2, (istate, ref_invariants)
    v = xsp2co1_elem_find_type4(gg._data)
    err = ""
    if istate == 0:
        ok = v <= 0
        if not ok: 
            err = "xsp2co1_elem_find_type4() succeeded but should not"
    else:
        ok = v > 0
        if not ok: 
            err = "xsp2co1_elem_find_type4() not successful"
        mv = MM0("c", v)**-1
        h = g * mv
        if ok and istate > 1: 
            ok = ok and  h.in_N_x0(), (g, h, ref_invariants)
            if not ok:
                err = "Could not conjugate element to N_x0"
            #print(g)
            #print("%-28s" % h, ref_invariants)
            if ok:
                h1 = conj_G_x0_to_Q_x0(g)
                ok = ok and  (g ** h1).in_Q_x0()
                if not ok:
                    err = "Could not conjugate element to Q_x0"
    if ok and istate == 3:
        _, a = gg.conjugate_involution(MM0)
        ok = ok and  g ** a == Z
        if not ok:
           err = "Could not conjugate element to centre of Q_x0" 
    if not ok:
        print("\nError in conjugating involution invariants")
        print(err)
        print("istate =", istate)
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


 



@pytest.mark.involution
def test_xx(verbose = 0):
    for i, (g, n) in enumerate(make_involution_samples()):
         if verbose:
             print("Test", i)
         do_test_involution_invariants(g, n, verbose=verbose)
