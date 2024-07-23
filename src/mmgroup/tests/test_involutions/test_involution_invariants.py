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
    time.sleep(0.1)
    from mmgroup.tests.test_involutions import involution_samples

INVOLUTION_SAMPLES = involution_samples.INVOLUTION_SAMPLES
#print(INVOLUTION_SAMPLES)





def make_involution_samples():
    # sample obtained from former error cases
    NEW_SAMPLES = [
        # The following sample led to a bug in MacOS
        (MM0("y_80fh*d_1h"), [(4, 4), (275, 43, 8, 0), (9, 4, 1, 4, 16)]),
    ]
    for h, invar in NEW_SAMPLES:
        yield  h, invar
    for invar, g in INVOLUTION_SAMPLES:
        g = MM0(g)
        g.in_G_x0()
        yield g, invar
        n_samples = 20 if invariant_status(invar) == 3 else 10
        last = g
        for i in range(n_samples): 
            t = MM0('r', 'G_x0')
            h = g**t
            h.in_G_x0()
            if h != last:
                yield  h, invar
            last = h



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
    lv = chk_qstate12(xsp2co1_elem_conj_G_x0_to_Q_x0(gg._data, a, 0))
    length, q = lv >> 25, lv & 0x1ffffff
    h = MM0('a', a[:length])
    gh = MM0(g) ** h
    assert gh.in_Q_x0()
    assert gh == MM0('q', q)
    return h

Z = MM0('x', 0x1000)

def do_test_involution_invariants(g, ref_invariants, verbose = 0):
    errors = 0
    error_text = None
    def check(condition, bit, text = None):
        nonlocal errors, error_text
        if not condition:
           errors |= 1 << bit 
           if not error_text:
               error_text = text
    g.reduce()
    if verbose:
        print("g =", g)
    gg = Xsp2_Co1(g)
    ref_ord, ref_chi, ref_involution_invariants = ref_invariants
    ref_chi = ref_chi[:4]
    invar, v1, v0  = gg._involution_invariants()
    ## invar[0] = int(invar[0]) & 0x3ffffff # this produces an error
    errors = 0
    check(ref_ord[0] == gg.order(), 0, "order of g")
    check(ref_ord[1] == XLeech2(gg**2).type, 1, "order of g**2")
    check(ref_chi == g.chi_G_x0(), 2, "characters of g")
    if int(invar[0]) & 0x8000000:
        check((gg**2).in_Q_x0(), 3, "Invariant in involution")
    check((int(invar[0]) >> 24) & 7 == ref_involution_invariants[1],
           4, "invar[0]")
    check((int(invar[1]) >> 24) & 7 == ref_involution_invariants[2],
           5, "invar[1]")
    check(gen_leech2_type(v0) == ref_involution_invariants[3],
           6, "leech2_type")
    check(invariant_count_type2(invar) == ref_involution_invariants[4],
           7, "type-2 vector count")
    invar = [int(x) for x in invar if int(x)]
    check(len(invar) == ref_involution_invariants[0],
           8, "length of invariants")
    for i, x in enumerate(invar[2:]):
        check(x & 0xff000000ff000000 == 0,
           9, "bits in invaraiant[%d]" % i)
    data =  [x for x in invar if (x >> 26) & 1 == 0] 
    pre_data =   [x for x in invar if (x >> 26) & 1 == 1]  
    preimages = [x >> 32 for x in data]
    images = [x & 0x1ffffff for x in data]
    g_images = [x for x in gg.xsp_conjugate(list(preimages))]
    g_images2 = [x for x in gg.xsp_conjugate(list(images))]
    zipp = zip(preimages, images, g_images, g_images2)
    if verbose:
       print("Involution invariants found (length = %d)" % len(invar))
       if len(pre_data) == 1:
            pre_l = [ int(pre_data[0]) >> 24, int(pre_data[0]) & 0xffffff]
            print("Prefix line", [hex(x) for x in pre_l])
       print(", ".join(["preimage", "image", "g_image", 
            "g_image2", "t"]))
    for i, (preim, im, g_im, g_im2) in enumerate(zipp):
        t = preim ^ im ^ g_im
        if verbose:
            print([hex(x) for x in (preim, im, g_im, g_im2, t)])
        check((preim ^ im ^ g_im) & 0xffffff == 0, 
            10, "(preimage ^ image ^ g_image)[%d]" % i)            
        check(g_im2 == im & 0xffffff,
            11, " g_image[%d]" % i) 
    # test conjugation
    istate = invariant_status(ref_invariants)
    v = xsp2co1_elem_find_type4(gg._data, 0)
    if not errors:
        if ref_chi[0] in (196883, 275, 4371):
            # The g is of type 1A, 2B or 2A in the monster
            check(istate >= 2,  12, "istate (= %s)" % istate)
    if not errors:   
        if istate == 0:
            check(v <= 0, 
                13, "xsp2co1_elem_find_type4(), succeeded but should not")
        else:
            check(v > 0, 
                14, "xsp2co1_elem_find_type4() not successful")
            mv = MM0("c", v)**-1
            h = g ** mv
            if errors == 0 and istate > 1: 
                check(h.in_N_x0(), 
                    15, "conjugating element to N_x0, h=" + str(h))
                if not errors:
                    h1 = conj_G_x0_to_Q_x0(g)
                    ok = (g ** h1).in_Q_x0()
                    check(ok, 16, "conjugating element to Q_x0")
    if not errors and istate == 3:
        _, a = gg.conjugate_involution(MM0)
        check(g ** a == Z, 
                17,  "conjugating element to centre of Q_x0")
    if errors:
        print("\nError in testing involution invariants")
        print("\nError status =", hex(errors))
        if error_text:
            print("Error in " + error_text)
        print("g =", g)
        print("g as an instance of Xsp2_Co1:")
        for i in range(26):
            print("  0x%016x" % gg.data[i])
        print("istate =", istate)
        print("Expected Invariants:", ref_invariants)
        print("Conjugtion status expected:", istate)
        print("Conjugation vector:", hex(v))
        print("Involution invariants obtained")
        for i, x in enumerate(invar): 
             print("%2d: 0x%014x" % (i, x))
        length = ref_invariants[2][0]
        if 8 <= length <= 9: 
            coset = ref_invariants[2][4] > 0
            _invar = np.array(invar + [0]*12, dtype= np.uint64)
            v1 =  xsp2co1_involution_find_type4(_invar, coset)
            print("Low level conjugation vector:", hex(v1))
        try:
            from mmgroup.clifford12 import xsp2co1_involution_error_pool
            error_pool = np.zeros(64, dtype = np.uint64)
            length = xsp2co1_involution_error_pool(
                error_pool, len(error_pool))
            if length:
                s = "Error pool from function xsp2co1_involution_invariants"
                print(s)
                for i, x in enumerate(error_pool[:length]): 
                    print("%2d: 0x%014x" % (i, x))
        except:
            pass
        err = "Error in involution invariants"
        raise ValueError(err)


 



@pytest.mark.involution
def test_xx(verbose = 0):
    for i, (g, n) in enumerate(make_involution_samples()):
         if verbose:
             print("Test", i)
         do_test_involution_invariants(g, n, verbose=verbose)

