from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample


import numpy as np
import pytest

from mmgroup import MM0, MMSpace, Cocode, MMV
from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_index_leech2_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
from mmgroup.mm_op import mm_aux_index_sparse_to_leech
from mmgroup.axes import Axis, get_sample_axes
from mmgroup.axes import beautify_axis
from mmgroup.axes import reduce_axis_G_x0

V = MMV(15)
G = Axis.group


START_AXIS =  ('I', 3, 2)

CENTRAL_INVOLUTIONS_OK = [
    G('x', x) for x in [0,  0x1000, 0x800, 0x1800]
]

def check_central_involution(axis):
    g2 = axis.central_involution() * axis.g_central
    cls_, h = g2.conjugate_involution_G_x0()
    g_test = (g2**-1 * g2**h)
    orbit = axis.axis_type()
    if orbit == '8B':
         ok = any(g_test == g for g in CENTRAL_INVOLUTIONS_OK)
    elif orbit in ['6A', '10A']:
         ok = any(g_test == g for g in CENTRAL_INVOLUTIONS_OK[:2])
    else:
         ok = g_test ==  CENTRAL_INVOLUTIONS_OK[0]
    if not ok:
        raise ValueError("BAD case: %s, %s, %s" % (orbit, cls_, g_test))







@pytest.mark.axes
def test_axes():
    sample_axes = get_sample_axes()
    for key in sample_axes:
        axis = sample_axes[key]
        v = V(*START_AXIS) * axis.g
        ref_axis = Axis(axis.g)
        assert v == axis.v15 == ref_axis.v15
        assert v.axis_type() == key
        for i in range(5):
             w = axis * G('r', 'G_x0')
             assert w.v15.axis_type() == key
             e = randint(1,2)
             we =  w.axis_type(e) 
             assert isinstance(we, str)
             if i < 2:
                 w1 = beautify_axis(key, w.g, check = 1)
                 for tag in "ABC":
                     if (w1[tag] != v[tag]).any():
                         S = "Orbit %s axes differ in tag %s!!!"
                         raise ValueError(S % (key, tag))
                 check_central_involution(w1)
             if i < 2:
                 g = reduce_axis_G_x0(w, check = False)
                 assert (w * g).v15 == axis.v15
