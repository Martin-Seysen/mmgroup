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
V = MMV(15)




START_AXIS =  ('I', 3, 2)



@pytest.mark.axes
def test_axes():
    from mmgroup.tests.axes import get_sample_axes
    from mmgroup.tests.axes import beautify_axis
    sample_axes = get_sample_axes()
    for key in sample_axes:
        axis = sample_axes[key]
        v = V(*START_AXIS) * axis.g
        assert v == axis.v15
        assert v.axis_type() == key
        for i in range(5):
             w = axis * MM0('r', 'G_x0')
             assert w.v15.axis_type() == key
             e = randint(1,2)
             we =  w.axis_type(e) 
             assert isinstance(we, str)
             if i < 3:
                 w1 = beautify_axis(key, w.g, check = 1)
                 for tag in "ABC":
                     if (w1[tag] != v[tag]).any():
                         S = "Orbit %s axes differ in tag %s!!!"
                         if key == '6F':
                             print(S % (key, tag))
                         else:
                             raise ValueError(S % (key, tag))

