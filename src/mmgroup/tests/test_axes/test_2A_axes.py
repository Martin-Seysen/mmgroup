from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample


import numpy as np
import pytest

from mmgroup import MM0, MMSpace, Cocode, MMV
from mmgroup.mm import mm_aux_index_extern_to_sparse
from mmgroup.mm import mm_aux_index_leech2_to_sparse
from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup.mm import mm_aux_index_sparse_to_leech
V = MMV(15)




START_AXIS =  ('I', 3, 2)

from mmgroup.tests.test_axes.get_sample_axes import import_sample_axes
sample_axes = import_sample_axes()

from mmgroup.tests.test_axes.test_import import AXES, BABY_AXES
# We have(V(START_AXIS) * MM(AXES[key]).axis_type() == key
# for all kes in the dictionar AXES 




@pytest.mark.axes
def test_axes():
    for key in AXES:
        v = V(*START_AXIS) * MM0(AXES[key])
        assert v.axis_type() == key
        for i in range(20):
             w = v * MM0('r', 'G_x0')
             assert w.axis_type() == key
             if (i > 5): continue
             e = randint(1,2)
             we =  w.axis_type(e) 
             assert isinstance(we, str)

