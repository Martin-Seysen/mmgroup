from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import pytest

from mmgroup.dev.mm_basics.mm_basics import MM_Basics
from mmgroup import mm

@pytest.mark.mm_op
def test_mm_aux_mmv_size():
    for i in range(2,9):
        p = 2**i - 1
        res = mm.mm_aux_mmv_size(p)
        ref = MM_Basics.sizes(p)["MMV_INTS"]
        assert res == ref, (res, ref, p)
    print("Test of function mm.mm_aux_mmv_size() passed")


