from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import numpy as np
from numbers import Integral
from random import randint, sample


import pytest

def check_order_vector_data():
    from mmgroup.mm_reduce import  mm_order_compare_vector
    from random import randint
    from mmgroup.structures.mm0_group import MM0
    from mmgroup.mm_space import MMV, MMVector, MMSpace
    #from mmgroup.dev.mm_reduce.order_vector import get_order_vector
    #from mmgroup.dev.mm_reduce.order_vector import order_vector_from_dict
    from mmgroup.dev.mm_reduce.order_vector import OrderVectorMod15
    MMV15 = MMV(15)
    ov = OrderVectorMod15('C')
    ov.check()
    ov_ref = ov.order_vector

    if (mm_order_compare_vector(ov_ref.data)):
       err = "Comparing with order vector failed" 
       raise ValueError(err) 
    positions = [(0,0), (72*2, 63), (72*2 + 759*4 - 1, 35)]
    for i in range(30, len(ov_ref.data) - 1, 537):
        bmax = 63 if 72*2 < i < 72*2 + 759*4 else 31
        positions.append( (i, randint(0,bmax)) ) 
    for i, sh in positions:
        a = np.copy(ov_ref.data)
        a[i] = int(a[i]) ^ (1 << sh)
        result = mm_order_compare_vector(a)  
        if (result != 1):
            err = "Comparing of bad vector with order vector failed" 
            raise ValueError(err) 
              





@pytest.mark.orders
def test_order_vector(verbose = 0):
    print("Testing order vector in C file")
    check_order_vector_data()
    #print("passed")



