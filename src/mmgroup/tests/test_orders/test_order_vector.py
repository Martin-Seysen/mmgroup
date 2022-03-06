from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import numpy as np
from numbers import Integral
from random import randint, sample


import pytest

from mmgroup.dev.mm_reduce.order_vector_tables import check_order_vector_data


@pytest.mark.orders
def test_order_vector(verbose = 0):
    print("Testing order vector in C file")
    check_order_vector_data()
    #print("passed")



