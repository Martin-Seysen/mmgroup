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
#from mmgroup.generators import gen_leech2_type
#from mmgroup.tests.test_involutions.make_involution_samples import invariant_count_type2
#from mmgroup.clifford12 import xsp2co1_elem_find_type4
#from mmgroup.clifford12 import xsp2co1_involution_find_type4
#from mmgroup.clifford12 import xsp2co1_elem_conj_G_x0_to_Q_x0
#from mmgroup.clifford12 import chk_qstate12



z = MM(('x', 0x1000))

def make_samples(n_samples = 5):
    for i in range(n_samples):
        g = MM.rand_mm(5)
        yield z**g




def do_test_involution(g, verbose = 0):
    h = g.conjugate_2B_involution(verbose = verbose)
    assert g**h == z




def do_test_2B_involution(n_tests = 3, verbose = 0):
    for i, g in enumerate(make_samples(n_tests)):
        if verbose:
             print("Test", i+1)
        h = g.conjugate_2B_involution(verbose = verbose)
        assert g**h == z


@pytest.mark.involution
def test_2B_involution():     
    do_test_2B_involution(n_tests = 3, verbose = 0)



@pytest.mark.involution
@pytest.mark.slow
def test_2B_involution_extensive():     
    do_test_2B_involution(n_tests = 10, verbose = 0)

