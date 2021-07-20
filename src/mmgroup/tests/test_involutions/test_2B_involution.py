from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral
import datetime
import time
from math import floor
from random import randint, shuffle, sample
from collections import defaultdict

import numpy as np
import pytest

from mmgroup import MM, AutPL, PLoop, Cocode, Xsp2_Co1
from mmgroup.structures.involutions import reduce_via_power


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


@pytest.mark.involution
@pytest.mark.slow
@pytest.mark.very_slow
def test_reduce_via_power(verbose = 1, ntests = 10):     
    start_time = datetime.datetime.now()
    header = "\nTest shortening element of the monster via powers"
    print(header)
    print("%d tests" % (ntests))
    print("started: ", start_time)

    for i in range(ntests):
        if verbose:
            print("\nTest", i+1)
        g = MM.rand_mm(18)
        g1 = reduce_via_power(g, ntrials=20, verbose = verbose)
        assert g1 == g
    end_time = datetime.datetime.now()
    diff_time = end_time - start_time
    t = diff_time.total_seconds() 
    print("started: ", start_time)
    print("finished:", end_time)
    print("time = %.3f s, per test: %.3f s" % (t, t/ntests))

