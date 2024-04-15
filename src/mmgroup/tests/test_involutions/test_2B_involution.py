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

from mmgroup import MM0, AutPL, PLoop, Cocode, Xsp2_Co1

def import_all():
    global reduce_via_power 
    from mmgroup.structures.involutions import reduce_via_power


z = MM0('x', 0x1000)

def make_2B_samples(n_samples = 5):
    for i in range(n_samples):
        g = MM0('r', 5)
        yield z**g




def do_test_involution(g, verbose = 0):
    itype, h = g.conjugate_involution(verbose = verbose)
    assert itype == 2
    assert g ** h == z




def do_test_2B_involution(n_tests = 3, verbose = 0):
    for i, g in enumerate(make_2B_samples(n_tests)):
        if verbose:
             print("Test", i+1)
        itype, h = g.conjugate_involution(verbose = verbose)
        assert itype == 2
        assert g ** h == z


@pytest.mark.involution
def test_2B_involution():
    import_all()     
    do_test_2B_involution(n_tests = 3, verbose = 0)



@pytest.mark.involution
@pytest.mark.slow
def test_2B_involution_extensive():     
    import_all()     
    do_test_2B_involution(n_tests = 10, verbose = 0)


@pytest.mark.involution
@pytest.mark.slow
@pytest.mark.very_slow
def test_reduce_via_power(verbose = 1, ntests = 10):     
    import_all()     
    start_time = datetime.datetime.now()
    header = "\nTest shortening element of the monster via powers"
    print(header)
    print("%d tests" % (ntests))
    print("started: ", start_time)

    for i in range(ntests):
        if verbose:
            print("\nTest", i+1)
        g = MM0('r', 18)
        g1 = reduce_via_power(g, ntrials=40, verbose = verbose)
        assert g1 == g
    end_time = datetime.datetime.now()
    diff_time = end_time - start_time
    t = diff_time.total_seconds() 
    print("started: ", start_time)
    print("finished:", end_time)
    print("time = %.3f s, per test: %.3f s" % (t, t/ntests))



Z_2A = MM0('d', Cocode([3, 2]))


def make_2A_samples(n_samples = 5):
    for i in range(n_samples):
        g = MM0('r', 3)
        yield Z_2A**g


def do_test_2A_involution(n_tests = 30, verbose = 1):
    for i, g in enumerate(make_2A_samples(n_tests)):
        if verbose:
            print("Test", i+1)
            print("g=", g)
        itype, h = g.conjugate_involution(ntrials = 40, verbose = verbose)
        assert itype == 1
        res = (g ** h)
        res.in_G_x0()
        if verbose:
            print("h=", h)
            print("g**h=", res)
        assert res in (Z_2A, MM0('x', 0x200) * Z_2A)


@pytest.mark.involution
def test_2A_involution():  
    import_all()     
    itype, h = MM0().conjugate_involution()
    assert (itype, h) == (0, MM0())  
    do_test_2A_involution(n_tests = 3, verbose = 0)


