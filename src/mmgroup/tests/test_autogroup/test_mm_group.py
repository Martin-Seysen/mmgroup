
import sys
import os
import numpy as np
from numbers import Integral
from random import randint



import pytest

from mmgroup.tests.groups.mgroup_n import MGroupN
from mmgroup.mm_group import MMGroup

ref_group = MGroupN(0)
group = MMGroup()

group.set_preimage(ref_group, tuple)
ref_group.set_preimage(group, tuple)


def str_group(g):
    return str(list(map(hex, g.data)))

####################################################################
# Test cases
####################################################################


def make_testcases():
    testdata = [
        (ref_group(("y", 0x1000),('p', 234579)), ("y", 0x1800)),
        (ref_group(("y", 0x1000),('t',2)), ("y", 0x1800)),
    ]
    for d in testdata:
        yield d
    for i in range(5):
        for g1 in "dpxyt":
            for g2 in "dpxyt":
                yield (g1,), (g2,)
    for i in range(150):
        g1 = ref_group.rand_word("dpxyt", 10)
        g2 = ref_group.rand_word("dpxyt", 10)
        yield g1, g2


####################################################################
# Test group multiplication
####################################################################

def ll_mul_test_case(g1, g2, verbose = 1):
    g1 = group(g1)
    g2 = group(g2)
    g1ref = ref_group(g1) 
    g2ref = ref_group(g2) 
    g12ref = g1ref * g2ref
    g12 = g1 * g2
    assert group(g12ref) == g12, (g12ref, g12)
    assert ref_group(g12) == g12ref, (g1, g2) 
    # next test mixed operation
    assert g12 == g1 * g2ref
    assert g12ref == g1ref * g2


@pytest.mark.auto_group
def test_ll_mul(verbose = 0):
    print("testing low-level multiplication in group N")
    for g1, g2 in make_testcases():
        ll_mul_test_case(g1, g2, verbose)
    print("passed")



####################################################################
# Test group division 
####################################################################

def ll_div_test_case(g1, g2, verbose = 1):
    g1 = group(g1)
    g2 = group(g2)
    g1ref = ref_group(g1) 
    g2ref = ref_group(g2) 
    g12ref = g1ref / g2ref
    g12 = g1 / g2
    assert group(g12ref) == g12, (g12ref, g12)
    assert ref_group(g12) == g12ref
    # next test mixed operation
    assert g12 == g1 / g2ref
    assert g12ref == g1ref / g2


@pytest.mark.auto_group
def test_ll_div(verbose = 0):
    print("testing low-level multiplication in group N")
    for g1, g2 in make_testcases():
        ll_div_test_case(g1, g2, verbose)
    print("passed")


