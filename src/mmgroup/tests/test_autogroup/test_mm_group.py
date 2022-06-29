
import sys
import os
import numpy as np
from numbers import Integral
from random import randint
from random import choice


import pytest

from mmgroup.tests.groups.mgroup_n import MGroupN
from mmgroup.structures.mm0_group import MM0Group, MM0, StdMM0Group



####################################################################
# Conversion between MM0 (i.e. the Monster group) and MGroupN 
# (i.e. an implementation of the subgroup N_0 of the monster).
####################################################################

ref_group = MGroupN()       # Implementation of subgroup N_0 of the monster
group = MM0Group()          # Implementation of the monster group
assert group is StdMM0Group # Check that class MMGroup is a singleton.




####################################################################
# Test cases
####################################################################

def rand_word(tags, length):   
    return [(choice(tags), 'r') for i in range(length)]


def make_testcases():
    testdata = [
        ([("t", 2), ("x",0)], []),
        ([("d", 2), ("d",2)], []),
        ([("p", 220249159)], []),
        ([("y", 0x1000),('p', 234579)], [("y", 0x1800)]),
        ([("y", 0x1000),('t',2)], [("y", 0x1800)]),
        ([("x", 0x800), ('t',1)],  [("y", 0x5a2), ("y", 0x1da2)]),  
    ]
    for d in testdata:
        yield d
    for i in range(5):
        for g1 in "dpxyt":
            for g2 in "dpxyt":
                yield [(g1,'r')], [(g2,'r')]
    for i in range(150):
        g1 = rand_word("dpxyt", 10)
        g2 = rand_word("dpxyt", 10)
        yield g1, g2


####################################################################
# Test group multiplication
####################################################################

def ll_mul_test_case(g1, g2, verbose = 0):
    g1, g2 = ref_group(g1), ref_group(g2)
    g1 = MM0(g1)
    g2 = MM0(g2)
    g1ref = ref_group(g1) 
    g2ref = ref_group(g2) 
    g12ref = g1ref * g2ref
    g12 = g1 * g2
    if verbose:
        print("g1 =", g1)
        print("g2 =", g2)
        print("g1 * g2 =", ref_group(g12))
        print("ref     =", g12ref)
    assert MM0(g12ref) == g12, (g1, g2, g12ref, g12)
    assert ref_group(g12) == g12ref, (g1, g2, ref_group(g12), g12ref)


@pytest.mark.auto_group
def test_ll_mul(verbose = 0):
    print("testing low-level multiplication in group N")
    for g1, g2 in make_testcases():
        ll_mul_test_case(g1, g2, verbose)
    print("passed")



####################################################################
# Test group division 
####################################################################

def ll_div_test_case(g1, g2, verbose = 0):
    g1 = MM0(g1)
    g2 = MM0(g2)
    g1ref = ref_group(g1) 
    g2ref = ref_group(g2) 
    g12ref = g1ref / g2ref
    g12 = g1 / g2
    assert MM0(g12ref) == g12, (g12ref, g12)
    assert ref_group(g12) == g12ref


@pytest.mark.auto_group
def test_ll_div(verbose = 0):
    print("testing low-level multiplication in group N")
    for g1, g2 in make_testcases():
        ll_div_test_case(g1, g2, verbose)
    print("passed")


