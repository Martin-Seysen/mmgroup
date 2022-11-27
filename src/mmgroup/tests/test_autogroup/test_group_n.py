
import sys
import os
import numpy as np
from numbers import Integral
from random import randint, choice



import pytest

from mmgroup.tests.groups.mgroup_n import MGroupN
from mmgroup.tests.groups.group_n import GroupN

ref_group = MGroupN()
group = GroupN()



def str_group(g):
    return str(list(map(hex, g.data)))

####################################################################
# Test cases
####################################################################


def make_rand_word(tags, length):   
    return [(choice(tags), 'r') for i in range(length)]

def make_testcases():
    for i in range(5):
        for g1 in "dpxyt":
            for g2 in "dpxyt":
                yield [(g1,'r')], [(g2,'r')]
    for i in range(150):
        g1 = make_rand_word("dpxyt", 10)
        g2 = make_rand_word("dpxyt", 10)
        yield g1, g2


####################################################################
# Test group multiplication
####################################################################

def ll_mul_test_case(g1, g2, verbose = 0):
    g1 = group(g1)
    g2 = group(g2)
    g1ref = ref_group(g1) 
    g2ref = ref_group(g2) 
    g12ref = g1ref * g2ref
    g12 = g1 * g2
    assert group(g12ref) == g12, (g12ref, g12)
    assert ref_group(g12) == g12ref
    g12mul = g1.copy().mul(g2)
    assert g12 == g12mul
    g12mul_atoms = g1.copy().mul_atoms(g2)
    assert g12 == g12mul_atoms
    # next test mixed operation
    #assert g12 == g1 * g2ref
    #assert g12ref == g1ref * g2
    g1_std = group('a', g1.std_word())
    assert g1_std == g1, (g1_std, g1)


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
    g1 = group(g1)
    g2 = group(g2)
    g1ref = ref_group(g1) 
    g2ref = ref_group(g2) 
    g12ref = g1ref / g2ref
    g12 = g1 / g2
    assert group(g12ref) == g12, (g12ref, g12)
    assert ref_group(g12) == g12ref
    g12div = g1.copy().div(g2)
    assert g12 == g12div
    g12div_atoms = g1.copy().div_atoms(g2)
    assert g12 == g12div_atoms
    # next test mixed operation
    #assert g12 == g1 / g2ref
    #assert g12ref == g1ref / g2


@pytest.mark.auto_group
def test_ll_div(verbose = 0):
    print("testing low-level multiplication in group N")
    for g1, g2 in make_testcases():
        ll_div_test_case(g1, g2, verbose)
    print("passed")


####################################################################
# Test coset decomposition
####################################################################

@pytest.mark.auto_group
def test_coset_decomp(verbose = 0):
    for i in  range(50):
        g = group(make_rand_word("dpxyt", 10))
        g_Nx0, g_t = g.right_coset_N_x0()
        assert g == g_Nx0 * g_t
        assert g_Nx0.data[0] == 0
        d = g_t.data 
        assert d[1] | d[2] | d[3] | d[4] == 0
        



