
import sys
import os
import numpy as np
from numbers import Integral
from random import randint



import pytest

from mmgroup.tests.groups.mgroup_n import MGroupN
from mmgroup.mm_group import MMGroup, MM


#######################################################################
# The following statement would be bad style,
# since ``group`` is a singleton.
#######################################################################
#group.set_preimage(ref_group, tuple)


####################################################################
# Conversion between MM (i.e. the Monster group) and MGroupN 
# (i.e. an implementation of the subgroup N_0 of the monster).
####################################################################

ref_group = MGroupN(0) # Implementation of subgroup N_0 of the monster
group = MMGroup()      # Implementation of the monster group
assert group is MM     # Check that class MMGroup is a singleton.

##################################################################
# We can convert from MM to MMGroup automatically.
ref_group.set_preimage(group, tuple)
# Here is the conversion function from MM to MMGroup
to_ref_group = ref_group  

##################################################################
# For automatic conversion from MMGroup to MM we would have to
# manipulate the instance MM of class MMGroup. This is bad style, 
# since MMGroup is a singleton. 
# So we write an explicit conversion function here.

def to_group(target_group, g):
    try:
        return target_group(g)
    except:
        return target_group(*(g.reduce(True).as_tuples()))

def to_MM(g):
    return to_group(MM, g)


####################################################################
# Test cases
####################################################################


def str_group(g):
    return str(list(map(hex, g.data)))

####################################################################
# Test cases
####################################################################


def make_testcases():
    testdata = [
        (ref_group(("y", 0x1000),('p', 234579)), ("y", 0x1800)),
        (ref_group(("y", 0x1000),('t',2)), ("y", 0x1800)),
        (ref_group(("x", 0x800)), ('t',1)),
        (ref_group(("y", 0x5a2)), ("y", 0x1da2)),  
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
    g1 = to_MM(g1)
    g2 = to_MM(g2)
    g1ref = to_ref_group(g1) 
    g2ref = to_ref_group(g2) 
    g12ref = g1ref * g2ref
    g12 = g1._mul(g2)
    assert to_MM(g12ref) == g12, (g12ref, g12)
    assert to_ref_group(g12) == g12ref, (g1, g2) 
    # next test mixed operation
    #assert g12 == g1 * g2ref # does not work with singleton MMGroup
    assert g12ref == g1ref * g2
    assert g12ref == g1 * g2ref


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
    g1 = to_MM(g1)
    g2 = to_MM(g2)
    g1ref = to_ref_group(g1) 
    g2ref = to_ref_group(g2) 
    g12ref = g1ref / g2ref
    g12 = g1._div(g2)
    assert to_MM(g12ref) == g12, (g12ref, g12)
    assert to_ref_group(g12) == g12ref
    # next test mixed operation
    #assert g12 == g1 / g2ref  # does not work with singleton MMGroup
    assert g12ref == g1ref / g2


@pytest.mark.auto_group
def test_ll_div(verbose = 0):
    print("testing low-level multiplication in group N")
    for g1, g2 in make_testcases():
        ll_div_test_case(g1, g2, verbose)
    print("passed")


