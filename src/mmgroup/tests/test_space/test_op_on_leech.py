
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import collections
import re
import warnings
from random import randint, shuffle
import pytest
import numpy as np

from mmgroup import MMSpace, MMGroup, MM0, MMV
from mmgroup.tests.spaces.leech_space import mul24_mmgroup


M = MM0
V = MMV(3)


def check_abs_equal(v, v_ref):
    if (v == v_ref).all():
        return 
    if (v == -v_ref).all():
        return 
    print("Leech lattice vector expected:")
    print(v_ref)
    print("Leech lattice vector obtained:")
    print(v)
    raise ValueError("Leech vectors are not equal or opposite")
        

def do_test_leech_variants(v, v_leech):
    index = MMSpace.index_to_linear(v)
    v1 = MMSpace.index_to_short(index)
    check_abs_equal(v1, v_leech)
    tuple_= MMSpace.index_to_tuple(v)
    v2 = MMSpace.index_to_short(*tuple_)
    check_abs_equal(v2, v_leech)



def one_test_leech_op(v, g, verbose = 0):
    """``v`` a vector in ``V``, ``g`` an element of ``M``""" 
    #assert v in V
    #assert g in M
    v_leech =  MMSpace.index_to_short(v)
    do_test_leech_variants(v, v_leech)
    vg = v * g
    vg_leech_ref = MMSpace.index_to_short(vg)
    if verbose:
        print("\ntest ", g)
        print(v_leech)
        print(vg_leech_ref)
    vg_leech = mul24_mmgroup(v_leech, g)
    if verbose:
        print(vg_leech)
    check_abs_equal(vg_leech, vg_leech_ref)


def leech_op_testdata(n_tests = 100):
    for i in range(n_tests):
        for vtag in "BCTTTXXX":
            v = V(vtag,'r')
            l = [x for x in  "dpyxl"]
            shuffle(l)
            g = M([(x, 'r') for x in l])
            yield v, g
    



@pytest.mark.space
def test_leech_op_xy(verbose = 0):
    print("Check vector against Leech lattice operation ...")
    for v,g in leech_op_testdata():
        one_test_leech_op(v, g)
    print("passed")






