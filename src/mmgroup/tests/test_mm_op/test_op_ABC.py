from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from random import randint
from importlib import import_module

from mmgroup import MM0, MMV, AutPL
from mmgroup.mm_space import characteristics 
from mmgroup.mm import mm_group_prepare_op_ABC


import pytest

OP_ABC = {}
for p in characteristics():
    mm = import_module('mmgroup.mm%d' % p)
    OP_ABC[p] = mm.op_word_ABC


def perm(src, dest, d = 0):
    return MM0(AutPL(d, zip(list(src),list(dest)), 0))


def op_word_ABC_testcases(ntests = 10):
    test_values = [
        (3, [('A',2,3)], perm([2,3],[0,1])),
        (3, [('A',2,3)], [('l', 1)]),
    ]
    g_once = [
        [],
        [('x','n'),('d', 0x837)],
        [('y','n')],
        [('t','n')],
        [('p','n')],
        [('l','n')],
    ] 
    g_more = [
        [('r','N_x0')],
        [('r','G_x0')],
        [('r','G_x0'), ('r','N_0')],
    ] 
    for p, vs, gs in test_values:
        yield MMV(p)(vs), MM0(gs)
    for g_values, n in zip([g_once, g_more], [1, ntests]):
        for i in range(n):
            for p in characteristics():
                v = MMV(p)('R')
                for g in g_values:
                    yield v, MM0(g)


def display_mm_group_prepare_op_ABC(g):
    a1 = np.zeros(12, dtype = np.uint32)
    len_a1 = mm_group_prepare_op_ABC(g.mmdata, len(g.mmdata), a1)
    g1 = MM0('a', a1[:len_a1])
    print("mm_group_prepare_op_ABC(g)=")
    print(g1)
    if g1 != g:
        print(" This difffers from g!!")
    return g1


def one_test_op_word_ABC(v, g, verbose = 0):
    p = v.p
    if verbose:
        print("Test function mm_op%d_word_ABC, g =" % p)
        print(" ", g)
    w_ref = v * g
    w = MMV(p)()
    res = OP_ABC[p](v.data, g.mmdata, len(g.mmdata), w.data)
    ok = res >= 0
    for tag in "ABC":
        ok = ok and not (w[tag] != w_ref[tag]).any()
        if not ok:
            print("Test function mm_op%d_word_ABC, g =" % p)
            print(" ", g)
            print("mm%d_op_word_ABC returns %d" % (p, res))
            print("Arrays differ in tag %s, diff =" % tag)
            delta = w[tag] - w_ref[tag]
            print(delta)
            display_mm_group_prepare_op_ABC(g)
            if res != 0:
                err = "Function mm_op%d_word_ABC returns %d" % (p, res)
            else:
                err = "Error in function mm_op%d_word_ABC" % p
            raise ValueError(err)
    

@pytest.mark.mm_op
def test_op_word_ABC(verbose = 0):
    print("Testing function mm_op<p>_word_ABC")
    for v, g in op_word_ABC_testcases():
        one_test_op_word_ABC(v, g, verbose = verbose)
    print("Test passed")





