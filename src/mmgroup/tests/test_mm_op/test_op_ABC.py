from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from random import randint
from importlib import import_module
import time
from functools import partial

from mmgroup import MM0, MMV, AutPL
from mmgroup.mm_space import characteristics 
from mmgroup.mm_op import mm_group_prepare_op_ABC
from mmgroup.mm_op import mm_op_word_ABC


import pytest




#######################################################################
# Test function mm_group_prepare_op_ABC
#######################################################################

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
    res = mm_op_word_ABC(p, v.data, g.mmdata, len(g.mmdata), w.data)
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


#######################################################################
# Benchmark for function mm_group_prepare_op_ABC
#######################################################################


def make_benchmark_sample(n):
    a = np.zeros((n,6), dtype = np.uint32)
    for i in range(n):
        m = MM0('c', 'r').mmdata
        np.copyto(a[i,:len(m)], m)
    return a

def benchmark_mm_op15_word_ABC(a):
    f = partial(mm_op_word_ABC, 15)
    v = MMV(15)('R').data
    w = MMV(15)().data
    t_start = time.process_time()
    for g in a:
        f(v, g, 6, w)
    t = time.process_time() - t_start
    return t / len(a)
    
def benchmark_mm_op15_word(a):
    from mmgroup.mm_op import mm_op_word as f
    v = MMV(15)('R').data
    w = MMV(15)().data
    t_start = time.process_time()
    for g in a:
        f(15, v, g, 6, 1, w)
    t = time.process_time() - t_start
    return t / len(a)


def benchmark_mm_op15_map_t(a):
    from collections import defaultdict
    from mmgroup.mm_reduce import mm_reduce_op_2A_axis_type
    d = defaultdict(int)
    axes = []
    for i in range(64):
        axes.append( (MMV(15)('I', 2, 3) * MM0('r', 5)).data )
    mask = len(axes) - 1
    assert mask & (mask + 1) == 0
    w = MMV(15)().data
    w1 = np.zeros(24*4, dtype = np.uint64)
    t_start = time.process_time()
    mode = 0x16
    for i, g in enumerate(a):
        axtypes = mm_reduce_op_2A_axis_type(axes[i & mask], g, 6, mode)
        d[(axtypes >> 8) & 0xff] += 1
        d[(axtypes >> 16) & 0xff] += 1
    t = time.process_time() - t_start
    #print(dict(d))
    return t / len(a)



def benchmark_empty(a):
    v = MMV(15)().data
    w = MMV(15)().data
    t_start = time.process_time()
    for g in a:
        pass
    t = time.process_time() - t_start
    return t / len(a)



@pytest.mark.slow
@pytest.mark.bench
@pytest.mark.mm_op
def test_benchmark_op_word():
    n = 2000
    a = make_benchmark_sample(n)
    t_0 =  benchmark_empty(a)
    t_f =  benchmark_mm_op15_word(a)
    t_fABC =  benchmark_mm_op15_word_ABC(a)
    t_ax =  benchmark_mm_op15_map_t(a)
    print("Average run times:")
    names = ["mm_op15_word", 
       "mm_reduce_op_2A_axis_type",
       "mm_op15_word_ABC",
       "empty loop",
    ]
    for name, value in zip(names, [t_f, t_ax, t_fABC, t_0]):
        print(" Function %-25s: %7.4f ms" % (name, 1000 * value))






    

