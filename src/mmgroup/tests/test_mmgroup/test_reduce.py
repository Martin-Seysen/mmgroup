
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import sample
import datetime
import time


import pytest

import numpy as np
from multiprocessing import Pool, TimeoutError, cpu_count

from mmgroup.mm_reduce import mm_reduce_M
from mmgroup import MM0, MM, MMV
from mmgroup.tests.test_axes.test_reduce_axis import g_complexity

#####################################################################################
# Auxiliary fuctions
#####################################################################################

MM_TAGS = dict(enumerate(" dpxyTl?"))
for i in [0]: MM_TAGS[i] = ""

def mm_pattern(g):
    """return string of tags in monster group element"""
    return "".join([MM_TAGS[(a >> 28) & 7] for a in g.mmdata])


#####################################################################################
# Test fast reduction in monster group with C function
#####################################################################################

reduce_mm_time = None


def reduce_mm_C(g, check = True, mode = 0):
    """The fastest reduction procedure for a monster element ``g``"""
    global reduce_mm_time
 
    g1 = np.zeros(256, dtype = np.uint32)
    t_start = time.perf_counter() 
    res = mm_reduce_M(g._data, g.length, mode, g1)
    reduce_mm_time = time.perf_counter() - t_start
    if (res < 0):
        err = "Reduction of element of monster failed"
        raise ValueError(err)
    length = res
    if check:
        w = mm_vector(15)
        work = mm_vector(15)
        mm_order_load_vector(w.data)
        mm_op15_word(w, g._data, len(g), 1, work)
        mm_op15_word(w, g1, length, -1, work)
        mm_order_load_vector(work.data)
        assert not mm_op15_compare(w, work)
         
    g._extend(length)
    g._data[:length] = g1[:length]
    g.length = length
    g.reduced = 0
    g.reduce()
    return g

def reduce_testcases_C():
    for quality in range(1,16):
        for i in range(2):
              yield  MM0('r', quality), i & 1  
    for i in range(4):
        yield  MM0('r', 16), i & 1


@pytest.mark.mmgroup
def test_reduce_mm_C(verbose = 0):
    for n, (g, mode) in enumerate(reduce_testcases_C()):
        g1 = reduce_mm_C(g.copy(), check = False, mode = mode)
        ok = g == g1
        if verbose:
            print("Test", n + 1)
        if verbose or not ok:
            print("g =", g)
            print("reduced:", g1)
            print("Time: %.3f ms" % (1000 * reduce_mm_time),
                 ", complexity;", g_complexity(g), ",", g_complexity(g1))
        if not ok:
            err = "Reduction of monster group element failed"
            raise ValueError(err)



#####################################################################################
# Test fast reduction in monster group
#####################################################################################

MIN_LEN_ALWAYS_UNREDUCED = 80

def reduce_testcases(ncases = 1000):
    for i in range(ncases):
        for complexity in range(15):
            yield MM0('r', complexity)
    g = MM0('r', 18)
    assert len(g.mmdata) > MIN_LEN_ALWAYS_UNREDUCED
    yield g

def one_test_reduce(g, verbose = 0):
    g1 = MM0(g)
    g2 = MM(g).reduce()
    if verbose:
        print(g1)
        print(g2)
    if len(g1.mmdata) > MIN_LEN_ALWAYS_UNREDUCED:
        assert len(g2.mmdata) < len(g1.mmdata)
    assert MM0(g2) == g1


# Support for multiprocessing         
POOL_MAGIC = 0x23f124ee
NPROCESSES =  max(1, cpu_count() - 1)
#NPROCESSES = 1

def single_test_reduce(ncases, verbose = 0):
     for i, g in enumerate(reduce_testcases(ncases)):
         if verbose:
              print("Test", i+1)
         one_test_reduce(g, verbose = verbose)
     return POOL_MAGIC 


# The final test programm
@pytest.mark.mmgroup 
def test_reduce(ncases = 10, verbose = 0):
    if verbose or NPROCESSES <= 1:
        single_test_reduce(ncases, verbose = verbose)
        return    
    with Pool(processes = NPROCESSES) as pool:
        num_cases = (ncases - 1) // (NPROCESSES)  + 1
        testvalues = [num_cases] * NPROCESSES
        results = pool.map(single_test_reduce, testvalues)
    pool.join()
    assert results ==  [POOL_MAGIC] * NPROCESSES


#####################################################################################
# Test multiplication in monster group
#####################################################################################


N_MUL_SAMPLES = 10

def make_mul_samples(n = N_MUL_SAMPLES):
    indices = list(range(n))
    glist = [ MM('r', 'M').reduce() for i in range(n) ]
    return indices, glist


@pytest.mark.mmgroup 
def test_mul(ncases = 20):
    vtest= MMV(127)('R')
    indices, glist = make_mul_samples()
    for n in range(ncases):
        i, j =  sample(indices, 2)
        vt1 = vtest * glist[i] * glist[j]
        glist[i] *= glist[j]
        vt2 = vtest * glist[i]
        assert vt1 == vt2

#####################################################################################
# Benchmark multiplication in monster group
#####################################################################################


def benchmark_mul(ncases = 20):
    indices, glist = make_mul_samples()
    index_pairs = [sample(indices, 2) for i in range(ncases)]
    #print(glist, "\n", index_pairs)
    glist[0] *= glist[1]
    t = []
    for i, j in index_pairs:
        t_start = time.process_time()
        glist[i] *= glist[j]
        t.append(time.process_time() - t_start)
        #print(mm_pattern(glist[i]))
    n = len(index_pairs) + 0.0
    mu = sum(t) / n
    var = sum([(x - mu)**2 for x in t]) / (n-1)
    return n, mu, var**0.5


@pytest.mark.bench 
@pytest.mark.mmgroup 
def test_benchmark_mul(ncases = 100):
    print("")
    for i in range(1):
        n, mu, sigma = benchmark_mul(ncases) 
        s = "Runtime of multiplication in class MM, %d tests: %.3f+-%.3f ms" 
        print(s % (n, 1000*mu, 1000*sigma))


