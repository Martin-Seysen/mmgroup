
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import sample
import datetime
import time
from collections import defaultdict

import pytest

import numpy as np
from multiprocessing import Pool, TimeoutError, cpu_count

from mmgroup import MM0, MM, MMV, Xsp2_Co1

def import_all():
    global mm_reduce_M 
    global mm_reduce_set_order_vector_mod15 
    global g_complexity 
    from mmgroup.mm_reduce import mm_reduce_M
    from mmgroup.mm_reduce import mm_reduce_set_order_vector_mod15
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
    for i in range(30):
          yield  MM0('r','N_0') * MM0('r','G_x0'), i & 1
    for quality in range(1,16):
        for i in range(2):
              yield  MM0('r', quality), i & 1  
    for i in range(4):
        yield  MM0('r', 16), i & 1


@pytest.mark.mmgroup
def test_reduce_mm_C(verbose = 0):
    import_all()
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
    import_all()
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
# Test fast reduction with function  mm_order_find_Gx0_via_v1_mod3
#####################################################################################

def G_x0_samples(n = 100):
    for i in range(n):
         yield  Xsp2_Co1('r', 'G_x0')

@pytest.mark.mmgroup 
def test_mm_order_find_Gx0_via_v1_mod3(verbose = 0):
    import_all()
    from mmgroup.mm_reduce import mm_order_load_vector_v1_mod3
    from mmgroup.mm_reduce import mm_order_find_Gx0_via_v1_mod3
    from mmgroup.mm_reduce import mm_order_compare_v1_mod3
    MMV3 = MMV(3) 
    v = MMV3()
    g_one = Xsp2_Co1()
    gi_buf = np.zeros(10, dtype = np.uint32)
    print("\nTesting function mm_order_find_Gx0_via_v1_mod3")
    for i, g in  enumerate(G_x0_samples()):
        if verbose:
            print("Test %d, g = %s" % (i+1, g))
        mm_order_load_vector_v1_mod3(v.data)
        v *= g
        length = mm_order_find_Gx0_via_v1_mod3(v.data, gi_buf)
        gi = Xsp2_Co1('a', gi_buf[:length])
        res = mm_order_compare_v1_mod3((v * gi).data)
        assert res == 0
        assert g * gi == g_one
        
         
        




#####################################################################################
# Test multiplication in monster group
#####################################################################################


N_MUL_SAMPLES = 30

_mul_samples = None

def get_mul_samples():
    global _mul_samples
    n = N_MUL_SAMPLES 
    if _mul_samples is None: 
        indices = list(range(n))
        glist = [ MM('r', 'M').reduce() for i in range(n) ]
        _mul_samples = indices, glist
    return _mul_samples



@pytest.mark.mmgroup 
def test_mul(ncases = 20):
    import_all()
    vtest= MMV(127)('R')
    indices, glist = get_mul_samples()
    for n in range(ncases):
        i, j =  sample(indices, 2)
        vt1 = vtest * glist[i] * glist[j]
        glist[i] *= glist[j]
        vt2 = vtest * glist[i]
        assert vt1 == vt2

#####################################################################################
# Benchmark multiplication in monster group
#####################################################################################


class MM_Pattern():
    def __init__(self):
        self.d = defaultdict(int)
        self.n = 0
    def add_pattern(self, s):
        for tag in s:
            self.d[tag] += 1
        self.n += 1
    def stat(self, tag):
        return self.d[tag] / self.n
    def stat_as_str(self):
        freq = ["%s:%.2f" % (tag,self.stat(tag)) for tag in "dxyplT"]
        return ", ".join(freq)  


def benchmark_mul(ncases = 20, verbose = 0, order_vector_mod15 = 0):
    indices, glist = get_mul_samples()
    index_pairs = [sample(indices, 2) for i in range(ncases)]
    glist[0] *= glist[1]
    glist[0].reduce()
    stat = MM_Pattern()
    t = []
    mm_reduce_set_order_vector_mod15(order_vector_mod15)
    for i, j in index_pairs:
        t_start = time.process_time()
        glist[i] *= glist[j]
        glist[i].reduce()
        t.append(time.process_time() - t_start)
        pattern = mm_pattern(glist[i])
        stat.add_pattern(pattern)
        if verbose:
            print(pattern)
    mm_reduce_set_order_vector_mod15(0)
    n = len(index_pairs) + 0.0
    mu = sum(t) / n
    var = sum([(x - mu)**2 for x in t]) / (n-1)
    return n, mu, var**0.5, stat.stat_as_str()


@pytest.mark.bench 
@pytest.mark.mmgroup 
def test_benchmark_mul(ncases = 100, verbose = 0, order_vector_mod15 = 0):
    import_all()
    print("")
    for i in range(1):
        n, mu, sigma, stat = benchmark_mul(ncases, verbose = verbose,
                               order_vector_mod15 = order_vector_mod15) 
        s = "Runtime of multiplication in class MM, %d tests: %.3f+-%.3f ms" 
        print(s % (n, 1000*mu, 1000*sigma))
        print("Average number of tags per reduced group element:")
        print(stat)



@pytest.mark.bench 
@pytest.mark.mmgroup 
def test_benchmark_mul_G_x0(ncases = 5000, verbose = 0):
    import_all()
    glist = []
    for i in range(100):
        glist.append(MM('r', 'G_x0'))
        glist[i].reduce()
    indices = range(len(glist))
    index_pairs = [sample(indices, 2) for i in range(ncases)]
    t_start = time.process_time()
    for i, j in index_pairs:
        glist[i] *= glist[j]
        glist[i].reduce()
    t = time.process_time() - t_start
    s = "\nRuntime of multiplication in subgroup G_x0 in class MM,"
    s += " %d tests: %.3f us" 
    print(s % (ncases, 1.0e6*t/ncases))




