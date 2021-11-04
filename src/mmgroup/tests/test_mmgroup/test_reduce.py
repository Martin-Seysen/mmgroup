
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import sample
import datetime
import time


import pytest

from mmgroup import MM0, MM, MMV
from multiprocessing import Pool, TimeoutError, cpu_count

#####################################################################################
# Test fast reduction in monster group
#####################################################################################

MIN_LEN_ALWAYS_UNREDUCED = 80

def reduce_testcases(ncases = 1000):
    for i in range(ncases):
        for complexity in range(15):
            yield MM0('r', complexity)
    g = MM0('r', 18)
    assert len(g) > MIN_LEN_ALWAYS_UNREDUCED
    yield g

def one_test_reduce(g, verbose = 0):
    g1 = MM0(g)
    g2 = MM(g).reduce()
    if verbose:
        print(g1)
        print(g2)
    if len(g1) > MIN_LEN_ALWAYS_UNREDUCED:
        assert len(g2) < len(g1)
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
    t_start = time.process_time()
    for i, j in index_pairs:
        glist[i] *= glist[j]
    return time.process_time() - t_start


@pytest.mark.bench 
@pytest.mark.mmgroup 
def test_benchmark_mul(ncases = 20):
    print("")
    for i in range(1):
        t_ms = 1000.0 * benchmark_mul(ncases) / ncases
        s = "Runtime for multiplication in class MM: %.3f ms" 
        print(s % t_ms)


