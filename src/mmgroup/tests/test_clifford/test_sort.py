
import time
from random import randint, sample
import numpy as np
import pytest

from mmgroup.clifford12 import bitvector32_heapsort
from mmgroup.clifford12 import bitvector64_heapsort
from mmgroup.clifford12 import bitvector32_sort
from mmgroup.clifford12 import bitvector64_sort

dtypes = {
  32 : np.uint32, 64 : np.uint64
}



def one_test_sort(sorter, a, verbose = 0):
    if verbose: 
        print("a=")
        print(a)
    sorted_a = np.copy(a)
    sorter(sorted_a, len(a))
    if verbose: 
        print("sorted")
        print(sorted_a)
    sorted_ref = np.sort(np.copy(a))
    if verbose: 
        #print("sorted ref")
        #print(sorted_ref)
        pass
    assert (sorted_a == sorted_ref).all(), (sorted_a, sorted_ref)



def sort_testdata(dtype = 32):
    d_t = dtypes[dtype]
    for i in range(1,40):
        yield   np.array( 0xfffffffe * np.random.rand(i), dtype = d_t)
  
    testdata = [ (4, 10), (10, 100), (100,10), (256, 20), (1299, 774), 
            (10, 10**9),
          (1000000, 2**dtype-1), 
    ]
    for length, maxval in testdata:
        yield   np.array( maxval * np.random.rand(length), dtype = d_t )


HEAP_SORTERS =  {32 : bitvector32_heapsort, 64 : bitvector64_heapsort}
SORTERS =  {32 : bitvector32_sort, 64 : bitvector64_sort}

def do_test_heapsort(dtype = 32, verbose = 0):
    sorter = HEAP_SORTERS[dtype]
    print("testing %d-bit heapsort...       " % dtype, end = "")
    for a  in sort_testdata(dtype): 
         one_test_sort(sorter, a, verbose = verbose)
    print("passed")

def do_test_sort(dtype = 32, verbose = 0):
    sorter = SORTERS[dtype]
    print("testing %d-bit bitvector sort... " % dtype, end = "")
    for a  in sort_testdata(dtype): 
         one_test_sort(sorter, a, verbose = verbose)
    print("passed")



@pytest.mark.qstate
def test_sorters(verbose = 0):
    for dtype in [32, 64]:
        do_test_heapsort(dtype, verbose)
        do_test_sort(dtype, verbose)



def do_test_benchmark_sort(n, nsamples = 100, n_repeat = 100, alg = 1):
    from mmgroup.clifford12 import bitmatrix32_test_sort
    if alg in (1,2):
        print("Runtime %ssort in micoroseconds, %d 32-bit integers, %s samples"
        % ("heap" if alg == 2 else "", n, nsamples))
    a = np.random.randint(0, 0xffffffff, n * nsamples, dtype=np.uint32)
    a = a.reshape((nsamples, n))
    #print(a)
    timings = []  #time of experiment in microseconds
    t_empty = bitmatrix32_test_sort(0, a[0], n_repeat)
    for i in range(nsamples):
        t1 = bitmatrix32_test_sort(alg, a[i], n_repeat)
        timings.append(1.0e6 * (t1 - t_empty) /  n_repeat)
    #print(timings)
    t_min, t_max = min(timings), max(timings)
    t_ave = sum(timings) / nsamples
    t_sigma = (sum((t - t_ave)**2 for t in timings) / (nsamples - 1)) ** 0.5
    print("ave: %.3f +- %.3f, min: %.3f, max: %.3f" % (
        t_ave, t_sigma, t_min, t_max))


@pytest.mark.bench
@pytest.mark.qstate
@pytest.mark.slow
def test_benchmark_sort(verbose = 0):
    print("")
    do_test_benchmark_sort(100, 100, 100, 1)
    do_test_benchmark_sort(100, 100, 100, 2)
    do_test_benchmark_sort(1000, 100, 100, 1)
    do_test_benchmark_sort(1000, 100, 100, 2)
    do_test_benchmark_sort(10000, 100, 10, 1)
    do_test_benchmark_sort(10000, 100, 10, 2)
    do_test_benchmark_sort(100000, 100, 1, 1)
    do_test_benchmark_sort(100000, 100, 1, 2)
    do_test_benchmark_sort(1000000, 10, 1, 1)
    do_test_benchmark_sort(1000000, 10, 1, 2)




