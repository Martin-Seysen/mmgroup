
import time
from random import randint, sample
import numpy as np
import pytest

from mmgroup.clifford12 import heapsort32, heapsort64
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


def do_test_heapsort(dtype = 32, verbose = 0):
    sorter = {32 : heapsort32, 64 : heapsort64}[dtype]
    print("testing %d-bit heapsort...       " % dtype, end = "")
    for a  in sort_testdata(dtype): 
         one_test_sort(sorter, a, verbose = verbose)
    print("passed")

def do_test_sort(dtype = 32, verbose = 0):
    sorter = {32 : bitvector32_sort, 64 : bitvector64_sort}[dtype]
    print("testing %d-bit bitvector sort... " % dtype, end = "")
    for a  in sort_testdata(dtype): 
         one_test_sort(sorter, a, verbose = verbose)
    print("passed")



@pytest.mark.qstate
def test_sorters(verbose = 0):
    for dtype in [32, 64]:
        do_test_heapsort(dtype, verbose)
        do_test_sort(dtype, verbose)


