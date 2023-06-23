from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, shuffle, sample
from functools import reduce
from operator import __or__
from multiprocessing import Pool
import time

import numpy as np
import pytest

from mmgroup import Xsp2_Co1, PLoop, AutPL, Cocode, MM0, MM, XLeech2
from mmgroup import MMV
from mmgroup.generators import gen_leech2_op_word
from mmgroup.clifford12 import xsp2co1_set_elem_word_scan
from mmgroup.clifford12 import xsp2co1_unit_elem
from mmgroup.clifford12 import xsp2co1_reduce_word
from mmgroup.clifford12 import xsp2co1_elem_read_mod3
from mmgroup.structures.construct_mm import iter_mm


#####################################################################
# Create test matrices
#####################################################################



def construct_xsp2co1_testcases():
    def a(string_data):
        data = [(x, 'n') for x in string_data]
        return np.array(list(iter_mm(None, data)), dtype = np.uint32)

    for tag0 in "pdxyl":
        yield a(tag0), 1
        for tag1 in "pdxyl":
            yield a(tag0 + tag1 + 't'), 2

    for i in range(5):
        data = Xsp2_Co1('r', 'G_x0').mmdata
        yield data, len(data)

def dirty_elem():
    return np.array([randint(0, 0xffffffffffffffff) for i in range(26)],
        dtype = np.uint64)

@pytest.mark.xsp2co1
def test_xsp2co1_set_elem_word_scan():
    for a, len_a in construct_xsp2co1_testcases():
        elem0 = dirty_elem()
        elem1 = dirty_elem()
        s0 =  xsp2co1_set_elem_word_scan(elem0, a, len(a), 0)
        assert s0 == len_a
        xsp2co1_unit_elem(elem1)
        s1 =  xsp2co1_set_elem_word_scan(elem1, a, len(a), 1)
        assert s1 == len_a
        assert s0 == s1
        xsp = Xsp2_Co1()
        np.copyto(xsp._data, elem0)
        assert Xsp2_Co1('a', a[:len_a]) == xsp



#####################################################################
# Test method as_Co1_bitmatrix of class Xsp2_Co1
#####################################################################

@pytest.mark.xsp2co1
def test_xsp2co1_as_Co1_bitmatrix():
    or_sum = 0
    for i in range(3):
         v  = XLeech2(randint(0, 0x1fffff))
         v_ord = v.ord
         or_sum |= v_ord
         g =  Xsp2_Co1('r', 'G_x0')
         v_bits = v.as_Leech2_bitvector()
         ref_bits = [(v_ord >> i) & 1 for i in range(24)]
         assert list(v_bits) == ref_bits
         w = v * g
         #print("v", v,  v_bits)
         #print("g", g.as_Co1_bitmatrix())
         w_bits = list((v_bits @ MM0(g).as_Co1_bitmatrix()) & 1)
         w_ord = w.ord
         w_ref_bits = [(w_ord >> i) & 1 for i in range(24)] 
         assert w_bits == w_ref_bits
    assert or_sum != 0


#####################################################################
# Benchmarks
#####################################################################


def benchmark_mul_xsp2co1(ncases = 20):
    glist = []
    for i in range(64):
        glist.append(Xsp2_Co1('r', 'G_x0'))
    indices = range(len(glist))
    index_pairs = [sample(indices, 2) for i in range(ncases)]
    t_start = time.process_time()
    for i, j in index_pairs:
        glist[i] *= glist[j]
    t = time.process_time() - t_start
    return ncases, t



@pytest.mark.bench 
@pytest.mark.xsp2co1
def test_benchmark_mul(ncases = 20000, verbose = 0):
    s = "Runtime of multiplication in class Xsp2_Co1, %d tests: %.5f ms" 
    print("")
    for i in range(1):
        n, t = benchmark_mul_xsp2co1(ncases) 
        print(s % (n, 1000*t/n))



def benchmark_mul_rep_4096(ncases = 20):
    glist = []
    for i in range(64):
        glist.append(Xsp2_Co1('r', 'G_x0').qs)
    indices = range(len(glist))
    index_pairs = [sample(indices, 2) for i in range(ncases)]
    t_start = time.process_time()
    for i, j in index_pairs:
        glist[i] @= glist[j]
    t = time.process_time() - t_start
    return ncases, t

@pytest.mark.bench 
@pytest.mark.xsp2co1
def test_benchmark_mul_rep_4096(ncases = 20000, verbose = 0):
    s = "Runtime of multiplication in rep 4096x, %d tests: %.5f ms" 
    print("")
    for i in range(1):
        n, t = benchmark_mul_rep_4096(ncases) 
        print(s % (n, 1000*t/n))







def benchmark_reduce_word(ncases = 20):
    glist = []
    for i in range(64):
        glist.append(MM0([(tag, 'n') for tag in "xydplplplp"]).mmdata)
    buffer = np.zeros(10, dtype = np.uint32)
    t_start = time.process_time()
    for i in range(ncases):
        xsp2co1_reduce_word(glist[i & 63], len(glist[i & 63]), buffer)
    t = time.process_time() - t_start
    return ncases, t


@pytest.mark.bench 
@pytest.mark.xsp2co1
def test_benchmark_reduce_word(ncases = 5000, verbose = 0):
    s = "Runtime of word reduction in G_x0, %d tests: %.5f ms" 
    print("")
    for i in range(1):
        n, t = benchmark_reduce_word(ncases) 
        print(s % (n, 1000*t/n))





def benchmark_mod2_op_xsp2co1(ncases = 20):
    glist = []
    for i in range(64):
        glist.append(Xsp2_Co1('r', 'G_x0'))
    t_start = time.process_time()
    for i in range(ncases):
        a = glist[i & 63].leech_mod2_op
    t = time.process_time() - t_start
    return ncases, t


@pytest.mark.bench 
@pytest.mark.xsp2co1
def test_benchmark_reduce_to_Co1(ncases = 5000, verbose = 0):
    s = "Runtime of reduction from G_x0 to Co_1, %d tests: %.5f ms" 
    print("")
    for i in range(1):
        n, t = benchmark_mod2_op_xsp2co1(ncases) 
        print(s % (n, 1000*t/n))



#####################################################################
# Test function xsp2co1_elem_read_mod3
#####################################################################

V3 = MMV(3)

def make_testcases_mod3():
    for i in range(50):
        v = V3('R')
        g = Xsp2_Co1('r', 'G_x0')
        row = randint(0, 0xfff)
        col = randint(0, 23) if i & 1 else 24
        yield v, g, row, col



def ref_xsp2co1_elem_read_mod3(v, g, row, col):
    v1 = v.copy() * g**-1
    r = v1['Z', row] if row < 2048 else  v1['Y', row - 2048]  
    if col < 24:
        return r[col]
    else:
        return (r[2] + 3 - r[3]) % 3

OFS_Z = 116416 >> 5

@pytest.mark.xsp2co1
def test_xsp2co1_elem_read_mod3(verbose = 0):
    for n, (v, g, row, col) in enumerate(make_testcases_mod3()):
        x_ref = ref_xsp2co1_elem_read_mod3(v, g, row, col)
        x = xsp2co1_elem_read_mod3(v.data[OFS_Z:], g._data, row, col)
        ok = x == x_ref
        if verbose or not ok:
            print("Test %d, row = %d, column =%d" % (n+1, row,col))
            print("  inverse of g =", g)
            vg = v * g
            r = vg['Z', row] if row < 2048 else  vg['Y', row - 2048]  
            print("  (v*g)[row] =", list(r))
            print("  (v*g)[row, col] =", x, ", expected:", x_ref)
            if not ok:
                raise ValueError("Test failed")
        

