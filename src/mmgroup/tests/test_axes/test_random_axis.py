"""Generate random 2A involutions in the Monster and their axes.

The main function random_axis() in this module generates a
random 2A involution in the Monster and its axis. It is
optimized for utmost speed. 

The idea is that a large number of random 2A involutions
is  generated; and that certain subsets of the Monster are
watermarked by scalar products of axes of 2A involutions
related to these subsets.  
"""

import sys
if __name__ == "__main__":
    sys.path.append("../../../")
import time

import numpy as np
import pytest

from mmgroup import MM, MM0, Xsp2_Co1
from mmgroup.tests.axes.axis import Axis, BabyAxis
from mmgroup.tests.axes.random_axis import RandomAxis, VV
from mmgroup.tests.axes.random_axis import RandomBabyAxis
from mmgroup.tests.axes.random_axis import rand_mm_element
from mmgroup.tests.axes.random_axis import rand_G_x0
from mmgroup.tests.axes.random_axis import rand_G_x0_baby
from mmgroup.tests.axes.random_axis import BETA



#############################################################
# Functional tests
#############################################################

def do_test_rand_G_x0():
    for i in range(10):
        a, n = rand_G_x0()
        a_comp, n_comp = rand_G_x0(num = n)
        assert n == n_comp, (n, n_comp)
        assert (a == a_comp).all(), (a, a_comp)
        _ = Xsp2_Co1('a', a) 
    t_start = time.process_time()
    N = 20000
    for i in range(N):
        rand_G_x0()
    t = (time.process_time() - t_start) / N
    print("Runtime of rand_G_x0: %.5f ms" % (1000 * t))


def do_test_rand_G_x0_baby():
    beta = Xsp2_Co1('d', BETA)
    import time
    for i in range(10):
        a, n = rand_G_x0_baby()
        a_comp, n_comp = rand_G_x0_baby(num = n)
        assert n == n_comp, (n, n_comp)
        assert (a == a_comp).all(), (a, a_comp)
        assert beta ** Xsp2_Co1('a', a) == beta
    t_start = time.process_time()
    N = 20000
    for i in range(N):
        rand_G_x0_baby()
    t = (time.process_time() - t_start) / N
    print("Runtime of rand_G_x0_baby: %.5f ms" % (1000 * t))




def do_test_rand_axis(random_class, axis_class):
    print("Test of %s" % random_class, end = "")
    for i in range(0):
        # Generate random axis and corresponding involution
        ax =  random_class()
        # Check recomputation of axis from its number 
        ax1 = random_class(ax.number)
        ax_v15 = ax.v15()
        assert ax_v15 == ax1.v15()
        # check (standard) axis derives from random axis
        axis = ax.axis()
        assert ax_v15 == axis.v15
        if i >= 2:
            continue
        involution = ax.g_axis()
        # check involution corresponding to axis
        axis = Axis('i', involution)
        assert ax_v15  == axis.in_space(VV)
    print(" passed")


@pytest.mark.axes
def test_random_axes():
    print()
    do_test_rand_G_x0()
    do_test_rand_G_x0_baby()
    do_test_rand_axis(RandomAxis, Axis)
    do_test_rand_axis(RandomBabyAxis, BabyAxis)


#############################################################
# Benchmarks
#############################################################

def bench_rand_axis():
    import time
    v = RandomAxis()
    N_BENCH_RAND_AXIS = 1000
    t_start = time.process_time()
    for i in range(N_BENCH_RAND_AXIS):
        axis = RandomAxis()
        v = axis.v15()
        g_axis = axis.g_axis(group = 'a')
        n = axis.num
    t = time.process_time() - t_start
    return t / N_BENCH_RAND_AXIS

def bench_mul_axis():
    import time
    N_BENCH_MUL_AXIS = 500
    axis = RandomAxis().v15()
    g = [MM0(MM('r')) for i in range(16)]   
    t_start = time.process_time()
    for i in range(N_BENCH_MUL_AXIS):
        axis *= g[i & 15]
    t = time.process_time() - t_start
    return t / N_BENCH_MUL_AXIS


def bench_scalprod():
    import time
    from mmgroup import mmv_scalprod
    N_BENCH_SCALPROD = 5000
    a = [RandomAxis().v15() for i in range(16)]
    x = 0
    t_start = time.process_time()
    for i in range(N_BENCH_SCALPROD):
        x += mmv_scalprod(a[i & 15], a[(i + 1) & 15])
    t = time.process_time() - t_start
    return t / N_BENCH_SCALPROD

@pytest.mark.bench 
@pytest.mark.slow 
@pytest.mark.axes
def test_bench_all():
    print("\nBenchmarking operations with random axes...")
    t = bench_rand_axis()
    print("Runtime for generating a random axis: %.3f ms" %
         (1000*t) ) 
    t = bench_mul_axis()
    print("Runtime for multiplying axis with group element: %.3f ms" %
         (1000*t) )
    t = bench_scalprod()
    print("Runtime for scalar product of axes: %.3f ms" %
         (1000*t) )


