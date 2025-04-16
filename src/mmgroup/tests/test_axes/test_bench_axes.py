import time

import numpy as np
import pytest


from mmgroup import Xsp2_Co1
from mmgroup.mm_reduce import mm_reduce_analyze_2A_axis


def bench_analyze_axis(axis_type, axis, ntests, chunksize = 1):
    n_axes = min(ntests, 32)
    out = np.zeros(1000, dtype = np.uint32)
    v = []
    for i in range(n_axes):
        v.append((axis.v15 *  Xsp2_Co1('r')).data)
    t = []
    n = (ntests // chunksize) + 2
    for n1 in range(n):
        t_start = time.perf_counter() 
        for n2 in range(chunksize):
            mm_reduce_analyze_2A_axis(v[n1 % n_axes], out) 
        t.append(time.perf_counter() - t_start)
    #print(t)
    m = sum(t)
    mu = m / n
    s = sum([x * x for x in t])
    sigma = (s - mu * m) / (n - 1)
    return mu / chunksize, sigma ** 0.5 / chunksize


@pytest.mark.axes
@pytest.mark.bench
@pytest.mark.slow
@pytest.mark.very_slow
def test_analyze_axes():
    M = "Axis type %3s: %7.3f +- %7.3f us"
    F = 1.0e6
    from mmgroup.tests.axes.axis import Axis
    ntests = 4000
    REP = Axis.representatives()
    print("\nRuntime for analyzing axis in mmgroup operation")
    for axis_type, axis in REP.items():
        mu, v = bench_analyze_axis(axis_type, axis, ntests)
        print(M % (axis_type, mu * F, v * F))
