
from random import randint, shuffle
import numpy as np

import datetime
import time
import pytest

from mmgroup import mat24
from mmgroup import MM, MM0, MMSpace, MMV, Cocode, AutPL, XLeech2

from mmgroup.tests.axes.griess import product_axis_2A # Preliminary!!


def import_all():
    global Axis, product_axis_2A, Griess
    from mmgroup.axes import Axis, product_axis_2A, Griess



def product_axis_2A_test_multipliers():
    #yield MM()
    for i in range(3):
        yield MM('r')


@pytest.mark.axes
def test_product_axis_2A(verbose = 0):
    #print("\nStart")
    import_all()
    INV1 = MM('d', [2,3])
    INV2 = MM('d', [1,2])
    INV3 = INV1 * INV2
    AX1 = Axis('i', INV1)
    AX2 = Axis('i', INV2)
    #print("Start loop")
    for i, g in enumerate(product_axis_2A_test_multipliers()):
        if verbose:
            print("Test", i)
        ax1 = AX1 * g
        ax2 = AX2 * g
        ax3 = product_axis_2A(ax1, ax2)
        inv3 = ax3.g_axis
        assert inv3 == INV3 ** g
        v3a = Griess(ax1, ax2, p=15, n=32)
        v3a.check()
        ref_v3a = 4 * (ax1.v15 + ax2.v15 - ax3.v15)
        ok = v3a == ref_v3a
    #print("End loop")




@pytest.mark.bench
@pytest.mark.slow
@pytest.mark.axes
def test_bench_product_axis_2A(ntests=5):
    import_all()
    INV1 = MM('d', [2,3])
    INV2 = MM('d', [1,2])
    AX1 = Axis('i', INV1)
    AX2 = Axis('i', INV2)
    cases = []
    for i in range(ntests):
        g = MM0('r')
        cases.append((AX1 * g, AX2 * g))
    t_start = time.perf_counter()
    for ax1, ax2 in cases:
        product_axis_2A(ax1, ax2)
    t = time.perf_counter() - t_start
    t1 = 1000 * t / ntests
    print("\nRuntime of function test_product_axis_2A: %.1f ms" % t1)


