import math
from random import randint, shuffle, sample
import numpy as np


import datetime
import time
import pytest

from mmgroup import mat24
from mmgroup import mmv_scalprod
from mmgroup import MM, MM0, MMVector, MMV, Cocode, AutPL, XLeech2

from mmgroup.tests.axes.griess import product_axis_2A # Preliminary!!


def import_all():
    global Axis, product_axis_2A, Griess, Griess_scalar
    from mmgroup.axes import Axis, product_axis_2A
    from mmgroup.axes import Griess, Griess_scalar



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
    last = None
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
        if last:
            b0, b1 = last
            prod_l = Griess((ax1, ax2), (b0, b1), p=15, n=32)
            prod_r = Griess((b0, b1), (ax1, ax2), p=15, n=32)
            assert prod_l == prod_r
        last = ax1, ax2
    #print("End loop")




@pytest.mark.slow
@pytest.mark.bench
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




def trilinear_testcases():
    def MV(*args):
        return MMVector(*args)
    PAR = {'p':3, 'n':128}
    def par(p, n, n_axes=None):
       return {'p':p, 'n':n, 'n_axes':  n_axes}

    DATA = [
       # A entry contains: arg1, arg2, arg3, (p,n,n_axes)
       (Axis('r'), MV(15,'R'), MV(15,'R'), (15,32,1)),
       (Axis('r'), MV(127,'R'), MV(127,'R'), (127,128,1)),
       (Axis('r'), MV(3,'R'), MV(3,'R'), (3,2,1)),
       (Axis('r'), MV(7,'R'), MV(7,'R'), (7,8,1)),
       (Axis('r'), MV(31,'R'), MV(31,'R'), (31,512,1)),
       (Axis('r'), MV(255,'R'), MV(255,'R'), (255,8,1)),
    ]
    for d in DATA:
        yield d[:-1] + (par(*d[-1]),)


def isqrt(x):
    w = math.sqrt(x)
    w0 = int(round(w))
    assert abs(w - w0) < 1.0e-6
    return w0

@pytest.mark.axes
def test_trilinear(verbose = 0):
    import_all()
    for i, (a, b, c, PAR) in enumerate(trilinear_testcases()):
        if verbose:
            print("\nTest", i+1, PAR)
        a1, b1, c1 = sample([a, b, c], 3)
        s1 = Griess_scalar(a1, b1, c1, **PAR)
        s2 = mmv_scalprod(c, Griess(a, b, **PAR))
        assert s1 == s2
        n_old = PAR['n']
        n_axes = PAR['n_axes']
        if n_axes is not None:
            factor = 2 ** randint(1,3)
            PAR['n'] = PAR['n'] * factor**2
            s12 = Griess_scalar(a1, b1, c1, **PAR)
            scale = factor ** n_axes
            p = PAR['p']
            assert s12 ==  (scale * s2) % p, (PAR, scale, n_old, n_axes)


@pytest.mark.axes
def test_axes_product_pairs():
    import_all()
    AX1 = Axis('i', MM('d', [2,3]))
    AX2 = Axis('i', MM('d', [1,2]))
    PAR_DAT = [
        (15, 8, 32), (127, 8, 8),  (255, 8, 128), (31, 32, 128)
    ]

    for  p, n1, n2 in PAR_DAT:
        assert n2 % n1 == 0
        r12, r34 = MM('r'), MM('r')
        ax1, ax2 = AX1 * r12, AX2 * r12
        ax3, ax4 = AX1 * r34, AX2 * r34
        v1 = Griess((ax1, ax2), (ax3, ax4), p=p, n=n1)
        v2 = Griess((ax4, ax3), (ax2, ax1), p=p, n=n2)
        ok = v2 == v1 * isqrt(n2/n1)**4
        assert ok


