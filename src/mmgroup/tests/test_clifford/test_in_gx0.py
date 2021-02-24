import numpy as np
import pytest

from mmgroup import MM
import mmgroup.mm_order
from mmgroup.mm_order import check_mm_in_g_x0

def in_gx0_testdata():
    testdata = [
      ([], True),
      ([('d', 2)], True),
      ([('d', 0x801)], True),
      ([('x', 0xA71)], True),
      ([('x', 0x1A71)], True),
      ([('y', 0xCD1)], True),
      ([('y', 0xCD1), ('x', 0x1B7f)], True),
      ([('p', 'r')], True),
      ([('l', '1')], True),
      ([(tag, 'n') for tag in "xydpl"*3], True),
      ([('t', 1)], False),
    ]
    for g, in_g_x0 in testdata:
       yield MM(*g), in_g_x0
    for n in range(5):
        for i in range(2,5):
            yield MM(*[(tag, 'n') for tag in "xydpl"*i]), True


@pytest.mark.mm
def test_in_gx0(verbose = 0):
    print("Testing function check_mm_in_g_x0()")
    for n, (g, in_g_x0) in enumerate(in_gx0_testdata()):
        g = MM(g)
        g1 = check_mm_in_g_x0(g.copy())
        if verbose:
            print("Test", n+1)
            print("g =", g)
            print("reduced", g1)
            if g1 is None:
                r = mmgroup.mm_order.err_in_g_x0
                print("Reason for non-membership:", r)
        if in_g_x0:
            assert g1 is not None
        elif in_g_x0 == False:
            assert g1 is None
        if g1:
            assert g == g1 
     
