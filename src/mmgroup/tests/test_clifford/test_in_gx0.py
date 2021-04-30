import numpy as np
import pytest

from mmgroup import MM
import mmgroup.mm_order
from mmgroup.mm_order import check_mm_in_g_x0

def in_gx0_testdata():
    """Test data for testing function mgroup.mm_order.check_mm_in_g_x0

    The function yields pairs ``(g, b)``. Here ``g`` an element of the
    monster group ``MM``, and ``b`` is True iff ``g`` is in the
    subgroup ``G_x0`` of the monster group.
    """
    # Yield some fixed test data
    testdata = [
      ([], True),
      ([('d', 2)], True),
      ([('d', 0x801)], True),
      ([('x', 0xA71)], True),
      ([('x', 0x1A71)], True),
      ([('y', 0xCD1)], True),
      ([('y', 0xCD1), ('x', 0x1B7f)], True),
      ([('p', 'r')], True),
      ([('l', 1)], True),
      ([(tag, 'n') for tag in "xydpl"*3], True),
      ([('t', 1)], False),
      (["M<y_4d1h*x_0b7fh>"], True),
    ]
    for g, in_g_x0 in testdata:
       yield MM(*g), in_g_x0
    # Yield elements ``g0 * t * g1``, where ``t`` is the triality 
    # element or its inverse, and ``g0 , g1`` are in ``G_x0``.
    for i in range(8):
        e = [1 + ((i >> j) & 1) for j in range(3)]
        g = [('p', 'r'), ('y', 'r'), ('l', e[0]), ('t', e[1]), ('l', e[2]),
             ('p', 'r'), ('y', 'r'), ('x', 'r'), ('d', 'r'), ]
        yield  MM(*g), False
    # Yield some less obvous examples. Let ``d`` be a random 2A 
    # involution and ``z`` be the 2B involution in the center of
    # ``G_x0``. Then ``g = z * d`` has even order ``o``, and we 
    # yield ``g ** (o/2)``, which is in ``G_x0``.
    z, d0 = MM(('x', 0x1000)), MM(('d', [2,3]))
    for n in range(5):
        g0 = MM(*[(tag, 'n') for tag in "xydplt"*2])
        d =  d0 ** g0
        g = z * d
        o = g.order()
        assert o & 1 == 0
        yield g ** (o >> 1), True
    # Yield some elements of ``G_x0``.
    for n in range(5):
        for i in range(2,5):
            yield MM(*[(tag, 'n') for tag in "xydpl"*i]), True
    



@pytest.mark.mm
def test_in_gx0(verbose = 0):
    print("Testing function check_mm_in_g_x0()")
    for n, (g, in_g_x0) in enumerate(in_gx0_testdata()):
        g = MM(g)
        if verbose:
            print("Test", n+1)
            print("g =", g)
        g1 = check_mm_in_g_x0(g.copy())
        if verbose:
            print("reduced", g1)
            if g1 is None:
                r = mmgroup.mm_order.err_in_g_x0
                print("Reason for non-membership:", r)
        if in_g_x0:
            assert g1 is not None
        else:
            assert g1 is None
        if g1:
            assert g == g1 
            print("in G_x0", g1)
            
     
