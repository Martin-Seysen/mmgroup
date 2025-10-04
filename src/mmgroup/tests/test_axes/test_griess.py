
from random import randint, shuffle
import numpy as np

import datetime
import time
import pytest

from mmgroup import mat24
from mmgroup import MM, MMSpace, MMV, Cocode, AutPL, XLeech2




def import_all():
    global Axis, product_axis_2A
    from mmgroup.axes import Axis, product_axis_2A



def product_axis_2A_test_multipliers():
    #yield MM()
    for i in range(50):
        yield MM('r')


@pytest.mark.mmm
@pytest.mark.axes
def test_product_axis_2A(verbose = 0):
    print("\nStart")
    import_all()
    INV1 = MM('d', [2,3])
    INV2 = MM('d', [1,2])
    INV3 = INV1 * INV2
    AX1 = Axis('i', INV1)
    AX2 = Axis('i', INV2)
    print("Start loop")
    for i, g in enumerate(product_axis_2A_test_multipliers()):
        ax1 = AX1 * g
        ax2 = AX2 * g
        ax3 = product_axis_2A(ax1, ax2)
        inv3 = ax3.g_axis
        assert inv3 == INV3 ** g
    print("End loop")
