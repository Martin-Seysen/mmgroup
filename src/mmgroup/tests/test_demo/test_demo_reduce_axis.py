r"""Test demo program

"""

from random import randint
import pytest

from mmgroup.demo.reduce_axis import reduce_axis


if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup.demo import Leech2, Mm, MmV15
from mmgroup.demo.reduce_axis import axis_orbit
from mmgroup.demo.reduce_axis import reduce_axis
from mmgroup.demo.reduce_feasible import reduce_feasible_axis

from mmgroup.tests.test_axes.test_import import AXES, BABY_AXES

V_START = MmV15("v+")
V_BABY_START = MmV15("v-")


def axis_testdata(NTESTS = 3):
    n_types = 0
    for axis_type, g_str in AXES.items():
        # Construct an axis v of the given axis type
        v = V_START * Mm(g_str)
        for i in range(NTESTS):
            yield v * Mm('r', 'G_x0'), axis_type
        n_types += 1
    assert n_types == 12


@pytest.mark.demo
def test_axis_type():
    for v, ref_ax_type in axis_testdata():
        ax_type =  axis_orbit(v)
        assert ax_type == ref_ax_type
        g = reduce_axis(v)
        assert v * g == V_START

  
def feasible_axis_testdata(NTESTS = 3):
    n_types = 0
    for axis_type, g_str in BABY_AXES.items():
        # Construct an axis v of the given axis type
        v = V_BABY_START * Mm(g_str)
        for i in range(NTESTS):
            yield v * Mm('r', 'H+'), axis_type
        n_types += 1
    assert n_types == 10


@pytest.mark.demo
def test_feasible_axis_type():
    for v, ref_ax_type in feasible_axis_testdata():
        g = reduce_feasible_axis(v)
        assert V_START * g == V_START
        assert v * g == V_BABY_START




