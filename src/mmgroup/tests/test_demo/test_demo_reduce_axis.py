r"""Test demo program

"""

from random import randint
import pytest



if __name__ == "__main__":
    sys.path.append("../../../")



def import_all():
    global V_START, V_BABY_START
    global reduce_axis
    global Mm, MmV15
    global axis_orbit
    global reduce_axis
    global reduce_feasible_axis
    from mmgroup.demo.reduce_axis import reduce_axis
    from mmgroup.demo import Mm, MmV15
    from mmgroup.demo.reduce_axis import axis_orbit
    from mmgroup.demo.reduce_axis import reduce_axis
    from mmgroup.demo.reduce_feasible import reduce_feasible_axis
    V_START = MmV15("v+")  
    V_BABY_START = MmV15("v-")


def axis_testdata(NTESTS = 3):
    from mmgroup.tests.axes.sample_axes import g_classes, g_strings
    n_types = 0
    for axis_type, g_str in zip(g_classes, g_strings):
        # Construct an axis v of the given axis type
        v = V_START * Mm(g_str)
        for i in range(NTESTS):
            yield v * Mm('r', 'G_x0'), axis_type
        n_types += 1
    assert n_types == 12


@pytest.mark.demo
def test_axis_type():
    import_all()
    for v, ref_ax_type in axis_testdata():
        ax_type =  axis_orbit(v)
        assert ax_type == ref_ax_type
        g = reduce_axis(v)
        assert v * g == V_START

  
def feasible_axis_testdata(NTESTS = 3):
    from mmgroup.tests.axes.baby_sample_axes import g_classes, g_strings
    n_types = 0
    for axis_type, g_str in zip(g_classes, g_strings):
        # Construct an axis v of the given axis type
        v = V_BABY_START * Mm(g_str)
        for i in range(NTESTS):
            yield v * Mm('r', 'B'), axis_type
        n_types += 1
    assert n_types == 10


@pytest.mark.demo
def test_feasible_axis_type():
    import_all()
    for v, ref_ax_type in feasible_axis_testdata():
        g = reduce_feasible_axis(v)
        assert V_START * g == V_START
        assert v * g == V_BABY_START




