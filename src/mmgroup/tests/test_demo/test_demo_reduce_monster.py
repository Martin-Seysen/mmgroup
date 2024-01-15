r"""Test demo program

"""

import pytest


if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup.demo import Mm, MmV15
from mmgroup.demo.reduce_axis import reduce_monster_element

V1 = MmV15("v1")


def monster_element_testdata(NTESTS = 3):
    for i in range(10):
        g = Mm('r') * Mm('r')
        yield g


@pytest.mark.demo
def test_monster_reduction():
    for g in monster_element_testdata():
        h = reduce_monster_element(g)
        assert h.count_triality_elements <= 7




       


