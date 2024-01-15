r"""Test demo program

"""

import pytest


if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup.demo import Mm, MmV15
from mmgroup.demo.reduce_monster import check_reduce_monster_element


def monster_element_testdata(NTESTS = 3):
    for i in range(10):
        g = Mm('r') * Mm('r')
        yield g


@pytest.mark.demo
def test_monster_reduction(NTESTS = 10):
    for i in range(NTESTS):
        check_reduce_monster_element()




       


