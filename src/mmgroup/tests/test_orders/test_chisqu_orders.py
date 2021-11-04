
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from multiprocessing import cpu_count

from mmgroup.tests.test_orders.check_monster_orders import check_chisqu_orders

import pytest


NTESTS = 600
NPROCESSES =  max(1, cpu_count() - 1)


#@pytest.mark.very_slow
@pytest.mark.slow
@pytest.mark.orders
@pytest.mark.user
def test_chisqu_orders(ntests = NTESTS, nprocesses = NPROCESSES):
    ok = check_chisqu_orders(ntests, nprocesses)
    if not ok:
        err = "Chisquare test of Monster group order failed"
        raise ValueError(err)


if __name__ == "__main__":
    check_chisqu_orders(NTESTS, NPROCESSES)
