

import sys
import os
from multiprocessing import cpu_count

sys.path.append("src")
from mmgroup.tests.test_mm.check_monster_orders import check_chisqu_orders

NTESTS = 1000
NPROCESSES =  max(1, cpu_count() - 1)

if __name__ == "__main__":
    check_chisqu_orders(NTESTS, NPROCESSES)
