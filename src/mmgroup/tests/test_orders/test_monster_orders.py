
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os

import pytest

if __name__ == "__main__":    
    #print(__file__)
    FILE_DIR = os.path.dirname(os.path.realpath(__file__))
    #print("FILE_DIR ==", FILE_DIR)
    TEST_DIR = os.path.dirname(FILE_DIR)
    #print("TEST_DIR ==", TEST_DIR)
    PACKAGE_DIR = os.path.dirname(TEST_DIR)
    #print("PACKAGE_DIR ==", PACKAGE_DIR)
    SRC_DIR = os.path.dirname(PACKAGE_DIR)
    print("SRC_DIR ==", SRC_DIR)
    sys.path.append(SRC_DIR)

from mmgroup.tests.test_orders.check_monster_orders import check_mm_orders

NTESTS = 100
BUILD_NTESTS = 20

@pytest.mark.orders 
def test_mm_orders():
    check_mm_orders(ntests = NTESTS, display = False)


@pytest.mark.build
@pytest.mark.slow     # We don't want this in the usual test
def test_mm_orders():
    check_mm_orders(ntests = BUILD_NTESTS, display = False)


if __name__ == "__main__":
    check_mm_orders(ntests = NTESTS, display = False)
    




