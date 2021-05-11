"""Configuration for development"""

import sys
import os



# directory names
DEV_DIR = os.path.dirname(os.path.realpath(__file__))
PACKAGE_DIR = os.path.dirname(DEV_DIR)
SRC_DIR = os.path.dirname(PACKAGE_DIR)
ROOT_DIR = os.path.dirname(SRC_DIR)

sys.path.insert(0, ROOT_DIR)
from config import INT_BITS, PRIMES, C_DIR, PXD_DIR
assert sys.path.pop(0) == ROOT_DIR



