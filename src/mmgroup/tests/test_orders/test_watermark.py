import sys
import os
import warnings

from random import randint
from collections import OrderedDict

import numpy as np
import pytest

from mmgroup import MM0, MMV
from mmgroup.mat24 import MAT24_ORDER
from mmgroup.mm_op import mm_op_word_tag_A
from mmgroup.mm_op import mm_op_watermark_A
from mmgroup.mm_op import mm_op_watermark_A_perm_num




#######################################################################
# Test leech3matrix_watermark and leech3matrix_watermark_perm_num
#######################################################################
  

MMV15 = MMV(15)
MM = MM0    
         
def one_test_watermark(verbose = 0):
    v0 = MMV15('R')
    v = v0.data[:48]
    w0 = np.zeros(24, dtype = np.uint32)
    result = mm_op_watermark_A(15, v, w0)
    if result < 0:
        return 0
    pi_num = randint(0, MAT24_ORDER-1)
    TAG_y, TAG_p = 0x40000000, 0x20000000
    y1, y2 = randint(0, 0xfff), randint(0, 0xfff)
    op = np.array([TAG_y + y1, TAG_p + pi_num, TAG_y + y2], 
         dtype = np.uint32)
    mm_op_word_tag_A(15, v, op, len(op), 1)
    pi_num_obt = mm_op_watermark_A_perm_num(15, w0, v)   
    assert pi_num_obt == pi_num 
    return 1


WATERMARK_TESTS = 10
WATERMARK_MIN_SUCCESS = 4

@pytest.mark.qstate
def test_watermark():
    success = 0
    for n in range(WATERMARK_TESTS):
        success += bool(one_test_watermark())
    if success < WATERMARK_MIN_SUCCESS:
        err = "%d of %s permutation watermark tests failed"
        raise ValueError(err, WATERMARK_TESTS-success, WATERMARK_TESTS)
