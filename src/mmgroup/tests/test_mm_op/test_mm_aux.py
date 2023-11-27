from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import pytest

from mmgroup.dev.mm_basics.mm_basics import MM_Basics
from mmgroup import mm_op
from mmgroup.mm_space import characteristics
from mmgroup.mm_space import MMSpace, MMVector
from random import randint



@pytest.mark.mm_op
def test_mm_aux_mmv_size():
    for i in range(2,9):
        p = 2**i - 1
        res = mm_op.mm_aux_mmv_size(p)
        ref = MM_Basics.sizes(p)["MMV_INTS"]
        assert res == ref, (res, ref, p)
    print("Test of function mm_op.mm_aux_mmv_size() passed")





def hash_testdata(ntests=3):
    for p in characteristics():
        MMV_INTS = MM_Basics.sizes(p)["MMV_INTS"]
        INT_FIELDS = MM_Basics.sizes(p)["INT_FIELDS"]
        FIELD_BITS = MM_Basics.sizes(p)["FIELD_BITS"]
        for i in range(ntests):
            v1 = MMVector(p, 'R')
            v2 = v1.copy()
            d1 = v1.data 
            d2 = v2.data 
            masks = []
            for j in range(64):
                field = randint(0, INT_FIELDS - 1)
                mask = p << (field * FIELD_BITS)
                masks.append(mask)

            for j in range(0, MMV_INTS, 3):
                mask = masks[i & 15]
                d1[j] = int(d1[j]) | mask
                d2[j] = int(d2[j]) & ~mask
            yield p, v1, v2

            v1 = MMVector(p)
            v2 = MMVector(p)
            d1 = v1.data 
            d2 = v2.data 
            for j in range(0, MMV_INTS, 33):
                mask = masks[j & 15]
                d1[j] = int(d1[j]) | mask
            if i == 0:
                yield p, v1, v2
            for j in range(100):
                x = randint(1, p-1)
                index = randint(0, MMV_INTS)
                mm_op.mm_aux_put_mmv(p, x, d1, index)
                mm_op.mm_aux_put_mmv(p, x, d2, index)
            yield p, v1, v2


@pytest.mark.mm_op
def test_mm_aux_hash(ntests = 5, verbose = 0):
    hashes = []
    print("Test hash function on vectors")
    for p, v1, v2 in hash_testdata(ntests):
        #assert sum(v1.data) == sum(v1.data)
        assert v1 == v2
        hashes.append(v1.hash())
        assert hashes[-1] == v2.hash()
    if verbose:
        s = "%d tests passed, %d different hash values"
        print(s % (len(hashes),len(set(hashes))))
    



        