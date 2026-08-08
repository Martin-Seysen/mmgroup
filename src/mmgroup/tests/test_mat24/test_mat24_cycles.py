from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint

from mmgroup import mat24


from mmgroup.bitfunctions import bit_mat_mul, v2
from mmgroup.dev.mat24.mat24theta import theta_to_basis_vector





#####################################################################
# Test operation of the group Mat24 on Golay code and cocode
#####################################################################




def cycle_structure(l):
    """Return the cycle structure of the permutation i -> l[i]."""
    n = len(l)
    visited = [False] * n
    result = {}

    for i in range(n):
        if not visited[i]:
            length = 0
            j = i
            while not visited[j]:
                visited[j] = True
                j = l[j]
                length += 1

            result[length] = result.get(length, 0) + 1

    return result


def Mat24_cycle_testcases():
    """Yields permutations in Mat24 for test function test_mat24_cycles.

    Permutation are represented as list of length 24.
    """
    testdata = [0, 244823040-1]
    for i in range(20):
        testdata.append(randint(1, 44823040))
    for t in testdata:
        yield mat24.m24num_to_perm(t)


@pytest.mark.mat24
def test_mat24_cycles():
    """test computation of cycle structure in the Mathieu group Mat24"""
    from math import lcm
    for p in Mat24_cycle_testcases():
        order, c = mat24.perm_cycles(p)
        d = {}
        for i in range(0, len(c), 2):
            d[c[i]] = c[i+1]
        d_ref = cycle_structure(p)
        if d != d_ref:
            print("Permutation", p)
            print("Cycles expected:", d_ref)
            print("Cycles obtained:", d)
            raise(ValueError, "Computation of cycle structure failed")
        order_ref = lcm(*d_ref.keys())
        if order != order_ref:
            print("Permutation", p)
            print("Order expected: %d, obtained: %d" % (order_ref, order))
            raise(ValueError, "Computation of corder failed")
