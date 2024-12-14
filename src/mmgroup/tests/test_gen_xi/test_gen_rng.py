from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from mmgroup.bitfunctions import bw24, lmap

import random
import numpy as np
from numpy import array, zeros, uint8, uint16, uint32

import pytest



from mmgroup.generators import gen_rng_int_dist
from mmgroup.generators import rand_get_seed
from mmgroup.tests.chisquare import chisquare

        

#######################################################################
# Test functions
#######################################################################

seed = rand_get_seed()


def sample_int_dist(length, prob):
    a = np.zeros(length, dtype = np.uint32)
    s = 1.0 / sum(prob)
    p = np.array(prob,  dtype = float) * s
    dist = np.copy(p)
    for i in range(1, len(dist)):
        dist[i] += dist[i-1]
    assert abs(1 - dist[-1]) < 1.0e-14
    gen_rng_int_dist(seed, len(dist) - 1, dist, a, length)
    return p, a

def test_rng_int_dist():
    LENGTH = 20000
    prob = [100, 1, 57, 39, 11, 323, 3]
    for i in range(5):
        expected, a = sample_int_dist(LENGTH, prob)
        found = np.bincount(a, minlength = len(prob))
        chisq, p = chisquare(found, expected)
        #print(p)
        if p < 0.95:
            return
    ERR = "Test of function gen_rng_int_dist has failed"
    raise ValueError(ERR)
    


