from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix
from mmgroup.structures.xs1_co1 import Xs12_Co1
from mmgroup.tests.spaces.clifford_space import Xs12_Co1



#####################################################################
# Convert matrix of type QStateMatrix to a complex matrix
#####################################################################

def create_test_vectors():
   vector_data = [ 
      (0x001, ("B", 2, 3))
   ]
   for x4096, x24 in vector_data:
        #print (x4096, x24)
        yield Xs12_Co1.unit(x4096, x24)



@pytest.mark.qstate1
def test_vector(verbose = 1):
    FORMAT_REDUCED = False
    for ntest, v in enumerate(create_test_vectors()):
        if verbose:
            print("TEST %s" % (ntest+1))
            print(v.as_tuples())
            print(v)
            v.dump()
