
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os

import numpy as np

import pytest

from mmgroup import MM0
from mmgroup.mm_crt_space import MMVectorCRT 


from mmgroup.tests.test_axes.get_sample_axes import import_sample_axes
from mmgroup.tests.test_axes.get_sample_axes import do_test_sample_axes
from mmgroup.tests.test_axes.get_baby_sample_axes import import_baby_sample_axes
from mmgroup.tests.test_axes.get_baby_sample_axes import do_test_baby_sample_axes
from mmgroup.tests.test_axes.beautify_axes import compute_beautifiers



NTESTS = 100


def display_A(A):
   for i in range(24):
      print("  ", end = "")
      for j in range(24):
          if i == j or A[i,j] != 0:
              print("%2d" % A[i,j], end = " ")
          else:
              print(" .", end = " ")
      print("")


@pytest.mark.axes 
@pytest.mark.slow 
def test_import(verbose = 0):
    baby_sample_axes = import_baby_sample_axes(verbose = verbose)
    do_test_baby_sample_axes(baby_sample_axes)
    sample_axes = import_sample_axes(verbose = verbose)
    do_test_sample_axes(sample_axes)
    beautfiers = compute_beautifiers(sample_axes.g_strings)
    if verbose:
        print("A matrices of axes")
    for i in range(12):
        g = MM0(sample_axes.g_strings[i])
        v =  MMVectorCRT(20, sample_axes.v_start)
        v *= g
        g1 = MM0(beautfiers[i])
        v *= g1
        Afloat = 128 * v["A"]
        A = np.array(Afloat, dtype = np.int32)
        assert (A == Afloat).all()
        if verbose:
            print("\nclass " + sample_axes.g_classes[i])
            display_A(A)

             




