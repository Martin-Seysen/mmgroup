
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os

import pytest



from mmgroup.tests.test_axes.get_sample_axes import import_sample_axes
from mmgroup.tests.test_axes.get_sample_axes import do_test_sample_axes
from mmgroup.tests.test_axes.get_baby_sample_axes import import_baby_sample_axes
from mmgroup.tests.test_axes.get_baby_sample_axes import do_test_baby_sample_axes

NTESTS = 100

@pytest.mark.axes 
@pytest.mark.slow 
def test_import():
    baby_sample_axes = import_baby_sample_axes(verbose = True)
    do_test_baby_sample_axes(baby_sample_axes)
    sample_axes = import_sample_axes(verbose = True)
    do_test_sample_axes(sample_axes)
    




