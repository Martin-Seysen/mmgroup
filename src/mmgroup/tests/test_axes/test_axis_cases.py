
from random import randint
import numpy as np
from collections import defaultdict

import datetime
import time
import pytest

from mmgroup import MM0, MMSpace, MMV
from mmgroup.mm15 import op_2A_axis_type as mm_op15_2A_axis_type

from mmgroup.tests.test_axes.test_import import AXES, BABY_AXES
from mmgroup.tests.test_axes.test_reduce_axis import short, span, radical
from mmgroup.tests.test_axes.test_reduce_axis import leech_type

V = MMV(15)


V_START_TUPLE = ("I", 3, 2)
V_START = V(*V_START_TUPLE)
#V_OPP = V_START * MM0('x', 0x200)

v_start = 0x200


#######################################################################################


def reduce_2B(v):
    return span(v, 4)

def reduce_4A(v):
    return [mm_op15_2A_axis_type(v.data) & 0xffffff]

def reduce_4BC(v):
    return radical(v, 1)

def reduce_6A(v):
    vt = mm_op15_2A_axis_type(v.data) & 0xffffff
    a = short(v, 5)
    return [v ^ vt for v in a]

def reduce_6C(v):
    return span(v, 3)

def reduce_6F(v):
    return radical(v, 7)

def reduce_8B(v):
    v2 = short(v, 1)
    v2_0 =  v2[0]
    return [x ^ v2_0 for x in v2]

def reduce_10A(v):
    v2 = short(v, 1)
    assert len(v2) == 100, hex(len(v2))
    v0 = short(v, 3)
    assert len(v0) == 1, hex(len(v0))
    return [x ^ v0[0] for x in v2]

def reduce_10B(v):
    return radical(v, 4)

def reduce_12C(v):
    return radical(v, 7)

reduce_cases = {
   "2A": None,
   "2B": reduce_2B,
   "4A": reduce_4A,
   "4B": reduce_4BC,
   "4C": reduce_4BC,
   "6A": reduce_6A,
   "6C": reduce_6C,
   "6F": reduce_6F,
   "8B": reduce_8B,
   "10A": reduce_10A,
   "10B": reduce_10B,
   "12C": reduce_12C,
}



#######################################################################################


reduce_targets = {
   "2A": None,    
   "2B": ["2A"],    
   "4A": ["2A"],    
   "4B": ["2B"],    
   "4C": ["2B"],    
   "6A": ["4A"],    
   "6C": ["4A"],    
   "6F": ["4C"],    
   "8B": ["4A"],    
   "10A": ["6A"],    
   "10B": ["4B", "4C"],    
   "12C": ["4B", "6A"],    
}

#######################################################################################

@pytest.mark.axes
def test_cases(verbose = 0):
    for axis_type, g_str in AXES.items():
        if verbose:
            print("Test reduction of axis type %s" % axis_type)
        # Construct an axis v of the given axi type
        v = V_START * MM0(g_str)
        target_axes = reduce_targets[axis_type]
        if target_axes is None:
            if verbose:
                print("  Reduction terminates here")
            continue
        nfound = 0
        dict_found = defaultdict(int)
        for w in reduce_cases[axis_type](v):
            if leech_type(w) != 4:
                continue
            nfound += 1
            v1 = v * MM0('c', w)**(-1)
            ok = False
            for e in [1,2]:
                t_axis_type = v1.axis_type(e)
                if t_axis_type in target_axes:
                   ok = True
                   dict_found[t_axis_type] += 1
            assert ok
        assert nfound > 0
        if verbose:
            if len(dict_found):
                print("  Reduced to axis types", dict(dict_found))

 


