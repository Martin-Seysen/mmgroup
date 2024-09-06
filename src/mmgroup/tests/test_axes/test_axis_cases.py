import sys
from random import randint
import numpy as np
from collections import defaultdict

import datetime
import time
import pytest

if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import MM0, MMSpace, MMV


def import_all():
    global mm_reduce_2A_axis_type
    global short, span, radical
    global leech_type
    global display_norm_A
    global get_sample_axes
    from mmgroup.mm_reduce import mm_reduce_2A_axis_type
    from mmgroup.tests.axes.get_sample_axes import get_sample_axes
    from mmgroup.tests.test_axes.test_reduce_axis import short, span, radical
    from mmgroup.tests.test_axes.test_reduce_axis import leech_type
    from mmgroup.tests.test_axes.test_import import display_norm_A


V = MMV(15)


V_START_TUPLE = ("I", 3, 2)
V_START = V(*V_START_TUPLE)
#V_OPP = V_START * MM0('x', 0x200)

v_start = 0x200


#######################################################################################



S_KER = """  The kernel of U contains a Leech lattice vector of v type %d."""

def reduce_2B(v, verbose):
    return span(v, 4, verbose)

def reduce_4A(v, verbose):
    if verbose:
        print(S_KER % ( 4))
    return [mm_reduce_2A_axis_type(v.data) & 0xffffff]

def reduce_4BC(v, verbose):
    return radical(v, 1, verbose)

def reduce_6A(v, verbose):
    if verbose:
        print(S_KER % (2))
    vt = mm_reduce_2A_axis_type(v.data) & 0xffffff
    a = short(v, 5, verbose)
    if verbose:
        s = "  Check vectors v + x for all x in that set of vectors"
        print(s)
    return [v ^ vt for v in a]

def reduce_6C(v, verbose):
    return span(v, 3, verbose)

def reduce_6F(v, verbose):
    return radical(v, 7, verbose)

def reduce_8B(v, verbose):
    v2 = short(v, 1, verbose)
    v2_0 =  v2[0]
    if verbose:
        s0 = "  Let v be any vector in that set"
        s1 = "  We add v to all other vectors iin that set"
        print(s0, "\n" + s1)
    return [x ^ v2_0 for x in v2]

def reduce_10A(v, verbose):
    v2 = short(v, 1, verbose)
    assert len(v2) == 100, hex(len(v2))
    v0 = short(v, 3, verbose)
    assert len(v0) == 1, hex(len(v0))
    if verbose:
        s0 = "  Check the sums x + y, with x in 1st and y in 2nd set."
        print(s0)
    return [x ^ v0[0] for x in v2]

def reduce_10B(v, verbose):
    return radical(v, 4, verbose)

def reduce_12C(v, verbose):
    return radical(v, 7, verbose)

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
    import_all()
    s = "For an axis let A be the symmetric matrix corresponding to part 300x"
    if verbose:
        print("\n%s\n" % s)
    for axis_type, axis in get_sample_axes().items():
        if verbose:
            print("\nTest reduction of axis type %s" % axis_type)
            text = display_norm_A(axis).split("\n")
            for s in text:
                if s: print("  " + s)
        # Construct an axis v of the given axi type
        v = axis.v15
        target_axes = reduce_targets[axis_type]
        if target_axes is None:
            if verbose:
                print("  Reduction terminates here")
            continue
        nfound = 0
        dict_found = defaultdict(int)
        for w in reduce_cases[axis_type](v, verbose):
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

 

if __name__ == "__main__":
    test_cases(verbose = 1)

