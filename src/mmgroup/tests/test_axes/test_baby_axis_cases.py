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
    global gen_leech2_reduce_type2_ortho
    global short, span, radical
    global leech_type
    global display_norm_A
    global eval_A_vstart
    global get_baby_sample_axes, get_sample_axes
    from mmgroup.mm_reduce import mm_reduce_2A_axis_type
    from mmgroup.generators import gen_leech2_reduce_type2_ortho
    from mmgroup.tests.test_axes.test_reduce_axis import short, span, radical
    from mmgroup.tests.test_axes.test_reduce_axis import leech_type
    from mmgroup.tests.test_axes.test_import import display_norm_A
    from mmgroup.tests.test_axes.test_reduce_axis import eval_A_vstart
    from mmgroup.tests.axes.get_baby_sample_axes import get_baby_sample_axes
    from mmgroup.tests.axes.get_sample_axes import get_sample_axes



V = MMV(15)


V_START_TUPLE = ("I", 3, 2)
V_START = V(*V_START_TUPLE)
V_OPP = V_START * MM0('x', 0x200)

v_start = 0x200


FINAL_2A0_AXES = [  V_OPP * MM0('t', e) for e in range(3) ]


#######################################################################################

S_KER = """  The kernel of U contains a Leech lattice vector of v type %d."""



def reduce_2A0(v, verbose = 0):
    vt = mm_reduce_2A_axis_type(v.data) & 0xffffff
    return [vt ^ v_start]

def reduce_2B(v, verbose = 0):
    return span(v, 4, verbose)

def reduce_4A(v, verbose = 0):
    if verbose:
        print(S_KER % (4))
    return [mm_reduce_2A_axis_type(v.data) & 0xffffff]

def reduce_4BC(v, verbose = 0):
    return radical(v, 1, verbose)

def reduce_6A(v, verbose = 0):
    if verbose:
        print(S_KER % (2))
    vt = mm_reduce_2A_axis_type(v.data) & 0xffffff
    a = short(v, 5, verbose)
    if verbose:
        s = "  Check vectors v + x for all x in that set of vectors"
        print(s)
    return [v ^ vt for v in a]

def reduce_6C(v, verbose = 0):
    return span(v, 3, verbose)



def reduce_10A(v, verbose = 0):
    v2 = short(v, 1, verbose)
    assert len(v2) == 100, hex(len(v2))
    v0 = short(v, 3, verbose)
    assert len(v0) == 1, hex(len(v0))
    if verbose:
        s0 = "  Check the sums x + y, with x in 1st and y in 2nd set."
        print(s0)
    return [x ^ v0[0] for x in v2]



reduce_cases = {
   "2A1": reduce_2A0,
   "2A0": reduce_2A0,
   "2B1": reduce_2B,
   "2B0": reduce_2B,
   "4A1": reduce_4A,
   "4B1": reduce_4BC,
   "4C1": reduce_4BC,
   "6A1": reduce_6A,
   "6C1": reduce_6C,
   "10A1": reduce_10A,
}



#######################################################################################


reduce_targets = {
   "2A1": None,    
   "2A0": ["2A1"],    
   "2B1": ["2A0"],    
   "2B0": ["2A0"],    
   "4A1": ["2A0"],    
   "4B1": ["2B0", "2B1"],    
   "4C1": ["2B0"],    
   "6A1": ["4A1"],    
   "6C1": ["4A1"],    
   "10A1": ["6A1"],    
}

#######################################################################################

TAU = MM0('t',1)

def baby_axis_type(v, e = 0):
    axis_type = v.axis_type(e)
    axis_type += str(int(eval_A_vstart((v * TAU**e).data) != 0))   
    return axis_type






@pytest.mark.axes
def test_cases(verbose = 0):
    import_all()
    if verbose: print("\n")
    r = np.zeros(10, dtype = np.uint32)
    for axis_type, axis in get_baby_sample_axes().items():
        # Construct an axis v of the given axi type
        v = axis.v15
        if verbose:
            print("\nTest reduction of axis type %s" % axis_type)
            for s in display_norm_A(axis).split("\n"):
                if s: print("  " + s)
            print("  A(v_start) = %d mod 15" % eval_A_vstart(v.data))
        target_axes = reduce_targets[axis_type]
        if target_axes is None:
            if verbose:
                print("  Reduction terminates here")
            continue
        nfound = 0
        dict_found = defaultdict(int)
        for w in reduce_cases[axis_type](v, verbose):
            if leech_type(w) != 4 or leech_type(w ^ v_start) != 2:
                continue
            w ^= v_start
            len_r = gen_leech2_reduce_type2_ortho(w, r)
            assert len_r >= 0
            v1 = v * MM0('a', r[:len_r])
            nfound += 1
            ok = False
            if axis_type == "2A0":
                e = 1 + (v1["C",3,2] == 2)
                v2 = v1 * MM0('t', e)
                ok = v2 == V_OPP
                dict_found = {"2A1":1}
            else:
                for e in [1,2]:
                    t_axis_type = baby_axis_type(v1, e)
                    #print(t_axis_type, target_axes)
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

