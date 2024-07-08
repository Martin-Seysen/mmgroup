r"""Investigate certain orbits of 2A axes of the baby monster group

The functions in this module are used to find the orbits of 2A axes
of the subgroup :math:`G2` (of structure 
:math:`2^{1+22}.\mbox{Co}_2`) of the monster group. 

"""

import sys
import os
import time
from math import floor, ceil
from random import randint, shuffle, sample
from collections import defaultdict, OrderedDict
from multiprocessing import Pool, freeze_support
import numpy as np
from operator import __or__
from functools import reduce

if __name__ == "__main__":
    sys.path.append("../../../")

from mmgroup import MM0, MMV, MMVector, Cocode, XLeech2, Parity, PLoop
from mmgroup import GcVector, AutPL
from mmgroup.bitfunctions import unnumpy


from mmgroup.tests.axes.get_sample_axes import G_CENTRAL
from mmgroup.tests.axes.get_sample_axes import G_AXIS_OPP
from mmgroup.tests.axes.get_sample_axes import V_AXIS_OPP
from mmgroup.tests.axes.get_sample_axes import g_axis, v_axis
from mmgroup.tests.axes.get_sample_axes import g_axis_opp
from mmgroup.tests.axes.get_sample_axes import v_axis_opp

from mmgroup.tests.axes.get_sample_axes import g_central, GVector
from mmgroup.tests.axes.get_sample_axes import next_generation_pool
from mmgroup.tests.axes.get_sample_axes import axis_type


########################################################################
########################################################################
# The group and the vector space to be used
########################################################################
########################################################################

G = MM0                     # The monster group
V15 = MMV(15)               # Its representation space

PROCESSES = 0


########################################################################
########################################################################
# Obtaining samples of transformed axes 
########################################################################
########################################################################


########################################################################
# Generate a random element of Co_2
########################################################################

FIXED_PAIR = [3,2]
FIXED_TUPLE = ("I", 3, 2)


PI22 = [0,1] + list(range(4,24))
PI7 = [2,3,0,1,4,5,8]

StdCocodeVector = Cocode(FIXED_PAIR)
EvenGCodeMask = None
for i in range(12):
    if (PLoop(1 << i) & StdCocodeVector).ord:
        assert EvenGCodeMask == None
        EvenGCodeMask = 0xfff ^ (1 << i)
assert EvenGCodeMask is not None




def rand_pi_m22():
    r = randint(0,1)
    img = [2+r, 3-r] + sample(PI22, 3)
    syn = GcVector(img).syndrome_list()
    compl = list(set(range(24)) - set(img + syn))
    img += sample(syn, 1) + sample(compl, 1)
    result =  AutPL('r', zip(PI7,img))
    return result

def rand_xy():
    return randint(0, 0xfff) & EvenGCodeMask



def rand_Co_2():
    m = G()
    for i in range(3):
        m *= G(rand_pi_m22())
        m *=  G('l', "n")
    for tag in "xy":
        m *= G(tag, rand_xy())
    assert v_axis * m == v_axis
    return m


########################################################################
# Next generation of axes
########################################################################



def spread(gv):
    g = rand_Co_2() * G('t','n')
    return gv * g



########################################################################
# Profiling an axis: count entries for tags B, C, T, X
########################################################################



def mark(gv):
    return (gv.baby_value_A(), gv.mark())


def score(gv):
    return gv.count_zeros_in_A()


########################################################################
# Collect one sample axis for each class
########################################################################

# Here will be a large directory samples found so far, if anybody wants 
# these samples. In that directory each key is a watermark, and each 
# value is a list of strings corresponding to group elements.
all_samples = defaultdict(list)


def _show_sample(i, sample):
    stage, gv, mark = sample
    print("Vector %d, stage = " % i, stage)
    print("g = ", gv.g)
    print("mark = ", mark)
    print("dihedral class:", baby_axis_type(gv))


def explore_axes(stages, n_spread, n_keep, verbose = 0):
    global all_samples
    gv0 = GVector(opp = True)
    gv_list = [gv0]
    m = mark(gv0)
    marks = set([m])
    sample_list = [(0, gv0, m)]
    if verbose:
        print("Start vector =", gv0.v)
        _show_sample(0, sample_list[0])
    for i in range(1, stages):
        gv_list,  new_samples = next_generation_pool(
           gv_list,  
           marks,
           f_spread = spread, 
           f_mark = mark, 
           f_score = score, 
           n_spread = n_spread, 
           n_keep = n_keep, 
           processes = PROCESSES,
           verbose = verbose,
        )
        for m in sorted(new_samples):
            if not m in marks: 
                sample_list.append((i, new_samples[m], m))
                if verbose:
                    _show_sample(len(sample_list)-1, sample_list[-1])
        marks |= new_samples.keys()
        for v in gv_list:
            m = v.mark()
            all_samples[m].append(str(v.g))
    num_samples = sum(len(x) for x in all_samples.values())
    print("Number of axes considered:", num_samples)
    #assert len(sample_list) == 10
    return sample_list



########################################################################
########################################################################
# Write list of axes to file "get_baby_sample_axes.py"
########################################################################
########################################################################


DIR = os.path.split(__file__)[0]
FILE = "baby_sample_axes"
PATH = os.path.join(DIR, FILE)



def baby_axis_type(gv):
    at = axis_type(gv.g, gv.g0)
    asub = int(gv.baby_value_A() != 0)
    return at + str(asub) 






f_text = """# This file has been generated automatically. Do not change!
# It contains samples of the 10 cosets of 2A axes othogonal to the
# standard 2A axis wrt 2^{1+23}.Co_2.
#

g_central = "%s"
g_axis = "%s"
v_start = "%s" 

g_strings = [
%s
]

g_stages = [
%s
]

g_classes = [
%s
]

g_marks = [
%s
]

"""



def sample_list_sort_key(sample_entry):
     stage, sample, _ = sample_entry 
     s_class = baby_axis_type(sample)
     return stage, int(s_class[:-2]), s_class[-2:-1], -int(s_class[-1])


def write_axes(sample_list, beautify = True, verbose = False):
    if beautify:
        from mmgroup.tests.axes.beautify_baby_axes import beautify_baby_axis
    sample_list = explore_axes(5, 40, 30, verbose = verbose)
    sample_list.sort(key = sample_list_sort_key)
    s_samples, s_stages, s_marks = "", "", "" 
    s_classes = ""
    for i, (stage, sample, mark) in enumerate(sample_list):
        g = sample.g.raw_str()
        s_samples += "\"" + g + "\",\n"
        s_stages += str(stage) + ", "
        s_marks += str(unnumpy(mark))  + ",\n"
        class_ = baby_axis_type(sample)
        #print("class", class_)
        s_classes += '"' + class_ + '", '
        try:
            if beautify:
                axis = beautify_baby_axis(class_, g, _vb)
                axis.rebase(reduce = True)          
                s_beautiful_samples += '"' + axis.g.raw_str() + '",\n'
        except:
            beautify = False
    if beautify:
        s_samples = s_beautiful_samples
       
    with open(PATH + ".py", "wt") as f:
        f.write(f_text % (G_CENTRAL,  
            G_AXIS_OPP, V_AXIS_OPP, 
          s_samples, s_stages, s_classes, s_marks))



########################################################################
########################################################################
# Import list of axes 
########################################################################
########################################################################



BABY_AXES= {}

def get_baby_axes_from_sample_axes(sample_axes):
    global BABY_AXES
    if len(BABY_AXES):
        return BABY_AXES
    from mmgroup.tests.axes.beautify_baby_axes import BabyAxis
    for i, g1 in enumerate(sample_axes.g_strings):
        g = MM0(g1)
        axis = BabyAxis(g)
        g_class = sample_axes.g_classes[i]
        axis.g_class = g_class
        axis.mark = sample_axes.g_marks[i]
        axis.stage = sample_axes.g_stages[i]
        axis.constant = True
        BABY_AXES[g_class] = axis
    return BABY_AXES



def compute_and_write_baby_axes(beautify = True, verbose = 0):
    v0 = GVector()
    sample_list = explore_axes(5, 40, 30, verbose = verbose)
    write_axes(sample_list, beautify, verbose)
    time.sleep(0.1)


def do_test_baby_sample_axes(baby_sample_axes):
    sax =  baby_sample_axes
    l = [G(sax.g_central), G(sax.g_axis), MMVector(7, sax.v_start)]
    lg =  [G(x) for x in sax.g_strings]
    assert len(lg) == 10


def get_baby_sample_axes(calculate = False, beautify = True, verbose = 0):
    if calculate:
        compute_and_write_baby_axes(beautify = True, verbose = verbose)
        time.sleep(0.1)
    try:
        from mmgroup.tests.axes import baby_sample_axes
    except ImportError: 
        write_axes(verbose)
        time.sleep(0.1)
        from mmgroup.tests.axes import baby_sample_axes
    do_test_baby_sample_axes(baby_sample_axes)
    return get_baby_axes_from_sample_axes(baby_sample_axes)



########################################################################
########################################################################
# Main program 
########################################################################
########################################################################





if __name__ == "__main__":
    if PROCESSES != 1: 
        freeze_support()
    axes = get_baby_sample_axes(calculate = 1, beautify = 1, verbose = 1)
    for orbit, axis in axes.items():
        assert orbit == axis.axis_type()
        assert v_axis * axis.g == v_axis
        # print(orbit, axis.axis_type())



