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




from mmgroup.tests.axes.axis import G # A Monster element is of type G
from mmgroup.tests.axes.axis import Axis, BabyAxis
from mmgroup.tests.axes.axis import G_CENTRAL, G_AXIS, G_AXIS_OPP
from mmgroup.tests.axes.axis import V_AXIS, V_AXIS_OPP

V15 = MMV(15)


from mmgroup.tests.axes.get_sample_axes import next_generation_pool

PROCESSES = 0


########################################################################
########################################################################
# Obtaining samples of transformed axes 
########################################################################
########################################################################


########################################################################
# Generate a random element of the subgroup 2^{2+22}.Co_2 of G_x0
########################################################################

def rand_Co_2():
    return G('r', 'B & G_x0')


########################################################################
# Next generation of axes
########################################################################


def spread(gv):
    g = rand_Co_2() * G('t','n')
    g_new = gv * g
    g_new.stage = gv.stage + 1
    return g_new



########################################################################
# Profiling an axis
########################################################################



def mark(gv):
    return  gv.axis_type(), gv.fixed_value(), gv.fixed_value('B')

def score(gv):
    return 576 - np.count_nonzero(gv['A'])


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
    gv0 = BabyAxis()
    gv0.stage = 1
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
            m = mark(v)
            all_samples[m].append(str(v.g))
    num_samples = sum(len(x) for x in all_samples.values())
    print("Number of axes considered:", num_samples)
    #assert len(sample_list) == 10
    return sample_list




########################################################################
########################################################################
# Tweak axes in the Monster so that they fit into the Baby Monster 
########################################################################
########################################################################

PI = AutPL(0, zip(range(16), list(range(8,16)) + list(range(8))))

TWEAKS = {
    '2A' : [ [('x', 0x200)], [('y', 0x200)] ],
    '2B' : [ [('p', PI)] ],
    '4C' : [ [('p', PI)] ],
}

def parent_samples():
    from mmgroup.tests.axes.get_sample_axes import get_sample_axes
    M1 = MM0(1)
    for parent_orbit, parent_axis in get_sample_axes().items():
        yield parent_orbit, parent_axis, M1
        if parent_orbit in TWEAKS:
            for g in TWEAKS[parent_orbit]:
                 yield parent_orbit, parent_axis, MM0(g)




def useful_baby_sample_axes():
    """Tweak axes to fit into the Baby Monster

    This function tries to tweak the (nice) standard 2A axes so that
    we obtain a sample axis for each class of the Baby Monster axes.
    The function returns a dictionary mapping the classes of 2A
    axes to nice sample axes.
    """
    ref_axes = {}
    for parent_orbit, parent_axis, g in parent_samples():
        #print("case", parent_orbit, g)
        try:
            baby_axis = BabyAxis(parent_axis * g)
        except:
            baby_axis = None
        if baby_axis:
            baby_orbit = baby_axis.axis_type()
            ref_axes[baby_orbit] = baby_axis
            #print("Represetative for baby axis orbit %s is axis %s * %s" %
            #    (baby_orbit, parent_orbit, g)
            #)
        #print("case", parent_orbit, "done")
    return ref_axes



########################################################################
########################################################################
# Write list of axes to file "get_baby_sample_axes.py"
########################################################################
########################################################################


DIR = os.path.split(__file__)[0]
FILE = "baby_sample_axes"
PATH = os.path.join(DIR, FILE)



def baby_axis_type(gv):
    #at = axis_type(gv.g, gv.g0)
    #asub = int(gv.baby_value_A() != 0)
    #return at + str(asub) 
    return gv.axis_type()





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


Qx0_text = """

Qx0_equations = [
%s
]

"""




def sample_list_sort_key(sample_entry):
     stage, sample, _ = sample_entry 
     s_class = baby_axis_type(sample)
     return stage, int(s_class[:-2]), s_class[-2:-1], -int(s_class[-1])


def write_axes(sample_list, beautify = True, verbose = False):
    if beautify:
        ref_baby_axes = useful_baby_sample_axes()
        from mmgroup.tests.axes.equations import compute_Qx0_equations_str
    sample_list = explore_axes(5, 40, 30, verbose = verbose)
    sample_list.sort(key = sample_list_sort_key)
    s_samples, s_stages, s_marks = "", "", ""
    s_beautiful_samples = ""
    s_classes = ""
    s_equations = ""
    for i, (stage, sample, mark) in enumerate(sample_list):
        g = sample.g.raw_str()
        s_samples += "\"" + g + "\",\n"
        s_stages += str(stage) + ", "
        s_marks +=   "%s,\n" % str(unnumpy(mark))
        class_ = baby_axis_type(sample)
        #print("class", class_)
        s_classes += '"' + class_ + '", '
        try:
            if beautify:
                #axis = beautify_baby_axis(class_, g)
                axis = ref_baby_axes[class_]
                assert isinstance(axis, BabyAxis)
                axis.rebase()
                s_beautiful_samples += '"' + axis.g.raw_str() + '",\n'
                s_equations += compute_Qx0_equations_str(axis)
        except:
            beautify = False
            raise
    if beautify:
        s_samples = s_beautiful_samples
       
    with open(PATH + ".py", "wt") as f:
        f.write(f_text % (G_CENTRAL,  
            G_AXIS_OPP, V_AXIS_OPP, 
          s_samples, s_stages, s_classes, s_marks))
        if beautify:
            f.write(Qx0_text % s_equations)



########################################################################
########################################################################
# Import list of axes 
########################################################################
########################################################################



BABY_AXES= OrderedDict()

def get_baby_axes_from_sample_axes(sample_axes):
    global BABY_AXES
    if len(BABY_AXES):
        return BABY_AXES
    for i, g1 in enumerate(sample_axes.g_strings):
        g = G(g1)
        axis = BabyAxis(g)
        g_class = sample_axes.g_classes[i]
        axis.g_class = g_class
        axis.mark = sample_axes.g_marks[i]
        axis.stage = sample_axes.g_stages[i]
        axis.Qx0_equations = [np.array(a, dtype = np.uint32)
            for a in sample_axes.Qx0_equations[i]
        ]
        axis.constant = True
        BABY_AXES[g_class] = axis
    return BABY_AXES



def compute_and_write_baby_axes(beautify = True, verbose = 0):
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
        assert axis.v_axis15 * axis.g == axis.v_axis15
        # print(orbit, axis.axis_type())



