r"""Investigate certain orbits of 2A axes of the monster group

The functions in this module are used to find the orbits of 2A axes
of the subgroup :math:`G_{x0}` (of structure 
:math:`2^{1+24}.\mbox{Co}_1`) of the monster group. According to
|Nor98| there are 12 such orbits.

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


from mmgroup import MM0, MMV, MMVector, Cocode, XLeech2
from mmgroup.bitfunctions import unnumpy


from mmgroup.tests.axes.axis import G, Axis, BabyAxis
from mmgroup.tests.axes.axis import G_CENTRAL, G_AXIS, G_AXIS_OPP
from mmgroup.tests.axes.axis import V_AXIS, V_AXIS_OPP

MM = MM0 = G
V15 = MMV(15)


try:
    import  mmgroup.tests.axes.sample_axes
    PROCESSES = 0
except:
    # No multiprocessing before all data have been calculated
    PROCESSES = 1

# Use extended watermarking if the following variable is set.
# If this is True then the calculation will take a long time!
EXTENDED_WATERMARK = False

########################################################################
########################################################################
# An abstract random walk through a part of a graph with profiling
########################################################################
########################################################################

def next_generation(
        obj_list,   # a list of objects
        marks,      # a set of watermarks already found
        f_spread,   # function, returns a random object from an object
        f_mark   ,  # function, returns a watermark of an object
        f_score,    # function, returns a score of an object
        n_spread,   # No of objects to generate from each object
        n_keep,     # No of objects to be kept for each watermark
        verbose=0
        ):
    """Find orbits of subgroup of a permutation group.

    Assume that a finite group ``G`` operates of a finite set ``V`` 
    of objects. Assume further that the objects in the set ``V`` can
    be watermarked in such a way that two such objects get the same
    watermark if they are in the same orbit of a subgroup ``H`` 
    of ``G``.

    Here we also assume that ``G`` and ``V`` are extremely large
    sets, but ``H`` has few orbits on ``V``. Then we have a chance
    to find all (or at least many) orbits of ``H`` of ``V`` as 
    follows.

    Let ``K`` be a subgroup of ``G`` containing ``H``. Let
    ``obj_list`` be a set of certain elements of ``V``. We multiply
    each object ``o`` in ``obj_list`` with a number of different
    random elements of ``K``. We watermark all products obtained
    that way. For each yet unknown watermark we keep a certain 
    number of objects with that watermark, and  we return the list
    ``new_obj_list`` of all objects (i.e. elements of ``V``) 
    obtained in that way.

    We also return a dictionary ``samples`` that maps each newly 
    found watermark to an element of ``V`` with that watermark.

    Repeating this function with output ``new_obj_list`` as new
    input ``obj_list`` returns more orbits of the group ``H``. 

    In a typical use case ``G`` is the monster group, ``V`` is
    the set of 2A axes in the 196884-dimensional representation
    of the monster, ``H`` is the subgroup :math:`G_{x0}` of the
    monster, and ``K = H T``, where ``T`` is the group generated
    by the triality element of the monster, see :cite:`Con85`. 
    Then the 12 orbits of ``G`` in ``V`` are known, see
    :cite:`Nor98`, and we can easily find them with this function.
    
    :param obj_list:

        List of elements of the set ``V`` to be modified. The 
        functions given by parameters  ``f_spread``, ``f_mark``,
        ``f_score`` must accept each entry of this list as a single
        parameter.

    :param marks:

        A set of watermarks already known. If new elements of the
        set ``V`` with a known watermark are generated then they
        are ignored.

    :param f_spread:

        This function takes an element ``v`` of ``V`` as input. It
        should apply a random element of the group ``K`` to ``v``
        and return the result as an element of ``V``. The function
        may return ``None`` if element ``v`` has no offspring.
        
    :param f_mark:

        This function takes an element ``v`` of ``V`` as input
        and returns the watermark of ``v``. The returned result
        must be immutable.

    :param f_score:

        This function takes an element ``v`` of ``V`` and returns a
        score for that element as a real number. If several elements
        with the same watermark are found then an element with the
        highest score is returned as a sample element for that
        watermark. Put ``f_score = 0`` to ignore the score.

    :param n_spread:   

        Number of calls to function ``f_spread`` applied to each
        entry in the input list ``obj_list``.

    :param n_keep:   

        Number of entries to be kept in the output list
        ``new_obj_list`` for each new watermark.

    :param: verbose:

        Verbose operation if True

    
    :return: A pair ``(new_obj_list,  samples)``.
             Here ``new_obj_list`` is the list of new elements of
             ``V`` that have been generated and kept. ``samples``
             is a dictionary that maps each new watermark found
             to an element of ``V`` with that watermark. 
 

    """
    if not f_score: f_score = lambda x : 0
    shuffle(obj_list)
    new_marks = defaultdict(list)
    samples = {}
    for i in range(n_spread):
        for gv in obj_list:
            gv_new = f_spread(gv)
            if gv_new is not None:
                m = f_mark(gv_new)
                if m not in marks:
                    new_marks[m].append(gv_new)
    if verbose:
        len_t = sum([len(x) for x in new_marks.values()]) 
        print("No of candidates tested: ", len_t)
    new_obj_list = []
    for m, obj_list in new_marks.items():
        new_list = sample(obj_list, min(n_keep, len(obj_list)))
        score_list = [f_score(x)  for x in new_list]
        best = new_list[score_list.index(max(score_list))]
        samples[m] = best      
        new_obj_list += new_list
    #marks = marks | new_marks.keys()
    if verbose:
        print("No of candidates kept: ", len(new_obj_list))
    return new_obj_list,  samples



########################################################################
# Multitasking version of the last function 
########################################################################


def _spread(obj_list, n_spread, f_spread, f_mark, f_score, marks):
    new_marks = defaultdict(list)
    for i in range(n_spread):
        for gv in obj_list:
            gv_new = f_spread(gv)
            if gv_new is not None:
                m = f_mark(gv_new)
                if m not in marks:
                    new_marks[m].append((f_score(gv_new), gv_new))
    return new_marks

def next_generation_pool(
        obj_list,   # a list of objects
        marks,      # a set of watermarks already found
        f_spread,   # function, returns a random object from an object
        f_mark,     # function, returns a watermark of an object
        f_score,    # function, returns a score of an object
        n_spread,   # No of objects to generate from each object
        n_keep,     # No of objects to be kept for each watermark
        processes=PROCESSES, # No of processes
        verbose=0
        ):
    """Multitasking version of function ``next_generation``

    Input parameters and return value are as in function
    ``next_generation``.

    The additional parameter ``processes`` specifies the 
    maximum number of processes to be started. ``processes = 0``
    (default) means to start as many processes as possible.
    """
    if processes == 1:
        return next_generation( obj_list, marks, f_spread,
             f_mark, f_score, n_spread,  n_keep,  verbose)
    if not f_score: f_score = lambda x : 0
    shuffle(obj_list)
    new_marks = defaultdict(list)
    samples = {}
    if not processes:
        processes = os.cpu_count()
    else:
        processes = max(1, min(processes, os.cpu_count()))
    n_rounds = ceil(n_spread / processes)

    arg_list = [(obj_list, n_rounds, f_spread, f_mark, f_score, marks)
         ] * processes
    with Pool(processes = processes) as pool:
        results = pool.starmap(_spread, arg_list, chunksize = 1)
    pool.join()
    for result in results:
        for m, value in result.items():
            new_marks[m] += value
    if verbose:
        len_t = sum([len(x) for x in new_marks.values()]) 
        print("No of candidates tested: ", len_t)
    new_obj_list = []
    for m, (obj_list) in new_marks.items():
        scores =  [x[0] for x in obj_list]
        best = max(scores)
        samples[m] = obj_list[scores.index(best)][1]      
        new_list = sample(obj_list, min(n_keep, len(obj_list)))
        new_obj_list += [x[1] for x in new_list]
    marks = marks | new_marks.keys()
    if verbose:
        print("No of candidates kept: ", len(new_obj_list))
    return new_obj_list,  samples



########################################################################
########################################################################
# Obtaining samples of transformed axes (there are 12 classes)
########################################################################
########################################################################


########################################################################
# Next generation of axes
########################################################################



def spread(gv):
    g = G([('r','G_x0'), ('t','n')])
    g_new = gv * g
    g_new.stage = gv.stage + 1
    return g_new



########################################################################
# Profiling an axis: count entries for tags B, C, T, X
########################################################################



def mark(gv):
    return gv.axis_type(), gv.v15.count_short(), gv.axis_type_info()

if not EXTENDED_WATERMARK:
    def mark(gv):
        return gv.axis_type(), gv.v15.count_short()

def score(gv):
    return 576 - np.count_nonzero(gv['A'])


########################################################################
# Collect one sample axis for each class
########################################################################

# Number of axis classes
NUM_AXIS_CLASSES = 12

# Here will be a large directory samples found so far, if anybody wants 
# these samples. In that directory each key is a watermark, and each 
# value is a list of strings corresponding to group elements.
all_samples = defaultdict(list)


def _show_sample(i, sample):
    stage, gv, mark = sample
    print("Vector %d, stage = " % i, stage)
    print("g = ", gv.g)
    print("mark = ", mark)
    print("dihedral class:", axis_type(gv.g))


def explore_axes(gv0, stages, f_spread, f_mark, f_score, n_spread, n_keep, verbose = 0):
    global all_samples
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
           f_spread = f_spread, 
           f_mark = f_mark, 
           f_score = f_score, 
           n_spread = n_spread, 
           n_keep = n_keep, 
           processes = 1,
           verbose = verbose,
        )
        for m in sorted(new_samples):
            if not m in marks: 
                sample_list.append((i, new_samples[m], m))
                if verbose:
                    _show_sample(len(sample_list)-1, sample_list[-1])
        marks |= new_samples.keys()
        for v in gv_list:
            m = f_mark(v)
            all_samples[m].append(str(v.g))
        if not EXTENDED_WATERMARK and len(sample_list) >= NUM_AXIS_CLASSES: 
            break
    num_samples = sum(len(x) for x in all_samples.values())
    print("Number of axes considered:", num_samples)
    assert len(sample_list) == NUM_AXIS_CLASSES
    mark_list = [mark[0] for stage, axis, mark in sample_list]
    assert len(set(mark_list)) == NUM_AXIS_CLASSES
    return sample_list




########################################################################
########################################################################
# Matching the names of the axis orbits with those in |Nor98|
########################################################################
########################################################################


"""Character of products of a 2A and a 2B involution

From |Nor98| we see that such a product must be in one of the 
following classes in the monster (in ATLAS notation):

2AB, 4ABC, 6ACF, 8B, 10AB, 12C

From the ATLAS we obtain the information about the characters of 
these classes stored in the following dictionary AXES_PROD_CLASSES.
In that dictionary we replace the character value of a class by 
None whenever the corresponding class powers up to class 2A.
"""

AXES_PROD_CLASSES = {
   2:  ("AB", [None, 275]),
   4:  ("ABC", [275, None, 19]),
   6:  ("ACF", [None, 14, -1]),
   8:  ("B", [11]),
  10:  ("AB", [None, 5]),
  12:  ("C", [None]),
}


AXES_CLASSES_DICT = OrderedDict()
for number, (letters, characters) in AXES_PROD_CLASSES.items():
    for letter, character in zip(letters, characters):
        AXES_CLASSES_DICT[(number, character)] = str(number) + letter


def compute_axis_type(axis, check = True):
    """Compute type of a 2A axis from character information

    Let ``g`` be a  2A-element  of the Monster group corresponding
    to the given ``axis``, and let ``x`` be the central involution
    of the subgroup :math:`G_{x0}` of the monster. We return the
    class of the ``h = g*x`` in ATLAS notation as a string. Here we
    assume that all possible classes are mentioned in |Nor98|   
    """
    o, chi, _ = axis.axis_type_info()
    ax_type = AXES_CLASSES_DICT[o, chi]
    expected_type = axis.axis_type()
    if check:
        assert ax_type == expected_type, (ax_type, expected_type)
    return ax_type




########################################################################
########################################################################
# Write list of axes to file "get_sample_axes.py"
########################################################################
########################################################################


DIR = os.path.split(__file__)[0]
FILE = "sample_axes"
PATH = os.path.join(DIR, FILE)


AXIS_GROUPS = {
   "2A": "2^(1+23).Co_2",
   "2B": "2^(2+8+16).O_8^+(2)",
   "4A": "2^(1+22).M_23",
   "4B": "(2^7 x 2^(1+8)).S_6(2)",
   "4C": "2^(1+14+5).A_8",
   "6A": "2^2.U_6(2).2",
   "6C": "2^(3+8).(3 x U_4(2)).2",
   "6F": "2^(1+8).A_9",
   "8B": "2.2^10.M_11",
  "10A": "2.HS.2",
  "10B": "2^(1+8).(A_5 x A_5).2",
  "12C": "2 x S_6(2)",
}

AXIS_POWERS = {
   "2A": "",
   "2B": "",
   "4A": "2B",
   "4B": "2A",
   "4C": "2B",
   "6A": "2A",
   "6C": "2B",
   "6F": "2B",
   "8B": "4A,2B",
  "10A": "2A",
  "10B": "2B",
  "12C": "6A,4B,2A",
}



f_text = """# This file has been generated automatically. Do not change!
# It contains samples of the 12 cosets of 2A axes wrt 2^{1+24}.Co_1.
#

g_central = "%s"
g_axis = "%s"
v_start = "%s" 
g_axis_opp = "%s"
v_start_opp = "%s" 

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


axes_text = """

powers= [
%s
]

groups = [
%s
]

"""


Qx0_text = """

Qx0_equations = [
%s
]

"""


def sample_list_sort_key(sample_entry):
     stage, axis, _ = sample_entry 
     s_class = axis.axis_type() 
     return stage, int(s_class[:-1]), s_class[-1:]


def beautify_12C(sample_list):
    from mmgroup.tests.axes.beautify_axes import beautify_axis
    from mmgroup.tests.axes.beautify_axes import Axis
    samples = {}
    for stage, axis, mark in sample_list:
        class_ = axis.axis_type()
        if class_ in ['4B', '12C']:
            samples[class_] = beautify_axis(axis.copy(), check = True)
    t, z = samples['12C'].g_axis, samples['12C'].g_central
    old_4B = samples['4B']
    new_4B = Axis('i', t * z * t * z * t)
    for tag in 'ABC':
        assert (old_4B[tag] == new_4B[tag]).all()
    samples['4B'] = new_4B
    return samples
       

def write_axes(sample_list, beautify = True, verbose = False):
    from mmgroup.tests.axes.beautify_axes import Axis
    if beautify:
        from mmgroup.tests.axes.beautify_axes import beautify_axis
        from mmgroup.tests.axes.equations import compute_Qx0_equations_str
    sample_list.sort(key = sample_list_sort_key)
    s_samples, s_stages, s_marks = "", "", ""
    s_classes = ""
    s_beautiful_samples = ""
    s_powers = ""
    s_groups = ""
    s_equations = ""
    axes = [x[1] for x in sample_list]
    _vb = min(verbose, 1)
    special_samples = beautify_12C(sample_list)
    for i, (stage, axis, mark) in enumerate(sample_list):
        g = axis.g.raw_str()
        s_samples += "\"" + g + "\",\n"
        s_stages += str(stage) + ", "
        s_marks += str(unnumpy(mark))  + ",\n"
        class_ = compute_axis_type(axis, check = True)
        s_classes += '"' + class_ + '", '
        try:
            if beautify:
                if class_ in special_samples:
                    axis = special_samples[class_]
                else:     
                    axis = beautify_axis(Axis(g), False, True, _vb)
                axis.rebase()
                s_beautiful_samples += '"' + axis.g.raw_str() + '",\n'
                s_equations += compute_Qx0_equations_str(axis)
        except:
            beautify = False
            raise
        s_powers += '"' + AXIS_POWERS[class_] + '",\n'
        s_groups += '"' + AXIS_GROUPS[class_] + '",\n'
    if beautify:
        s_samples = s_beautiful_samples

    with open(PATH + ".py", "wt") as f:
        f.write(f_text % (G_CENTRAL, G_AXIS, V_AXIS,
            G_AXIS_OPP, V_AXIS_OPP,
            s_samples, s_stages, s_classes, s_marks))
        if beautify:
            f.write(Qx0_text % s_equations)
        f.write(axes_text % (s_powers, s_groups))
        

########################################################################
########################################################################
# Import list of axes 
########################################################################
########################################################################



AXES = OrderedDict()

def get_axes_from_sample_axes(sample_axes):
    global AXES
    if len(AXES):
        return AXES
    for i, g1 in enumerate(sample_axes.g_strings):
        g = MM0(g1)
        axis = Axis(g)
        g_class = sample_axes.g_classes[i]
        axis.g_class = g_class
        axis.mark = sample_axes.g_marks[i]
        axis.auto_group = sample_axes.groups[i]
        axis.powers = sample_axes.powers[i]
        axis.stage = sample_axes.g_stages[i]
        axis.Qx0_equations = [np.array(a, dtype = np.uint32)
            for a in sample_axes.Qx0_equations[i]
        ]
        axis.constant = True
        AXES[g_class] = axis
    return AXES


def compute_and_write_axes(beautify = True, verbose = 0):
    from mmgroup.tests.axes.beautify_axes import Axis
    v0 = Axis()
    v0.stage = 0
    sample_list = explore_axes(v0, 5, spread, mark, score, 
            100, 50, verbose = verbose)
    write_axes(sample_list, beautify, verbose)
    time.sleep(0.1)


def do_test_sample_axes(sample_axes):
    sax =  sample_axes
    l = [G(sax.g_central), G(sax.g_axis), MMVector(7, sax.v_start)]
    lg =  [G(x) for x in sax.g_strings]
    assert len(lg) == NUM_AXIS_CLASSES


def get_sample_axes(calculate = False, beautify = True, verbose = 0):
    if calculate:
        compute_and_write_axes(beautify, verbose)
    try:
        from mmgroup.tests.axes import sample_axes
    except ImportError: 
        compute_and_write_axes(beautify, verbose)
        from mmgroup.tests.axes import sample_axes
    do_test_sample_axes(sample_axes)
    return get_axes_from_sample_axes(sample_axes)



########################################################################
########################################################################
# Main program 
########################################################################
########################################################################





if __name__ == "__main__":
    if PROCESSES != 1: 
        freeze_support()
    _ = get_sample_axes(calculate = False, beautify = True, verbose = True)

