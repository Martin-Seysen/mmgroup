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


from mmgroup import MM0, MMV, MMVector, Cocode, XLeech2, Parity, PLoop


########################################################################
########################################################################
# Start axis vectors and involutions in the monster group
########################################################################
########################################################################

# We give the standard axes and involutions as strings. So we can use
# it in any suitable constructor of the for the monster or its rep.

# Standard axes v^+ of 2A involution x_\beta, \beta = Cocode([2,3])
V_AXIS = "A_2_2 - A_3_2 + A_3_3 - 2*B_3_2" 
# Opposite axis v^- of of 2A involution x_{-1} * x_\beta  
V_AXIS_OPP = "A_2_2 - A_3_2 + A_3_3 + 2*B_3_2" 

# 2A involution x_\beta corresponding to standard axis v^+
G_AXIS = "d_200h"
# 2A involution x_{-1} x_\beta corresponding to opposite axis v^-
G_AXIS_OPP = "x_1000h * d_200h"

# Central involution in the subgroup G_x0 of the monster
G_CENTRAL = "x_1000h"

# Group element mapping v^+ to v^-
G_MAP_STD_OPP = "x_200h"

########################################################################
########################################################################
# The group and the vector space to be used
########################################################################
########################################################################

G = MM0                     # The monster group
V15 = MMV(15)               # Its representation space

PROCESSES = 0

def check_std_axes():
    """Some consistency checks for standard axes and involutions"""
    g_map_std_opp = G(G_MAP_STD_OPP)
    g_axis, g_axis_opp = G(G_AXIS), G(G_AXIS_OPP)
    assert g_axis ** g_map_std_opp == g_axis_opp
    v_axis, v_axis_opp = V15(V_AXIS), V15(V_AXIS_OPP)
    assert v_axis * g_map_std_opp == v_axis_opp
    
check_std_axes()


########################################################################
########################################################################
# The group elements and vectors and the vector space to be used
########################################################################
########################################################################

# The central involution in the subgroup ``G_x0``-
g_central = G(G_CENTRAL)  

# The standard 2A element in the monster group
g_axis = G(G_AXIS)

# The opposite standard 2A element in the monster group
g_axis_opp = G(G_AXIS_OPP)

# The standard 2A axis vector
v_axis = V15(V_AXIS)

# Opposite ofstandard 2A axis vector
v_axis_opp = V15(V_AXIS_OPP)




########################################################################
########################################################################
# Class for storing a pair (g, v), g in MM, v an axis
########################################################################
########################################################################



class GVector:
    r"""Models a 2A axis in the monster group

    The constructor takes three arguments ``opp``. If this is
    ``False`` (default) the instance is intialized with the
    standard "A axis ``v^+``. Otherwise it is initialized with
    the opposite axis ``v^-``. See |Seysen22| for the description
    of the 2A axes ``v^+``and ``v^-``.

    Attributes ``g0`` and ``v0``of an instance ""w"" of this class 
    contains that involution and its axis as an instance of class 
    ``mmgroup.MM0`` and ``mmgroup.MMVector``, respectively. Here 
    attribute ``v0`` is stored as a vector of integers modulo  15. 
    These two attributes are never changed.

    After construction, that instance ``w`` contains another attribute 
    ``g``, which is initialized to to the neutral element of the 
    monster group, and also an attribute ``v``, which is initialized
    to ``v0``.  Then instance ``w`` can be right multiplied with an
    arbitrary element ``g1`` of the monster group. If such a 
    multiplication is performed the both attributes, ``g`` and ``v0``
    are right multiplied with ``g1``. Note that all elements of the
    monster are implemented as instances of class ``mmgroup.MM0``.

    A consequence of that synchronous multiplication is that
    ``g**(-1) * g0 * g`` is a 2A involution which has axis ``v``
    for any instance of this class with attributes ``g``, ``g0``,
    and ``v``.
    """
    g_central = g_central
    __slots__ = "v", "g", "stage", "g0", "input"
    def __init__(self, opp = False):
        if opp is not None:
            g, v = self.get_gv_start(opp)
            self.g = G(1)
            self.g0 = G(g)
            self.v = V15(v)
            self.stage = 0
            self.input = opp

    @staticmethod 
    def get_gv_start(opp):
        return (G_AXIS_OPP, V_AXIS_OPP) if opp else (G_AXIS, V_AXIS)

    @property
    def v0(self):
        return V15(self.get_gv_start(self.input)[1])
        
    def __mul__(self, g1):
        g1 = G(g1)
        gv_new = self.__class__(None)
        gv_new.g = self.g * g1
        gv_new.v = self.v * g1
        gv_new.g0 = self.g0
        gv_new.stage = self.stage+1       
        gv_new.input = self.input       
        return gv_new

    def g_axis(self):
        r"""Return 2A involution corresponding to current axis

        The current 2A axis of an instance ``w`` of this class is
        given by ``w.v``. The function returns the 2A involution
        ``w.g**(-1) * w.g0 * w.g`` corresponding th that axis.
        """
        return self.g0 ** self.g

    def dihedral(self, axis = None):
        r"""Return product of current 2A and central involution

        For an instance ``w`` of this class let ``a`` be the
        2A involution ``w.g_axis()`` and let ``z`` be the central 
        involution in the subgroup :math:`G_{x0}`, which is a 2B
        involution. The function returns the product ``a * z``

        By |Nor98| there are 12 classes of pairs of involutions of 
        type 2A and 2B in the monster. The class of such e pair is 
        determined by the class of product ``a * z`` in the monster.
        Note that function ``axis_type`` can compute the class of
        ``a * z`` in the monster group.
        """
        axis = g_axis if axis is None else axis
        return self.g_axis() * self.g_central

    def mark(self):
        r"""Watermark the orbit of the current 2A axis

        For an instance ``w`` of this class let ``v = w.v`` be its
        current 2A axis vector. We want to watermark ``v`` in such
        a way that all 2A axis vectors in the same orbit under the
        group :math:`G_{x0}` obtain the same mark.

        Since :math:`G_{x0}` operates monomially on a 98280-dimensional
        subspace of the real vector space ``V`` containing the 2A
        axes, it is natural to count the absolute values of the 98280
        coordinates of ``V`` corresponding to that subspace.

        Since coordinates in ``V`` are stored modulo 15, these
        absolute values correspond to the integers 0,...,7. We return
        an 8-tupel ``t``, where ``t[i]`` contains the number of
        coordinates in that subspace with absolute value ``i``. 

        It turns out that this kind of watermarking allows us to
        distinguish between all 12 orbits of 2A axes. The orbits
        of 2A axes under :math:`G_{x0}` are described in |Nor98|.      
        """
        return self.v.count_short()

    def baby_value_A(self):
        r"""Yet to be documented!!!"""
        return self.v.eval_A(XLeech2(Cocode([2,3])))
    

    def count_zeros_in_A(self):
        r"""Count number of zeros in part 'A' of the sxis

        For an instance ``w`` of this class let ``v = w.v`` be its
        current 2A axis vector. The function returns the number of
        zero entries in the 'A' part of vector 'v', which can be
        considered as a symmetric 24 times 24 matrix.
        """
        return 576 - np.count_nonzero(self.v['A'])


    def score_A(self):
        A = self.v['A']
        score = num_zeros = 576 - np.count_nonzero(self.v['A'])
        DIAG = np.diagonal(A)
        diag_bins = sorted([x for x in np.bincount(DIAG) if x])
        axis_type = self.axis_type()
        # We shall use the axis type for finding the 'nicest' axis
        # of a given type only!
        #print(self.axis_type())
        bonus = 0
        if axis_type == "10A" and max(diag_bins) >= 22:
            # Special treatment for case 10A: Prefer solution with
            # large No of off-diagonal elements of same absolute value
            d = defaultdict(list)
            for i, x in enumerate(DIAG):
                d[x].append(i) 
            for x in d.values():
                if len(x) >= 22:
                    I = x
            B0 = A[I,:][:,I]
            B = np.where(B0 > 7, 15 - B0 , B0)
            max_occur = max(np.bincount(B.ravel()))
            if max_occur >= 462:
                return 1000 + max_occur
        elif axis_type == "10B":
            bonus = diag_bins == [4,20]
            # Special treatment for case 10B: Prefer solution with
            # block [4, 24]
            #print("YEEEAHHH", diag_bins)
        elif axis_type == "12C":
            #print(diag_bins)
            bonus =   1000 * (diag_bins == [6, 8, 10])          
        return score + bonus
        

    def axis_type(self):
        return self.v.axis_type(e = 0)

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
        processes=0, # No of processes
        verbose=0
        ):
    """Multitasking version of function ``next_generation``

    Input parameters and return value are as in function
    ``next_generation``.

    The additional parameter ``processes`` specifies the 
    maximum number of processes to be started. ``processes = 0``
    (default) means to start as many processes as possible.

    Not yet functional. Use function ``next_generation`` instead!
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
# Compute axis type
########################################################################
########################################################################


def axis_char(g0, axis = None):
    """Compute type of a 2A axis

    Let ``axis`` be a  2A-element  of the monster group, let ``g0`` 
    be an arbitrary element of the monster group. Put
    ``g = g0**(-1) * axis * g``. Let ``x`` be the central involution
    of the subgroup :math:`G_{x0}` of the monster. ``axis``
    defaults to the element ``MM0("d", Cocode([2,3]))``
    
    We want to classify the pair ``(g, x)`` of involutions. Put
    ``h = g*x``. Then ``h`` is of even order ``o``. We try to
    compute the character ``chi`` of ``h`` in the 196883-dimensional
    representation of the monster. This can be done with our
    capabilities if ``h**(o/2)`` is a '2B' involution; otherwise
    we put ``chi = None``.

    We return the pair ``(g, o)``. All computations are done in
    the class ``mmgroup.MM0`` modelling the monster group, since
    the class ``mmgroup.MM`` modelling the monster group (provided 
    for the user) uses the results of such computations.

    It turns out that the result of this function contains
    sufficient information for obtaining the class of ``h``
    in the monster group as given in |Nor98|.
    """
    g0 = MM0(g0)
    axis = g_axis if axis is None else MM0(axis)
    g = g0**(-1) * axis * g0 * g_central
    o, h = g.half_order()
    if h.chi_G_x0()[0] != 275:
        return o, None  # The prod is a 2A involution
    # Here prod is a 2B involution
    itype, m = h.conjugate_involution()
    assert itype == 2
    chi = (g**m).chi_G_x0()[0]
    return o, chi
    


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


AXES_CLASSES_DICT = {}
for number, (letters, characters) in AXES_PROD_CLASSES.items():
    for letter, character in zip(letters, characters):
        AXES_CLASSES_DICT[(number, character)] = str(number) + letter


def axis_type(g0, axis = None):
    """Compute type of a 2A axis

    Let ``axis`` be a  2A-element  of the monster group, let ``g0`` 
    be an arbitrary element of the monster group. Put
    ``g = g0**(-1) * axis * g``. Let ``x`` be the central involution
    of the subgroup :math:`G_{x0}` of the monster. ``axis``
    defaults to the element ``MM0("d", Cocode([2,3]))``
    
    We return class of ``h = g*x`` in ATLAS notation as a string.
    """
    return AXES_CLASSES_DICT[axis_char(g0, axis)]



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
    return gv * g



########################################################################
# Profiling an axis: count entries for tags B, C, T, X
########################################################################



def mark(gv):
    return gv.mark()


def score(gv):
    return gv.score_A()
    #return gv.count_zeros_in_A()


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
        if len(sample_list) >= 12: 
            break
    num_samples = sum(len(x) for x in all_samples.values())
    print("Number of axes considered:", num_samples)
    assert len(sample_list) >= 12
    return sample_list



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

g_beautifiers = [
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


def sample_list_sort_key(sample_entry):
     stage, sample, _ = sample_entry 
     s_class = axis_type(sample.g) 
     return stage, int(s_class[:-1]), s_class[-1:]


def write_axes(sample_list, verbose = False):
    from mmgroup.tests.test_axes.beautify_axes import compute_beautifiers
    sample_list.sort(key = sample_list_sort_key)
    s_samples, s_stages, s_marks = "", "", "" 
    s_classes = ""
    s_beautifiers = ""
    s_powers = ""
    s_groups = ""
    g_strings = [x[1].g for x in sample_list]
    g_classes = [axis_type(x[1].g) for x in sample_list]
    beautifiers = compute_beautifiers(g_strings, g_classes)
    for i, (stage, sample, mark) in enumerate(sample_list):
        s_samples += "\"" + sample.g.raw_str() + "\",\n"
        s_stages += str(stage) + ", "
        s_marks += str(mark)  + ",\n"
        class_ = axis_type(sample.g)
        s_classes += '"' + class_ + '", '
        s_beautifiers += '"' + beautifiers[i] + '",\n'
        s_powers += '"' + AXIS_POWERS[class_] + '",\n'
        s_groups += '"' + AXIS_GROUPS[class_] + '",\n'
       
    with open(PATH + ".py", "wt") as f:
        f.write(f_text % (G_CENTRAL, G_AXIS, V_AXIS, 
            G_AXIS_OPP, V_AXIS_OPP,
          s_samples, s_stages, s_classes, s_marks, s_beautifiers))

        f.write(axes_text % (s_powers, s_groups)) 
        

########################################################################
########################################################################
# Import list of axes 
########################################################################
########################################################################


def compute_and_write_axes(verbose = 0):
    v0 = GVector()
    sample_list = explore_axes(v0, 5, spread, mark, score, 
            100, 50, verbose = verbose)
    write_axes(sample_list, verbose)
    time.sleep(0.1)


def do_test_sample_axes(sample_axes):
    sax =  sample_axes
    l = [G(sax.g_central), G(sax.g_axis), MMVector(7, sax.v_start)]
    lg =  [G(x) for x in sax.g_strings]
    assert len(lg) == 12


def import_sample_axes(calculate = False, verbose = 0):
    if calculate:
        compute_and_write_axes(verbose)
    try:
        from mmgroup.tests.test_axes import sample_axes
    except ImportError: 
        compute_and_write_axes(verbose)
        from mmgroup.tests.test_axes import sample_axes
    do_test_sample_axes(sample_axes)
    return sample_axes



########################################################################
########################################################################
# Main program 
########################################################################
########################################################################





if __name__ == "__main__":
    if PROCESSES != 1: 
        freeze_support()
    import_sample_axes(calculate = True, verbose = True)
