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

from mmgroup.tests.test_axes.beautify_axes import compute_beautifiers

from mmgroup.tests.test_axes.get_sample_axes import G, V15 
from mmgroup.tests.test_axes.get_sample_axes import V_AXIS_OPP 
from mmgroup.tests.test_axes.get_sample_axes import GVector
from mmgroup.tests.test_axes.get_sample_axes import next_generation_pool
from mmgroup.tests.test_axes.get_sample_axes import axis_type
from mmgroup.tests.test_axes.get_sample_axes import _show_sample
from mmgroup.tests.test_axes.get_sample_axes import explore_axes
from mmgroup.tests.test_axes.get_sample_axes import write_axes
from mmgroup.tests.test_axes.get_sample_axes import do_test_sample_axes

from mmgroup.tests.test_axes.get_baby_sample_axes import rand_Co_2




########################################################################
########################################################################
# 
########################################################################
########################################################################


PROCESSES = 0



########################################################################
########################################################################
# The group elements and vectors and the vector space to be used
########################################################################
########################################################################

# The central involution in the subgroup ``G_x0``-
g_central = G("x", 0x1000)  

# The standard 2A element in the monste froup
g_axis = G("d", Cocode([2,3]))

# Tuple describing the standard 2A axis vector
v_start_tuple = ("I", 3, 2) 

# The standars 2A axis vector
v_axis = V15(*v_start_tuple)


########################################################################
########################################################################
# Class for storing a pair (g, v), g in MM, v an axis
########################################################################
########################################################################

TYPE_SCORES = {
   '2A': 11, '2B': 10, '4A': 9, '4B': 38, '4C': 7, '6A': 6, '6C': 5,
   '8B': 4, '10A': 3, '6F': 2, '10B': 1, '12C': 0
}

class GVector1(GVector):
    r"""Models a 2A axis in the monster group

    Creating an instance of this class is equivalent to creating
    an instance of class ``GVector`` with default parameters.

    Here in contrast to class ``GVector`` we also keep tack to 
    the image of the oppsite axis ``v^-`` of the standard axis
    ``v^+``.
    """
    v_opp = None  
    def __init__(self, opp = False):
        assert opp in (False, None)
        super(GVector1, self).__init__(False)
        self.v_opp = V15(V_AXIS_OPP)

    def __mul__(self, g1):
        gv_new = super(GVector1, self).__mul__(g1) 
        gv_new.v_opp = self.v_opp * g1    
        return gv_new

    def score_A(self):
        score = super(GVector1, self).score_A() 
        if self.v_opp:
            opp_type = self.v_opp.axis_type()
            opp_score = TYPE_SCORES[self.v_opp.axis_type()]  
            #if opp_type == "4A" and self.v.axis_type() == "8B":
            #     opp_score += 2
            score += 1000000 * opp_score      
        return score

    def reduce_g(self):
        from mmgroup import MM
        g = MM(self.g).reduce()
        self.g = G(g)
        





########################################################################
########################################################################
# Obtaining samples of transformed axes (there are 12 classes)
########################################################################
########################################################################


########################################################################
# Next generation of axes
########################################################################


def calc_conjugate_to_opp(i0, i1):
    coc = Cocode([i0,i1]) 
    g = G('d', coc)
    v = V15('I', i0, i1)
    coc_value = coc.ord
    x = PLoop(coc_value & -coc_value)
    return G('x', x)

Co_2_OPP = calc_conjugate_to_opp(2, 3)
   


    

def spread(gv):
    ax_type = gv.axis_type()
    if ax_type in ['6C', '8B', '6F', '10A', '10B', '12C']:
         return None
    if gv.stage < 3 :
        if ax_type == "4A" and randint(0,99) < 20:
            gv_new = gv * G([('r','G_x0'), ('t','n')])
            if gv_new.axis_type() == "8B":
                return gv_new
        g = (rand_Co_2() * G('t','n')) ** Co_2_OPP 
    else:
        g = G([('r','G_x0'), ('t','n')])
    return gv * g



########################################################################
# Profiling an axis: count entries for tags B, C, T, X
########################################################################



XLeech2_Opp = -XLeech2(Cocode([2,3]))

def mark(gv):
    v = gv.v
    ax_type = v.axis_type()
    value = v.eval_A(XLeech2_Opp) if ax_type in ['2A', '2B'] else 0
    return ax_type, value


def score(gv):
    return gv.score_A()
    #return gv.count_zeros_in_A()




    

########################################################################
########################################################################
# Import list of axes 
########################################################################
########################################################################



def filter_sample_list(sample_list):
    out_list = []
    axis_types = set()
    for sample in sample_list:
        gv = sample[1]
        gv.reduce_g()
        ax_type = gv.axis_type()
        if not ax_type in axis_types:
            axis_types.add(ax_type)
            out_list.append(sample)
    return out_list 
        


def compute_and_write_axes(verbose = 0):
    v0 = GVector1()
    sample_list = explore_axes(v0, 5, spread, mark, score, 
            100, 50, verbose = verbose)
    sample_list = filter_sample_list(sample_list)
    write_axes(sample_list, verbose)
    time.sleep(0.1)



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
    sample_axes = import_sample_axes(calculate = True, verbose = True)
