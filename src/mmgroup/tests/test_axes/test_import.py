
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os

import numpy as np

import pytest


if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import MM0, MMV
from mmgroup.mm_crt_space import MMVectorCRT 


from mmgroup.tests.test_axes.get_nice_sample_axes import import_sample_axes
from mmgroup.tests.test_axes.get_sample_axes import do_test_sample_axes
from mmgroup.tests.test_axes.get_baby_sample_axes import import_baby_sample_axes
from mmgroup.tests.test_axes.get_baby_sample_axes import do_test_baby_sample_axes
from mmgroup.tests.test_axes.beautify_axes import compute_beautifiers, beautify
from mmgroup.tests.test_axes.beautify_axes import adjacent_blocks
from mmgroup.tests.test_axes.beautify_axes import block_eigenvalues
from mmgroup.generators import gen_leech3to2_short, gen_leech3to2_type4
from mmgroup.mm15 import op_eval_A as mm_op15_eval_A
from mmgroup.mm15 import op_norm_A as mm_op15_norm_A
from mmgroup.mm15 import op_eval_A_rank_mod3 as mm_op15_eval_A_rank_mod3


baby_sample_axes = import_baby_sample_axes()
sample_axes = import_sample_axes()

AXES = dict(zip(sample_axes.g_classes, sample_axes.g_strings))
BABY_AXES = dict(zip(baby_sample_axes.g_classes, baby_sample_axes.g_strings))


def display_A(A):
   """Display a 24 times 24 matrix A"""   
   fmt = [4 if max(abs(A[i])) > 99 else 2 for i in range(24)] 
   for i in range(24):
      print(" ", end = "")
      for j in range(24):
          l = fmt[j]
          if i == j or A[i,j] != 0:
              print("%*d" % (fmt[j], A[i,j]), end = " ")
          else:
              print("%*s" % (fmt[j], "."), end = " ")
      print("")


norms_A_mod15 = None

def axis_index(axis_type):
    sample_axes = import_sample_axes()
    NUM_AXES = len(sample_axes.g_strings)
    for i, s in enumerate(sample_axes.g_classes):
        if s == axis_type:
            return i
    err = "No axis of type %s found"
    raise ValueError(err % str(axis_type))


def norm_A_mod15(i):
    """Return norm(A) mod 15 for entry i in the list in sample_axes.py

    The lists in file sample_axes.py refer to a collection of 12-axes 
    of different types.
   
    """
    global norms_A_mod15
    if isinstance(i, str): i = axis_index(i)
    try:
        return norms_A_mod15[i]
    except:
        norms_A_mod15 = []
        sample_axes = import_sample_axes()
        V15 = MMV(15)
        v_start = V15(sample_axes.v_start)
        for g in sample_axes.g_strings:
             v = v_start * MM0(g)
             norms_A_mod15.append(mm_op15_norm_A(v.data))
        return norms_A_mod15[i]




def display_norm_A(i):
    """Display additional information of a 2A axis

    
    """
    sample_axes = import_sample_axes()
    
    if isinstance(i, str): i = axis_index(i)
    norm_A = norm_A_mod15(i)
    s = "norm(A) = %d (mod 15)" % (norm_A)
    common = [s for j, s in enumerate(sample_axes.g_classes) 
                if  norm_A_mod15(j) == norm_A and j != i]
    if len(common) == 0:
        return s + "\n"
    s += ", same norm for class " + " and ".join(map(str,common))
    d = 2 if norm_A == 4 else 0
    V15 = MMV(15)
    v = V15(sample_axes.v_start) * MM0(sample_axes.g_strings[i])
    r = mm_op15_eval_A_rank_mod3(v.data, d)    
    rank = r >> 48
    v3 = r & 0xffffffffffff
    s_A = "(A - %d * 1)" % d if d else "A"
    s += "\nMatrix U = %s (mod 3) has rank %d" % (s_A, rank)
    if (rank != 23):
         return s + "\n"
    v2 = gen_leech3to2_short(v3)
    if v2:
         f = "Kernel of U is spanned by a short vector v with A(v) ="
         a2 = mm_op15_eval_A(v.data, v2)
         s += "\n%s %d (mod 15)" % (f, a2)
    if gen_leech3to2_type4(v3):
         f = "Kernel of U is spanned by a vector of type 4"
         s += "\n" + f
    return s + "\n"
    





HEADER = r"""Classes of 2A axes

We display information about the classes of 2A axes of the monster 
group. Each class corresponds to an orbit of 2A axes under the group 
G = 2^{1+24}.Co_1. A 2A axis corresponds to a 2A involution in the
monster group. The class of a 2A involution is named by the name of
the class of z*g in the monster group. Here z is the central 
involution of the group G_x0. Class names of the monster group are 
given in ATLAS notation.

A 2A axis v is a vector in the 196884-dimensional representation of 
the monster group. The 'A' part of v is the projection into a subpace
of that representation that corresponds to the 299-dimensional 
represententation of G_x0. That 'A' part has a natural interpretation
as a symmetric 24 times 24 matrix operating on Leech lattice.

For each class of 2A axes we give a sample matrix A operating on
the Leech lattice in standard cooordinates. The matrices displayed
below and also their eigenvalues should be divided by 256.

It turns out that for distinguishing between the classes of 2A axes
it suffices to know the the corresponding matrix A modulo 15. We use 
the norm (i.e. the sum of the squared entries) of the matrix A for
distinguishing between classes. In cases where this is not sufficient
we also use the rank of matrix U = (A - d*1) modulo 3. Here d is an 
integer depending on the norm of A, and '1' is the unit matrix.
If the kernel of U has dimension 1 and is spanned by a short Leech
lattice vector then we may also use A(v) modulo 15 for distinguishing
between classes.

""" 

@pytest.mark.axes 
@pytest.mark.slow 
def test_2A_axes_classes(verbose = 0):
    """Test computation of samples classes of 2A axes.

    If this function is called with parameter ``verbose == 1``
    then function displays the information relevant for 
    distinguishing between classes of 2A axes.
    """
    do_test_baby_sample_axes(baby_sample_axes)
    do_test_sample_axes(sample_axes)
    if verbose:
        print(HEADER)

    NUM_AXES = len(sample_axes.g_strings)
    for i in range(NUM_AXES): # sorted_axes(sample_axes):
        g = MM0(sample_axes.g_strings[i])
        v =  MMVectorCRT(20, sample_axes.v_start)
        v *= g
        g1 = MM0(sample_axes.g_beautifiers[i])
        v *= g1
        Afloat = 256 * v["A"]
        A = np.array(Afloat, dtype = np.int32)
        assert (A == Afloat).all()
        if verbose:
            class_ = sample_axes.g_classes[i]
            print("\nClass " + class_, end = ", ")
            print("stage = " + str(sample_axes.g_stages[i]), end = "")
            powers = sample_axes.powers[i]
            s_powers = ", powers: " + powers if len(powers) else ""
            print(s_powers)
            print("Automorphism group:", sample_axes.groups[i])
            print(display_norm_A(i), end = "")
            print("Eigenvalues", block_eigenvalues(A))
            display_A(A)

             
if  __name__ == "__main__":
    test_2A_axes_classes(verbose = 1)





