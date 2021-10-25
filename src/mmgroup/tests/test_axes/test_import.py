
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os

import numpy as np

import pytest

from mmgroup import MM0, MMV
from mmgroup.mm_crt_space import MMVectorCRT 


from mmgroup.tests.test_axes.get_sample_axes import import_sample_axes
from mmgroup.tests.test_axes.get_sample_axes import do_test_sample_axes
from mmgroup.tests.test_axes.get_baby_sample_axes import import_baby_sample_axes
from mmgroup.tests.test_axes.get_baby_sample_axes import do_test_baby_sample_axes
from mmgroup.tests.test_axes.beautify_axes import compute_beautifiers
from mmgroup.tests.test_axes.beautify_axes import adjacent_blocks
from mmgroup.tests.test_axes.beautify_axes import block_eigenvalues
from mmgroup.clifford12 import leech_matrix_norm_A, leech3matrix_kernel
from mmgroup.clifford12 import leech2matrix_eval_A
from mmgroup.generators import gen_leech3to2_short, gen_leech3to2_type4

NTESTS = 100


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

def norm_A_mod15(i):
    """Return norm(A) mod 15 for entry i in the list in sample_axes.py

    The lists in file sample_axes.py refer to a collection of 12-axes 
    of different types.
   
    """
    global norms_A_mod15
    try:
        return norms_A_mod15[i]
    except:
        norms_A_mod15 = []
        sample_axes = import_sample_axes()
        V15 = MMV(15)
        v_start = V15(sample_axes.v_start)
        for g in sample_axes.g_strings:
             v = v_start * MM0(g)
             norms_A_mod15.append(leech_matrix_norm_A(15, v.data))
        return norms_A_mod15[i]




def display_norm_A(i):
    """Display additional infomration of a 2A axis

    
    """
    sample_axes = import_sample_axes()
    
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
    r = leech3matrix_kernel(15, v.data, d)    
    rank = r >> 48
    v3 = r & 0xffffffffffff
    s_A = "(A - %d * 1)" % d if d else "A"
    s += "\nMatrix U = %s (mod 3) has rank %d" % (s_A, rank)
    if (rank != 23):
         return s + "\n"
    v2 = gen_leech3to2_short(v3)
    if v2:
         f = "Kernel of U is spanned by a short vector v with A(v) ="
         a2 = leech2matrix_eval_A(15, v.data, v2)
         s += "\n%s %d (mod 15)" % (f, a2)
    if gen_leech3to2_type4(v3):
         f = "Kernel of U is spanned by a vector of type 4"
         s += "\n" + f
    return s + "\n"
    



def sorted_axes(sample_axes):
    """Return sorted list of axes"""
    NUM_AXES = len(sample_axes.g_strings)
    def axes_order(i):
        cl = (sample_axes.g_classes)[i]
        st = (sample_axes.g_stages)[i]
        return st, int(cl[:-1]), cl, i 
    axes_data = [axes_order(i) for i in range(NUM_AXES)] 
    axes_data.sort()
    sorted_axes_list = [x[-1] for x in axes_data]
    assert set(sorted_axes_list) == set(range(NUM_AXES))
    return sorted_axes_list 



HEADER = r"""Classes of 2A axes

We display information about the classes of 2A axes of the monster 
group. Each class corresponds to an orbit of 2A axes under the group 
G = 2^{1+24}.Co_1. A 2A axis corresponds to a 2A involution in the
monster group. The class of a 2A involution is named by the name of
the class of z*g in the monster group. Here z is the central 
involution of the group G. Class names of the monster group are 
given in ATLAS notation.

A 2A axis v is a vector in the 196884-dimensional representation of 
the monster group. The 'A' part of v is the projection into a subpace
of that representation that corresponds to the 299-dimensional 
represententation of G. That 'A' part has a natural interpretation
as a symmetric 24 times 24 matrix operating on Leech lattice.

For each class of 2A axes we give a sample matrix A operating on
the Leech lattice in standard cooordinates. The matrices displayed
below and also their eigenvalues should be divided by 128.

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
def test_2A_axes_classes(verbose = 1):
    """Test computation of samples classes of 2A axes.

    If this function is called with parameter ``verbose == 1``
    then function displays the information relevant for 
    distinguishing between classe of 2A axes.
    """
    baby_sample_axes = import_baby_sample_axes(verbose = verbose)
    do_test_baby_sample_axes(baby_sample_axes)
    sample_axes = import_sample_axes(verbose = verbose)
    do_test_sample_axes(sample_axes)
    beautfiers = compute_beautifiers(sample_axes.g_strings)
    if verbose:
        print(HEADER)

    for i in sorted_axes(sample_axes):
        g = MM0(sample_axes.g_strings[i])
        v =  MMVectorCRT(20, sample_axes.v_start)
        v *= g
        g1 = MM0(beautfiers[i])
        v *= g1
        Afloat = 128 * v["A"]
        A = np.array(Afloat, dtype = np.int32)
        assert (A == Afloat).all()
        if verbose:
            print("\nClass " + sample_axes.g_classes[i], end = ", ")
            print("stage = " + str(sample_axes.g_stages[i]))
            print(display_norm_A(i), end = "")
            print("Eigenvalues", block_eigenvalues(A))
            display_A(A)

             




