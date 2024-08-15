
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os

import numpy as np

import pytest


if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import MM0, MMV
from mmgroup.generators import gen_leech3to2_short, gen_leech3to2_type4
from mmgroup.mm_op import mm_op_eval_A
from mmgroup.mm_op import mm_op_norm_A
from mmgroup.mm_op import mm_op_eval_A_rank_mod3
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup.tests.axes.get_baby_sample_axes import get_baby_sample_axes
from mmgroup.tests.axes.get_sample_axes import get_sample_axes




baby_sample_axes = get_baby_sample_axes()
sample_axes = get_sample_axes()


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




def display_norm_A(axis):
    """Display additional information of a 2A axis

    
    """
    sample_axes = get_sample_axes()
    
    norm_A = axis.norm_A_mod15
    g_class = axis.g_class
    s = "norm(A) = %d (mod 15)" % (norm_A)
    common = [other.g_class for other in sample_axes.values() 
                if  other.norm_A_mod15 == norm_A and 
                other.g_class != g_class]
    if len(common) == 0:
        return s + "\n"
    s += ", same norm for class " + " and ".join(map(str,common))
    d = 2 if norm_A == 4 else 0
    v = axis.v15
    r = mm_op_eval_A_rank_mod3(15, v.data, d)    
    rank = r >> 48
    v3 = r & 0xffffffffffff
    s_A = "(A - %d * 1)" % d if d else "A"
    s += "\nMatrix U = %s (mod 3) has rank %d" % (s_A, rank)
    if (rank != 23):
         return s + "\n"
    v2 = gen_leech3to2_short(v3)
    if v2:
         f = "Kernel of U is spanned by a short vector v with A(v) ="
         a2 = mm_op_eval_A(15, v.data, v2)
         s += "\n%s %d (mod 15)" % (f, a2)
    if gen_leech3to2_type4(v3):
         f = "Kernel of U is spanned by a vector of type 4"
         s += "\n" + f
    return s + "\n"
    

#################################################################
# Eigenvalues of matrix A
#################################################################


def purge_diag(diag, power = 1):
    """Convert list of real numbers to a dictionary

    Given a list 'diag' of rea numbers, the function returns a
    dictionary mapping each number occuring in that list to it
    multiplicity. The function makes reasonable assumptions
    regarding equality of approximately equal entries, so that
    we obtain a dictionary from a list of eigenvalues of
    a 24 times 24 matrix acting on the Leech lattice.
    """
    EPS = 1.0e-8
    data = {}
    non_ints = []
    for x0 in diag:
        x = float((x0**power).real) 
        if abs(x - round(x)) < EPS:
            i = int(round(x))
            if i in data:
                data[i] += 1
            else:
                data[i] = 1
        else: 
            done = False           
            for d in non_ints:
                if abs(d - x) < EPS:
                    data[d] += 1
                    done = True
            if not done:
                data[x] = 1
                non_ints.append(x)
    return data


def format_eigen_values(A):
   """Return description of eigenvalues of a matrix ``A``

   The description of the eigenvalues is returned as a string.
   """
   eigen =  purge_diag(np.linalg.eigvals(A))
   eigenv = sorted(eigen.items(), key = lambda x: x[0])
   data = [("%d^%d" if isinstance(k, int) else "%.2f^%d") % 
       (k,v) for k, v in eigenv ]
   return "(" + ", ".join(data) + ")"
 


#################################################################
# Main test
#################################################################



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

@pytest.mark.slow 
def test_2A_axes_classes(verbose = 0):
    """Test computation of samples classes of 2A axes.

    If this function is called with parameter ``verbose == 1``
    then function displays the information relevant for 
    distinguishing between classes of 2A axes.
    """
    if verbose:
        print(HEADER)

    assert len(sample_axes) == 12
    for orbit, axis in sample_axes.items():
        g = axis.g
        v =  axis.in_space(MMVectorCRT, 20)
        #g1 = MM0(sample_axes.g_beautifiers[i])
        #v *= g1
        Afloat = 256 * v["A"]
        A = np.array(Afloat, dtype = np.int32)
        assert (A == Afloat).all()
        if verbose:
            class_ = axis.g_class
            print("\nClass " + class_, end = ", ")
            print("stage = " + str(axis.stage), end = "")
            powers = axis.powers
            s_powers = ", powers: " + powers if len(powers) else ""
            print(s_powers)
            print("Automorphism group:", axis.auto_group)
            print(display_norm_A(axis), end = "")
            print("Eigenvalues", format_eigen_values(A))
            display_A(A)

             
if  __name__ == "__main__":
    test_2A_axes_classes(verbose = 1)





