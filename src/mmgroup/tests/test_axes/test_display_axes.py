import sys
import time
from collections import OrderedDict, defaultdict
from random import sample

import numpy as np
import pytest

if __name__ == "__main__":
    sys.path.append("../../../")

from mmgroup import MM, MM0, AutPL, Cocode, GCode, MMSpace, MMV, GcVector
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import bitmatrix64_solve_equation
from mmgroup import mat24
from mmgroup.generators import gen_leech3to2_short, gen_leech3to2_type4
from mmgroup.mm_op import mm_op_eval_A
from mmgroup.mm_op import mm_op_norm_A
from mmgroup.mm_op import mm_op_eval_A_rank_mod3



from mmgroup.tests.test_axes.test_import import display_A
from mmgroup.tests.test_axes.test_import import display_norm_A
from mmgroup.tests.test_axes.test_import import format_eigen_values
from mmgroup.tests.axes.get_baby_sample_axes import get_baby_sample_axes
from mmgroup.tests.axes.get_sample_axes import get_sample_axes



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




#######################################################################
# Display central involtion of dihedral group of axis
#######################################################################


def x_equations(axis, t = 0):
    #print(".", end = "")
    v = axis.v15 * MM0('t', t)
    data = v["E", 300:300+98280]
    short = [MMSpace.index_to_short_mod2(i+300) 
        for i, x in enumerate(data) if x != 0]
    matrix= np.zeros(24, dtype = np.uint64)
    rows = 0
    for sh in short:
        new_rows = leech2matrix_add_eqn(matrix, rows, 24, sh)
        rows += new_rows
    #print(".", end = "")
    return matrix[:rows]

def gc(i):
    g = GCode(i & 0xfff)
    if len(g.bit_list) > 12: g = ~g
    return g.bit_list

def cocode(i, gc_ref):
    syn0 = 0
    if gc_ref:
       gc_list =  gc(gc_ref)
       if len(gc_list): syn0 = gc_list[0]
    c = Cocode(i)
    return c.syndrome_list(syn0)    

def analyze_xy(g): 
    t = g.as_tuples()
    pi_ref = None
    for tag, i in t:
       if tag in  "xy":
           print(" ", tag, gc(i)) 
           pi_ref = i
       if tag in  "d":
           print(" ", "d", cocode(i, pi_ref))
       if tag in  "p":
           print(" ", "p", AutPL(0,i).perm)





#######################################################################
# Diplay information for axes
#######################################################################


@pytest.mark.axes 
@pytest.mark.slow 
def test_display_axes(verbose = 0):
    if verbose: 
        print(HEADER)
    else:
        print("\nClasses of 2A axes (relative to G_x0):\n")
    for cl, axis in  get_sample_axes().items():
        print("Class:", axis.g_class, ", stage =", axis.stage, 
                ", powers:", axis.powers)
        print("Automorphism group:", axis.auto_group)
        print(display_norm_A(axis), end = "")
        if verbose:
            print("Eigenvalues of 256 * A part of axis v:",
                format_eigen_values(axis.in_space(MMVectorCRT, 20)['A'] * 256))
        if verbose:
            z_i = axis.central_involution()
            orb, c = z_i.conjugate_involution_G_x0()
            print("central involution: ", z_i)
            if c: 
                print("Character of central involution z_i:",  z_i.chi_G_x0(),
                ", class:", orb)
                #print("Product of z_i and standard involution", z_i * z_i**c)
            print("Axis type of axis v * MM('t',2):", axis.axis_type(2))
            #print("Dim intersection with Q_x0:", 24 - len(x_equations(axis,0)))
            #print("Dim intersection with Q_y0:", 24 - len(x_equations(axis,1)))
            #print("Dim intersection with Q_z0:", 24 - len(x_equations(axis,2)))
            #opp_axis_type = axis.axis_type()
            #print("Image of opposite axis: ", opp_axis_type) 
            print("A part of axis v:")
            axis.display_sym('A')
            axis.display_sym(2)
        print("")



if  __name__ == "__main__":
    test_display_axes(verbose = 1)


