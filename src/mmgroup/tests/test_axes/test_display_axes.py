import sys
import time
from collections import OrderedDict, defaultdict
from random import sample

import numpy as np
import pytest

if __name__ == "__main__":
    sys.path.append("../../../")

from mmgroup.tests.test_axes import sample_axes 
from mmgroup.tests.test_axes.beautify_axes import Axis
from mmgroup.tests.test_axes.beautify_axes import make_blocks_positive
from mmgroup.tests.test_axes.beautify_axes import change_signs_A
from mmgroup import MM0, AutPL, Cocode, GCode, MMSpace, MMV, GcVector
from mmgroup.mm_crt_space import MMVectorCRT 
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import bitmatrix64_solve_equation
from mmgroup import mat24
from mmgroup.generators import gen_leech3to2_short, gen_leech3to2_type4
from mmgroup.mm15 import op_eval_A as mm_op15_eval_A
from mmgroup.mm15 import op_norm_A as mm_op15_norm_A
from mmgroup.mm15 import op_eval_A_rank_mod3 as mm_op15_eval_A_rank_mod3

from mmgroup.tests.test_axes.get_sample_axes import G, V15
from mmgroup.tests.test_axes.get_sample_axes import g_central
from mmgroup.tests.test_axes.get_sample_axes import V_AXIS
from mmgroup.tests.test_axes.get_sample_axes import V_AXIS_OPP
from mmgroup.tests.test_axes.get_sample_axes import g_axis, g_axis_opp

v_axis = MMVectorCRT(16, V_AXIS)
v_axis_opp = MMVectorCRT(16, V_AXIS_OPP)
v_axis15 = V15(V_AXIS)
v_axis_opp15 = V15(V_AXIS_OPP)


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
# Display an integer matrix A
#######################################################################


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


#######################################################################
# Display central involtion of dihedral group of axis
#######################################################################


def x_equations(axis, t = 0):
    #print(".", end = "")
    v = v_axis15 * axis.g_all() * MM0('t', t)
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
# Display eigenvalues of matrix A
#######################################################################



def purge_diag(diag, power = 1):
    EPS = 1.0e-8
    data = {}
    non_ints = []
    for x0 in diag:
        x = (x0**power).real
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

def eigen(A):
    eigen = purge_diag(np.linalg.eigvals(A))
    values = sorted(eigen.keys())
    l = []
    for v in values:
        s = str(v) if v == int(v) else "%.3f" % v
        l.append(s + '^' + str(eigen[v]))
    return ", ".join(l)



#######################################################################
# Read axes from file sample_axes.py
#######################################################################



AXES= {}

def get_axes():
    global AXES
    if len(AXES):
        return AXES
    for i, g1 in enumerate(sample_axes.g_strings):
        g2 = sample_axes.g_beautifiers[i]
        g = MM0(g1) * MM0(g2)
        g_class = sample_axes.g_classes[i]
        axis =  Axis(g, g_class)
        axis.mark = sample_axes.g_marks[i]
        axis.group = sample_axes.groups[i]
        axis.powers = sample_axes.powers[i]
        axis.stage = sample_axes.g_stages[i]
        axis.v15 =  v_axis15 * MM0(sample_axes.g_strings[i])
        axis.norm15 = mm_op15_norm_A(axis.v15.data)
        AXES[g_class] = axis
    return AXES


#######################################################################
# Display information of axes mod 15
#######################################################################

#norms_A_mod15 = [axis.norm15 for axis in AXES]

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
        for g in sample_axes.g_strings:
             v = v_axis * MM0(g)
             norms_A_mod15.append(mm_op15_norm_A(v.data))
        return norms_A_mod15[i]




def display_norm_A(axis):
    """Display additional information of a 2A axis

    
    """
    #sample_axes = import_sample_axes()
    
    norm_A = axis.norm15
    s = "norm(A) = %d (mod 15)" % (norm_A)
    common = [a.g_class for a in get_axes().values()
       if  a.norm15 == norm_A and a.g_class != axis.g_class]
    if len(common) == 0:
        return s + "\n"
    s += ", same norm for class " + " and ".join(map(str,common))
    d = 2 if norm_A == 4 else 0
    V15 = MMV(15)
    v = axis.v15
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
    for cl in get_axes():
        axis = get_axes()[cl]
        print("Class:", axis.g_class, ", stage =", axis.stage, 
                ", powers:", axis.powers)
        print("Automorphism group:", axis.group)
        print(display_norm_A(axis), end = "")
        print("Eigenvalues of A part of axis v:", eigen(axis.A()))
        if verbose:
            n, c = (g_central * axis.reflection()).half_order()
            print("Half order:", n, ", central involution: ", c)
            if c: 
                print("Character of central involution:",  c.chi_G_x0())
                (analyze_xy(c))
            print("Axis type of axis v * MM('t',1):", axis.axis_type(1))
            print("Dim intersection with Q_x0:", 24 - len(x_equations(axis,0)))
            print("Dim intersection with Q_y0:", 24 - len(x_equations(axis,1)))
            #print("Dim intersection with Q_z0:", 24 - len(x_equations(axis,2)))
            opp_axis_type = (v_axis_opp15 * MM0(axis.g)).axis_type()
            print("Image of opposite axis: ", opp_axis_type) 
            print("A part of axis v:")
            display_A(axis.A())
            print("\nA part of axis v * MM('t', 1):")
            display_A(axis.Ax())
            print("A part of opposite axis v:")
            display_A(axis.A_opp())
        print("")



if  __name__ == "__main__":
    test_display_axes(verbose = 1)


