import sys
import numpy as np

from mmgroup import MM0, mat24, MMSpace 
from mmgroup.mm_op import mm_op_compare_abs
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import leech2matrix_prep_eqn
from mmgroup.clifford12 import leech2matrix_solve_eqn
from mmgroup.mm_op import mm_aux_mmv_extract_sparse_signs
from mmgroup.mm_op import mm_aux_mmv_extract_sparse




#######################################################################
# subroutines
#######################################################################

def equal_upto_sign(axis1, axis2):
    v1 = axis1.v15
    v2 = axis2.v15
    eq =  mm_op_compare_abs(15, v1.data, v2.data) == 0
    for part in "ABC":
        eq = eq and  (v1[part] == v2[part]).all()
        pass
    return eq


def abs_mod15(a):
    return np.where(a > 7, 15 - a , a)


def display_axis_equal_upto_sign(ax, orbit):
    rep = AXES[orbit].copy()
    eq = equal_upto_sign(rep, ax)
    s = "" if eq else "not "
    print("%s Axes are %sequal up to sign" % (orbit, s))
    if (not eq):
            for e in "ABC":
                if (ax[e] == rep[e]).all():
                     continue
                s = "Axes differ at part %s" % e
                ax.display_sym(e, rep, text = s)
            for e in "TXY":
                 if not (abs_mod15(ax[e]) == abs_mod15(rep[e])).all():
                     print("Axes differ at part", e)
                     for i in range(759):
                         diff = abs_mod15(ax[e, i]) - abs_mod15(rep[e, i])
                         if (diff != 0).any():
                              #print("Axes differ at part", e, ", row", i,
                              #         hex(mat24.octad_to_vect(i)))
                              #print(diff)
                              #break
                              pass   




#######################################################################
# equations
#######################################################################


def augment_v_data(v, data):
    a = np.array(data, dtype = np.uint32)
    mm_aux_mmv_extract_sparse(v.p, v.data, a, len(a))
    return a





def eqn_x(tag, i, j):
    #from mmgroup.mm_op import  mm_aux_index_sparse_to_leech2
    #x_index = MMSpace.index_to_sparse(tag, i, j)
    #v2 = mm_aux_index_sparse_to_leech2(x_index)
    v = MMSpace.index_to_short_mod2(tag, i, j)
    #assert v == v2
    assert v > 0 
    return ((v & 0xfff) << 12) + ((v >> 12) & 0xfff) 


 

def display_equations(text, eqn):
    if text:
        print(text)
    al = 0
    for i, a in enumerate(eqn):
        print("%2d" % i, "".join([str((a >> i) & 1) for i in range(24)]))
        al |= a
    print("AL", "".join([str((al >> i) & 1) for i in range(24)]))


def x_equations(ref_axis, baby = False):
    v = ref_axis.v15.copy()
    nrows = 0
    solve_t = np.zeros(24, dtype = np.uint64)
    sp_values = []
    TAGS = "BCTX"
    for tag in TAGS:
        a = v[tag]
        for i, row in enumerate(a):
            if (row == 0).all():
                continue
            for j, entry in enumerate(row):
                if entry != 0:
                    eqn = eqn_x(tag, i, j)
                    n = leech2matrix_add_eqn(solve_t, nrows, 24, eqn)
                    if n:
                        sp = MMSpace.index_to_sparse(tag, i, j)
                        sp_values.append(sp)
                        nrows += n
    if baby:
        eqn = eqn_x('B', 2, 3)
        nrows += leech2matrix_add_eqn(solve_t, nrows, 24, eqn)
    sp_values = augment_v_data(v, sp_values)
    #print("SOLVE_X", [hex(x) for x in sp_values])
    lsp = len(sp_values)
    ok = mm_aux_mmv_extract_sparse_signs(15, v.data, sp_values, lsp) == 0
    assert ok
    equations = np.zeros(nrows, dtype = np.uint32)
    ok = leech2matrix_prep_eqn(solve_t, nrows, 24, equations) == 0
    assert ok
    #print("Ax_x", ref_axis.g_class, nrows)
    #display_equations("x equations", equations)
    return sp_values, equations


   


def solve_x_equations(axis, sp_values, equations):
    w = mm_aux_mmv_extract_sparse_signs(15, axis.v15.data, 
        sp_values, len(sp_values))
    if w < 0:
        raise ValueError("Equation has no solution_")
    x = leech2matrix_solve_eqn(equations, len(equations), w)
    return MM0('q', x & 0xffffff)



#######################################################################




def central_equations(ref_axis):
    v = ref_axis.v15
    nrows = 0
    solve_t = np.zeros(24, dtype = np.uint64)
    sp_values = []
    map_eqn = {'X': 1, 'Z':2, 'Y':3} 
    for tag in "XZY":
        if nrows >= 2:
            break
        a = v[tag]
        for i, row in enumerate(a):
            if (row != 0).any():
                j = np.flatnonzero(row)[0]
                assert v[tag, i, j] == row[j]
                assert v[tag, i, j] != 0
                n = leech2matrix_add_eqn(solve_t, nrows, 2, map_eqn[tag])
                if n:
                   sp = MMSpace.index_to_sparse(tag, i, j) + row[j]
                   sp_values.append(sp)
                   nrows += 1
                   break
    sp_values = np.array(sp_values, dtype = np.uint32)
    equations = np.zeros(nrows, dtype = np.uint32)
    assert leech2matrix_prep_eqn(solve_t, nrows, 2, equations) == 0
    #print("Ax", ref_axis.g_class, nrows)
    return sp_values, equations



def solve_central_equations(axis, sp_values, equations):
    w = mm_aux_mmv_extract_sparse_signs(15, axis.v15.data, 
        sp_values, len(sp_values))
    if w < 0:        
        back = np.get_printoptions()
        np.set_printoptions(formatter={'int':hex})
        print("Trying to solve equations for values:") 
        print(sp_values)
        print(hex(w))
        np.set_printoptions(**back) 
        raise ValueError("Equation has no solution")
    x = leech2matrix_solve_eqn(equations, len(equations), w)
    return MM0('x', (x & 3) << 11)


######################################################################



EQUATION_ORBITS = [
     '4A', '4B', '4C', '6A', '6C', '8B', '6F', '10A', '10B', '12C'
]


BABY_EQUATION_ORBITS = [
   '2B1', '2B0', '4A1', '4B1', '4C1', '6A1', '6C1', '10A1'
]




def compute_Qx0_equations(axis):
    baby = axis.axis_class == "baby"
    x_equ = x_equations(axis, baby)
    central_equ = central_equations(axis)
    return x_equ + central_equ

def compute_Qx0_equations_str(axis):
    equ = compute_Qx0_equations(axis)
    s = "(\n"
    for a in equ:
        s += "[" + ", ".join([hex(x) for x in a]) + "],\n"
    return s + "),\n"



def solve_Qx0_equations(orbit, axis):
    if isinstance(orbit, str):
        if orbit in EQUATION_ORBITS:
            from mmgroup.tests.axes.get_sample_axes import get_sample_axes
            equ = get_sample_axes()[orbit].Qx0_equations
        elif orbit in BABY_EQUATION_ORBITS:
            from mmgroup.tests.axes.get_baby_sample_axes import get_baby_sample_axes
            equ = get_baby_sample_axes()[orbit].Qx0_equations
        else:
            return axis.group()
    else:
        equ = orbit
    g1 = solve_x_equations(axis, equ[0], equ[1])
    ax = axis * g1
    g2 = solve_central_equations(ax, equ[2], equ[3])
    return axis.group(g1 * g2)






