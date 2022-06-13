"""Test the C function xsp2co1_elem_conjugate_involution_Gx0

Given a involution  :math:`g` in the subgroup :math:`G_{x0}` of the 
monster, the C function ``xsp2co1_elem_conjugate_involution_Gx0`` in
file ``xspco1_traces.c`` computes a representative of the class of 
:math:`g` in  :math:`G_{x0}`. That function also computes a number as 
an indication for the class, as give by the C function
``xsp2co1_elem_involution_class`` in the same file. For the 
numbering of the classes of involutions in :math:`G_{x0}` we also
also refer to module ``mmgroup.tests.test_involutions.test_xp2_trace``.
"""

import sys
import os


if __name__ == "__main__":
    sys.path.append("../../../")

from random import randint
from collections import OrderedDict
import numpy as np
import pytest


from mmgroup import MM0, Xsp2_Co1, XLeech2, mat24
from mmgroup.generators import mm_group_n_mul_word_scan
from mmgroup.generators import mm_group_n_conj_word_scan
from mmgroup.generators import mm_group_n_reduce_element
from mmgroup.generators import mm_group_n_mul_element
from mmgroup.generators import mm_group_n_inv_element
from mmgroup.generators import mm_group_n_to_word
from mmgroup.generators import mm_group_invert_word
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import gen_leech2_reduce_n
from mmgroup.generators import gen_leech2_op_word
from mmgroup.clifford12 import xsp2co1_elem_find_type4
from mmgroup.clifford12 import xsp2co1_elem_to_N0
from mmgroup.clifford12 import xsp2co1_conjugate_elem
from mmgroup.clifford12 import xsp2co1_elem_involution_class
from mmgroup.clifford12 import xsp2co1_elem_conjugate_involution_Gx0
from mmgroup.clifford12 import xsp2co1_map_involution_class_Gx0
from mmgroup.tests.test_involutions.test_involution_invariants import INVOLUTION_SAMPLES

#print("Hello world")


#######################################################################
# Auxiliary class for function xsp2co1_elem_conjugate_involution_Gx0
#######################################################################


class N_x0_involution:
    """Model the natural mapping :math:`G_{x0} \rightarrow N_{x0}`


    """  
    IND_Y = 1
    IND_X = 2
    
    def __init__(self, elem):
        self.gn0 =  np.zeros(5, dtype = np.uint32)
        self.tf_n0 =  np.zeros(5, dtype = np.uint32)
        if (xsp2co1_elem_to_N0(elem, self.gn0)) < 0:
            raise ValueError("xsp2co1_elem_from_N0 failed")
        if (self.gn0[0] | self.gn0[4]):
            raise ValueError("Bad N_x0 element")
        mm_group_n_mul_element(self.gn0, self.gn0, self.tf_n0)
        if mm_group_n_reduce_element(self.tf_n0):
            raise ValueError("Element is not an involution")
        mm_group_n_reduce_element(self.gn0)

    def get_xy(self, index):
        return self.gn0[index] & 0x7ff

    def get_q(self):
        assert self.gn0[0] | self.gn0[1] | self.gn0[4]  == 0
        q = self.gn0[2]
        return (q << 12) ^ mat24.ploop_theta(q) ^ self.gn0[3]

    def transform(self, a):
        mm_group_n_conj_word_scan(self.gn0,  a, len(a))
        mm_group_n_mul_word_scan(self.tf_n0, a, len(a))
        mm_group_n_reduce_element(self.gn0)

    def out(self, a):
        inv = np.zeros(5, dtype = np.uint32)
        mm_group_n_inv_element(self.tf_n0, inv)
        length = mm_group_n_to_word(inv, a)
        mm_group_invert_word(a, length)
        return length

    def display_transformed(self):
        a = np.zeros(5, dtype = np.uint32)
        length = mm_group_n_to_word(self.gn0, a)
        print(MM0('a', a[:length]).reduce())
    
  
#########################################################################
# Python implementation of function xsp2co1_elem_conjugate_involution_Gx0
#########################################################################



def xsp2co1_elem_conjugate_involution_in_Gx0_py(elem, guide, a):
    elem1 = np.zeros(26, dtype = np.uint64)
    v4 = xsp2co1_elem_find_type4(elem, guide)
    if (v4 < 0):
        raise ValueError("xsp2co1_elem_find_type4 failed")
    len_a = gen_leech2_reduce_type4(v4, a);
    assert len_a >= 0
    np.copyto(elem1, elem)
    if (xsp2co1_conjugate_elem(elem1, a, len_a)) < 0:
        raise ValueError("xsp2co1_conjugate_elem failed")  
    invol = N_x0_involution(elem1)

    if invol.get_xy(invol.IND_Y):
        b = np.zeros(3, dtype = np.uint32)
        if invol.get_xy(invol.IND_X):
            b[0] = 0x10000800
            invol.transform(b[:1])
        b[0] = 0x50000002   # (this is t**2)
        invol.transform(b[:1])
        vy = invol.get_q()
        gen_leech2_reduce_n(vy, b)
        invol.transform(b[:3])
        b[0] = 0x50000001   # (this is t**1)
        invol.transform(b[:1])
        len_a2 = invol.out(a[len_a:])
        return len_a + len_a2
    else:
        assert len_a == 0
        vx = invol.get_q()
        t = gen_leech2_type(vx)
        if t == 0:
            if  (gen_leech2_type(guide)) == 4:
                len_a = gen_leech2_reduce_type4(guide, a);
                assert len_a >= 0
                return len_a
            else:
                return 0
        elif t == 4:
            len_a = gen_leech2_reduce_type4(vx, a);
            assert len_a >= 0
            vx = gen_leech2_op_word(vx, a, len_a)
            if vx & 0x1000000:
                 a[len_a] = 0x90000800
                 len_a += 1
        elif t == 2:
            len_a = 0
            if  (gen_leech2_type(guide) == 4 and 
                gen_leech2_type(guide ^ vx)) == 2:
                len_a = gen_leech2_reduce_type4(guide, a);
                vx =  gen_leech2_op_word(vx, a, len_a)       
            len_a2 = gen_leech2_reduce_type2(vx, a[len_a:]) 
            assert len_a2 >= 0
            vx =  gen_leech2_op_word(vx, a[len_a:], len_a2)       
            len_a += len_a2
            if vx & 0x1000000:
                a[len_a] = 0xB0000200
                len_a += 1
        else:
            raise ValueError("Not an involution")
        assert len_a  <= 10
        return len_a

#######################################################################
# Pythonic version of function xsp2co1_elem_conjugate_involution_Gx0
#######################################################################

 
                
def conjugate_involution_in_Gx0(elem, guide = 0):
    elem1 = Xsp2_Co1(elem)
    a = np.zeros(10, dtype = np.uint32)
    len_a = xsp2co1_elem_conjugate_involution_in_Gx0_py(
        elem._data, guide, a)
    a = a[:len_a]
    if (xsp2co1_conjugate_elem(elem1._data, a, len_a)) < 0:
        print(elem1, [hex(x) for x in a])
        raise ValueError("xsp2co1_conjugate_elem failed")  
    return MM0('a', list(a)), elem1


#######################################################################
# Wrapper for C function xsp2co1_elem_conjugate_involution_Gx0
#######################################################################

def xsp2co1_elem_conjugate_involution_Gx0_C(elem, guide = 0):
    a = np.zeros(10, dtype = np.uint32)
    data = Xsp2_Co1(elem)._data
    res = xsp2co1_elem_conjugate_involution_Gx0(data, guide, a)
    assert res >= 0
    iclass, length = res >> 8, res & 0xff
    assert 0 <= length <= 10
    return iclass, MM0('a', a[:length])


#######################################################################
# Dictionary of rpresentatives of classes of involutions in G_x0
#######################################################################


STD_REP = None

def get_std_rep():
    global STD_REP
    neutral = Xsp2_Co1()
    if STD_REP is not None:
        return STD_REP
    std_rep_dict = OrderedDict()
    for itype, gs in INVOLUTION_SAMPLES:
        g = Xsp2_Co1(gs)
        if g**2 == neutral:
            iclass = xsp2co1_elem_involution_class(g._data)
            c, g1 = conjugate_involution_in_Gx0(Xsp2_Co1(g))
            assert g ** Xsp2_Co1(c) == g1
            std_rep_dict[iclass] = g1 
    STD_REP = std_rep_dict
    neg = MM0(STD_REP[0x1121] *  STD_REP[0x1122]).reduce()
    assert neg == MM0('x', 0x1000), neg
    return STD_REP

   

#######################################################################
# Test data for tesing function xsp2co1_elem_conjugate_involution_Gx0
#######################################################################


def conjugate_testdata(n_samples = 10):
    unit, neg = Xsp2_Co1(), Xsp2_Co1('x', 0x1000)
    Omega = XLeech2(0x800000)
    std_rep = get_std_rep()
    for itype, gs in INVOLUTION_SAMPLES:
        iclass = xsp2co1_elem_involution_class(Xsp2_Co1(gs)._data)
        if itype[0][0] <= 2:
            g =  Xsp2_Co1(gs)
            std_g = std_rep[iclass]
            n = 4 if g in [unit, neg] else n_samples
            for i in range(n):
                transform = Xsp2_Co1('r', 'G_x0') 
                if i & 1: # and  iclass != 0x21:
                    img_omega = Omega * transform
                    guide = img_omega.ord & 0xffffff
                    yield std_g ** transform, guide  
                else:
                    yield g ** transform, 0  
          


#######################################################################
# Test the C function xsp2co1_elem_conjugate_involution_Gx0
#######################################################################


@pytest.mark.involution
def test_std_rep(ntests = 50, verbose = 0):
    Omega = XLeech2(0x800000)
    print("Test conjugating involutions in G_x0 to standard form")
    std_rep = get_std_rep()
    for g, guide in conjugate_testdata(ntests):
        iclass = xsp2co1_elem_involution_class(g._data)
        if verbose:
            print("Involution class %s:\n g =" % hex(iclass), g)
            if guide:
                print(" guide =", hex(guide))
        c, g1 = conjugate_involution_in_Gx0(g, guide)
        if verbose:
            print(" t**-1 =", (MM0(c)**(-1)).reduce())
            print(" g**t= ", g1)
        assert g1 == std_rep[iclass], (g1, std_rep[iclass])
        if guide:
            Omega1 = XLeech2(guide) * c
            assert (Omega1.ord ^ Omega.ord) & 0xffffff == 0
        c_class, cc = xsp2co1_elem_conjugate_involution_Gx0_C(g, guide)
        if (cc != c or c_class != iclass):
            print("Involution class %s:\n g =" % hex(iclass), g)
            if c_class != iclass:
                print("Class computed by C function", hex(c_class))
            if guide:
                print(" guide =", hex(guide))
            print(" t**-1 =", (MM0(c)**(-1)).reduce())
            print(" array", [hex(x) for x in c.mmdata])
            print(" g**t= ", g1)
            print(" Result obtained by C function")
            print(" t**-1 =", (MM0(cc)**(-1)).reduce())
            print(" array", [hex(x) for x in cc.mmdata])
            err = "Error in C function"
            raise ValueError(err)
        a = np.zeros(2, np.uint32)
        len_a =  xsp2co1_map_involution_class_Gx0(iclass, a)
        g2 = Xsp2_Co1('a', a[:len_a])
        assert g1 == g2, (hex(iclass), g1, g2)
    print("passed")


#######################################################################
# Display representatives of classes of involutions in G_x0
#######################################################################

def vector_to_bitlist(x):
    w, l =  mat24.vect_to_bit_list(x)
    return l[:w]

def display_q_xy(g):
    g_n0 = N_x0_involution(Xsp2_Co1(g)._data)
    try:
        print("  in N_x0:", hex(g_n0.get_q()))
    except (ValueError, AssertionError):
        b = np.array([0x50000002], dtype = np.uint32)
        g_n0.transform(b)
        print("  in N_y0:",  hex(g_n0.get_q()))
       

       
def display_g(g):
    point = 0
    def display_code(x, prefix):
        sign = (x >> 12) & 1
        gc = mat24.gcode_to_vect(x)
        cpl = mat24.bw24(gc) > 12
        if cpl:
           gc =  mat24.gcode_to_vect(x ^ 0x800)
        l = vector_to_bitlist(gc)
        if len(l) == 0:
            s = "Omega" if cpl else "1"
        else:
            point = l[0]
            s = str(l)
            if cpl: s = '~' + s
        if sign:
            s = '-' + s
        return '  ' + prefix + "_" + s
        
    print(g, "=")

    g = MM0(g).reduce()
    if g == MM0():
        print("  1")
        return
    for x in g.mmdata:
        tag = x >> 28
        if tag == 1:
            print("  d_" + str(mat24.cocode_to_bit_list(x, point)))
        elif tag == 3:
            print(display_code(x, 'x'))
        elif tag == 4:
            print(display_code(x, 'y'))
        else:
            raise ValueError("Cannot display group element")
    display_q_xy(g)
        
def display_std_rep():
    print("\nRepresentatives of classes of involutions in G_x0")
    for iclass, g in get_std_rep().items():
        print(hex(iclass), end = ": ")
        display_g(g)



#######################################################################
# Display a table mapping numbers of classes to representatives in G_x0
#######################################################################

def display_involution_map():
    s = """
// Here is a mapping from the numbers of the involution classes
// to the elements of ``G_x0``. Images are given as alements of 
// the Monster. This mapping has been computed by function 
// ``display_involution_map`` in module ``test_involution_G_x0.py.
"""
    print(s)
    classes, data = [], []
    for iclass, g in get_std_rep().items():
        classes.append(iclass)
        g_data = list((MM0(g).reduce()).mmdata)
        assert len(g_data) <= 2, g_data
        while len(g_data) < 2:
            g_data.append(0)
        data.append(g_data)
    print("static uint16_t _MAP_INVOLUTION_KEYS[] = {")
    for (i, c) in enumerate(classes):
            separator = ", " if i < len(data) - 1 else ""
            print("0x%04x%s" % (c, separator))
    print("};")
    print("static uint32_t _MAP_INVOLUTION_VALUES[][2] = {")
    for i, (d1, d2) in enumerate(data):
            separator = ", " if i < len(data) - 1 else ""
            print("{%s, %s}%s" % (hex(d1), hex(d2), separator))
    print("};")
    

#######################################################################
# Execution as a stand-alone program
#######################################################################




if __name__ == "__main__":
    test_std_rep(100, 1)
    display_std_rep()
    display_involution_map()

