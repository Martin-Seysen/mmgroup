import numpy as np
import pytest

from mmgroup import MM

from mmgroup.mat24 import ploop_theta
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import mm_group_invert_word
from mmgroup.clifford12 import chk_qstate12
from mmgroup.clifford12 import uint64_parity
from mmgroup.clifford12 import leech2matrix_solve_eqn
from mmgroup.clifford12 import xsp2co1_check_word_g_x0
from mmgroup.clifford12 import xsp2co1_set_elem_word
from mmgroup.clifford12 import xsp2co1_elem_to_word

from mmgroup import MM0
from mmgroup.mm_op import mm_vector

from mmgroup.mm_op import mm_aux_mmv_extract_sparse_signs
from mmgroup.mm_op import mm_vector
from mmgroup.mm_op import mm_op_copy
from mmgroup.mm_op import mm_op_compare
from mmgroup.mm_op import mm_op_word
from mmgroup.mm_op import mm_op_word_tag_A        
from mmgroup.mm_op import mm_op_norm_A
from mmgroup.mm_op import mm_op_eval_A_rank_mod3 
from mmgroup.mm_op import mm_op_watermark_A_perm_num
from mmgroup.dev.mm_reduce.order_vector import OrderVectorMod15
from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord
from mmgroup.mm_space import MMVector

###########################################################################
# Python implementation of check_mm_in_g_x0() using low-level functions
###########################################################################

# Here we implement the python version of function check_mm_in_g_x0()
# that has been used in the development phase.


try:
    import mm_reduce
    ov = OrderVectorMod15('C')
except (ImportError, ModuleNotFoundError):
    ov = OrderVectorMod15('py')



err_in_g_x0_py = 0 
err_in_g_x0_c = 0 





def find_in_Q_x0(w):
    global err_in_g_x0_py
    w_x = mm_aux_mmv_extract_sparse_signs(15, w, ov.TAGS_X, 24)
    if w_x < 0:
        err_in_g_x0_py = 7
        return None
    x = leech2matrix_solve_eqn(ov.SOLVE_X, 24, w_x)
    w_sign = ((x >> 12) & 0x7ff) ^ (x & 0x800)
    aa = np.array([ov.TAG_SIGN ^ (w_sign << 14)], dtype = np.uint32)
    sign = mm_aux_mmv_extract_sparse_signs(15, w, aa, 1)
    if sign < 0:
        err_in_g_x0_py = 8
        return None
    x &= 0xffffff
    sign ^= uint64_parity(x & (x >> 12) & 0x7ff)
    x ^=  (sign & 1) << 24
    x ^= ploop_theta(x >> 12)
    #print("final x =", hex(x))
    return x


def find_in_G_x0(w):
    global err_in_g_x0_py
    g1i = np.zeros(11, dtype = np.uint32)

    if mm_op_norm_A(15, w) != ov.NORM_A: # ORDER_TAGS[OFS_NORM_A]:
        err_in_g_x0_py = 1
        return None        
    w3 = mm_op_eval_A_rank_mod3(15, w, 0)
    w3 &= 0xffffffffffff
    if w3 == 0: 
        err_in_g_x0_py = 2
        return None
    w_type4 = gen_leech3to2_type4(w3)
    if w_type4 == 0: 
        err_in_g_x0_py = 3
        return None
    wA = np.array(w[:2*24], copy = True)
    len_g1 = gen_leech2_reduce_type4(w_type4, g1i)
    assert 0 <= len_g1 <= 6 
    res = mm_op_word_tag_A(15, wA, g1i, len_g1, 1)
    assert res == 0
    perm_num = mm_op_watermark_A_perm_num(15, ov.WATERMARK_PERM, wA)
    if perm_num < 0: 
        err_in_g_x0_py = 4
        return None
    if perm_num > 0:
        g1i[len_g1] = 0xA0000000 + perm_num 
        res = mm_op_word_tag_A(15, wA, g1i[len_g1:], 1, 1)
        assert res  == 0
        len_g1 += 1
    v_y = mm_aux_mmv_extract_sparse_signs(15, wA, ov.TAGS_Y, 11)
    if v_y < 0:
        err_in_g_x0_py = 5
        return None
    y = leech2matrix_solve_eqn(ov.SOLVE_Y, 11, v_y)
    if y > 0:
        g1i[len_g1] = 0xC0000000 + y
        res = mm_op_word_tag_A(15, wA, g1i[len_g1:], 1, 1)
        assert res  == 0
        len_g1 += 1
    if (wA != ov.order_vector.data[:2*24]).all():
        err_in_g_x0_py = 6
        return None
    #print("g1i", g1i[:len_g1])
    return g1i[:len_g1]


def py_order_check_in_g_x0(g, mode = 0):
    """Check if ``g`` is in the subgroup ``G_x0`` of the monster
   
    If ``g`` is in the subgroup ``G_x0`` of the monster then the
    function changes the word representing ``g`` to a (uniquely 
    defined) word in the generators of the subgroup ``G_x0`` and
    returns ``g``.  
    
    Otherwise the function does not change ``g`` and returns ``None``.
    
    ``g`` must be an instance of class 
    ``mmgroup.mm_group.MMGroupWord``, or of class ``AbstractMMGroup``.
    """
    global err_in_g_x0_py
    err_in_g_x0_py = 0
    work = mm_vector(15)
    w = mm_vector(15)
    v = ov.order_vector.data
    if isinstance(g, AbstractMMGroupWord):
        mm_op_copy(15, v, w)
        res = mm_op_word(15, w, g.mmdata, len(g.mmdata), 1, work)
        assert res == 0
    elif isinstance(g, MMVector):
        assert g.p == 15
        mm_op_copy(15, g.data, w)
    else:
        E = "g must be in the Monster or in its representation mod 15"
        raise TypeError(E)  

    g1i = find_in_G_x0(w)
    if g1i is None:
        return None
    res = mm_op_word(15, w, g1i, len(g1i), 1, work)
    assert res == 0

    x = np.zeros(0, dtype = np.uint32) if mode & 4 else find_in_Q_x0(w)
    if x == None:
        return None

    g2i = np.array([0x90000000 + (x & 0xfff), 
        0xB0000000 + ((x >> 12) & 0x1fff)], dtype = np.uint32)
    res = mm_op_word(15, w, g2i, 2, 1, work)  
    assert res == 0   
    g1i = np.append(g1i, g2i)
 
    if mm_op_compare(15, v, w):
        #print("vW", v, "\n",  w)
        err_in_g_x0_py = 9
        return None

    if (mode & 1) == 0:
        mm_group_invert_word(g1i, len(g1i))
    return MM0('a', g1i)
    


def py_order_element_Gx0(g, o = 119):
    o = max(1, min(o, 119))
    assert isinstance(g, AbstractMMGroupWord)
    status = xsp2co1_check_word_g_x0(g.mmdata, len(g.mmdata));
    if status == 0:
        if len(g.mmdata) == 0:
            return 1, MM0()
        else:
            h = np.zeros(10, dtype = np.uint32)
            elem = np.zeros(26, dtype = np.uint64)
            res = xsp2co1_set_elem_word(elem, g.mmdata, len(g.mmdata))
            assert res >= 0, res
            res = xsp2co1_elem_to_word(elem, h)
            assert res >= 0, res
            return 1, MM0('a', h[:res])

    w = ov.order_vector.copy()
    for i in range(o): 
        w *= g
        h = py_order_check_in_g_x0(w)
        if h:
            return i + 1, h
  
    return 0, None
  

def py_order_element_M(g, o = 119):
    e1, h = py_order_element_Gx0(g, o)
    assert 0 <= e1 <= 119
    if e1:
        e2 =  xsp2co1_order_word(h.mmdata, len(h.mmdata))
        assert e2 > 0
        return e1 * e2
    return 0 




