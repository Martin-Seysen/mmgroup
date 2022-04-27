import time
import numpy as np

from mmgroup import structures
from mmgroup.mat24 import vect_to_cocode
from mmgroup.mat24 import ploop_theta
from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup.mm import mm_vector
from mmgroup.mm import mm_aux_mmv_extract_sparse_signs
from mmgroup.structures.mm0_group import MM0Group, MM0
from mmgroup.mm_group import MMGroup, MM
from mmgroup.mm_space import MMSpace, MMV
from mmgroup.generators import mm_group_check_word_n
from mmgroup.generators import mm_group_words_equ
from mmgroup.generators import mm_group_n_mul_element
from mmgroup.generators import mm_group_n_reduce_element 
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.clifford12 import chk_qstate12
from mmgroup.clifford12 import uint64_parity
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import leech2matrix_solve_eqn
from mmgroup.clifford12 import bitmatrix64_t
from mmgroup.clifford12 import xsp2co1_half_order_word
from mmgroup.clifford12 import xsp2co1_power_word
from mmgroup.mm15 import op_eval_A_rank_mod3 as mm_op15_eval_A_rank_mod3 

from mmgroup.mm15 import op_copy as mm_op15_copy
from mmgroup.mm15 import op_compare as mm_op15_compare
from mmgroup.mm15 import op_word as mm_op15_word
from mmgroup.mm15 import op_norm_A as mm_op15_norm_A
from mmgroup.mm15 import op_watermark_A as mm_op15_watermark_A
from mmgroup.mm_reduce import mm_order_element_M
from mmgroup.mm_reduce import mm_order_element_Gx0
from mmgroup.mm_reduce import mm_order_load_vector
from mmgroup.mm_reduce import mm_order_load_tag_data
from mmgroup.mm_reduce import mm_reduce_M



MMV3 = MMV(3)
MMV15 = MMV(15)
MM = MM0  #  TODO: Fixme




def get_order_vector(*args, **kwds):
    v = mm_vector(15)
    mm_order_load_vector(v.data)
    return v

def get_order_tag_vector(*args, **kwds):
    a = np.zeros(97, dtype = np.uint32)
    res = mm_order_load_tag_data(0, a, len(a))
    assert res >= 0
    return a[:res]


###########################################################################
# Check equality of two elements of the monster
###########################################################################
 

def check_mm_equal(g1, g2, mode = 0):
    """Return ``g1 == g2`` for elements ``g1, g2`` of the monster.

    If ``mode == 0`` (default) we first try to check equality inside 
    in the subgroup ``N_0`` of the monster, which may be considerbly 
    faster. 

    If ``mode != 0`` or this is not possible we check if 
    ``v * g1 * g2**(-1) == v`` holds for the *ORDER_VECTOR* ``v``. 

    We just check the data in ``g1`` and ``g2``, ingnoring
    ``g1.group`` and ``g2.group``.
    """
    assert isinstance(g1, (MM, MM0))
    assert isinstance(g2, (MM, MM0))
    g3 = np.zeros(2 * (g1.length + g2.length) + 1, dtype = np.uint32)
    status = mm_group_words_equ(g1._data, g1.length,
        g2._data, g2.length, g3)
    if status < 2:
        return not status
    v = mm_vector(15)
    mm_order_load_vector(v.data)
    work = mm_vector(15)
    mm_op15_word(v, g3, status - 2, 1, work)
    mm_order_load_vector(work.data)
    return not mm_op15_compare(v, work)



###########################################################################
# Computing the order of an element of the monster
###########################################################################






def check_mm_order(g, max_order = 119):
    """Return order of monster group element ``g``.

    ``g`` must be an instance of class ``MM``.  The function
    returns the order of ``g`` if ``max_order`` is set to its default
    value. 

    Computing the order of an element of the monster is time consuming; 
    and in some cases we are interested in small orders only.   
    If ``max_order``is given then the function may return 0 if the
    order of ``g`` is greater than ``max_order``.
    """
    assert isinstance(g, (MM0, MM))
    g.reduce()
    o = mm_order_element_M(g._data, g.length, max_order)
    return  chk_qstate12(o)

def check_mm_half_order(g, max_order = 119):
    """Return (halved) order of monster group element ``g``.

    ``g`` must be an instance of class ``MM``.  The function
    returns a pair ``(o, h)`` where ``o`` is the order of ``g``, and
    ``h = g**(o/2)`` for an even ``o``. We put ``h = None`` if
    ``o`` is odd.

    Parameter ``max_order`` is as in function ``check_mm_order``. 

    If ``h`` is in the subgroup :math:`G_{x0}`` then ``h`` is 
    returned as a word in the generators of that subgroup.
    """
    assert isinstance(g, (MM0, MM))
    g.reduce()
    h = np.zeros(10, dtype = np.uint32)
    o1 = mm_order_element_Gx0(g._data, g.length, h, max_order)
    chk_qstate12(o1)
    if o1 == 0:
        return 0, None
    h = h[:o1 & 0xff]
    o1 >>= 8
    h2 = np.zeros(10, dtype = np.uint32)
    o2 = xsp2co1_half_order_word(h, len(h), h2)
    h2 = h2[:o2 & 0xff]
    o2 >>= 8
    o = o1 * o2
    if (o2 & 1) == 0:
        return o, g.group('a', h2).reduce()
    if (o & 1):
        return o, None
    if o2 == 1:
        return o, g ** (o1 >> 1)
    # compute q, r, such that the result is  (g**o1)**q  * g**r
    q, r = divmod(o >> 1, o1) 
    w = g.group('a', h)
    return o, (w**q * g**r).reduce()

###########################################################################
# Check if an element of the monster is in the subgroup G_x0
###########################################################################
 




def check_mm_in_g_x0(g):
    """Check if ``g`` is in the subgroup ``G_x0`` of the monster
   
    If ``g`` is in the subgroup ``G_x0`` of the monster then the
    function changes the word representing ``g`` to a (uniquely 
    defined) word in the generators of the subgroup ``G_x0`` and
    returns ``g``.  
    
    Otherwise the function does not change ``g`` and returns ``None``.
    
    ``g`` must be an instance of class 
    ``mmgroup.mm_group.MM``.
    """
    g1 = np.zeros(10, dtype = np.uint32)
 
    res = chk_qstate12(mm_order_element_Gx0(g._data,  g.length, g1, 1))
    #print("RES", hex(res))
    if ((res >> 8) != 1):
        return None 
    length = res & 0xff
    assert length <= 10
    g._extend(10)
    g._data[:length] = g1[:length]
    g.length = length
    g.reduced = 0
    g.reduce()
    return g





