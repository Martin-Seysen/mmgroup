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
from mmgroup.mm_reduce import mm_order_store_vector
from mmgroup.mm_reduce import mm_order_element_Gx0
from mmgroup.mm_reduce import mm_order_load_vector
from mmgroup.mm_reduce import mm_order_load_tag_vector
from mmgroup.mm_reduce import mm_reduce_M



MMV3 = MMV(3)
MMV15 = MMV(15)
MM = MM0  #  TODO: Fixme

ORDER_VECTOR = None
ORDER_TAGS = None


OFS_NORM_A = 0
OFS_DIAG_VA = 1
OFS_WATERMARK_PERM = 2
OFS_TAGS_Y = 26
OFS_SOLVE_Y = 37
OFS_TAGS_X = 48
OFS_SOLVE_X = 72
OFS_TAG_SIGN = 96


#######################################################################
# Find a vector stabilizing by an element of order n
#######################################################################

def stabilizer_vector(v, g, n):
    """Compute a vector stabilized by an element of the monster

    Le ``g`` be an element of the monster group of order ``n`` and 
    ``v`` a vector in a represention of the monster. We return the
    vector ``sum(v * g**i for i  in range(n))`` which is stabilized
    by ``g``. We always return ``None`` if that sum is 0 or a 
    multiple of the 1 element in the representation space. The 
    last condition is checked with a fast crude check only.  
    """
    vg = v.copy()
    w = v.copy()
    for i in range(1, n):
        vg *= g 
        w += vg
    assert v == vg * g
    if (w['B'] == 0).all():
        return None
    return w


#######################################################################
# Assemble a test vector mod 15 from the input data
#######################################################################



def make_order_vector(s_g71, s_v71, s_gA, diag, s_g94, s_v94):
    v71 = MMV15(10, s_v71)
    #print("v71 =", v71)
    g71 = MM(s_g71)
    w71 = stabilizer_vector(v71, g71, 71)
    assert w71 is not None
    w71 *= MM(s_gA)
    v94 = 6 * MMV15(s_v94)
    g94 = MM(s_g94)
    w94 = stabilizer_vector(v94 - v94 * g94, g94**2, 47)
    assert w94 is not None
    w = w71 + w94
    v3 = mm_op15_eval_A_rank_mod3(w.data, diag) & 0xffffffffffff
    assert v3 != 0
    v_type4 = gen_leech3to2_type4(v3)
    assert v_type4 == 0x800000
    w.reduce()
    return w



def map_y(y_index):
    i, j = (y_index >> 14) & 0x1f, (y_index >> 8) & 0x1f
    vect = (1 << i) + (1 << j)
    gc = vect_to_cocode(vect)
    assert 0 <= gc < 0x800
    return gc 
    
    
def map_x(x_index):
    v2 = mm_aux_index_sparse_to_leech2(x_index) 
    return ((v2 & 0xfff) << 12) | ((v2 >> 12) & 0xfff)    

ORDER_VECTOR_PRESENT = False

def compute_order_vector(recompute = False, verbose = 0):
    global  ORDER_VECTOR_PRESENT
    if ORDER_VECTOR_PRESENT and not recompute:
        return  
    try:
        from mmgroup.structures import order_vector_data
        assert not recompute
    except (ImportError, AssertionError):
        from mmgroup.structures import find_order_vector
        result = find_order_vector.find_order_vector(verbose)
        find_order_vector.write_order_vector(result)
        from mmgroup.structures import order_vector_data
        del find_order_vector
    from mmgroup.structures.order_vector_data import S_G71, S_V71
    from mmgroup.structures.order_vector_data import S_GA, DIAG_VA
    from mmgroup.structures.order_vector_data import S_G94, S_V94
    ORDER_VECTOR =  make_order_vector(
        S_G71, S_V71, S_GA, DIAG_VA, S_G94, S_V94
    )
    assert ORDER_VECTOR is not None
    OV = ORDER_VECTOR.data
    TAGS_Y = np.array(order_vector_data.TAGS_Y, dtype = np.uint32) 
    TAGS_X = np.array(order_vector_data.TAGS_X, dtype = np.uint32)
    TAG_SIGN =  np.array(order_vector_data.TAG_SIGN, dtype = np.uint32) 
    WATERMARK_PERM = np.zeros(24, dtype = np.uint32)
    ok = mm_op15_watermark_A(OV, WATERMARK_PERM)
    assert ok >= 0

    SOLVE_YT = np.zeros(11, dtype = np.uint64)
    assert len(TAGS_Y) == 11
    nrows = 0
    for y in TAGS_Y:
        eqn = map_y(y)
        nrows += leech2matrix_add_eqn(SOLVE_YT, nrows, 11, eqn)
        #print(i, hex(y), hex(eqn), nrows)
    assert nrows == 11, nrows
    SOLVE_Y = list(bitmatrix64_t(SOLVE_YT, 11))
    assert len(SOLVE_Y) == 11
    assert mm_aux_mmv_extract_sparse_signs(15, OV, TAGS_Y, 11) == 0
    
    SOLVE_XT = np.zeros(24, dtype = np.uint64)
    assert len(TAGS_X) == 24
    nrows = 0
    for i, x in enumerate(TAGS_X):
        eqn =  map_x(x)  
        nrows += leech2matrix_add_eqn(SOLVE_XT, nrows, 24, eqn)
        #print("SOLVE_XT", i, hex(x), hex(eqn), nrows)
    assert nrows == 24, nrows
    SOLVE_X = list(bitmatrix64_t(SOLVE_XT, 24))
    
    # Concatenate computed lists to the global numpy array 'ORDER_TAGS'
    ORDER_TAGS = np.array(sum(map(list, [
        [mm_op15_norm_A(OV), order_vector_data.DIAG_VA], 
        WATERMARK_PERM, TAGS_Y, SOLVE_Y, TAGS_X, SOLVE_X, TAG_SIGN
    ]), []), dtype = np.uint32)
    assert len(ORDER_TAGS) == 97, len(ORDER_TAGS)
    t0 = mm_aux_mmv_extract_sparse_signs(
        15, OV, ORDER_TAGS[OFS_TAGS_X:], 24)
    assert t0 == 0
    mm_order_store_vector(ORDER_TAGS, ORDER_VECTOR.data)
    del ORDER_VECTOR
    ORDER_VECTOR_PRESENT = True

def get_order_vector(recompute = False, verbose = 0):
    compute_order_vector(recompute, verbose)
    v = mm_vector(15)
    mm_order_load_vector(v.data)
    return v

def get_order_tag_vector(recompute = False, verbose = 0):
    compute_order_vector(recompute, verbose)
    tags = np.zeros(97, dtype = np.uint32)
    mm_order_load_tag_vector(tags)
    return tags


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
    if not ORDER_VECTOR_PRESENT:
        compute_order_vector()
    v = mm_vector(15)
    mm_order_load_vector(v.data)
    work = mm_vector(15)
    mm_op15_word(v, g3, status - 2, 1, work)
    mm_order_load_vector(work.data)
    return not mm_op15_compare(v, work)



###########################################################################
# Computing the order of an element of the monster
###########################################################################


def check_mm_order_old(g, max_order = 119, mode = 0):
    """Return order of monster group element ``g``.

    if ``order(g) < max_order`` return ``order(g)``; else return ``0``.

    If mode is ``0`` (default) we first check if ``g`` is in the 
    subgroup ``N_0 of`` the monster. If this is the case the we check 
    the order of ``g``  by calculating in ``N_0``.

    Othewise we compute the minimum ``i`` such that
    ``v * g**i == v`` for the *order vector* ``v`. 
    """
    assert isinstance(g, (MM0, MM))
    g.reduce()
    if mode == 0:
        n0 = np.zeros(5, dtype = np.uint32)
        status = mm_group_check_word_n(g._data, g.length, n0)
        if status == 0:
            return 1
        if status == 1:
            n1 = np.copy(n0)
            for i in range(2, max_order+1):
                mm_group_n_mul_element(n1, n0, n1)
                if not mm_group_n_reduce_element(n1):
                    return i
            return 0

    if not ORDER_VECTOR_PRESENT:
        compute_order_vector()
    v = mm_vector(15)
    w = mm_vector(15)
    work = mm_vector(15)
    mm_order_load_vector(v.data)
    mm_order_load_vector(w.data)
    for i in range(1, max_order+1):
        mm_op15_word(w, g._data, g.length, 1, work)
        if not mm_op15_compare(v, w):
            return i
    return 0




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
    if not ORDER_VECTOR_PRESENT:
        compute_order_vector()
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
    if not ORDER_VECTOR_PRESENT:
        compute_order_vector()
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
    if not ORDER_VECTOR_PRESENT:
        compute_order_vector()
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


###########################################################################
# Fast reduction of an element of the monster 
###########################################################################


reduce_mm_time = None

def reduce_mm(g, check = True):
    """The fastest reduction procedure for a monster element ``g``"""
    global reduce_mm_time
    if not ORDER_VECTOR_PRESENT:
        compute_order_vector()
    g1 = np.zeros(256, dtype = np.uint32)
    t_start = time.perf_counter() 
    res = mm_reduce_M(g._data, g.length, g1)
    reduce_mm_time = time.perf_counter() - t_start
    if (res < 0):
        err = "Reduction of element of monster failed"
        raise ValueError(err)
    length = res
    if check:
        w = mm_vector(15)
        work = mm_vector(15)
        mm_order_load_vector(w.data)
        mm_op15_word(w, g._data, len(g), 1, work)
        mm_op15_word(w, g1, length, -1, work)
        mm_order_load_vector(work.data)
        assert not mm_op15_compare(w, work)
         
    g._extend(length)
    g._data[:length] = g1[:length]
    g.length = length
    g.reduced = 0
    g.reduce()
    return g


###########################################################################
# Main program (for testing)
###########################################################################


if __name__ == "__main__":
    compute_order_vector(recompute = 0, verbose = 1)
        


