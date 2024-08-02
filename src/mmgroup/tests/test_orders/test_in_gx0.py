import numpy as np
import pytest

from mmgroup import MM

from mmgroup.mat24 import ploop_theta
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.clifford12 import chk_qstate12
from mmgroup.clifford12 import uint64_parity
from mmgroup.clifford12 import leech2matrix_solve_eqn


from mmgroup import MM0Group, MM0, MM
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
from mmgroup.dev.mm_reduce.py_mm_order import py_order_check_in_g_x0



###########################################################################
# Python implementation of check_mm_in_g_x0() using low-level functions
###########################################################################

# Here we implement the python version of function check_mm_in_g_x0()
# that has been used in the development phase.



def import_all():
    global mm_order_check_in_Gx0
    global OrderVectorMod15, ov
    from mmgroup.mm_reduce import mm_order_check_in_Gx0
    from mmgroup.dev.mm_reduce.order_vector import OrderVectorMod15
    ov = OrderVectorMod15('C')




def py_check_mm_in_g_x0(g, mode = 0):
    """Check if ``g`` is in the subgroup ``G_x0`` of the monster
   
    If ``g`` is in the subgroup ``G_x0`` of the monster then the
    function changes the word representing ``g`` to a (uniquely 
    defined) word in the generators of the subgroup ``G_x0`` and
    returns ``g``.  
    
    Otherwise the function does not change ``g`` and returns ``None``.
    
    ``g`` must be an instance of class 
    ``mmgroup.mm_group.MMGroupWord``.
    """

    g1 = py_order_check_in_g_x0(g, mode)
    if g1 is None:
        return None

    g._extend(11)
    g.length = length = len(g1.mmdata)
    g._data[:length] = g1.mmdata
    g.reduced = 0
    return g

    
###########################################################################
# Function to be tested (based on C function mm_order_check_in_Gx0)
###########################################################################

def check_mm_in_g_x0(g, mode = 0):
    global err_in_g_x0_c
    assert isinstance(g, (MM0, MM))
    w = mm_vector(15)
    work = mm_vector(15, 2).ravel()
    mode = (mode | 8) & ~2
    v = ov.order_vector.data
    mm_op_copy(15, v, w)
    res = mm_op_word(15, w, g.mmdata, len(g.mmdata), 1, work)
    assert res == 0
    h = np.zeros(11, dtype = np.uint32)
    res = mm_order_check_in_Gx0(w, h, mode, work)
    if res < 0:
        s = "mm_order_check_in_Gx0 failed with error %d"
        raise ValueError(s % res)
    if res >= 0x100:
         err_in_g_x0_c = res
         return None
    return MM0('a', h[:res])


###########################################################################
# Test data for function  check_mm_in_g_x0() 
###########################################################################


def mm_order(g):
    w, order = ov.order_vector.copy().data, 0
    work = mm_vector(15)
    while order < 120:
        mm_op_word(15, w, g.mmdata, len(g.mmdata), 1, work)
        order +=1
        if mm_op_compare(15, w, ov.order_vector.data) == 0:
            return order
    raise ValueError("Order computation failed")
    


def in_gx0_testdata():
    """Test data for testing function mgroup.mm_order.check_mm_in_g_x0

    The function yields pairs ``(g, b)``. Here ``g`` an element of the
    monster group ``MM``, and ``b`` is True iff ``g`` is in the
    subgroup ``G_x0`` of the monster group.
    """
    # Yield some fixed test data
    testdata = [
      ([], True),
      ([('d', 2)], True),
      ([('d', 0x801)], True),
      ([('x', 0xA71)], True),
      ([('x', 0x1A71)], True),
      ([('y', 0xCD1)], True),
      ([('y', 0xCD1), ('x', 0x1B7f)], True),
      ([('p', 'r')], True),
      ([('l', 1)], True),
      ([(tag, 'n') for tag in "xydpl"*3], True),
      ([('t', 1)], False),
      ("M0<y_4d1h*x_0b7fh>", True),
    ]
    for g, in_g_x0 in testdata:
       yield MM0(g), in_g_x0
    # Yield elements ``g0 * t * g1``, where ``t`` is the triality 
    # element or its inverse, and ``g0 , g1`` are in ``G_x0``.
    for i in range(8):
        e = [1 + ((i >> j) & 1) for j in range(3)]
        g = [('p', 'r'), ('y', 'r'), ('l', e[0]), ('t', e[1]), ('l', e[2]),
             ('p', 'r'), ('y', 'r'), ('x', 'r'), ('d', 'r'), ]
        yield  MM0(g), False
    # Yield some less obvious examples. Let ``d`` be a random 2A 
    # involution and ``z`` be the 2B involution in the center of
    # ``G_x0``. Then ``g = z * d`` has even order ``o``, and we 
    # yield ``g ** (o/2)``, which is in ``G_x0``.
    z, d0 = MM0('x', 0x1000), MM0('d', [2,3])
    for n in range(5):
        g0 = MM0([(tag, 'n') for tag in "xydplt"*2])
        d =  d0 ** g0
        g = z * d
        o = mm_order(g)
        assert o & 1 == 0
        yield g ** (o >> 1), True
    # Yield some elements of ``G_x0``.
    for n in range(10):
        for i in range(2,5):
            yield MM0([(tag, 'n') for tag in "xydpl"*i]), True
    

###########################################################################
# Test function  check_mm_in_g_x0() 
###########################################################################


@pytest.mark.orders
def test_in_gx0(verbose = 0):
    import_all()
    print("Testing function check_mm_in_g_x0()")
    for n, (g, in_g_x0) in enumerate(in_gx0_testdata()):
        g = MM0(g)
        if verbose:
            print("Test", n+1)
            print("g =", g)
        g1 = check_mm_in_g_x0(g)
        if verbose:
            print("reduced", g1)
            if g1 is None:
                r = err_in_g_x0_c
                print("Reason for non-membership (.c):", r)
        g1_py = py_check_mm_in_g_x0(g)
        if verbose:
            print("reduced (py)", g1_py)
            if g1_py is None:
                r = err_in_g_x0_py
                print("Reason for non-membership (py):", r)
        if in_g_x0:
            assert g1_py is not None
        else:
            assert g1_py is None
        if g1_py:
            assert g == g1_py 
        assert g1 == g1_py
        #assert g.in_G_x0() == in_g_x0
        pass








