"""We implement soe simple cases of the Griess algebra.
 
"""

# If import fails then we want a quick failure
import mmgroup
from mmgroup.tests.axes import sample_axes

import numpy as np

import mmgroup
from mmgroup import MM0, MM, MMV, MMVector, characteristics
try:
   from mmgroup import mat24
except (ImportError, ModuleNotFoundError):
   from mmgroup.dev.mat24.mat24_ref import Mat24
   mat24 = Mat24
from mmgroup.mm_op import mm_op_mul_std_axis, mm_op_word
from mmgroup.tests.axes.axis import Axis, BabyAxis
from mmgroup.tests.axes.axis import G_AXIS, V_AXIS

characteristics()
VV = MMV(15)

#############################################################
# Random elements in G_x0 and in G_x0 \cap 2.B
#############################################################


def _fast_mul_v15(v, g, e = 1):
    """Change v to v * g**e inplace

    Here v is a vector in the representation of the Monster mod 15,
    g is an element of the Monster, and e is an integer. 
    """
    work = VV()
    mm_op_word(15, v.data, g.mmdata, len(g.mmdata), e, work.data)



def product_axis_2A(ax1, ax2):
    """Return axis of product of two involutions

    Let ``ax1``, ``ax2`` be the axes of two 2A involutions ``t1``
    and ``t2`` given  as instances of class ``Axis``. If the product
    ``t`` of ``t1`` and ``t2`` is another 2A involution, then the
    function returns the Axis of ``t``.

    If ``check`` is True (default) then the function fails
    if ``t`` is not a 2A involution. Otherwise the function
    returns garbage in that case.

    Implementation:

    We use the formula in :cite:`Con85`, Table 3 for the case that
    the product ``t`` is of type 2A. This requires the Griess algebra
    product ``ax1 * ax2``, which is implemented for the case that
    ``ax1`` is the standard 2A axis. So we first transform ``ax1``
    to the standard axis.
    """
    assert isinstance(ax1, Axis)
    assert isinstance(ax2, Axis)
    if ax1.scalprod15(ax2, True) != 1:
        ERR = "Product of 2A involutions of axes must be in class 2A"
        raise ValueError(ERR)
    """ 
    # This is a standard implementation
    ax2red = ax2 * ax1.g ** -1
    v2red = ax2red.v15
    v_prod = ax1.v_axis15 + v2red
    mm_op_mul_std_axis(15, v2red.data)
    v_prod -= v2red >> 2
    ax3 = Axis(v_prod * ax1.g) 
    return ax3
    """
    # This is a faster implementation using low-level functions
    v2red = ax2.v15.copy()
    _fast_mul_v15(v2red, ax1.g, -1)
    v_prod = v2red + ax1.v_axis15
    mm_op_mul_std_axis(15, v2red.data)
    v_prod -= v2red >> 2
    _fast_mul_v15(v_prod, ax1.g)
    return Axis(v_prod)

