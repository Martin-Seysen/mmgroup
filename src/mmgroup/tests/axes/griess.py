"""We implement soe simple cases of the Griess algebra.
 
"""

from numbers import Integral

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
from mmgroup.generators import mm_group_invert_word 
from mmgroup.mm_op import mm_op_mul_std_axis, mm_op_word
from mmgroup.tests.axes.axis import Axis
from mmgroup.tests.axes.axis import G_AXIS, V_AXIS
from mmgroup.mm_reduce import mm_reduce_vector_vp

characteristics()
VV = MMV(15)

#############################################################
# Random elements in G_x0 and in G_x0 \cap 2.B
#############################################################


def _fast_mul_v15(v15, g, e = 1):
    """Change v to v * g**e inplace

    Here v15 is a vector in the representation of the Monster mod 15,
    g is an element of the Monster, and e is an integer.

    g may also be a numpy array of generators of the Monster.
    """
    work = VV()
    if not isinstance(g, np.ndarray):
        g = g.mmdata
    mm_op_word(15, v15.data, g, len(g), e, work.data)


def _fast_reduce_axis_v15(v15, e = 1):
    """Group element data transforming axis to standard axis

    Here v15 is a vector in the representation of the Monster mod 15,
    that must be a 2A axis. The function computes a group element g
    transforming axis v15 to the standard axis and returns g as a
    numpy array of generators.

    In case e = -1 then the inverse of g is computed instead. 
    """
    v0 = np.zeros(1, dtype = np.uint32)
    v = v15.copy()
    w = MMV(15)(0)
    g = np.zeros(256, dtype = np.uint32)
    l_g = mm_reduce_vector_vp(v0, v.data, 1, g, w.data)
    assert 0 <= l_g < 128, hex(l_g)
    assert abs(e) == 1
    if e == -1:
        mm_group_invert_word(g, l_g)
    return g[:l_g]

def product_axis_2A(ax1, ax2):
    """Return axis of product of two involutions

    Let ``ax1``, ``ax2`` be the axes of two 2A involutions ``t1``
    and ``t2``, given  as instances of class ``Axis``. If the product
    ``t`` of ``t1`` and ``t2`` is another 2A involution, then the
    function returns the Axis of ``t``.

    The function fails if ``t`` is not a 2A involution.

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
    g_data = ax1._fast_g_mmdata()
    v2red = ax2.v15.copy()
    _fast_mul_v15(v2red, g_data, -1)
    v_prod = v2red + ax1.v_axis15
    mm_op_mul_std_axis(15, v2red.data)
    v_prod -= v2red >> 2
    _fast_mul_v15(v_prod, g_data)
    g_out_data = _fast_reduce_axis_v15(v_prod, -1)
    return Axis('a', g_out_data)


def adjust_vector_p(vector, p):
    if p == vector.p:
        return vector
    elif vector.p % p == 0:
        return vector % p
    else:
        ERR = "Incompatible moduli of representation vectors"
        raise ValueError(ERR)


def auto_modulus(a):
    if isinstance(a, MMVector):
        return a.p
    if isinstance(a, tuple) and len(a) == 2:
        p0, p1 = auto_modulus(a[0]), auto_modulus(a[1])
        if p0 == p1 or p1 is None:
            return p0
        elif p0 == None:
            return p1
        ERR = "Incompatible moduli of representation vectors"
        raise ValueError(ERR)
    return None


def mul_axis_vector(axis, vector):
    g_data = axis._fast_g_mmdata()
    res = vector.copy()
    _fast_mul_v15(res, g_data, -1)
    mm_op_mul_std_axis(15, res.data)
    _fast_mul_v15(res, g_data)
    return res



def mul_ax_ax(ax1, ax2):
    p = ax1.p
    assert p == ax2.p
    if ax1.tag == "ax":
        if ax2.tag == "ax":
            v1, v2 = ax1.value, ax2.value
            assert isinstance(v1, Axis), type(v1)
            assert isinstance(v2, Axis)
            if v1.scalprod15(v2, None) == 0:
                return GriessIntermediate(0, p)
            res = ax2.as_vector()
            if res.tag == '1':
                assert res.value == 0
                return res
            g_data = v1._fast_g_mmdata()
            v2res = ax2.as_vector().value
            _fast_mul_v15(v2res, g_data, -1)
            mm_op_mul_std_axis(15, v2res.data)
            _fast_mul_v15(v2res, g_data)
            result = GriessIntermediate(v2res, p)
            result.n_axes =  ax1.n_axes + ax2.n_axes
            result.n_axes =  ax1.factor * ax2.factor
            return result
    ERR = "Griess algebra product is to complicated"
    raise ValueError(ERR)

class GriessIntermediate:
    TYPE_DICT = {
        "v": MMVector, "ax": Axis, "1": Integral 
    }
    COMPLEXITY = {
        "1": 0, "ax":1, "v":2
    }
    def __init__(self, a, p):
        self.p = p
        if isinstance(a, Axis):
            self.set_axis(a)
        elif isinstance(a, MMVector):
            self.set_vector(a)
        elif isinstance(a, Integral):
            self.set_int(a)
        elif isinstance(a, GriessIntermediate):
            self.set_mul_int(a)
        elif isinstance(a, tuple):
            if len(a) != 2:
                ERR = "Griess algebra product must have 2 factors"
                raise ValueError(ERR)
            a, b = a
            if not isinstance(a, GriessIntermediate):
                a = GriessIntermediate(a, p)
                a.check()
            if not isinstance(b, GriessIntermediate):
                b = GriessIntermediate(b, p)
                b.check()
            ca, cb = a.complexity(), b.complexity()
            if cb > ca:
                a, b = b, a
            if a.tag == "1":
                self.set_mul_int(b, a.factor)
            elif cb > 1:
                if ca > 1:
                    error_complicated()
                vb = b.as_vector_out()
                v = self.mul_ax_vect(a, vb)
                self.set_vector(v, a.factor * b.factor,
                    a.n_axes + b.n_axes)
            elif a.tag == b.tag == "ax":
                self.set_mul_ax_ax(a, b)
            else:
                error_complicated()
        else:
            ERR = "Cannot convert %s object to representation vector"
            raise ValueError(ERR % type(a))

    def complexity(self):
        return self.COMPLEXITY[self.tag]

    def error_complicated():
        ERR = "Expression to complicated for Griess algebra"
        raise ValueError(ERR)

    def check(self):
        assert self.p & 1
        assert isinstance(self.factor, Integral)
        assert self.tag in self.TYPE_DICT, (self.tag)
        assert isinstance(self.value,  self.TYPE_DICT[self.tag])

    def set_int(self, value):
        self.tag = "1"
        self.factor = int(value)
        self.value = 1
        self.n_axes = 0
        self.check()

    def set_mul_int(self, other, factor=1):
        assert self.p == other.p
        if factor * other.factor == 0:
            self.set_int(0)
        else:
            other.check()
            self.tag = other.tag
            self.factor = factor * other.factor
            self.value = other.value
            self.n_axes = other.n_axes
            self.check()

    def set_axis(self, axis, factor = 1):
        assert isinstance(axis, Axis)
        self.tag = "ax"
        self.factor = int(factor)
        self.value = axis
        self.n_axes = 1
        self.check()

    def set_vector(self, vector, factor = 1, n_axes = 0):
        assert isinstance(vector, MMVector)
        self.tag = "v"
        self.factor = int(factor)
        self.value = adjust_vector_p(vector, self.p)
        self.n_axes = 0
        self.check()

    def set_mul_ax_ax(self, ax1, ax2):
        factor = ax1.factor * ax2.factor
        n_axes = ax1.n_axes + ax2.n_axes
        assert ax1.tag ==  ax2.tag == "ax"
        a1, a2 = ax1.value, ax2.value
        scalprod = a1.scalprod15(a2, None)
        if scalprod == 0:
            self.set_int(0)
        else:
            v = mul_axis_vector(a1, ax2.as_vector())
            self.set_vector(v, factor, n_axes)

    def set_mul_ax_vect(self, ax1, v):
        factor = ax1.factor * ax2.factor
        n_axes = ax1.n_axes + ax2.n_axes
        if ax1.tag == 'ax':
            v = mul_axis_vector(ax1.value, ax2.as_vector())
            self.set_vector(v, factor, n_axes)
        else:
            error_complicated()

    def as_vector(self, factor = 1):
        self.check()
        p = self.p
        factor = self.factor
        if self.tag == "ax":
            if p in [3, 15]:
                v = (self.value.v15).copy() % p
            else:
                v = self.value.in_space(MMV, p)
        elif self.tag == "v":
            v = self.value
        elif self.tag == "1":
            if factor == 0:
                return MMVector(p, 0)
            v = MMVector(p, 'U' )
            factor *= (p + 1) >> 2
        else:
            ERR = "Griess algebra product is to complicated"
            raise ValueError(ERR)
        return v * factor

    @staticmethod
    def n2_ax_to_factor(n2_ax, n_axes, p=0):
        ok = n2_ax > 0 and (n2_ax & (n2_ax -1)) == 0
        if ok:
            log_n2_ax =  n2_ax.bit_length() - 1
            ok = ok and log_n2_ax & 1
        if not ok:
            ERR = "Squared norm of an axis must be an odd power of 2"
            raise ValueError(ERR)
        shift = (log_n2_ax - 5) >> 1
        exp = shift * n_axes
        return pow(2, exp, p) if p else pow(2, exp)


    def as_vector_out(self, n2_ax = 8):
        if self.factor == 0:
            return MMV(self.p)
        factor = self.n2_ax_to_factor(n2_ax, self.n_axes, self.p)
        return self.as_vector(factor)



def Griess(a, b, **kwds):
    p = getattr(kwds, 'p', None)
    if p is None:
        p = auto_modulus((a,b))
        if p is None:
            p = 15
    result = GriessIntermediate((a, b), p)
    return result.as_vector_out(getattr(kwds, 'n', 32))
