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
    """Adjust a vector to modulus ``p``, if possible"""
    if p == vector.p:
        return vector
    elif vector.p % p == 0:
        return vector % p
    else:
        ERR = "Incompatible moduli of representation vectors"
        raise ValueError(ERR)


def auto_modulus(a):
    """Extract modulus from Griess algebra product expression

    Here ``a`` is an entry.
    An entry is either a tuple of entries or an atom.
    A atom is anything that is not a tuple.

    A atom ``a```with is an instance of class  ``MMVector`` has
    modulus ``MMVector.p``. Any other atom has moduli ``None``,
    meaning 'undefined'.

    The function searches recursively for a defined modulus ``p``.
    It returns ``p`` if this is the only defined modulus found,
    and ``None`` if no defined modulus has been found. It raises
    ValueError if different defined moduli have been found.
    """
    if isinstance(a, MMVector):
        return a.p
    if isinstance(a, tuple):
        p = None
        for x in a:
            p1 = auto_modulus(x)
            if p1 is None:
                p = p1
            elif p1 != p:
                ERR = "Incompatible moduli of representation vectors"
                raise ValueError(ERR)
    return None


def mul_axis_vector(axis, vector):
    """Return Griess product of a 2A axis by a vector

    ``axis`` must be an instance of class Axis, and ``vector``
    must be an instance of class ``MMVector``. The result is
    returned as an instance of class ``MMVector``.
    """
    g_data = axis._fast_g_mmdata()
    res = vector.copy()
    _fast_mul_v15(res, g_data, -1)
    mm_op_mul_std_axis(15, res.data)
    _fast_mul_v15(res, g_data)
    return res



class GriessIntermediate:
    TYPE_DICT = {
        "v": MMVector, "ax": Axis, "1": Integral 
    }
    COMPLEXITY = {
        "1": 0, "ax":1, "v":3
    }
    PAIR_COMPLEXITY = {
        "2A": 2, None:4
    }

    def __init__(self, a, p, lazy=False):
        self.p = p
        if isinstance(a, Axis):
            self.set_axis(a)
        elif isinstance(a, MMVector):
            self.set_vector(a)
        elif isinstance(a, Integral):
            self.set_int(a)
        elif isinstance(a, GriessIntermediate):
            self.set_instance(a)
        elif isinstance(a, tuple):
            if len(a) != 2:
                ERR = "Griess algebra product must have 2 factors"
                raise ValueError(ERR)
            a, b = a
            if not isinstance(a, GriessIntermediate):
                a = GriessIntermediate(a, p, lazy)
                a.check()
            if not isinstance(b, GriessIntermediate):
                b = GriessIntermediate(b, p, lazy)
                b.check()
            ca, cb = a.complexity(), b.complexity()
            if cb > ca:
                a, b = b, a
            if a.tag == "1":
                self.set_instance(b, a.factor)
            elif cb > 1:
                if ca > 2:
                    self.set_pair_lazy(a, b, lazy)
                if b.to_vector() == 0:
                    self.set_int(0)
                else:
                    v = self.set_mul_axis_vector(a, b)
                    self.set_vector(v, a.factor * b.factor,
                        a.n_axes + b.n_axes)
            elif a.tag == b.tag == "ax":
                self.set_mul_ax_ax(a, b)
            else:
                self.set_pair_lazy(a, b, lazy)
        else:
            ERR = "Cannot convert %s object to representation vector"
            raise ValueError(ERR % type(a))

    def complexity(self):
        if self.tag == "pair":
            return self.HINT_COMPLEXITY[self.hint]
        return self.COMPLEXITY[self.tag]

    @staticmethod
    def error_complicated():
        ERR = "Expression to complicated for Griess algebra"
        raise ValueError(ERR)

    def check(self):
        assert self.p & 1
        assert isinstance(self.factor, Integral)
        if self.tag == "pair":
            assert self.hint in self.PAIR_COMPLEXITY
            assert isinstance(self.value[0], GriessIntermediate)
            assert isinstance(self.value[1], GriessIntermediate)
        else:
            assert self.tag in self.TYPE_DICT, (self.tag)
            assert isinstance(self.value,  self.TYPE_DICT[self.tag])

    def set_int(self, value):
        self.tag = "1"
        self.factor = int(value)
        self.value = 1
        self.n_axes = 0
        self.check()
        return self.factor

    def set_instance(self, node, factor=1):
        assert isinstance(node, GriessIntermediate)
        assert self.p == node.p
        if factor * node.factor == 0:
            self.set_int(0)
        else:
            node.check()
            self.tag = node.tag
            self.factor = factor * node.factor
            self.value = node.value
            self.n_axes = node.n_axes
            self.check()

    def set_axis(self, axis, factor = 1):
        assert isinstance(axis, Axis)
        self.tag = "ax"
        self.factor = int(factor)
        self.value = axis
        self.n_axes = 1
        self.check()

    def set_vector(self, vector, factor = 1, n_axes = 0):
        self.n_axes = n_axes
        if vector == 0 or factor == 0:
            set_int(self, 0)
        else:
            self.tag = "v"
            self.factor = int(factor)
            self.value = adjust_vector_p(vector, self.p)
            self.check()

    def set_pair(self, node1, node2, hint = None):
        """Set this node to a the pair ``(node1, node2)``.

        This means that the node should compute the Griess algebra
        product  ``node1`` * ``node2``.

        Caution: This function modifies ``node1`` and ``node2``!

        ``node1`` and ``node2`` will have members ``factor = 1``
        and ``n_axes = 0``.

        An option ``hint`` for the evaluating the pair
        ``(node1, node2)`` can be given.

        Valid hints are:

        "2A" : This indicates a product of two 2A axes, such that the
        product of the corresponding involutions is in class 2A.

        None :  Nothing is known about the entries of a pair.
        """
        assert isinstance(node1, GriessIntermediate)
        assert isinstance(node2, GriessIntermediate)
        assert hint in self.PAIR_COMPLEXITY
        assert self.p == node1.p == node2.p
        self.tag = "pair"
        self.factor = node1.factor * node2.factor
        self.n_axes = node1.n_axes + node2.n_axes
        if self.factor == 0:
            return set_int(self, value)
        node1.factor = node2.factor = 1
        node1.n_axes = node2.n_axes = 0
        self.value = (node1, node2)
        self.hint = hint

    def set_pair_if_lazy(self, node1, node2, lazy):
        """Set this node to a the pair ``(node1, node2)``

        if ``lazy`` is True and fail otherwise.
        """
        if not lazy:
            self.error_complicated()
        set_pair(self, node1, node2)

    def set_mul_ax_ax(self, ax1, ax2):
        assert ax1.tag ==  ax2.tag == "ax"
        factor = ax1.factor * ax2.factor
        n_axes = ax1.n_axes + ax2.n_axes
        a1, a2 = ax1.value, ax2.value
        scalprod = a1.scalprod15(a2, None)
        if scalprod == 0 :
            self.set_int(0)
        elif scalprod == 1:
           self.set_pair(ax1, ax2, hint = '2A')
        else:
            if ax2.to_vector() == 0:
                self.set_int(0)
            else:
                v = mul_axis_vector(a1, ax2.value)
                self.set_vector(v, factor, n_axes)

    @staticmethod
    def set_mul_axis_vector(ax, v):
        """Return ``ax * v`` for an node ``ax`` and a vector ``v``.

        ``ax`` must of type ``GriessIntermediate``. Here ``ax.factor``
        is ignored.
        """
        if ax.tag == 'ax':
            return mul_axis_vector(ax.value, v)
        if ax.tag == 'pair':
            if ax.hint == "2A":
                ax1, ax2 = ax.value
                ax1, ax2 = ax1.value, ax2.value
                ax3 = product_axis_2A(ax1, ax2)
                v1, v2, v3 = [a * v for a in [ax1, ax2, ax3]]
                return 4 * (v1 + v2 - v3)
            else:
                self.error_complicated()
        else:
            self.error_complicated()

    def axis_to_vector(self, axis):
        """Convert an axis to a vector"""
        p = self.p
        if p in [3, 15]:
            return (axis.v15).copy() % p
        else:
            return axis.in_space(MMV, p)

    def pair_to_vector(self):
        assert self.tag == "pair"
        a, b = self.value[0], self.value[1]
        if self.hint == "2A":
            assert a.tag == b.tag == "ax"
            v = mul_axis_vector(a.value, self.axis_to_vector(b.value))
            self.set_vector(v, self.factor, self.n_axes)
        else:
            a, b = self.value
            res = GriessIntermediate((a,b), self.p, lazy=False)
            res.to_vector()
            self.set_instance(res, self.factor * res.factor,
                self.n_axes * res.factor)
        return self


    def to_vector(self):
        """Change intermedate node to a vector

        The function changes the node so that it is either a
        vector (with tag "v") or zero.

        It either returns the modified node itself or zero.
        """
        self.check()
        factor = self.factor
        if factor == 0:
            return self.set_int(0)
        if self.tag == "v":
            return self
        p = self.p
        if self.tag == "ax":
            v = self.axis_to_vector(self.value)
        elif self.tag == "1":
            v = MMVector(p, 'U' )
            factor *= (p + 1) >> 2
        elif self.tag == "pair":
            return self.pair_to_vector()
        else:
            ERR = "Griess algebra product is to complicated"
            raise ValueError(ERR)
        self.set_vector(v, factor)
        return self

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


    def vector_out(self, n2_ax = 8):
        if self.to_vector() == 0:
            return MMV(self.p)
        factor = self.n2_ax_to_factor(n2_ax, self.n_axes, self.p)
        return (factor * self.factor) * self.value



def Griess(a, b, **kwds):
    p = getattr(kwds, 'p', None)
    if p is None:
        p = auto_modulus((a,b))
        if p is None:
            p = 15
    result = GriessIntermediate((a, b), p)
    return result.vector_out(getattr(kwds, 'n', 32))
