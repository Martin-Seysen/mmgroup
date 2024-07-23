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

from numbers import Integral
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
from mmgroup.structures.xsp2_co1 import DICT_INVOLUTION_G_x0
from mmgroup.tests.test_involutions.test_involution_invariants import INVOLUTION_SAMPLES

#print("Hello world")


def hexx(x):
    return hex(x) if isinstance(x, Integral) else x


#######################################################################
# Auxiliary class for function xsp2co1_elem_conjugate_involution_Gx0
#######################################################################


class N_x0_involution:
    r"""Model the mapping :math:`G_{x0} \rightarrow N_{x0}`

    An instance of this class is initialized with an instance
    :math:`g` of the group :math:`G_{x0}`. Actually, :math:`g`
    must be in the subgroup  :math:`N_{x0}` of  :math:`G_{x0}`;
    otherwise the function fails.

    Method ``transform`` of this class conjugates the element
    :math:`g` with an element :math:`a` of the supergroup
    :math:`N_{0}` of  :math:`N_{x0}`. All subsequent 
    calculations are done in the group :math:`N_{0}`.

    The (conjugated) element :math:`g` is kept in attribute
    ``gn0`` and the accumulated transformation :math:`t`
    conjugating that element is kept in attribute ``tf_n0``. 
    Both elements are stored as alements of the group 
    :math:`N_{0}`, as described in section   
    **Implementing the group N_0** in
    **The C interface of the mmgroup project**.

    Method ``get_q`` returns the conjugated element :math:`g`
    if it is in the subgroup :math:`Q_{x0}`; otherwise the
    method fails.

    Method ``out`` returns the accumulated transformation
    :math:`t` as a word in the generors of  :math:`N_{0}`.
    """  
    IND_Y = 1
    IND_X = 2
    
    def __init__(self, elem):
        r"""Initialize instance with an element of :math:`G_{x0}`

        Here the element :math:`g` of the group :math:`G_{x0}` must
        be given by parameter ``elem`` in **G_x0 representation** 
        as described in section
        **C functions dealing with the subgroup G_x0 of the monster**
        in  **The C interface of the mmgroup project**. So parameter
        ``elem`` should be a ``numpy`` array of legnth 26 and 
        dtype ``numpy.uint64``.

        The function succeeds if  :math:`g` is of order at most 2, 
        and :math:`g` is in the subgroup  :math:`N_{x0}` of
        :math:`G_{x0}`, and if the image  of :math:`g` in the 
        Mathieu  group :math:`M_{24}` (which is a factor group of  
        :math:`N_{x0}`) is one.  Otherwise the function fails.
        """ 
        self.gn0 =  np.zeros(5, dtype = np.uint32)
        self.tf_n0 =  np.zeros(5, dtype = np.uint32)
        # Convert g to an element of N_0, store it in self.gn0.
        # Fail if g is not in N_x0.
        if (xsp2co1_elem_to_N0(elem, self.gn0)) < 0:
            raise ValueError("xsp2co1_elem_from_N0 failed")
        # Fail if g is not in N_x0 or its image in M_24 is not 1
        if (self.gn0[0] | self.gn0[4]):
            raise ValueError("Bad N_x0 element")
        # Put t = g**2
        mm_group_n_mul_element(self.gn0, self.gn0, self.tf_n0)
        # Fail if g**2 != 1. As a side effect, this puts t = 1.
        if mm_group_n_reduce_element(self.tf_n0):
            raise ValueError("Element is not an involution")
        mm_group_n_reduce_element(self.gn0)

    def get_xy(self, index):
        r"""Return x-part or y-part of the internal element :math:`g`

        Depending on parameter ``index`` the method returns the
        x-part or y-part of the (conjugated) element :math:`g` 
        stored inside the instance.  That part is an element of 
        the Parker loop. It is returned modulo its sign and modulo
        :math:`\Omega` as an 11-bit integer.

        We return the x-part if ``index = self.IND_X`` and
        the  y-part if ``index = self.IND_Y``.
        """
        return self.gn0[index] & 0x7ff

    def get_q(self):
        r"""Return internal element :math:`g`

        The method ``get_q`` returns the (conjugated) element 
        :math:`g` if it is in the subgroup :math:`Q_{x0}`; 
        otherwise the method fails.

        In case of success the conjugated element :math:`g` is
        returned as an integer in **Leech lattice encoding**.
        """
        assert self.gn0[0] | self.gn0[1] | self.gn0[4]  == 0
        q = self.gn0[2]
        return (q << 12) ^ mat24.ploop_theta(q) ^ self.gn0[3]

    def transform(self, a):
        r"""Conjugate the element :math:`g` 

        The function conjugates the element :math:`g` stored
        inside this instance with a transformation  :math:`a`.

        Here  :math:`a` must be given as a word of generators of
        :math:`N_{0}` in the parameter ``a``. That parameter must 
        be a ``numpy`` array of dtype = ``np.uint32``.
        """
        mm_group_n_conj_word_scan(self.gn0,  a, len(a))
        mm_group_n_mul_word_scan(self.tf_n0, a, len(a))
        mm_group_n_reduce_element(self.gn0)

    def out(self, a):
        r"""Return the accumulated transformation :math:`t`

        The method writes the accumulated transformation
        :math:`t` (that conjugates the element :math:`g`) into 
        the array ``a`` as a word of generators of  :math:`N_{0}`.
        The method returns the length of the data in the array
        ``a``. Here ``a`` must be a ``numpy`` array of
        length at least 5 and of ``dtype = np.uint32``.

        The result ``a`` is presented such that the inverse of
        ``a`` has a simple form.
        """
        inv = np.zeros(5, dtype = np.uint32)
        mm_group_n_inv_element(self.tf_n0, inv)
        length = mm_group_n_to_word(inv, a)
        mm_group_invert_word(a, length)
        return length

    def display_transformed(self):
        r"""Display the conjugated element :math:`g`"""
        a = np.zeros(5, dtype = np.uint32)
        length = mm_group_n_to_word(self.gn0, a)
        print(MM0('a', a[:length]).reduce())
    
  
#########################################################################
# Python implementation of function xsp2co1_elem_conjugate_involution_Gx0
#########################################################################


def conj_involution_Gx0_type2(vx, guide, a):
    r"""Auxiliary function for xsp2co1_elem_conjugate_involution_in_Gx0_py

    This function deals with a specific case of the function
    ``xsp2co1_elem_conjugate_involution_in_Gx0_py``.
    """
    # Here ``vx`` has type 2 in Leech lattice mod 2.
    # The ``guide`` is useful if ``guide`` is of type 4
    # and ``guide ^ vx`` is of type 2.
    len_a = 0
    if  (gen_leech2_type(guide) == 4 and 
        gen_leech2_type(guide ^ vx)) == 2:
        # If the ``guide`` useful, transform it to the
        # standard frame in the Leech lattice mod 2.
        len_a = gen_leech2_reduce_type4(guide, a);
        vx =  gen_leech2_op_word(vx, a, len_a) 
    # Transform ``vx`` to the standard cocode word
    # :math:`x_\beta`, where :math:`\beta` is the
    # cocode word [2,3].
    len_a2 = gen_leech2_reduce_type2(vx, a[len_a:]) 
    assert len_a2 >= 0
    vx =  gen_leech2_op_word(vx, a[len_a:], len_a2)       
    len_a += len_a2
    # Correct a negative sign of ``vx`` by transforming with
    # an :math:`x_d, \langle d, \beta \rangle = 1`.
    if vx & 0x1000000:
        a[len_a] = 0xB0000200
        len_a += 1
    return len_a

def xsp2co1_elem_conjugate_involution_in_Gx0_py(elem, guide, a):
    r""" Map an involution in :math:G_{x0}` to a standard form.

    Here the element :math:`g` of the group :math:`G_{x0}` must
    be given by parameter ``elem`` in **G_x0 representation** 
    as described in section
    **C functions dealing with the subgroup G_x0 of the monster**
    in  **The C interface of the mmgroup project**. So parameter
    ``elem`` should be a ``numpy`` array of legnth 26 and 
    dtype ``numpy.uint64``.

    The function succeeds if  :math:`g` is of order at most 2, 
    Then it computes an element :math:`a` in :math:`G_{x0}` such
    that :math:`h = a^{-1} g a` is a (fixed) representative of 
    the class of :math:`g` in the group :math:`G_{x0}`. The chosen
    representatives of the involution classes in :math:`G_{x0}` 
    are described in the documentation of the C function
    ``xsp2co1_elem_conjugate_involution_Gx0`` in file
    ``xsp2co1_traces.c``.

    The element :math:`a` is stored in the array ``a`` as a word 
    of generators of  the monster group. In case of success the
    function returns  ``len(a)``, where ``len(a)`` is the length  
    of the data in the array ``a`` The function  raises
    ``ValueError`` in case of failure, e.g. if :math:`g` has
    order greater than 2. The array ``a`` must have length at
    least :math:`10`.

    Parameter ``guide`` should usually be zero. Otherwise ``guide`` 
    should be a type-4 vector :math:`v_4` in the Leech lattice mod 2
    (in **Leech lattice encoding**) such that the two  conditions 
    :math:`h = a^{-1} g a` and :math:`v_4 \cdot a = \Omega` can 
    both be achieved. Then we compute an element :math:`a` 
    satisfying these two conditions. Otherwise parameter ``guide`` 
    is ignored. 
    """
    # Find a frame ``v4`` in the Leech lattice mod 2 such that
    # a transformation mapping ``v4`` to the standard frame
    # :math:`\Omega` also changes ``g`` to a simpler form.
    # Therefore we use the C function 
    # ``xsp2co1_elem_find_type4`` which also deals with ``guide``.
    v4 = xsp2co1_elem_find_type4(elem, guide)
    if (v4 < 0):
        # Abort if ``xsp2co1_elem_find_type4`` has failed
        raise ValueError("xsp2co1_elem_find_type4 failed")
    # Store a transformation ``a_1`` that maps ``v4`` to 
    # :math:`\Omega` in the part ``a[:len_a]`` of the output 
    # array ``a``.
    len_a = gen_leech2_reduce_type4(v4, a);
    assert len_a >= 0
    # Put ``elem1 = elem ** a_1``
    elem1 = np.zeros(26, dtype = np.uint64)
    np.copyto(elem1, elem)
    if (xsp2co1_conjugate_elem(elem1, a, len_a)) < 0:
        raise ValueError("xsp2co1_conjugate_elem failed")
    # Convert ``elem1`` to an element of the group :\math:`N_{x0}`
    # and store the result in the instance ``invol`` of class
    # ``N_x0_involution``. Abort if this conversion fails.
    invol = N_x0_involution(elem1)

    # Now ``invol`` is equal to :math:`x_d x_\delta, y_d x_\delta`,
    # or  :math:`x_d y_d x_\delta`, where :math:`d` is in the
    # Parker loop and :math:`delta` is in the Golay cocode.
    # The function fails if this is not the case for ``invol``.
    if invol.get_xy(invol.IND_Y):
        # Here ``invol`` is of shape :math:`y_d x_\delta` or
        # :math:`x_d y_d x_\delta`.
        b = np.zeros(3, dtype = np.uint32)
        if invol.get_xy(invol.IND_X):
            # If ``invol`` has shape :math:`x_d y_d x_\delta` then
            # transform ``invol`` with :math:`x_\epsilon`, where
            # math:`\epsilon` is an odd cocode word. This changes
            # the shape of ``invol`` to :math:`y_d x_\delta`
            b[0] = 0x10000800
            invol.transform(b[:1])
        # We temporarily transfrom ``invol`` with the triality
        # element :math:`\tau^2`. This changes the shape of
        # ``invol`` to  :math:`x_d x_\delta.`
        b[0] = 0x50000002   # (this is t**2)
        invol.transform(b[:1])
        # So we may convert ``invol`` to an element ``vy``
        # of :math:`Q_{x0}`.
        vy = invol.get_q()
        # Compute a transformation ``b`` in ``N_x0`` that reduces
        # ``vy`` modulo the group :math:`N_{x0} to a standard form.
        # Therefore we use the C function ``gen_leech2_reduce_n``.
        # Note that this function ignores the sign of `vy`` and that
        # it computes an even element ``b`` of  :math:`N_{x0}.
        gen_leech2_reduce_n(vy, b)
        invol.transform(b[:3])
        # Reverse the transfromation of ``invol`` with :math:`\tau^2`. 
        # This changes the shape of ``invol`` back
        # to  :math:`x_d x_\delta`.
        b[0] = 0x50000001   # (this is t**1)
        invol.transform(b[:1])
        # Append the transformation accumulated in the object
        # ``invol`` to the transformation ``a_1`` and return result.
        len_a2 = invol.out(a[len_a:])
        return len_a + len_a2
    else:
        # Here ``invol`` is of shape :math:`x_d x_\delta`. The
        # transformation ``a_1``computed by ``xsp2co1_conjugate_elem``
        # should be trivial here.
        if guide == 0:
            assert len_a == 0
        # So we may obtain ``invol`` as an element ``vx``.
        vx = invol.get_q()
        # Let ``t`` be the type of ``vx`` in the Leech lattice mod 2.
        t = gen_leech2_type(vx)
        if t == 0:
            # ``vx`` has type 0 in the Leech lattice mod 2
            if  (gen_leech2_type(guide)) == 4:
                # If a suitable ``guide`` is given as a vector
                # of the type 4 in the Leech lattice mod 2 then
                # transform that ``guide`` to the standard frame.
                len_a = gen_leech2_reduce_type4(guide, a);
                assert len_a >= 0
                return len_a
            else:
                # Otherwise perform no action
                return 0
        elif t == 4:
            # ``vx`` has type 4 in the Leech lattice mod 2.
            # Transform ``vx`` to the standard frame
            len_a = gen_leech2_reduce_type4(vx, a);
            assert len_a >= 0
            vx = gen_leech2_op_word(vx, a, len_a)
            # If ``vx`` transforms to :math:`-\Omega` then correct 
            # the sign by transforming with an odd cocode element.
            if vx & 0x1000000:
                 a[len_a] = 0x90000800
                 len_a += 1
        elif t == 2:
            # ``vx`` has type 2 in the Leech lattice mod 2. Then
            # we execute function ``conj_involution_Gx0_type2``.
            len_a = conj_involution_Gx0_type2(vx, guide, a)
        else:
            # Any other type of ``vx`` is illegal
            raise ValueError("Not an involution")
        # The length of the output word my be atmost ten.
        assert len_a  <= 10
        return len_a

#######################################################################
# Pythonic version of function xsp2co1_elem_conjugate_involution_Gx0
#######################################################################

 
                
def conjugate_involution_in_Gx0(elem, guide = 0):
    r"""Pythonic function ``xsp2co1_elem_conjugate_involution_in_Gx0_py``

    Input is as in function 
    ``xsp2co1_elem_conjugate_involution_in_Gx0_py``, but here
    parameter ``elem`` is anything that can be converted to an 
    instance of class ``Xsp2_Co1``; and the function returns a pair
    ``(t, elem1)`` Here ``t`` is the transformation computed by 
    function ``xsp2co1_elem_conjugate_involution_in_Gx0_py``
    returned as an instance of class ``MM0``; and ``elem1`` is
    ``elem ** t`` returned as an instance of class ``Xsp2_Co1``. 
    """
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


REVERSE_DICT_INVOLUTION_G_x0 = {}
for key, value in DICT_INVOLUTION_G_x0.items():
    REVERSE_DICT_INVOLUTION_G_x0[value] = key




def xsp2co1_elem_conjugate_involution_Gx0_C(elem, guide = 0):
    r"""Wrapper for the corresponding C function

    Input is as in function ``conjugate_involution_in_Gx0``; but 
    instead of a transformation ``t`` the function  returns a 
    pair ``(iclass, t)``. Here ``iclass`` is the class of the 
    involution ``elem`` in :math:`G_{x0}`, as returned by 
    the C function ``xsp2co1_elem_involution_class``.
    """
    elem =  Xsp2_Co1(elem)
    iclass, g = elem.conjugate_involution_G_x0(guide, MM0) 
    hex_iclass = REVERSE_DICT_INVOLUTION_G_x0[iclass]
    return hex_iclass, g


#######################################################################
# Dictionary of rpresentatives of classes of involutions in G_x0
#######################################################################


STD_REP = None

def get_std_rep():
    r"""Return dictionary for involution classes in :math:`G_{x0}`

    The function returns a dictionary that maps the name of a class 
    of an involution in :math:`G_{x0}` to the representative chosen 
    for that class.

    The names of the classes are given as integers as returned by
    the C function ``xsp2co1_elem_involution_class``.

    Representatives of the classes are obtained by applying function
    ``conjugate_involution_in_Gx0`` to all involution samples stored 
    in the list ``INVOLUTION_SAMPLES`` in the python module
    ``mmgroup.tests.test_involutions.involution_samples".
    """
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
    r"""Test data for function ``xsp2co1_elem_conjugate_involution_Gx0``

    The function yields pairs ``(g, guide)``, where ``g`` is an 
    element of :math:`G_{x0}` returned as an instance of class 
    ``Xsp2_Co1``. The value ``guide`` is either 0 or an integer 
    that can be used as parameter ``guide`` in function 
    ``xsp2co1_elem_conjugate_involution_Gx0``.

    Test data are created by conjugating all involutions in
    the dictionary ``INVOLUTION_SAMPLES`` with random elements
    of :math:`G_{x0}`.

    If a useful nonzero ``guide`` is to be computed then the
    involutions are taken from the dictionary returned by 
    function ``get_std_rep()`` is taken instead.
    """ 
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
                if i & 1:
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
    r"""Test C function ``xsp2co1_elem_conjugate_involution_Gx0``
  
    This function tests that C function with test data taken 
    from the generator function ``conjugate_testdata``.

    It checks the C implementation against the corresponding
    python implementation ``conjugate_involution_in_Gx0``.
    It verifies that the transformed input element of 
    :math:`G_{x0}` is equal to the representative of the class 
    of the input obtained from the dictionary ``get_std_rep()``.
    """
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
            print("Involution class %s:\n g =" % hexx(iclass), g)
            if c_class != iclass:
                print("Class computed by C function", hexx(c_class))
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
    """Return list of the bis set in the integer x & 0xffffff"""
    w, l =  mat24.vect_to_bit_list(x)
    return l[:w]

def display_q_xy(g):
    """Auxiliary function for function ``display_g`` """
    g_n0 = N_x0_involution(Xsp2_Co1(g)._data)
    try:
        print("  in N_x0:", hex(g_n0.get_q()))
    except (ValueError, AssertionError):
        b = np.array([0x50000002], dtype = np.uint32)
        g_n0.transform(b)
        print("  in N_y0:",  hex(g_n0.get_q()))
       

       
def display_g(g):
    r"""Display an element of ``N_{x0}`` in readable form

    Here input ``g`` should be an element of  ``N_{x0}``
    that can be converted to an instance of class ``MM0``.
    """
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

