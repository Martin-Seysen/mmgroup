r"""This module implements Norton's generators of the Monster.

Norton :cite:`Nor02` has given a presentation of the Monster group
that greatly simplifies a mapping from the *projecive plane*
presentation of the Monster to the representation of the Monster
in :cite:`Gri82` and :cite:`Con85`. Our implemention of the 
Monster is based on the representation in :cite:`Con85`. So we may
use the presentation in :cite:`Nor02` to construct a homomorphism
from the  *projecive plane* representation of the Monster into our
implementation of the Monster that can be computed in practice.

"""
import os
import sys
from numbers import Integral
from random import randint, sample, choice

import numpy as np

from mmgroup import XLeech2, Cocode, PLoop, MM
from mmgroup.clifford12 import xsp2co1_elem_from_mapping
from mmgroup.bimm import inc_p3
from mmgroup.bimm.inc_p3 import P3_incidences
from mmgroup.bimm.inc_p3 import P3_point_set_type
from mmgroup.bimm.inc_p3 import AutP3
from mmgroup.bimm.inc_p3 import invert_perm_P3
from mmgroup.bimm.inc_p3 import mul_perm_P3
from mmgroup.bimm.inc_p3 import check_perm_P3

#####################################################################
# Constructing 'Points' and 'Stars' in the Monster
#####################################################################


precomputation_pending = True

# The following dictionary maps some of the points of :math:`P3`
# to the MOG as in the diagram on page 80 in :cite:`Nor02`.
DICT_POINT_MOG = {
    2:13,  6:17,  7:21,  
    8:14, 11:18,  4:22,
   12:15, 10:19,  5:23,  
}

# The following dictionary maps some pairs of points of :math:`P3`
# to the MOG as in the diagram on page 80 in :cite:`Nor02`.
# Entry `i:j` means that the pair :math:`P_0 P_i` is mapped to
# entry ``j`` of the  MOG.
DICT_POINT_MOG_COLUMN = {1:0, 3:4, 9:8}





def make_P(x = 0, delta = 0):
    r"""Create an element of the group :math:`2^{1+24}`

    Given parameters ``d`` and ``delta``,  the functions
    returns the element  :math:`x_d x_{\delta}` of 
    the group :math:`Q_{x0}` of structure :math:`2^{1+24}`
    as an instance of class ``XLeech2``.
    """
    return XLeech2(PLoop(x), Cocode(delta))
    #return MM([('x', PLoop(x)), ('d', Cocode(delta))])



# Dictionary of images of 'points' :math:`P_0 P_i` in :math:`Q_{x0}`
P0_DICT = {}

def compute_P0(x):
    r"""Return image of 'point' :math:`P_0 P_i` in :math:`Q_{x0}`

    Norton cite:`Nor02` defines a mapping of the 'points' :math:`P_i`
    of :math:`P3` into the subgroup :math:`Q_{x0}` of structure 
    :math:`2^{1+24}` of the Monster. More precisely, he defines
    a mapping from the products of an even number of 'points' into
    :math:`Q_{x0}` . Note that all 'points'  commute; so if suffices
    to store the products :math:`P_0 P_i, 0 \leq i < 13` in a
    dictionary.

    According to cite:`Nor02` these products are mapped to the
    Leech lattice modulo 2 as follows:

    We have :math:`P_0^2 = 1`.  In case :math:`i = 1,3,9` we have
    :math:`P_0 P_i = x_{o(i)}`, where :math:`o(i)` is the octad
    with entries set in row  0 of the MOG and in column
    ``DICT_POINT_MOG_COLUMN[i]/4`` of the MOG, excluding the 
    intersection of row 0 with that column.

    For the other cases of :math:`i` we have 
    :math:`P_0 P_i = x_{\delta(i)}`, where :math:`x_{\delta(i)}`
    is the Golay cocode word of weight 1 with a single entry at
    position ``DICT_POINT_MOG[i]``.

    This defines a mapping from the pairs :math:`P_0 P_i` to
    :math:`Q_{x0}` up to sign. Here we simple map all these pairs
    to elements of :math:`Q_{x0}` that by definition are considerd 
    as *positive* in our construction of :math:`Q_{x0}`.    
    """
    if x == 0: return XLeech2()
    if x in (1, 3, 9):
        c = DICT_POINT_MOG_COLUMN[x]
        octad = [0, 4, 8, 12, 16, 20, c+1, c+2, c+3]
        return make_P(octad)
    return make_P(0, [DICT_POINT_MOG[x]])


def PointP3(x):
    r"""Map product of 'points' of :math:`P3` into subgroup of Monster

    Let parameter ``x`` be a list of 'points' in :math:`P3` 
    of even length. Here an entry of that list may be anything
    that is accepted by the constructor of class ``P3_node``.

    The function compute the image of the product of the 
    entries of that list in the subgroup :math:`Q_{x0}` of the
    Monster under the mapping defined in function ``compute_P0``.
    The function returns that image as an instance of 
    class ``XLeech2``.
    """
    if precomputation_pending:
        precompute_all()
    assert len(x) & 1 == 0
    p = XLeech2() 
    for x_i in  x:
       p *= P0_DICT[x_i % 13]
    return p



# :cite:`Nor02` defines so-called 'stars', see function ``StarP3``.
# The following dictionary maps most of the 'stars' to entries
# of the MOG as in the diagram on page 80 in :cite:`Nor02`.
DICT_LINE_MOG = {
    1:12,  3:16,  9:20,
   12: 1, 11: 5,  7: 9,
    8: 2,  6: 6,  5:10,
    2: 3, 10: 7,  4:11, 
}

# Dictionary of images of 'stars' :math:`P_i^*` in :math:`Q_{x0}`
PSTAR_DICT = {}



def compute_StarP3(i, check = False):
    r"""Return image of 'star' :math:`P_i^*` in :math:`Q_{x0}`

    Norton cite:`Nor02` defines a mapping of the so-called 'stars' 
    :math:`P_i^*` to the subgroup :math:`Q_{x0}` of structure 
    :math:`2^{1+24}` of the Monster.  See also function ``StarP3`` 
    for a brief description of the 'stars'. On input ``i`` the 
    function returns the image of the 'star' :math:`P_i^*` as an 
    instance of class ``XLeech2``.

    According to cite:`Nor02` these 'stars' commute and they are 
    mapped to the Leech lattice modulo 2 as follows:

    We have :math:`P_0^* = x_{\Omega} x_{\omega}`, where 
    :math:`\Omega` is the (positive) element of the Parker loop
    with image :math:`(1,\ldots,1)` in the Golay code, and 
    :math:`\omega` is the cocode word of weight 4 corresponding
    to an (arbitrary) column of the MOG. 

    In case :math:`i = 1,3,9` we have 
    :math:`P_i^* = x_{\omega} x_{\delta(i)}`, where 
    :math:`x_{\delta(i)}` is the Golay cocode word of weight 2
    containing the entries in positions :math:`0, 4, 8` different
    from the value ``DICT_POINT_MOG_COLUMN[i]``.

    In the other cases let :math:`\delta(i)` be the cocode word 
    of weight 2 containing  the entries  at positions
    ``DICT_POINT_MOG[j]`` for all  ``j != i`` such that
    ``DICT_POINT_MOG[j]`` is in the same column as
    ``DICT_POINT_MOG[i]``. Let :math:`o(i)` be the octad
    containing the entries at positions :math:`0, 4, 8`,
    ``DICT_POINT_MOG[i]``, and the positions 
    ``DICT_LINE_MOG[l]`` for all lines ``l`` incident with
    point :math:`i`.

    This defines a mapping from the stars :math:`P_i^*` to
    :math:`Q_{x0}` up to sign. Here we simple map all these 'stars'
    to elements of :math:`Q_{x0}` that by definition are considerd 
    as *positive* in our construction of :math:`Q_{x0}` .    
    """
    if i == 0:
        return make_P(~PLoop(), [0,1,2,3])
    if i in (1,3,9):
        return make_P(0,  [1,2,3,4,8, DICT_POINT_MOG_COLUMN[i]])
    lines = [DICT_LINE_MOG[y.ord % 13] for y in P3_incidences(i)]
    octad = [0, 4, 8, DICT_POINT_MOG[i]] + lines
    if check:
         from mmgroup.mat24 import syndrome
         assert syndrome(sum(1 << k for k in set(octad))) == 0
    point =  DICT_POINT_MOG[i]
    col = point & (-4)
    cocode = [j for j in range(col+1, col+4) if j != point]
    return make_P(octad, cocode)




def StarP3(x):
    r"""Map product of 'stars' of :math:`P3` into subgroup of Monster

    :cite:`Nor02` defines so-called 'stars' :math:`P_i^*`as words in 
    the generators of ``P2`` (modulo the defining relations of the 
    Bimonster). These 'stars' are numbered from 0 to 12 in th same 
    way as the 'points'. :cite:`Nor02` also defines a mapping from 
    these 'stars' into the subgroup :math:`Q_{x0}` of the Monster.
     
    Let parameter ``x`` be a list integers encoding a sequence of 
    'stars'`. The function computes the image of the product of these 
    'stars'`. It function returns that image as an instance of 
    class ``XLeech2``.

    An integer 'x' is interpreted as a list of lenght 1.
    """
    if precomputation_pending:
        precompute_all()
    if isinstance(x, Integral):
        return PSTAR_DICT[x % 13]
    p = XLeech2() 
    for x_i in  x:
        p *= PSTAR_DICT[x_i % 13]
    return p 





def precompute_all():
    r"""Perform all precomputations required for this module"""
    # Do precompution on demand only. Otherwise Sphinx will fail.
    global precomputation_pending
    if precomputation_pending:
        global P0_DICT, PSTAR_DICT
        for x in range(13):
            P0_DICT[x] = compute_P0(x)
        for x in range(13):
            PSTAR_DICT[x] = compute_StarP3(x)
        precomputation_pending = False



    

#####################################################################
# class Precomputed_AutP3 
#####################################################################






def MM_from_perm(perm, verbose = 0):
    r"""Map an element of ``AutP3`` into the Monster group

    Here parameter ``perm`` is a list of 13 integers
    ``0 <= perm[i] < 13`` such that the mapping ``i -> perm[i]`` is
    an automorphism ``a`` of points of the projective plane ``P3``.

    The function maps automorphsm ``a`` into the subgroup ``G_x0``
    of the Monster. Therefore it uses the operation of ``a`` on the
    points and stars for constructing the iamge ``g`` of ``a`` in 
    ``G_x0``.  Here ``a`` operates (by conjugation) on the 'points' 
    and 'stars', which are elements of the extraspecial subgroup 
    ``Q_x0`` of structure :math:`2^{1+24}` of ``G_x0``.

    Note that this operation of ``a`` determines the element ``g``
    of the subgroup ``G_x0`` of the Monster up to sign only. We use
    the C function ``xsp2co1_elem_from_mapping``  for constructing
    one of the elements ``+g`` or ``-g``. At most one of these two
    elements may have odd order, and we compute the element of odd
    order if such an element exists. Thus, if ``a`` has odd order
    then the returned element ``g`` is the correct image of ``a``.

    More precisely, we return a triple ``(g, order, special)``,
    where ``g`` is the image of ``a``  in ``G_x0`` given as an
    instance of class ``MM``.  Here ``order`` is the order of ``g``,
    which is odd if possible. The part ``special`` of the result
    is ``True`` if the function can distinguish between ``g`` and 
    ``-g`` by computing the orders and certain characters of these 
    two elements in ``G_x0``.
    
    Remark:

    If ``MM_from_perm(a1)`` returns ``(g1, order, True)``, and the 
    automorphism ``a2`` of ``AutP3`` is conjugate to ``a1`` then
    ``MM_from_perm(a2)`` returns ``(g2, order, True)``; and it is
    guaranteed that ``g2`` is conjugate to ``g1``, but not to ``-g1``
    in ``G_x0``. So once we know the 'right' image ``+g1`` or ``-g1``
    of ``a1``, we can also deduce the  'right' image ``+g2`` or
    ``-g2`` of ``a2``.

    The group ``AutP3`` is isomorphic to :math:`L_3(3)` in ATLAS
    notation. From the ATLAS we see that for each even order the
    group :math:`L_3(3)` has at most on conjugacy class with that
    order. So if ``MM_from_perm(a)`` returns ``(g, order, True)``
    for an even ``order``, and we have found the 'right' image
    ``g`` or ``-g`` of ``a``, then we can easily find the 'right'
    images of all other elements of ``AutP3`` of the same order.
   
    If cannot find the 'right' image ``+-g`` of an element ``a`` of
    ``AutP3`` that way, then we may always decompose ``a`` into a 
    (random) product of elements of odd order.
    """
    a_src = np.zeros(24, dtype = np.uint32)
    a_dest = np.zeros(24, dtype = np.uint32)
    a = np.zeros(10, dtype = np.uint32)
    pi0 = perm[0]
    for i in range(12):
       d_src = [0,i+1]
       d_dest = [pi0, perm[i+1]]
       a_src[i] = PointP3(d_src).ord
       a_src[i + 12] = StarP3(d_src).ord
       a_dest[i] = PointP3(d_dest).ord
       a_dest[i + 12] = StarP3(d_dest).ord
    res = xsp2co1_elem_from_mapping(a_src, a_dest, a)
    if res < 0:
        err = "xsp2co1_elem_from_mapping returns %d"
        raise ValueError(err % res)
    g = MM('a', a[:res & 0xff])
    order = (res >> 8) & 0xff
    special = (res >> 16) & 1
    if verbose and not special:
        print("Ambiguous element of order %d found" % order)
    return g, order, special





class Precomputed_AutP3:
    r"""Auxiliary class for storing images of class ``AutP3`` in class ``MM``

    This class contains class methods only. Its purpose is to compute
    a (fixed) mapping of the automorphism group of the projectve plane
    ``P3`` into the subgroup :math:`G_{x0}` of 
    structure :math:`2^{1+24}.\mabox{Co}_1` of the Monster. 

    Method ``as_MM`` of ths class computes this mapping. Some of
    the images of this mapping will be remembered for reuse.
    """
    # Store tranversal of a subgroup, see method ``_split_transveral``
    MAX_IND = 2*13*13  # length of the following array ``transversal``
    transversal = np.zeros((MAX_IND, 27), dtype = np.uint8)
    # Enter the neutral element into cls.transversal[1]
    transversal[1] = list(range(13)) * 2  + [1]

    # Record the orders of elements of the group ``AutPL`` for which
    # the sign conflicts (in ``G_x0``) can be resolved
    good_orders = {    # dict of 'good' orders where elements can be
                       # distinguished from their negatives
        1:1, 3:1, 13:1
    }
    bad_orders = {}    # dict of 'bad' orders where elements cannot
                       # be distinguished from their negatives
    # Record also some statistices
    n_splits = 0       # No of calls to method _split_into_good_orders
    n_split_trials = 0 # No of trials in method _split_into_good_orders

    # Store the images of the elements in ``cls.transversal``
    # We have to store at most 13*12 + 48 images in array data_MM
    data_MM = np.zeros((1+192, 10), dtype = np.uint32)  
    # cls.data_MM[ind_MM[i]] stores image of cls.transversal[i]  
    # cls.ind_MM[i] = 0  means that the image has not yet been computed
    ind_MM = np.zeros(MAX_IND, dtype = np.uint8)
    # cls.num_MM is the number of entries of cls.data_MM used.
    # Here we deliberately waste entry cls.data_MM[0], and we store
    # the image of the neutral element in cls.data_MM[1].
    ind_MM[1], num_MM = 1, 2


    @classmethod
    def _split_transveral(cls, g):
        r"""Split inctance ``g`` of class ``AutP3`` into two factors

        The method splits the element ``g`` of ``AutP3`` into a product
        ``f1 * f2`` where ``f1`` fixes the points 0 and 1, and ``f2`` 
        is in a transversal of the group fixing these two points. 

        Here the tansversal has size 156 and the subgroup fixing points
        0 and 1 has size 36. The functions in this class will store the
        156 + 36 elements of ``AutP3`` obtained that way, and also the
        images of these elements in the Monster, in internal arrays.

        The function return a pair  ``(h1, h2)``, where ``h1, h2`` are 
        indices pointing to the array ``cls.transversal``, such that
        ``cls.transversal[h1, :13]`` and ``cls.transversal[h1, :13]``
        are the permutations of the points ``P_i`` corresponding 
        to ``f1`` and ``f2``.

        We now describe the structure of ``cls.transversal``

        For ``0 <= i, j < 13`` entry ``13*i+j`` stores an element
        of the transversal that maps the points ``0`` and ``1`` to
        the points ``i`` and ``j``. 
 
        For ``0 <= i, j < 13`` entry ``169+13*i+j`` stores the element
        fixing ponts ``0`` and ``1`` that maps the points ``2`` 
        and ``5`` to the points ``i`` and ``j``. 

        Each entry ``e`` of ``cls.transversal`` is an array of 27 
        unsigned 8-bit integers, and it stores an element 
        ``a`` of ``AutP3`` as follows:

        ``e[i]`` is the image of the point ``P_i``, and ``e[13+i]``
        is the preimage of ``P_i`` under the automporphism ``a``.
        Entry ``e[26]`` if 1 if ``e`` is in use and 0 otherwise.

        Remarks:

        The entry with index 1 (corresponding to the identity in
        ``AutP3``) is precomputed and always marked as used. 

        The transversal is computed on demand; here we take the
        element of a coset that is passed first to this function.
        
        The preimages of the points ``P_i`` are computed for the
        entries containing the transversal only. 
        """
        p = g.perm
        # We must check the input ``g``, since a bad input will
        # permanently spoil the computed  transversal
        check_perm_P3(p)
        h2 = 13*p[0] + p[1]
        t2 = cls.transversal[h2]
        if t2[26] == 0:
            t2[:13] = p
            t2[13:26] = invert_perm_P3(t2[:13])
            t2[26] = 1
            return 1, h2
        else:
            f1 = mul_perm_P3(p, t2[13:26]) 
            h1 = 169 + 13 * f1[2] + f1[5]
            t1 =  cls.transversal[h1]
            if t1[26] == 0:
                t1[0:13] = f1
                t1[26] = 1
            return h1, h2


    @classmethod
    def _split_into_good_orders(cls, h):
        """Split instance ``h`` of class ``AutP3`` into two factors

        In some cases the sign of an image of an element ``h`` of
        the group ``AutP3`` in :math:`G_{x0}` cannot be determined.
        
        Here we always split ``h`` into a product  ``h1 * h2``, 
        choosing one of the factors at random, so that the signs of 
        both, ``h1`` and ``h2``, can be determined. Then we return
        ``h1, h2`` as a pair of instances of class ``AutP3``.

        The question whether the sign of an element of ``AutP3`` can
        be determined depends on the order of the element only. Here
        we assume that all keys of dictionary ``cls.good_orders``
        are orders for which this is possible.
        """
        cls.n_splits += 1
        while (1):
            h1 = AutP3('r')
            h2 = h1**(-1) * h
            cls.n_split_trials += 1
            if (h1.order() in cls.good_orders  and
                h2.order() in cls.good_orders):
                    return h1, h2

    
    @classmethod
    def compute_image_in_MM(cls, h):
        r"""Map the element  ```h`` of ``AutP3`` into the Monster group
    
        The function maps the instance ``h`` of class ``AutP3`` to an
        element :math:`g` of the  subgroup 
        :math:`G_{x0} = 2^{1.24}.\mbox{Co}_1` of the 
        Monster.  It returns a numpy array ``a`` such that the result 
        is equal to ``MM('a', a)``.

        If the sign of the image of  ``h`` in :math:`G_{x0}` cannot be
        determined directly then the function splits ``h`` into
        a product of two factors such that the signs of these factors
        can be determined.
        """
        # Let ``order`` be the order of ``h``
        order = h.order()
        # If we know how embed an element of that order into ``G_x0``,
        # then just do so. Store result in dictionary and return it.
        if order in cls.good_orders:
            g, _, _1 = MM_from_perm(h.perm)
            if cls.good_orders[order] == 0:
                g *= MM('x', 0x1000)
            return g.mmdata
        # Else decompose ``g`` as a (random) a product ``g = h1 h2``.
        # such that we can deal with elements of the orders of ``h1``
        # and `h2`. Compute the images of ``h1`` and ``h2`` as in the
        # previous case. Store product of these images in dictionary.
        h1, h2 = cls._split_into_good_orders(h)
        g =  MM_from_perm(h1.perm)[0] * MM_from_perm(h2.perm)[0]
        g_data = g.mmdata 
        # Try to deal with the order of input ``perm``. If we already
        # know that we cannot do so, we simply return the image.
        if order in cls.bad_orders:
            return g_data
        # If we can distinguish the image of ``perm`` from its
        # negative (via character theory) then remember the sign of
        # the character of the correct image for the order of 
        # ``perm``.  If we cannot do so, mark that order as bad.
        g1, _, special = MM_from_perm(h.perm)
        if not special:
            cls.bad_orders[order] = True
        else:
            if g1 == g:
                cls.good_orders[order] = 1       
            else:
                assert g1 == g *  MM('x', 0x1000)
                cls.good_orders[order] = 0
        return g_data


    @classmethod
    def map_to_MM(cls, t):
        r"""Store an element of AutP3 as an element of the Monster group
    
        Given the integer ``t``, let  ``h``  be the instance 
        ``h = AutP3('p', cls.tranversal[t, :13])   of class ``AutP3``
        The function maps the instance ``h`` to an  element :math:`g` 
        of the  subgroup :math:`G_{x0} = 2^{1.24}.\mbox{Co}_1` of the 
        Monster.  It returns a numpy array ``a`` such that the result 
        math:`g` is equal to ``MM('a', a)``.

        The function remembers all such arrays ``a`` already computed in
        the array ``cls.data_MM``.
        """
        # Return image of element in ``G_x0`` if already computed
        index = cls.ind_MM[t]
        if index:
            return cls.data_MM[index]

        # Do global precomputations if not yet done
        if precomputation_pending:
            precompute_all()
        # Let ``h`` be the instance of class ``AutP3`` given by ``t``
        h = AutP3('p', cls.transversal[t, :13])
        # Comput the image ``g_data`` of ``h``
        g_data = cls.compute_image_in_MM(h)
        # Store that image in the  array ``cls.data_MM`` and return it
        cls.ind_MM[t] = index = cls.num_MM
        cls.data_MM[index, :len(g_data)] = g_data
        cls.num_MM += 1
        return g_data       
 
    @classmethod
    def as_MM(cls, h):
        """Map instance ``h`` of class ``AutP3`` into the Monster"""
        # Split ``h`` into two factors
        h1, h2 = Precomputed_AutP3._split_transveral(h)
        # Map the product of automorphisms corresponding to the
        # entries ``cls.transversal[h1]`` and  ``cls.transversal[h2]``
        # into the Monster group.     
        return MM('a', np.hstack(
            (cls.map_to_MM(h1),  cls.map_to_MM(h2))
        ))
  

      

def AutP3_MM(h):    
    r"""Embed  the element ``h`` of AutP3 into the Monster

    The function computes the image of the automorphism ``h`` of the
    projective plane ``P3`` in the Monster group, and returns the
    result as an instance of class ``MM``.

    Parameter ``h`` must be an instance of class ``AutP3``. We use
    a fixed embedding of the automorphsm group of ``P3`` into 
    the Monster.
    """
    return Precomputed_AutP3.as_MM(h)
       
           







#####################################################################
# Norton's generators and relations for the Monster
#####################################################################



def comm(a,b):
    r"""Return commutator of group elements ``a`` and ``b``."""
    return a**(-1) * a**b


def Norton_generators_stuv(check = True):
    r"""Auxiliary function for function ``Norton_generators``

    The function returns a tuple of the first four values to be
    returned by that function. Parameter ``check`` is as in that
    function.
    """
    if precomputation_pending:
        precompute_all()
    MM1 = MM()
    # define genertors s, t, u
    s_AutP3 = AutP3(zip([1,2,5,9,8,7], [2,5,9,8,7,1]))
    s = AutP3_MM(s_AutP3)
    t_AutP3 = AutP3(zip([0,12,3,1,2,4], [12,3,0,2,4,1]))
    t = AutP3_MM(t_AutP3)
    u_AutP3 = s_AutP3 * t_AutP3 * s_AutP3**2 * t_AutP3**2
    u = s * t * s**2 * t**2
    v = MM(PointP3(range(1,13)))

    if check:
        # A consistency check of generator u, just for fun
        assert u == AutP3_MM(u_AutP3)
        # Check that u acts as  z -> z-1,  as stated in [Nor02]
        assert u_AutP3.point_map() == [12] + list(range(12))

        # Check the relations (1) in [Nor02]
        assert s**6 == t**3 == (s*t)**4 == (s**3 *t)**3 == MM1
        assert comm(s**2, (t * s**2 * t)**2) == MM1

        # Check the relations (2) in [Nor02]
        assert comm(v, u * t**(-1)) == comm(v, u**3 * s * u **(-2)) == MM1
        assert v**2  == comm(v, v**t) == (v*u)**13 == MM1
    return s, t, u, v


# Norton's generator 'x' has the form  MM('t', 1)**EXP_X in our 
# notation, where EXP_X is 1 or 2. We could quickly check the
# relations in these generators for both cases of EXP_X, and we 
# have set EXP_X to the correct value manually.
EXP_X = 1

# List of Norton_generators (if already computed)
NORTON_GENERATORS = None


def Norton_generators(check = False):
    r"""Return the images of Norton's generators in the Monster

    Norton :cite:`Nor02` defines a presention of the Monster
    with generators :math:`(s,t,u,v,x)` and relations.

    The function returns the images of the presentation of
    these generators in the Monster under a (fixed)
    homomorphism. The result is returned as a quintuple of
    instances of class ``MM``, corresponding to the images of
    these generators in the Monster.

    If parameter ``check`` is ``True`` then we also check the
    relations in this presentation.
    """
    global NORTON_GENERATORS
    if NORTON_GENERATORS and not check:
        return NORTON_GENERATORS
    s, t, u, v = Norton_generators_stuv(check)
    x = MM('t', EXP_X)
    if check:
        MM1 = MM()
        # Check the relations (4) in [Nor02]
        assert comm(v*x, u * t**(-1)) == comm(v*x, s**(u**3)) == MM1

        # Check the relations (5) in [Nor02]
        assert x**3 == (v**u * v * x)**2 == (x**(-1) * x**s)**2 == MM1

        # Check the relations (7) in [Nor02]
        assert (x * (t**-1))**12 == MM1

        # Check the relations (8) in [Nor02]
        assert (x**(u**6) * s)**6 * (s * u * x**(-1) * u**(-1))**6 == s

        # Check the relations (9) in [Nor02]
        assert ((x * v**(u**4) * v**(u**10))**3 * u)**13 == MM1

    NORTON_GENERATORS = s, t, u, v, x
    return NORTON_GENERATORS



