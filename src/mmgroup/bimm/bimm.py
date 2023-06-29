r"""This module implements the Bimonster.

Class ``BiMM`` in this module implements an element of the Bimonster
:math:`\mathbb{M} \wr 2` as described in the documentation of this
application. Let :math:`\mbox{IncP3}` be the Coxeter group as in
that documentation. Function ``P3_BiMM`` maps a word of generators
of :math:`\mbox{IncP3}` into the Bimonster. The generators of
:math:`\mbox{IncP3}` correspond to the points and lines  of the
projective plane :math:`\mbox{P3}` over the field :math:`\mathbb{F}_3`.
A point or a line  of :math:`\mbox{P3}` is implemented as an
instance of class ``P3_node`` in module ``inc_p3``.

There is also a natural mapping from the automorphism group of
:math:`\mbox{P3}` into the Bimonster compatible with the mapping
from the Coxeter group into the Bimonster. Function
``AutP3_BiMM`` computes that mapping. An automorphism  of
:math:`\mbox{P3}` is implemented as an instance of class ``AutP3``
in module ``inc_p3``. For background see  :cite:`Nor02`.

"""
import os
import sys
from numbers import Integral
from random import randint, sample, choice
from sys import getrefcount

import numpy as np

from mmgroup import  MM
from mmgroup.structures.abstract_group import AbstractGroup
from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.abstract_group import singleton

from mmgroup.bimm import inc_p3
from mmgroup.bimm.inc_p3 import p3_list, P3_node
from mmgroup.bimm import p3_to_mm
from mmgroup.bimm.p3_to_mm import PointP3, AutP3, AutP3_MM
from mmgroup.bimm.p3_to_mm import Norton_generators


#####################################################################
# class BiMM
#####################################################################


def gcd(a, b):
    r"""Return greatest common divisor of ``a`` and ``b``"""
    assert a > 0 and b > 0
    while b > 0: a, b = b, a % b
    return a

ERR_BIMM1 = "Single argument for class BiMM must be 'r' or instance of class BiMM"

class BiMM(AbstractGroupWord):
    r"""This class models an element of the Bimonster.

    The Bimonster is the group  :math:`\mathbb{M} \wr 2` of
    structure  :math:`(\mathbb{M} \times \mathbb{M}).2`, where
    :math:`\mathbb{M}`  is the Monster group.

    :param m1:
     
       An element :math:`m_1` of the Monster that embeds into the
       Bimonster as :math:`(m_1, 1)`. Default is the neutral
       element :math:`1` of the Monster.

    :type m1:

       Instance of class ``MM``

    :param m2:
     
       An element :math:`m_2` of the Monster that embeds into the
       Bimonster as :math:`(1, m_2)`.  Default is the neutral
       element :math:`1` of the Monster.

    :type m2:

       Instance of class ``MM``

    :param e:

       An optional exponent of the involution :math:`\alpha`,
       default is 0.
       Conjugation of an element of 
       :math:`\mathbb{M} \times \mathbb{M}` by  :math:`\alpha`
       means swapping the two factors in the direct product.

    :type e:

       integer

    :return: 

       The element :math:`(m_1, m_2) \cdot \alpha^e` of the
       Bimonster
           
    :rtype:  

       Instance of class ``BiMM``

    Arguments ``m1`` or ``m2`` may be anything that is accepted by
    the constructor of class ``MM`` as a single argument.

    Alternatively, ``m1`` may be an instance of class ``BiMM``;
    then arguments ``m2`` and ``e`` must be dropped.

    Let ``g1`` and ``g2`` be instances of class ``BiMM`` representing 
    elements of the Bimonster group.
    ``g1 * g2``  means group multiplication, and ``g1 ** n`` means
    exponentiation of ``g1`` with the integer ``n``. ``g1 ** (-1)`` 
    is the inverse of ``g``. ``g1 / g2`` means ``g1 * g2 ** (-1)``.
    We have ``1 * g1 == g1 * 1 == g1`` and ``1 / g1 == g1 ** (-1)``.

    ``g1 ** g2`` means ``g2**(-1) * g1 * g2``.   


    """
    __slots__ = "m1", "m2", "alpha"
    group_name = "BiMM"
    group = None    # will be set to StdBiMMGroup later

    def __init__(self, *args):
        if len(args) == 0:
            self.m1 = MM()
            self.m2 = MM()
            self.alpha = 0
        elif len(args) == 1:
            m = args[0]
            if isinstance(m, BiMM):
                self.m1 = m.m1.copy()
                self.m2 = m.m2.copy()
                self.alpha = m.alpha
            elif m == 1:
                self.m1 = MM()
                self.m2 = MM()
                self.alpha = 0
            elif m == 'r':
                self.m1 = MM('r')
                self.m2 = MM('r')
                self.alpha = randint(0,1)
            else:
                raise TypeError(ERR_BIMM1) 
        else:
            self.m1, self.m2 = MM(args[0]),  MM(args[1])
            self.alpha = args[2] & 1 if len(args) > 2 else 0
            
    def orders(self):
        r"""Return orders(!) of element of the group BiMM"""
        self.reduce()
        if self.alpha & 1:
             a, par = self * self, 2
        else:
             a, par = self, 1
        return a.m1.order(), a.m2.order(), par

    def order(self):
        r"""Return the order of the element of the Bimonster"""
        o1, o2, s = self.orders()
        return s * o1 * o2 // gcd(o1, o2)


    def reduce(self):
        r"""Reduce the element of the Bimonster""" 
        self.m1.reduce()
        self.m2.reduce()
        self.alpha &= 1

    def decompose(self):
        r"""Decompose the element of the Bimonster 

        The function returns a triple ``(m1, m2, e)`` such
        that the element of the Bimonster is equal to
        :math:`(m_1, m_2) \cdot \alpha^e`. Here ``m1`` and ``m2``
        are instances of class ``MM`` representing the elements
        :math:`m_1` and :math:`m_2` of the Monster, and ``e`` is
        equal to 0 or 1.
        """
        self.reduce()
        return self.m1, self.m2, self.alpha

    def __hash__(self):
        m1, m2, alpha = self.decompose()
        return hash((hash(m1), hash(m2), alpha))

 
@singleton
class BiMMGroup(AbstractGroup):
    """Auxilary class for class ``BiMM`` 

    This makes the methods in class ``AbstractGroup`` available to
    instancs of class ``BiMM``.
    """
    word_type = BiMM   
    conversions = {}

    def __init__(self):
        super(BiMMGroup, self).__init__()


    def atom(self, tag = None, data = None):
        err = "Class BiMMGroup has no attribute 'atom'"
        raise AttributeError(err)

   

    @staticmethod
    def _imul(g1, g2):
        # Beware of imumutability:
        # If actually only one copy of g1 exits then we have g1,
        # a call to g1.__mul__, to cls. _imul, and to getrefcount
        if getrefcount(g1) > 4:
            g1 = BiMM(g1.m1, g1.m2, g1.alpha)
        if g1.alpha & 1:
            m1, m2 = g1.m1 * g2.m2, g1.m2 * g2.m1
        else:
            m1, m2 = g1.m1 * g2.m1, g1.m2 * g2.m2
        return BiMM(m1, m2, (g1.alpha ^ g2.alpha) & 1) 

    @staticmethod
    def _invert(g1):
        if g1.alpha & 1:
            m1, m2 = g1.m2, g1.m1
        else:
            m1, m2 = g1.m1, g1.m2
        return BiMM(m1**(-1), m2**(-1), g1.alpha & 1)

    @staticmethod
    def copy_word(g1):
        return BiMM(g1.m1, g1.m2, g1.alpha)

    @staticmethod
    def _equal_words(g1, g2):
        g1.reduce()
        g2.reduce()
        return g1.m1 == g2.m1 and g1.m2 == g2.m2 and g1.alpha == g2.alpha

    def str_word(self, g):
        r"""Convert group atom g to a string

        """
        s1, s2 = g.m1.raw_str(), g.m2.raw_str()
        al = ', *' if g.alpha & 1 else ''
        return "BiMM<%s, %s%s>" % (s1, s2, al)


StdBiMMGroup = BiMMGroup()   # This is the only instance of AutP3Group

BiMM.group = StdBiMMGroup

 

#####################################################################
# Enter points and lines into the Bimonster
#####################################################################

# This will be set to ``False``  after finishing the precomputation    
precomputation_pending = True



# ``ALPHA`` a will be set to the involution (in the Bimonster)
# swapping the two copies of the Monster.
ALPHA = None #  This will be set to BiMM([],[], 1)

# PL_DATA will be an array of shape (26, 2, MAXLEN_PL_DATA).
# Then the node ``P3_node(i)`` of the projective plane ``P3`` will be
# mapped to ``BiMM(MM('a', PL_DATA[i,0]), MM('a', PL_DATA[i,1])) * ALPHA 
PL_DATA = None
MAXLEN_PL_DATA = 0

def precompute_points_lines_list():
    """Perform the precomputations required for this module.

    The function sets the values ``PL_DATA`` and ``ALPHA`` as
    defined above. We map the point ``P_i`` to 
    ``(P_0 * P_1 * ... * P_12 * P_i)`` * ALPHA. Here the points in
    that product are given by function  ``PointP3`` in module
    ``inc_p3``. We map the line ``L_i`` to 
    ``ALPHA * BiMM(L_i, L_i ** (-1))``, with ``L_i = (v*x) ** (u**i)``.
    Here ``u, v, x`` are as returned by function ``Norton_generators``
    in module ``p3_to_mm``.    
    """
    global PL_DATA, MAXLEN_PL_DATA, ALPHA
    ALPHA = BiMM([],[], 1)
    PL_VALUES = []
    s, t, u, v, x =  Norton_generators()
    for i in range(13):
        P_i = MM(PointP3(list(range(13)) + [i]))
        bm = BiMM(P_i, P_i, 1)
        PL_VALUES.append(bm)
    for i in range(13):
        L_i = (v*x) ** (u**i)
        bm = ALPHA * BiMM(L_i, L_i ** (-1))
        PL_VALUES.append(bm)
    data = [[list(bm.m1.mmdata),list(bm.m2.mmdata)] for bm in PL_VALUES]
    MAXLEN_PL_DATA = max(max(len(m1), len(m2)) for m1, m2 in data)
    #print(maxlen)
    PL_DATA = np.zeros((26, 2, MAXLEN_PL_DATA), dtype = np.uint32)
    for i in range(26):
        for j in range(2):
            d = data[i][j]
            PL_DATA[i, j, :len(d)] = d

 


def P3_BiMM(pl = []):
    r"""Map a word of generators in :math:`\mbox{IncP3}` into the Bimonster

    :param pl:

       List of generators in :math:`\mbox{IncP3}`. Each entry in the 
       list should be an instance of class ``P3_node``. Such an entry
       may also be an integer or a string accepted by the constructor
       of class ``P3_node``.

    :type pl:

        List containing integers, strings, or instances of 
        class ``P3_node``

    :return: 

       The image of the word of the given generators in the  Bimonster
           
    :rtype:  

       Instance of class ``BiMM``

     
    An integer ``pl`` or an instance ``pl`` of class ``P3_node`` is
    considered as a word of length 1.

    ``pl`` may also be a string of alphanumric identifiers separated
    by commas. This is interpreted as a sequence of generators,
    where the names of the generators are interpreted as in the 
    constructor of class ``P3_node``.
    """
    if precomputation_pending:
        precompute_all()
    if isinstance(pl, (Integral, P3_node)):
        pl = [P3_node(pl).ord]
    else:
        pl = p3_list(pl)
    a = np.zeros((2, len(pl), MAXLEN_PL_DATA), dtype = np.uint32)
    len_pl = len(pl)
    for i in range(0, len_pl & -2, 2):
        k = pl[i]
        a[0,i] = PL_DATA[k,0]
        a[1,i] = PL_DATA[k,1]
        k = pl[i+1]
        a[0,i+1] = PL_DATA[k,1]
        a[1,i+1] = PL_DATA[k,0]
    al = len_pl & 1
    if al:
        k = pl[-1]
        a[0,-1] = PL_DATA[k,0]
        a[1,-1] = PL_DATA[k,1]
    return BiMM(MM('a', a[0].ravel()), MM('a', a[1].ravel()), al)
   


#####################################################################
# Enter automorphisms of ``P3`` into the Bimonster
#####################################################################

def AutP3_BiMM(g):
    r"""Map an automorphism of :math:`\mbox{P3}` into the Bimonster

    :param g:

       Automorphism of :math:`\mbox{P3}`

    :type g:
  
       Instance of class  ``AutP3``

    :return: 

       The image of the automorphism of :math:`\mbox{P3}` in the  
       Bimonster
           
    :rtype:  

       Instance of class ``BiMM``

    Parameters ``g`` may be anything that is accepted by
    the constructor of class ``AutP3`` as a single argument.
    
    """
    g = AutP3_MM(AutP3(g))
    return BiMM(g, g.copy())


#####################################################################
# Precomputation
#####################################################################

def precompute_all():
    # Do precompution on demand only. Otherwise Sphinx will fail.
    global precomputation_pending
    precompute_points_lines_list()
    precomputation_pending = False

