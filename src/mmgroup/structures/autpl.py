r"""We deal with a certain group of automorphisms of the Parker loop.

A *standard automorphism* is an automorphism of the Parker loop 
:math:`\mathcal{P}`  which maps to an automorphism of the Golay 
code :math:`\mathcal{C}` when  factoring out the elements 
:math:`\{1,-1\}` of :math:`\mathcal{P}`. Let
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` be the 
group of standard automorphisms of  :math:`\mathcal{P}`. We 
embed :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`
into the Monster :math:`{\mathbb{M}}` as in :cite:`Con85`.
In ATLAS notation, see :cite:`Atlas`, that group has structure

.. math::
      {{\rm Aut}_{{\rm St}} \mathcal{P}} =
         2^{12} . M_{24}\,  .
         

This means that :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` has an
elementary Abelian normal subgroup of order  :math:`2^{12}` with
factor group  :math:`M_{24}`. That group extension does not
split. Here the group :math:`2^{12}` is isomorphic
to the Golay cocode :math:`\mathcal{C}^*`  and  :math:`M_{24}`
is the Mathieu group acting on :math:`\mathbb{F}_2^{24}`
by permuting the basis vectors. :math:`M_{24}` is the
automorphism group of the Golay code  :math:`\mathcal{C}` .

The automorphisms in the subgroup :math:`\mathcal{C}^*` of
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` are called 
*diagonal* automorphisms. A diagonal automorphism  
:math:`d \in \mathcal{C}^*` maps an element :math:`a` of  
:math:`\mathcal{P}` to :math:`(-1)^s \cdot a` with 
:math:`s = |a \, \& \, d|`. 
Here :math:`\&` means the bitwise ``and`` operation of bit 
vectors and :math:`|.|` means the bit weight of a bit vector.


Since the extension :math:`2^{12} . M_{24}` does not split,
there is no canonical embedding of :math:`M_{24}` into 
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`. For practical 
purposes we need such an embedding. We choose a basis
:math:`b_0,\ldots,b_{11}` of the Golay code :math:`\mathcal{C}`,
see section :ref:`basis-golay-label`. Let :math:`a_i` be
the positive element of the Parker loop :math:`\mathcal{P}`
corresponding  to  :math:`b_i`. For any 
:math:`\tilde{g} \in M_{24}` we let :math:`g` be the unique 
element of :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` with
:math:`g \mathcal{C}^* = \tilde{g}` 
that maps all elements :math:`a_0,\ldots,a_{11}` of
:math:`\mathcal{P}` to positive elements of :math:`\mathcal{P}`.
We call :math:`g` the *standard representative* of 
:math:`\tilde{g}` in :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`. 


Elements of :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` are 
modelled as instances of class |AutPl|. If ``d`` is an instance
of class |Cocode| representing an element of  :math:`\mathcal{C}^*`
then ``AutPL(d)`` is the corresponding diagonal automorphism.

A permutation :math:`p \in M_{24}` can be represented as a
list ``l_p = [p(0),...,p(23)]``. If ``l_p`` is such a list
then ``AutPL(l_p)`` is the standard representative of ``p``
in  :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`. Here an error
occurs if the mapping from ``i`` to ``p(i)`` is not in
:math:`M_{24}`. A permutation in :math:`M_{24}` can also be 
given as a mapping from a subset of :math:`\{0,\ldots,23\}`
to  :math:`\{0,\ldots,23\}`, provided that this mapping
extends to a unique element of :math:`M_{24}`. For details,
see class  |AutPl|. 


The group :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` operates
naturally on :math:`\mathcal{P}` by right multiplication.
It also operates on  :math:`\mathbb{F}_2^{24}`, 
:math:`\mathcal{C}`,  and :math:`\mathcal{C}^*` by right 
multiplication. Here the basis vectors of 
:math:`\mathbb{F}_2^{24}` are permuted by :math:`M_{24}`;
so the kernel of that operation is the subgroup
:math:`\mathcal{C}^*` of 
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`.

We order the elements  :math:`p` of :math:`M_{24}` in lexicographic
order based on the lists ``l_p = [p(0),...,p(23)]``. We number the
``244823040`` elements of :math:`M_{24}` in that order from ``0`` 
to ``244823039``.  Thus the number ``0`` is assigned to the neutral
element of  :math:`M_{24}`, and the number ``244823039`` is 
assigned to the permutation ``p`` with ``l_p`` =


    :math:`\textstyle [ 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 
    11, 10, 9, 8, 0, 1, 2, 3, 4, 5, 6, 7 ]` .

These numbers are used for addressing the corresponding standard 
representatives in the subgroup 
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` of the
monster :math:`\mathbb{M}`. 

One way to create an element of :math:`M_{24}` is to
map an *umbral heptad* to another umbral heptad. Here an
umbral heptad is a set of ``7`` elements of 
:math:`\{0,\ldots,23\}` such that precisely ``6`` of these
elements are contained in an octad, see :cite:`CS99`, 
Chapter 10. Then the ``6`` elements in one octad must map 
to the ``6`` elements in the other octad.

The set ``[0,1,2,3,4,5,8]`` is an umbral heptad with ``8``
not contained in the octad. We may create the standard 
representative (in :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`)
of a random element of :math:`M_{24}` as follows:


.. code-block:: python

  from random import sample
  from mmgroup import Octad, AutPL
  # Create a random octad
  o = Octad("r")
  # Select 6 random elements from that octad
  hexad = sample(o.bit_list, 6)
  # Select one random element not in that octad
  extra = sample((~o).bit_list, 1)
  # Create umbral heptad from these random elements
  heptad = hexad + extra
  # Create a mapping from [0,1,2,3,4,5,8] to that heptad 
  hmap = zip([0,1,2,3,4,5,8], heptad)
  # Create a permutation in the Mathieu group from that mapping
  aut = AutPL(0, hmap) 
  # 'aut' is the standard representative of the permutation
  # Show 'aut' as a permutation in form of a list 
  print(aut.perm)     
  # Show the automorphism 'aut' with its permutation number 
  print(aut)     

**Remark**

Class |AutPl| is a subclass of class 
``mmgroup.structures.AbstractGroupWord`` which implements
elements of an abstract group and the operations on such elements.
For an instance ``g`` of |AutPl|, the object ``g.group`` is
an instance of a subclass of  ``mmgroup.structures.AbstractGroup``.
That subclass models the group
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`.


.. _aut_ploop_rand_label:


Generating random elements of certain subgroups of :math:`M_{24}`
..................................................................

We support the generation of random elements of certain subgroups
of :math:`M_{24}`. Here any such subgroup stabilizes a subset
(or a sets of subsets) of the set
:math:`\bar{\Omega} = \{0,\ldots,23\}` on which :math:`M_{24}` acts.
Furthermore, any such group is named by a flag (which is a character),
as indicated in the table below. Intersections of these subgroups
are described by combining the characters. The character ``r``
describes the whole group  :math:`M_{24}`. A string of characters
describing such a subgroup should start with ``r``. Whitespace is
ignored.

Table of subgroups:

.. math::
      \begin{array}{|c|c|c|c|}
      \hline 
      \mbox{Flag} & \mbox{Mnemonic} &
         \mbox{Subgroup stabilizing the set} &
         \mbox{Structure}  \\
      \hline
      \mbox{'r'} & \mbox{(random)} & \mbox{the whole set } \bar{\Omega} &
         M_{24} \\ 
      \mbox{'2'} & \mbox{2-set} & \{2,3\} &
         M_{22}:2 \\ 
      \mbox{'o'} & \mbox{octad}  &  \{0,\ldots,7\} &
         2^4:A_8  \\
      \mbox{'t'} & \mbox{trio}   & 
         {\{ \{8i,\ldots,8i+7\} \mid i < 3 \}} &
         2^6:(S_3 \times L_3(2))   \\
      \mbox{'s'} & \mbox{sextet} & 
         {\{ \{4i,\ldots,4i+3\} \mid i < 6 \}} &
         2^6:3.S_6 \\
      \mbox{'l'} & \mbox{(line)} & \{ \{2i, 2i+1\} \mid 4 \leq i < 12 \} &
       2^{1+6}:L_3(2)  \\
      \mbox{'3'} & \mbox{3-set} & \{1, 2,3\}   &
         L_3(4):S_3 \\
      \hline 
      \end{array}

For mathematical background we refer to section :ref:`subgroup-mat24-label`.


.. table:: Examples of strings describing subgoups of :math:`M_{24}`
      :widths: 25 75

      ========== ========================================================
      String     Action
      ========== ========================================================
      ``'r'``    Generate a random element of :math:`M_{24}`

      ``'r 2o'`` Generate a random element of :math:`M_{24}` stabilizing
                 :math:`\{0,\ldots,7\}` and :math:`\{2,3\}`
      ========== ========================================================

"""

import re

from functools import reduce
from operator import __or__
from numbers import Integral, Number
from random import randint

import numpy as np

from mmgroup.structures.abstract_group import singleton
from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.abstract_group import AbstractGroup
from mmgroup.structures.parse_atoms import AtomDict, ihex     
from mmgroup.structures.parse_atoms import eval_atom_expression        





ERR_TYPE = "Unsupported operand types for %s: '%s' and '%s'"
ERR_RAND = "Illegal string for constricting type %s element" 

ERR_RND_S = "String describing a random permutation must begin with 'r'"
ERR_RND_CH = "Illegal character '%s' for describing a random permutation"



#######################################################################
# Import derived classed
#######################################################################

import_pending = True
Cocode = None


def complete_import():
    """Internal function of this module

    If you need any of the objects declared above, do the following:

    if import_pending:
        complete_import()
    """
    global import_pending, mat24, MAT24_ORDER 
    global Cocode
    from mmgroup import mat24, MAT24_ORDER
    from mmgroup.structures.cocode import Cocode
    import_pending = False


#######################################################################
# Generating a random permutation number form a string
#######################################################################




def rand_perm_num(s):
    """Return number of a random permutation depending on string 's'"""
    if import_pending:
       complete_import()
    mode = 0
    if s[:1] != "r":
        raise ValueError(ERR_RND_S)
    for c in s[1:]:
        if not c.isspace():
            try:
                mode |= mat24.MAT24_RAND[c]
            except KeyError:
                raise ValueError(ERR_RND_CH % c)
    rand = randint(0, MAT24_ORDER - 1)
    return mat24.m24num_rand_local(mode, rand)



#######################################################################
# Function implementing the constructor for class AutPL
#######################################################################



ERR_AUTPL_P1 = "AutPL expects at most 1 parameter of type %s" 
ERR_AUTPL_TYPE = "Cannot construct AutPL instance from type %s" 


autpl_conversions = {
}

def autpl_from_obj(d = 0, p = 0, unique = 1):
    """Try to convert tuple (d, p) to a Parker loop automorphism.

    Parameter ``d`` decribes a element of the Golay cocode as in the
    constructor of class ``Cocode``. It may be of type ``int``, ``str``
    or an instance of class ``Cocode``. Pt defaults to ``0``.

    Alternatively, ``d`` may be an instance of class ``AutPL``. In this 
    case, ``p`` must be set to its default value.


    Parameter ``p`` describes a element of the Mathieu group ``Mat24``.
    It defaults to the neutral element of ``Mat24``. It may be
    
     * An integer encoding the number of a permutation in ``Mat24``.

     * A list encoding a permutation in the Mathieu group ``Mat24``.

     * A zip object or a dictionary encodes such a permutation as
       a mapping. That mapping must contain sufficiently many 
       entries to be unique.  

     * The string 'r' encodes a random permutation in ``Mat24``.

    If parameter ``unique`` is ``True`` (default), then the parameter
    ``p`` must either be ``r`` or describe a unique permutation in
    ``Mat24``. Otherwise lowest possible permutation (in 
    lexicographic order) is taken.

    The function returns a pair ``(cocode, permutation_number)``.
    """
    if import_pending:
        complete_import()
    if isinstance(d, Integral):
        d_out = int(d)
    elif isinstance(d, Cocode):
        d_out = d.value
    elif isinstance(d, AutPL):
        d_out, p_out =  d._cocode, d._perm_num
        if p:
            raise TypeError(ERR_AUTPL_P1 % type(d))
        return d_out, p_out
    elif isinstance(d, str):
        d_out = Cocode(d).value
    else: 
        try:
            f = autpl_conversions[type(d)]
        except KeyError:
            raise TypeError(ERR_AUTPL_TYPE % type(d))
        if p:
            raise TypeError(ERR_AUTPL_P1 % type(d))
        d_out, p_out =  f(d)
        return d_out, p_out
            
    if isinstance(p, str):
        return d_out, rand_perm_num(p)
    if isinstance(p, Integral):
        if 0 <= p <  MAT24_ORDER:
            return d_out, int(p)
        err = "Bad number for permutation in Mathieu group  Mat24"
        raise ValueError(err) 
    if isinstance(p, (list, range)):
        if mat24.perm_check(p):
            err = "Permutation is not in Mathieu group M_24"
            raise ValueError(err)
        return d_out, mat24.perm_to_m24num(p)
    if isinstance(p, zip):
        p = dict(p)
    if isinstance(p, dict):
        h1, h2 = [list(t) for t in zip(*p.items())]
        res, perm = mat24.perm_from_map(h1, h2)
        if res == 1 or (res > 1 and not unique):
            return d_out, mat24.perm_to_m24num(perm) 
        if res < 1:
            err = "Permutation is not in Mathieu group M_24"
        else:                 
            err = "Permutation in Mathieu group M_24 is not unique"
        raise ValueError(err)






#######################################################################
# Class AutPL
#######################################################################


            



class AutPL(AbstractGroupWord):
    r"""This class models a standard automorphism of the Parker loop.

    :param d:

      This parameter describes an element of the Golay cocode. If ``d``
      is an instance of class |AutPL| then this describes an
      automporphism of the Parker loop; in that case parameter ``p``
      must be set to its default value. Legal types of parameter ``d`` 
      see  below. The parameter defaults to zero word of the cocode.
       
    :param p:

      This parameter describes a permutation in the Mathieu group 
      ``M_24`` or, more precisely, the standard representative of 
      that permutation in the automporphism group of the Parker 
      loop. Legal types of that parameter see below.
      The parameter defaults to neutral element of ``Mat24``.

    :param unique:

       If this is ``True`` (default) and parameter ``p`` is not ``r``
       then parameter ``p`` must describe a unique permutation in
       the Mathieu group ``M_24``. Otherwise we take the least 
       possible permutation (in lexicographic order) that agrees
       with  parameter ``p``.
 
    :return: A standard Parker loop automorphism
    :rtype:  an instance of class |AutPL|

    :raise:
        * TypeError if ``type(value)`` is not in the table given below.
        * ValueError if the set of arguments cannot be converted to an
          instance of class  |AutPL|.


    .. table:: Legal types for parameter ``d`` in constructor of class ``AutPL``
      :widths: 25 75

      ===================== ==================================================
      type                  Evaluates to
      ===================== ==================================================
      ``int``               The number of an element of the Golay cocode.

      class |Cocode|        This represents an element of the Golay cocode.

      ``str``               Create random element depending on ``str``
                              | ``'r'``: Create an arbitrary cocode element
                              | ``'e'``: Create an even cocode element
                              | ``'o'``: Create an odd cocode element

                            Any other string is illegal.

      class |AutPL|         A deep copy of the given automorphism 
                            in |AutPL| is returned. Then parmeter ``p``
                            must be set to its default value.

      ===================== ==================================================


    .. table:: Legal types for parameter ``p`` in constructor of class ``AutPL``
      :widths: 25 25 50

      +-------------------+--------------------------------------------------+
      |type               |Evaluates to                                      |
      +===================+==================================================+
      |``int``            |Here the integer is the number of a permutation   |
      |                   |in the Mathieu group ``M_24``.                    |
      +-------------------+--------------------------------------------------+
      |``list`` of ``int``|A list ``l_p`` of ``24`` of disjoint integers     |
      |                   |``0 <= i < 24`` is interpreted as a permutation   |
      |                   |in ``M_24`` that maps ``i`` to ``l_p[i]``.        |
      +-------------------+--------------------------------------------------+
      |``dict``           |A dictionary specifies a mapping from a subset    |
      |                   |of the integers ``0 <= i < 24`` to integers       |
      |                   |``0 <= dict[i] < 24``. This mapping must          |
      |                   |extend to a permutation in ``M_24``.              |
      |                   |If parameter  ``unique`` is ``True`` (default)    |
      |                   |then that permutation must be unique in ``M_24``. |
      +-------------------+--------------------------------------------------+
      |``zip`` object     |``zip(x,y)`` is equivalent to ``dict(zip(x,y))``. |
      +-------------------+--------------------------------------------------+
      |``str``            |Create random element depending on the string     |
      |                   +---------------+----------------------------------+
      |                   |``'r'``        | Create random element of ``M_24``|
      |                   +---------------+----------------------------------+
      |                   |``'r <flags>'``| Create random element of subgroup|
      |                   |               | of ``M_24``, depending on        |
      |                   |               | ``<flags>``, see explanation     |
      |                   |               | above.                           |
      +-------------------+---------------+----------------------------------+
                         


    Let ``a`` be an instance of class |GcVector|, |GCode|, |Cocode|,
    or |PLoop|, and let ``g1`` , ``g2`` be instances of class |AutPL|. 
    Then  ``a * g1`` is the result of the natural operation of ``g1`` 
    on  ``a``, which belongs to the same class as ``a``.

    ``g1 * g2``  means group multiplication, and ``g1 ** n`` means
    exponentiation of ``g1`` with the integer ``n``. ``g1 ** (-1)`` 
    is the inverse of ``g``. ``g1 / g2`` means ``g1 * g2 ** (-1)``.

    ``g1 ** g2`` means ``g2**(-1) * g1 * g2``. 
   
    """
    __slots__ = "_cocode", "_perm_num", "_perm", "rep"
    ERR_UNDEF = "Parker loop automorphism is not defined by input data"
    _perm_ = list(range(24))  # neutral permutation
    group = None       # will be set to StdAutPlGroup later
    
    def __init__(self, d = 0, p = 0, unique = 1):
        if import_pending:
            complete_import()
        self._cocode, self._perm_num = autpl_from_obj(d, p, unique)
        self._compute_from_numbers()

    def _compute_from_numbers(self):
        self._perm = mat24.m24num_to_perm(self._perm_num)
        self.rep = mat24.perm_to_autpl(self._cocode, self._perm)
            
    def _compute_from_rep(self):
        self._perm =  mat24.autpl_to_perm(self.rep)
        self._perm_num = mat24.perm_to_m24num(self._perm)
        self._cocode = mat24.autpl_to_cocode(self.rep)

    def check(self):
        """Check automorphism for consistency via 'assert' statements

        ``a.check()`` returns ``a``.
        """
        assert self._perm ==  mat24.autpl_to_perm(self.rep)
        assert self._perm_num == mat24.perm_to_m24num(self._perm)
        assert 0 <= self._perm_num < mat24.MAT24_ORDER
        assert self._cocode == mat24.autpl_to_cocode(self.rep)
        return self

    @property    
    def parity(self):
        """Return parity of an automorphism
 
        This is ``s`` if ``PLoopOmega`` maps to ``(-1)**s * PLoopOmega``
        for ``s == 0`` or ``s == 1``.
        """
        return (self._cocode >> 11) & 1

    @property    
    def cocode(self):
        """Return number of cocode element in the automorphism.

        For an instance ``g`` of |AutPL| we have

          ``g == AutPL(Cocode(g.cocode)) * AutPL(g.perm)``.
        """
        return self._cocode

    @property        
    def perm_num(self):
        """Return the number of the permutation of the automorphism

        The Mathieu group ``M_24`` has ``244823040`` elements 
        numbered from ``0`` to ``244823039`` in lexicographic
        order. Thus the neutral element has number ``0``. 
        """
        return self._perm_num

    @property    
    def perm(self):
        """Return permutation of automorphism as a list ``p_l``.

        Then the permutation maps ``i`` to ``p_l[i]`` for
        ``i = 0,...,23``.
        """
        return self._perm

    @property    
    def mmdata(self):
        return np.array([
            0x10000000 + (self._cocode & 0xfff),
            0x20000000 + self._perm_num
        ], dtype = np.uint32)

    def as_tuples(self):
        return [('d', self._cocode), ('p', self._perm_num)]
       



#######################################################################
# Class AutPlGroup
#######################################################################


def autpl_element_from_obj(g, t):
    res =  AutPL()
    res._cocode, res._perm_num = autpl_from_obj(t)
    res._compute_from_numbers()
    return res

@singleton
class AutPlGroup(AbstractGroup):
    word_type = AutPL              # type of an element (=word) in the group
    is_mmgroup = True
    #atom_parser = {}               # see method parse()
    #rand_masks  = {"r":(0xfff,0), "e":(0x7ff,0), "o":(0x7ff,0x800)} 
    conversions = {
      #  list: autpl_element_from_obj,
      #  zip: autpl_element_from_obj,
      #  dict: autpl_element_from_obj,
      #  AutPL: autpl_element_from_obj,
    }
    #FRAME =  AUTPL_FRAME
    ERR_PERM_NUM = "Illegal permutation number for Mathieu group M_24"

    def __init__(self):
        super(AutPlGroup, self).__init__()
        self.atom_parser = AtomDict(self.atom)


    def __call__(*args, **kwds):
        err = "Class AutPlGroup object is not callable"
        raise TypeError(err)

    def atom(self, tag = None, data = None):
        err = "Class AutPlGroup has no attribute 'atom'"
        raise AttributeError(err)

    def _imul(self, g1, g2):
        g1.rep = mat24.mul_autpl(g1.rep, g2.rep)
        g1._compute_from_rep()
        return g1    

    def _invert(self, g1):
        res = AutPL()
        res.rep = mat24.inv_autpl(g1.rep)
        res._compute_from_rep()
        return res

    def copy_word(self, g1):
        res = AutPL()
        res._cocode, res._perm_num = g1._cocode, g1._perm_num
        res._perm, res.rep =  g1._perm[:], g1.rep[:]
        return res

    def _equal_words(self, g1, g2):
        return g1._cocode == g2._cocode and g1._perm_num == g2._perm_num

        

    def str_word(self, g):
        """Convert group atom g to a string

        """
        s = []
        if g._cocode:
            s.append("d_" + ihex(g._cocode))
        if g._perm_num:
            s.append("p_%d" % g._perm_num)
        s = "*".join(s) if len(s) else "d_0"
        return "AutPL<%s>" % s



StdAutPlGroup = AutPlGroup()   # This is the only instance of AutPlGroup

AutPL.group = StdAutPlGroup

