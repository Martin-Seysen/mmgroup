r"""We deal with a certain group of automorphisms of the Parker loop.

A *standard automorphism* is an automorphism of the Parker loop 
:math:`\mathcal{P}`  which maps to an automorphism of the Golay 
code :math:`\mathcal{C}` when  factoring out the elements 
:math:`\{1,-1\}` of :math:`\mathcal{P}`. Let
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` be the 
group of standard automorphisms of  :math:`\mathcal{P}`. We 
embed :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`
into :math:`{\mathbb{M}}` as in :cite:`Con85`. 
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
:math:`b_1,\ldots,b_{11}` of the Golay code :math:`\mathcal{C}`,
see section :ref:`basis-golay-label`. Let :math:`a_i` be
the positive element of the Parker loop :math:`\mathcal{P}`
corresponding  to  :math:`b_i`. For any 
:math:`\tilde{g} \in M_{24}` we let :math:`g` be the unique 
element of :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` with
:math:`g \mathcal{C}^* = \tilde{g}` 
that maps all elements :math:`a_1,\ldots,a_{11}` of
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

The ``244823040`` elements of :math:`M_{24}` are numbered 
from ``0`` to ``244823039``, with the number ``0`` assigned to 
the identity. These numbers are used for addressing
the corresponding standard representatives in the subgroup 
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` of the
monster :math:`\mathbb{M}`. This numbering system is
not part of the public interface, but the methods of
class |AutPL| provide access to this numbering.

The easiest way to create an element of :math:`M_{24}` is to
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
  aut = AutPL(hmap) 
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

"""

import re

from functools import reduce
from operator import __or__
from numbers import Integral, Number
from random import randint

from mmgroup.structures.auto_group import AbstractGroupWord
from mmgroup.structures.auto_group import AbstractGroup
from mmgroup.structures.parse_atoms import AtomDict, ihex     

try:
    # Try importing the fast C function
    from mmgroup import mat24 
except: # (ImportError, ModuleNotFoundError):
    # Use the slow python function if the C function is not available
    from mmgroup.dev.mat24.mat24_ref import  Mat24
    mat24 = Mat24






ERR_TYPE = "unsupported operand types for %s: '%s' and '%s'"
ERR_RAND = "Illegal string for constricting type %s element" 



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
    global Cocode
    from mmgroup.structures.cocode import Cocode
    import_pending = False



####################################################################


def diagonal_from_obj(t):
    """Convert object 't' to a diagnal Parker loop automorphism.

    Ojects of the following types are recognized:
    - None and 0 encode the identity automorphism
    - An instance of class Cocode encodes a diagonal automorphism.
    - An element of the group AutPLoop encodes itself. 

    The function returns the cocode value
    """
    if isinstance(t, AutPL) and t.perm_num == 0:
        return t.cocode
    if import_pending:
        complete_import()
    if isinstance(t, Cocode):
        return t.value
    else:
        return Cocode(t).cocode
        



def autpl_from_obj(t):
    """Try to convert the object 't' to a Parker loop automorphism.

    Ojects of the following types are recognized:
    - None and 0 encode the identity automorphism
    - A list encodes a permutation in the Mathieu group M_24.
    - A zip object or a dictionary encodes such a permutation as
      a mapping. That mapping must contain sufficiently many 
      entries to be unique.  
    - An instance of class Cocode encodes a diagonal automorphism.
    - An element of the group AutPLoop encodes itself. 

    The function returns a pair (cocode, pemrutation_number).
    """
    global Cocode
    if t is None or t == 0 or t == 1:
        return 0, 0
    if isinstance(t, list):
        if mat24.perm_check(t):
            raise ValueError("Permutation is not in Mathieu group M_24")
        return 0, mat24.perm_to_m24num(t)
    if isinstance(t, zip):
        t = dict(t)
    if isinstance(t, dict):
        h1, h2 = [list(t) for t in zip(*t.items())]
        res, perm = mat24.perm_from_map(h1, h2)
        if res == 1:
            return 0, mat24.perm_to_m24num(perm) 
        if res < 1:
            err = "Permutation is not in Mathieu group M_24"
        else:
            err = "Permutation in Mathieu group M_24 is not unique"
        raise ValueError(err)
    if isinstance(t, AutPL):
        return t._cocode, t._perm_num
    if isinstance(t, tuple) and len(t) and isinstance(t[0], str):
        return autpl_from_tag(*t)
    if import_pending:
        complete_import()
    if isinstance(t, Cocode):
        return t.value, 0
    if isinstance(t, str):
        if len(t) == 1 and t.islower():
            return Cocode(t).cocode, randint(0, mat24.MAT24_ORDER-1)
        else:
            a = AutPlGroup.parse(AutPL.group, t)
            if isinstance(a, AutPL):
                return a.cocode, a.perm_num
            elif a == 1:
                return 0, 0
    err = "Cannot convert type %s object to Parker loop automorphism"
    raise TypeError(err % type(t))
        


def autpl_from_tag(tag = None, data = None):
    if tag == "p":
        if isinstance(data, Integral):
            if not 0 <= data < mat24.MAT24_ORDER:
                raise ValueError(self.ERR_PERM_NUM)
            return 0, data
        elif isinstance(data, str) or data is None:
            return  0, randint(0, mat24.MAT24_ORDER - 1)
        else:
            return  autpl_from_obj(data)
    elif tag == "d":
        return diagonal_from_obj(data), 0
    else:
        ERR_TAG = "Bad tag '%s' for Parker loop automorphism"
        raise ValueError(ERR_TAG % tag)


            



class AutPL(AbstractGroupWord):
    r"""This class models a standard automorphism of the Parker loop.

    :param \*data:

      A variable number of arguments; each argument describes a
      Parker loop automorphism. These automorphisms are
      multiplied.  

    :type \*data: see table below for legal types

    :param \*\*kwds: for internal use only.


    :return: A standard Parker loop automorphism
    :rtype:  an instance of class |AutPL|

    :raise:
        * TypeError if ``type(value)`` is not in the table given below.
        * ValueError if an argument cannot be converted to an
          instance of class  |AutPL|.

    Depending on its type each parameter in **\*data** is  
    interpreted as follows:

    .. table:: Legal types for constructor of class ``AutPL``
      :widths: 25 75

      ===================== ==================================================
      type                  Evaluates to
      ===================== ==================================================
      tuple (``'p', n``)    Here ``n`` is the number of a permutation in 
                            the Mathieu group ``M_24``. Then the
                            standard representative of the corresponding
                            permutation in ``M_24`` is returned. 
                            
                            If ``n`` is a list, or a zip object or a
                            dictionary, this is also interpreted as a 
                            permutation. If ``n`` is the string 
                            ``'r'`` then a random permutation in ``M_24`` 
                            is generated.

      ``list`` of ``int``   A list ``l_p`` of disjoint integers
                            ``0 <= i < 24`` is interpreted as a permutation
                            in ``M_24`` that maps ``i`` to ``l_p[i]``.

      ``dict``              A dictionary specifies a mapping from a subset
                            of the integers ``0 <= i < 24`` to integers
                            ``0 <= dict[i] < 24``. This mapping must
                            extend to a unique permutation in ``M_24``. 
                            The standard representative of that 
                            permutation is returned.

      ``zip`` object        ``zip(x,y)`` is equivalent to 
                            ``dict(zip(x,y))``.

      class |Cocode|        This represents the *diagonal* automorphism of 
                            the Parker loop given by the |Cocode| element.

      tuple (``'d', n``)    Equivalent to ``Cocode(n)``. 

      class |AutPL|         A deep copy of the given
                            automorphism in |AutPL| is returned.

      ``str``               Create random automorphism depending on ``str``
                              | ``'r'``: Create an arbitrary automorphism
                              | ``'e'``: Create an even automorphism
                              | ``'o'``: Create an odd automorphism

                            For an instance ``g`` of |AutPL| we have
                            ``AutPL(str(g)) == g``. This is helpful for
                            rereading printed instances of |AutPL|.
      ===================== ==================================================

    Let ``a`` be an instance of class |GcVector|, |GCode|, |Cocode|,
    or |PLoop|, and let ``g1`` , ``g2`` be instances of class |AutPL|. 
    Then  ``a * g1`` is the result of the natural operation of ``g1`` 
    on  ``a``, which belongs to the same class as ``a``.

    ``g1 * g2``  means group multiplication, and ``g1 ** n`` means
    exponentiation of ``g1`` with the integer ``n``. ``g1 ** (-1)`` 
    is the inverse of ``g``. ``g1 / g2`` means ``g1 * g2 ** (-1)``.
    We have ``1 * g1 == g1 * 1 == g1`` and ``1 / g1 == g1 ** (-1)``.

    ``g1 ** g2`` means ``g2**(-1) * g1 * g2``.    
    """
    __slots__ = "_cocode", "_perm_num", "_perm", "rep"
    ERR_UNDEF = "Parker loop automorphism is not defined by input data"
    _perm_ = list(range(24))  # neutral permutation
    _rep_ = mat24.perm_to_autpl(0, _perm_)
    group = None       # will be set to StdAutPlGroup later
    
    def __init__(self, *data, **kwds):
        if len(data) == 0:
            self._cocode = 0
            self._perm_num = 0
            self._perm = self._perm_
            self.rep = self._rep_
        else:
            self._cocode, self._perm_num = autpl_from_obj(data[0])
            self._compute_from_numbers()
            for d in data[1:]:
                cocode, perm_num = autpl_from_obj(d) 
                perm = mat24.m24num_to_perm(perm_num)  
                rep = mat24.perm_to_autpl(cocode, perm)
                self.rep = mat24.mul_autpl(self.rep, rep)
            if len(data) > 1:
                self._compute_from_rep()

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
        numbered from ``0`` to ``244823039``. The identity has
        number ``0``. The numbering system is not part of the 
        public interface; but with this property the numbers of
        individual permutations may be obtained.
        """
        return self._perm_num

    @property    
    def perm(self):
        """Return permutation of automorphism as a list ``p_l``.

        Then the permutation maps ``i`` to ``p_l[i]`` for
        ``i = 0,...,23``.
        """
        return self._perm


def autpl_element_from_obj(g, t):
    res =  AutPL()
    res._cocode, res._perm_num = autpl_from_obj(t)
    res._compute_from_numbers()
    return res


class AutPlGroup(AbstractGroup):
    word_type = AutPL              # type of an element (=word) in the group
    atom_parser = {}               # see method parse()
    rand_masks  = {"r":(0xfff,0), "e":(0x7ff,0), "o":(0x7ff,0x800)} 
    conversions = {
        list: autpl_element_from_obj,
        zip: autpl_element_from_obj,
        dict: autpl_element_from_obj,
        AutPL: autpl_element_from_obj,
    }
    FRAME =  re.compile(r"^(?:AutPL)?<(.*)>$") # see method parse()
    ERR_PERM_NUM = "Illegal permutation number for Mathieu group M_24"

    def __init__(self):
        super(AutPlGroup, self).__init__()
        self.atom_parser = AtomDict(self.atom)

    def atom(self, tag = None, data = None):
        if isinstance(tag, str):
            cocode, perm_num = autpl_from_tag(tag , data)
        else:
            cocode, perm_num = autpl_from_obj(tag)
        res = AutPL()
        res._cocode = cocode
        res._perm_num = perm_num 
        res._compute_from_numbers()
        return res        

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

    def as_tuples(self, g):
        return [('d', g._cocode), ('p', g._perm_num)]
        

    def str_word(self, g):
        """Convert group atom g to a string

        """
        s = []
        if g._cocode:
            s.append("d_" + ihex(g._cocode))
        if g._perm_num:
            s.append("p_%d" % g._perm_num)
        s = "*".join(s) if len(s) else "1"
        return "AutPL<%s>" % s



StdAutPlGroup = AutPlGroup()   # This is the only instance of AutPlGroup

AutPL.group = StdAutPlGroup

