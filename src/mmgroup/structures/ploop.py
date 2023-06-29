r"""We deal with the Parker loop, which is a non-associative Moufang loop.

The Parker loop is written multiplicatively its elements
can be represented as pairs in
the set :math:`\mathcal{P} = \mathcal{C} \times \mathbb{F}_2`, see
e.g. :cite:`Asc86`, :cite:`Con85`.
We also write ``1`` for the neutral element 
:math:`(0,0) \in  \mathcal{C} \times \mathbb{F}_2 = \mathcal{P}` of
the Parker loop and ``-1`` for its negative :math:`(0,1)`. 
Let :math:`\Omega` = :math:`(\tilde{\Omega}, 0)`, where
:math:`\tilde{\Omega}` is the Golay code word ``1,...,1``. Then the 
center :math:`Z(\mathcal{P})` of :math:`\mathcal{P}` is 
:math:`\{1, -1, \Omega, -\Omega\}`.  An element :math:`(g,s)`,
:math:`g \in \mathcal{C}` is called *positive* if :math:`s = 0` 
and *negative* if :math:`s = 1`. 

Multiplication of two Parker loop elements 
:math:`(g_1,s_1), (g_2,s_2)` is given by a cocycle
:math:`\theta: \mathcal{C} \times \mathcal{C} \rightarrow \mathbb{F}_2`
as follows:

 :math:`(g_1,s_1) \cdot (g_2,s_2) = (g_1 + g_2, s_1 + s_2 + \theta(g_1, g_2))`.

There are several possibilities for the cocycle :math:`\theta`. We
choose a cocycle :math:`\theta` satisfying the requirements given 
by :cite:`Seysen20`, see :ref:`basis-golay-label` for details.


We represent elements of :math:`\mathcal{P}` as instances of 
class |PLoop|. For Golay code word  ``g`` given as an instance of 
class |GCode|, the value ``PLoop(g)``  is the positive element 
``(g,0)`` of the Parker loop, and  ``-PLoop(g)`` is the corresponding 
negative element ``(g,1)``. The constant  ``PLoopOne`` is equal
to the neutral element ``(0,0)`` and the constant ``PLoopOmega``
is equal to the central element ``~PLoopOne``. Here the ``~``
means bitwise complement of the corresponding Golay code element
without changing the sign of a Parker loop element. Thus 
``PLoopOmega`` is the positive Parker loop element corresponding
to the Golay code word containing ``24`` one bits. So the center 
:math:`Z(\mathcal{P})` of :math:`\mathcal{P}` is:

  ``[PLoopOne, -PLoopOne, PLoopOmega, -PLoopOmega]`` .

Let ``g1, g2`` be instances of 
class  |GCode|. For positive elements of the Parker loop the 
multiplication formula given above can be coded as follows: 

  ``PLoop(g1) * PLoop(g2) == (-1)**g1.theta(g2) * PLoop(g1 + g2)`` .

The elements of the Parker loop are numbered from ``0`` to 
``0x1fff``. Here the numbers of the  positive Parker loop
elements correspond to the numbers of the Golay code words.
The numbering system for the Parker loop is also used for 
indexing certain elements of the monster group :math:`\mathbb{M}`.

A transversal of :math:`Z(\mathcal{P})` in :math:`\mathcal{P}` is 
used for indexing certain basis vectors of the representation 
:math:`\rho` of :math:`\mathbb{M}`. Property ``ord`` of an 
instance ``a`` of class |Ploop| returns the number of that element 
of the Parker loop. 
For indexing vectors in  :math:`\rho` it is important to 
know that ``(-a).ord =  a.ord ^ 0x1000`` and 
``(~a).ord = a.ord ^ 0x800``. Thus the Parker loop elements
``a`` with ``0 <= a.ord < 0x800`` are a transversal of 
:math:`Z(\mathcal{P})` in :math:`\mathcal{P}`.
"""




from functools import reduce
from operator import __xor__
from numbers import Integral, Number
from random import randint



from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.parity import Parity
from mmgroup.structures.parse_atoms import ihex


from mmgroup.structures.gcode import GCode, GcVector
from mmgroup.structures.cocode import Cocode


ERR_RAND = "Illegal string for constricting type %s element" 

ERR_DIV4 = "%s object may be divided by 2 or 4 only"
ERR_DIV2 = "%s object may be divided by 2 only"


#######################################################################
# Import derived classed
#######################################################################

import_pending = True

def complete_import():
    """Internal function of this module

    If you need any of the objects declared above, do the following:

    if import_pending:
        complete_import()
    """
    global import_pending
    global mat24, AutPL, AutPlGroup, XLeech2
    from mmgroup import mat24
    from mmgroup.structures.autpl import AutPL, AutPlGroup
    from mmgroup.structures.xleech2 import XLeech2
    import_pending = False



#######################################################################
# Auxiliary functions
#######################################################################



def as_vector24(data):
    """Convert list of bit positions to vector in ``GF(2)**24``"""
    try:
        vector = reduce(__xor__, (1 << x for x in data), 0)
    except:
        err = "List of integers >= 0 required for a 24-bit vector"
        raise TypeError(err)
    if vector >= 0x1000000:
        raise ValueError("List entry too large for a 24-bit vector")
    return vector


def str_vector24(v):
    """Return the 24-bit vector ``v`` as a string of binary digits

    Here ``v`` must be an integer with the coefficient of the
    ``i``-th bit coded in the bit with valence ``2**i``.
    Coefficients  ``i >= 24`` are ignored.
    """
    l = ["01"[(v >> i) & 1] for i in range(24)]
    l4 = ("".join(l[i:i+4]) for i in range(0, 24, 4))
    return " ".join(l4)

def str_basis(text, letter, basis):
   s = text + "\n"
   for i, v in enumerate(basis):
       s += "  %3s: %s\n" % (letter + str(i), str_vector24(v)) 
   return s






#######################################################################
# Subclasss PLoop of class GCode
#######################################################################





class PLoop(GCode):
    """This class models an element of the Parker Loop.

    The Parker loop is a non-associative loop with  :math:`2^{1+12}` 
    elements. Each element can be interpreted as a signed Golay
    code word, see class |GCode|. Accordingly, class |PLoop| is
    a subclass of class |GCode|.

    The loop operation is written
    multiplicatively. The neutral element of the loop (corresponding
    to the positive zero word of the Golay code) is ``PLoopOne``, and 
    its negative is ``-PLoopOne``.  ``PLoopOmega`` is the (positive) 
    element corresponding  to the all-one vector of the Golay code. 
        
    The :math:`2^{1+12}` Parker loop elements are numbered from
    ``0`` to ``0x1fff``. Elements ``0`` to ``0xfff`` are are
    positive and corresponding to the Golay code words with the 
    same number. The element with number  ``0x1000 ^ i`` is the
    negative of the element with number  ``i``.

    :param value:

      This is the value of the Parker loop element to be returned. 

    :type value: see table below for legal types

    :return: A Parker loop element
    :rtype:  an instance of class |PLoop|

    :raise:
        * TypeError if ``type(value)`` is not in the table given above.
        * ValueError if ``value`` cannot be converted to an
          instance of class  |PLoop|.


    Depending on its type parameter **value** is  interpreted as follows:

    .. table:: Legal types for constructor of class ``PLoop``
      :widths: 20 80

      ===================== ================================================
      type                  Evaluates to
      ===================== ================================================
      ``int``               Here the code word with number ``value`` is
                            returned.  ``0 <= value < 0x2000`` must hold.
                            We have ``PLoop(0) == PLoopOne``,
                            ``PLoop(0x1000) == -PLoopOne``, and
                            ``PLoop(0x800) == ~PLoopOne == PLoopOmega``.
  
      ``list`` of ``int``   Such a list is converted to a Golay code word,
                            see class |GCode|, and the corresponding 
                            (positive) Parker loop element is returned.

      class |GCode|         The corresponding 
                            (positive) Parker loop element is returned. 

      class |PLoop|         A deep copy of parameter ``value`` is returned.

      class |GcVector|      This is converted to a Golay code word,
                            see class |GCode|, and the corresponding 
                            (positive) Parker loop element is returned.

      class |XLeech2|       The Parker loop part of the |XLeech2| object  
                            is  returned. 

      ``str``               Create random element depending on the string
                             | ``'r'``: Create arbitrary Parker loop element

      ===================== ================================================



    **Standard operations**

    Let ``a`` be a Parker loop element. 
    The multiplication operator ``*`` implements the non-associative
    loop operation. Division by a Parker loop element means multiplication
    by its inverse, and exponentiation means repeated multiplication,
    with ``a**(-1)`` the loop inverse of  ``a``,  as usual. 

    The   ``~``  operator maps the Parker loop element  ``a`` to 
    ``a * PLoopOmega``, leaving the sign of ``a`` unchanged.

    Multiplication with the integer ``1`` or ``-1`` means the 
    multiplication with the neutral element ``PLoopOne`` or 
    with its negative ``-PLoopOne``, respectively.
    
    If any of the operations  ``+``, ``-`` or ``&`` is applied to a Parker
    loop element, that element is converted to a Golay code word of 
    type |GCode|, ignoring the sign. This conversion also takes place if
    a Parker loop element is multiplied or divided by an integer
    different from ``1`` or ``-1``.
    

    **Standard functions**

    ``len(a)`` is a shorthand for ``len(GCode(a))``.
    Here function ``GCode(a)`` converts  ``a`` to
    a Golay code word of type |GCode|. 

    
    ``abs(a)`` returns the element in the set ``{a, -a}`` which is
    positive.
    """
    __slots__ = "value", "bit_list_"
    parity_class = GCode

    def __init__(self, value = 0):
        if isinstance(value, Integral):
            self.value = value & 0x1fff
            return
        if import_pending:
            complete_import()
        if isinstance(value, PLoop):
            self.value = value.value & 0x1fff
        elif isinstance(value, GCode):
            self.value = value.value & 0xfff
        elif isinstance(value, str):
            self.value = randint(0, 0x1fff)
        elif isinstance(value, GcVector):
            vector = value.value
            vector ^= mat24.syndrome(vector) 
            self.value = mat24.vect_to_gcode(vector)
        elif isinstance(value, XLeech2):
            self.value = (value.value >> 12) & 0x1fff
        else:
            vector = as_vector24(value)
            vector ^= mat24.syndrome(vector) 
            self.value = mat24.vect_to_gcode(vector)

    def __mul__(self, other):
        if isinstance(other, PLoop):
            return PLoop(mat24.mul_ploop(self.value, other.value))
        elif isinstance(other, Integral):
            if abs(other) == 1:
                return PLoop(self.value ^ ((other & 2) << 11))
            return GCode(self.value & -(other & 1))
        elif isinstance(other, AutPL):
            return PLoop(mat24.op_ploop_autpl(self.value, other.rep))
        else:           
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Integral):
            if abs(other) == 1:
                return PLoop(self.value ^ ((other & 2) << 11))
            return GCode(self.value & -(other & 1))
        else:           
            return NotImplemented
                              
    def __truediv__(self, other):
        if isinstance(other, PLoop):
            return PLoop(mat24.mul_ploop(self.value, 
                mat24.pow_ploop(other.value, 3)))
        elif isinstance(other, Integral):
            if abs(other) == 1:
                return PLoop(self.value ^ ((other & 2) << 11))
            elif abs(other) == 2:
                return Parity(0)
            elif abs(other) == 4:
                return Parity(mat24.gcode_weight(self.value) & 1)
            else:
                raise TypeError(ERR_DIV4 % type(self))
        else:           
            return NotImplemented
            
    def __rtruediv__(self, other):
        if isinstance(other, Integral):
            if abs(other) == 1:
                inv = mat24.pow_ploop(self.value, 3)
                return PLoop(inv ^ ((other & 2) << 11))
            else:
                raise TypeError("Can only divide unit by Parker loop element")
        else:           
            return NotImplemented
        
        
    def __pow__(self, other):
        if isinstance(other, Integral):            
            return PLoop(mat24.pow_ploop(self.value, other & 3))
        elif isinstance(other, PLoop):
            comm = mat24.ploop_comm(self.value, other.value) << 12
            return PLoop(self.value ^ comm)
        elif isinstance(other, Parity):
            if mat24.gcode_weight(self.value) & 1:
                raise ValueError("Parker Loop element has order > 2")
            return PLoop(mat24.pow_ploop(self.value, other.value))
        else:
            return NotImplemented

    def __eq__(self, other):
        return (isinstance(other, PLoop)  
            and (self.value ^ other.value) & 0x1fff == 0)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __pos__(self):
        return self
        
    def __neg__(self):
        return  PLoop(self.value ^ 0x1000)

    def __invert__(self):
        return  PLoop(self.value ^ 0x800)

    @property
    def ord(self):
        """Return the number of the Parker loop element.

        We have ``0 <= i < 0x2000`` for the returned number ``i``.
        """
        return self.value & 0x1fff


    @property
    def sign(self):
        """Return the sign of the Parker loop element.

        This is ``1`` for a positive and ``-1`` for a negative element.
        """
        return 1 - ((self.value >> 11) & 2)


    def split(self):
        """Split sign and Omega from Parker loop element ``a``. 

        :return:  a triple  ``(es, eo, v)`` with 
                  ``a = (-1)**e1 * PLoopOmega**eo * v`` .

        Here  ``v`` is the unique (positive) element of the Parker  
        loop of type |PLoop| with ``0 <= v.ord < 0x800`` satisfying
        that equation. ``es`` and ``eo`` are ``0`` or ``1``.
        """
        v = self.value
        return (v >> 12) & 1, (v >> 11) & 1, PLoop(v & 0x7ff)
 
    def split_octad(self):
        """Split Parker loop element ``a`` into central element and octad

        :return: a triple ``(es, eo, o)`` with  
                 ``a = (-1)**e1 * PLoopOmega**eo * o``.

        Here ``o`` is either the neutral Parker loop element
        ``PLoopOne`` or Parker loop element corresponding to a
        positive octad.  ``es`` and ``eo`` are ``0`` or ``1``.
        An *octad* is a Golay code word 
        (of type |GCode|)  of weight ``8``.

        :raise:
          * Raise ValueError if this is not possible.
        """
        v = self.value
        e1, eo, w = (v >> 12) & 1, 0, mat24.gcode_weight(v)
        if mat24.gcode_weight(v) > 3:
            v, eo, w = v ^ 0x800, 1, 6 - w
        if w <= 2:
            return e1, eo, PLoop(v & 0xfff)
        raise ValueError("Cannot convert Golay code word to octad")
     
    def str(self):
        return "<PLoop_%s>" % ihex(self.value & 0x1fff, 3)
    __repr__  = str


def PLoopZ(e1 = 0, eo = 0):
    """Return a specific central element of the Parker loop. 

    :param e1: Sign, the sign of the element will be ``(-1)**e1``
    :type e1: int
    :param eo: Exponent for ``PLoopOmega``, default is  ``0``
    :type eo: int
    :return:   The central element ``(-1)**e1  * PLoopOmega**eo``
               of the Parker loop.
    :return type:  an instance of class |PLoop|.

    Here ``PLoopOmega`` is the positive element of the Parker loop
    corresponding to the Golay code word ``1,...,1``.

    ``e1`` and ``eo`` may also be parity objects of type |Parity|.
    """
    e1 = e1.value if isinstance(e1, Parity) else e1
    eo = eo.value if isinstance(eo, Parity) else eo
    return PLoop(((e1 & 1) << 12) | ((eo & 1) << 11))


PLoopOne = PLoop(0)
PLoopOmega = ~PLoopOne    



