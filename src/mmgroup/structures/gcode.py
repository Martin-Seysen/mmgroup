r"""We deal with the Golay code and its cocode.

Let :math:`\tilde{\Omega}` be the set of integers
:math:`\{ i \in \mathbb{Z} \mid 0 \leq i < 24\}` of size 
:math:`24` and construct the vector space
:math:`\mathcal{V} = \mathbb{F}_2^{24}` as 
:math:`\prod_{i\in \tilde{\Omega}} \mathbb{F}_2`. In this 
document elements of :math:`\mathcal{V}` are called *bit vectors*.
We represent bit vectors as instances of class |GcVector|.

We identify  the power set of :math:`\tilde{\Omega}` with the
vector space :math:`\mathcal{V} = \mathbb{F}_2^{24}` by mapping 
each subset of :math:`\tilde{\Omega}` to its characteristic 
function, which is a bit vector in :math:`\mathcal{V}`.  So we may 
talk about union and intersection of bit vectors in a natural way 
and use the python operators ``|`` and ``&`` for these operations 
as usual. We use the  ``+`` operator for vector addition of bit 
vectors. Bit vectors are numbered from ``0`` to ``0xffffff``, with 
bit ``i`` of that number corresponding to the ``i``-th unit vector.

The Golay code :math:`\mathcal{C}` is a :math:`12`-dimensional 
linear subspace of :math:`\mathcal{V}` such that a nonzero 
vector in :math:`\mathcal{C}` has weight at least  :math:`8`.
This characterizes the Golay code up to permutation. In this
package we chose a specific basis of the Golay code that
satisfies the requirements in :cite:`Seysen20`, see
:ref:`basis-golay-label` for details. We represent Golay
code words as instances of class |GCode|. 
Golay code words of weight :math:`8` and :math:`12` are called 
*octads* and *dodecads*, respectively. 

We stress that  class |GCode| is not a subclass of class |GcVector|.
Instances of these two classes are numbered in a different way.
Golay code words are numbered from ``0`` to ``0xfff`` with bit
``i`` of that number corresponding to the ``i``-th basis vector
of the Golay code. See :ref:`basis-golay-label` for the chosen 
basis of the Golay code. 

The Golay cocode :math:`\mathcal{C}^*` is the ``12``-dimensional 
quotient space :math:`\mathbb{F}_2^{24} / \mathcal{C}`. For any
:math:`\delta \in \mathcal{C}^*` the weight of :math:`\delta`
is the weight of its lightest representative in
:math:`\mathbb{F}_2^{24}`.  That weight is at most ``4``.
There is a unique lightest representative  it the weight is
less than ``4``, and there are six mutually disjoint  
lightest representatives of an element of weight  ``4``.
A partition of  :math:`\mathbb{F}_2^{24}` given by such a set of 
disjoint representatives is called a *tetrad*. We represent 
elements of the cocode of the Golay code as instances of class 
|Cocode|. Cocode elements are numbered from ``0`` to ``0xfff`` 
with bit ``i`` of that number corresponding to the ``i``-th basis 
vector of the cocode. Our basis of the cocode is the reciprocal 
basis of the basis of the Golay code, see :ref:`basis-golay-label` 
for details.

The *scalar product* of a Golay code word and a cocode element is 
the parity of the intersection of the Golay code word and any 
representative of the cocode element. If ``g`` is an instance of 
class |GCode| and ``d`` is an instance of class |Cocode| then
``g & d`` is an instance of class |Parity|. Here class |Parity|
models an integer modulo ``2``, which is also considered as the
parity of a bit vector. Standard operations, e.g. ``(-1)**p``,
work as expected for an instance ``p`` of class |Parity|.

The standard function ``len(x)`` returns the (minimal) weight 
``|x|`` of an instance ``x`` of class |GcVector|, |GCode| or 
|Cocode|. For instances ``g1, g2, g3`` of class  |GCode|, the
expression ``g1 & g2`` is an instance of class  |Cocode| and
``g1 & g2 & g3`` an instance of class |Parity|; these values 
have a natural interpretations as bitwise intersections.  
``g1/4`` is a shorthand for ``Parity(len(g1)/4)``; and
``(g1 & g2)/2`` is a shorthand for ``Parity(len(g1 & g2)/2)``. 
These two expressions are called the *power map* and the
*commutator*, respectively, in  :cite:`Asc86`;
and ``g1 & g2 & g3`` is called the *associator* in 
:cite:`Asc86`.

The numbering of the Golay code words and of the cocode elements 
is important for indexing some generators of 
:math:`\mathbb{M}` and also some basis vectors in the
representations  :math:`\rho_p` of :math:`\mathbb{M}`. 
"""

from functools import reduce
from operator import __xor__
from numbers import Integral, Number
from random import randint




from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.parity import Parity
from mmgroup.structures.parse_atoms import ihex



ERR_RAND = "Illegal string for constructing type %s element" 

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
    global BASIS, mat24
    global import_pending, Cocode, PLoop, PLoopIntersection
    global AutPL, AutPlGroup, XLeech2
    from mmgroup import mat24, GCODE_BASIS
    from mmgroup.structures.cocode import Cocode
    from mmgroup.structures.cocode import PLoopIntersection
    from mmgroup.structures.ploop import PLoop
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
# Class GCode 
#######################################################################



class GCode():
    """This class models a code word of the Golay code.

    The Golay code is a binary ``[24, 12, 8]`` code. So there are 
    :math:`2^{12}` code words, each word is :math:`24` bit long, and
    the  Hamming distance between two different code words is at least 
    :math:`8`.

    The :math:`2^{12}` Golay code words are numbered from
    ``0`` to ``0xfff``. Bit ``i`` in that number refers to the
    coefficient of the  ``i``-th basis vector of the Golay code,
    which may be :math:`0` or :math:`1`.


    :param value:

      This is the value of the Golay code word to be returned. 

    :type value: see table below for legal types

    :return: A Golay code vector 
    :rtype:  an instance of class |GCode|

    :raise:
        * TypeError if ``type(value)`` is not in the table given below.
        * ValueError if ``value`` cannot be converted to an
          instance of class  |GCode|.


    Depending on its type parameter **value** is  interpreted as follows:

    .. table:: Legal types for constructor of class ``GCode``
      :widths: 20 80

      ====================== ================================================
      type                   Evaluates to
      ====================== ================================================
      ``int``                Here the code word with number ``value`` is
                             returned.  ``0 <= value < 0x1000`` must hold.
  
      ``list`` of ``int``    A list of integers  ``0 <= i < 24`` is 
                             interpreted as a list of bit positions
                             to be set in a bit vector. The Golay code word
                             corresponding to that bit vector is returned.
                             Up to ``3`` erroneous bits are corrected.

      class |GCode|          A deep copy of parameter ``value`` is returned. 
 
      class |PLoop|          The Parker loop element ``value`` is converted
                             to a Golay code word (by dropping its sign)
                             and returned.

      class |GcVector|       The bit vector ``value`` of type |GcVector| is 
                             converted to a Golay code word and returned. 
                             Up to ``3``  erroneous bits
                             in that vector are corrected.

      class |XLeech2|        The Golay code part of the |XLeech2| object  
                             is converted to a Golay code word and returned. 

       ``str``               Create random element depending on the string
                              | ``'r'``: Create arbitrary Golay code word
      ====================== ================================================


    **Standard operations**
    
    Addition of Golay code words and scalar multiplication (with scalars 
    of type  ``int`` or |Parity|) are done as usual.

    The bitwise and operator ``&`` means bitwise *and* of two code words
    of type  |GCode|. The result of such an operation on two Golay code
    words is a Golay cocode word of type |Cocode|.

    The bitwise ``&`` operator of a Golay code word and a Golay cocode
    word of type |Cocode| is not well defined as a bit vector, but it
    has a well-defined parity. So in this case the result of the 
    ``&`` operator is a parity object of type |Parity|.

    For the bitwise ``&`` operator bit vectors of type |GcVector|, see
    class |GcVector|.

    The   ``~``  operator means bitwise complement  as usual; it
    returns the complemented code word, which is also of type |GCode|.

    **Standard functions**

    ``len(g)`` returns the bit weight of the Golay code word ``g``, i.e.
    the number of bits set in the word ``g``. 


    **Special functions**

    Let ``g1``,  ``g2``, and ``g3`` be Golay code vectors of 
    type |GCode|.
  
    The expression ``g1/4`` is a shorthand for  ``Parity(len(g1)/4)``.
    It returns   ``len(g1)/4``  modulo ``2`` as a parity
    object of type  |Parity|.  
    
    ``GcVector(g1 & g2)`` returns the bitwise intersection of 
    ``g1`` and ``g2`` as a bit vector of type  |GcVector|.

    ``(g1 & g2)/2`` is a shorthand for 
    ``Parity(len(GcVector(g1 & g2))/2)``. It returns the
    halved bit weight of ``g1 & g2`` modulo ``2`` as a parity
    object of type  |Parity|.      
    This value is also called the *commutator* of ``g1`` and ``g2``. 

    ``g1 & g2 & g3`` returns the parity of the bitwise intersection 
    of   ``g1``,  ``g2``, and ``g3`` as a parity object
    of type  |Parity|. This value is also called the
    *associator* of ``g1``, ``g2``, and ``g3``.

    ``abs(g1)`` is a shorthand for ``abs(PLoop(g1))``, see 
    class   |PLoop| for details.

    """
    __slots__ = "value", "bit_list_"
    parity = 0
    cocode = 0
    sign = 1
    parity_class = True

    def __init__(self, value):
        """Constructor"""
        if import_pending:
            complete_import()
        if isinstance(value, Integral):
            self.value = value & 0x0fff
        elif isinstance(value, GCode):
            self.value = value.value & 0xfff
        elif isinstance(value, GcVector):
            vector = value.value
            vector ^= mat24.syndrome(vector) 
            self.value = mat24.vect_to_gcode(vector)
        elif isinstance(value, XLeech2):
            self.value = (value.value >> 12) & 0x1fff
        elif isinstance(value, str):
           self.value = randint(0, 0xfff)
        else:
            vector = as_vector24(value)
            vector ^= mat24.syndrome(vector) 
            self.value = mat24.vect_to_gcode(vector)       


    def __eq__(self, other):
        return (isinstance(other, GCode)  
            and ((self.value ^ other.value) & 0xfff) == 0)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __pos__(self):
        return self
        
    __neg__ = __pos__

    def __and__(self, other):
        if import_pending:
            complete_import()
        if isinstance(other, GCode):
            return PLoopIntersection(self, other)
        elif  isinstance(other, Cocode):
            return Parity(mat24.scalar_prod(self.value, other.value))
        elif isinstance(value, GcVector):
            return  GcVector(mat24.gcode_to_vect(self.value) & other.value)
        else:
            return NotImplemented
            

    def __add__(self, other):
        if (isinstance(other, GCode)):
            return GCode((self.value ^ other.value) & 0xfff)
        elif isinstance(other, GcVector):
            return  GcVector(mat24.gcode_to_vect(self.value) ^ other.value)
        elif other == 0:
            return GCode(self)
        else:
             return NotImplemented


    __radd__ = __add__
    __sub__ = __add__
    __rsub__ = __add__


    def __mul__(self, other):
        if import_pending:
            complete_import()
        if isinstance(other, Integral):
            mask = -(other & 1)
            return GCode(self.value & mask & 0xfff)
        elif isinstance(other, AutPL):
            return GCode(mat24.op_gcode_perm(self.value, other.perm))
        else:
             return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Integral):
            mask = -(other & 1)
            return GCode(self.value & mask & 0xfff)
        else:
            return NotImplemented

    def __len__(self):
        if import_pending:
            complete_import()
        return mat24.gcode_weight(self.value) << 2

    def __abs__(self):
        return PLoop(self.value & 0xfff)

    def __invert__(self):
        return  GCode(self.value ^ 0x800)

    def __truediv__(self, other):
        if import_pending:
            complete_import()
        if isinstance(other, Integral):
            other = abs(other)
            if other == 1:
                return self
            elif other == 2:
                return Parity(0)
            elif other == 4:
                return Parity(mat24.gcode_weight(self.value) & 1)
            else:
                raise ValueError(ERR_DIV4 % type(self))
        else:
             return NotImplemented

    def __mod__(self, other):
        return Parity(self.parity) % other


    @property
    def ord(self):
        """Return the number of the Golay code word.

        We have ``0 <= i < 0x1000`` for the returned number ``i``.
        """
        return self.value & 0xfff


    @property
    def octad(self):
        """Return number of octad of the Golay code word as an octad.
        
        Complements of octads are also accepted as octads.   We have 
        ``0 <= o < 759`` for the returned number ``o``.

       :raise:
        * ValueError if the Golay code word is not an octad.


        """
        if import_pending:
            complete_import()
        return mat24.gcode_to_octad(self.value & 0xfff, 0)

    @property
    def gcode(self):
        """Return the number of the Golay code word.

        Same as property ``ord``.
        """
        return self.value & 0xfff
        
    @property    
    def vector(self):
        """Return bit vector corresponding to the Golay code word as an int.

        We have ``0 <= v < 0x1000000`` for the returned number ``v``.
        Bit ``i`` of the number ``v`` is the  ``i``-th coordinate
        of the corresponding bit vector.
        """
        if import_pending:
            complete_import()
        return mat24.gcode_to_vect(self.value)
        
    @property    
    def bits(self):
        """Return the Golay code word as a ``24``-bit vector.

        The function returns a list with  ``24`` entries equal to
        ``0`` or ``1``.
        """
        try:
            return self.bit_vector
        except:
            if import_pending:
                complete_import()
            v = mat24.gcode_to_vect(self.value)
            self.bit_vector = ([(v >> i) & 1 for i in range(24)]) 
            return self.bit_vector

    @property    
    def bit_list(self):
        """Return the ordered list of positions of bits being set"""
        try:
            return self.bit_list_
        except AttributeError:
            if import_pending:
                complete_import()
            self.bit_list_ = mat24.gcode_to_bit_list(self.value)
            return self.bit_list_
             

    def split(self):
        """Split the Golay code word ``g``. 

        :return: A triple ``(0, eo, v)`` with  ``g  =  Omega * eo + v``.

        Here ``Omega`` is the Golay code word containing  one bits only,
        ``eo`` is ``0`` or ``1``, and ``v`` is a Golay code word 
        of type |GCode| with ``0 <= v1.ord < 0x800``.

        This method is for compatibility with the corresponding method in
        class |PLoop|.
        """
        return 0, (v >> 11) & 1, GCode(v & 0x7ff)
 
    def split_octad(self):
        """Split an octad from the Golay code word ``g``. 

        :return: A triple ``(0, eo, o)`` with  ``g  =  Omega * eo + o``.

        Here ``Omega`` is the Golay code word containing  one bits only,
        ``eo`` is ``0`` or ``1``, and ``v`` is a Golay code word 
        of type |GCode| which is either the zero code word or an
        octad, i.e. a code word of weight ``8``.

        This method is for compatibility with the corresponding method 
        in class |PLoop|.

        :raise:
          * ValueError if this is not possible.
        """
        if import_pending:
            complete_import()
        v = self.value
        eo, w = 0, mat24.gcode_weight(v)
        if mat24.gcode_weight(v) > 3:
            v, eo, w = v ^ 0x800, 1, 6 - w
        if w <= 2:
            return 0, eo, GCode(v & 0xfff)
        raise ValueError("Cannot convert Golay code word to octad")
        
    def theta(self, g2 = None):
        """Return cocycle of Golay code words.

        The cocycle ``theta`` maps a pair of Golay code words 
        to an integer modulo ``2``. 
        It is linear in its second argument, so it
        may also be considered as a mapping from the Golay code
        to the cocode of the Golay code.


        :param g2:  ``None`` (default) or another Golay code word of
                    type |GCode|.
        :returns:
            * If ``g2`` is a code word of type |GCode|, we return
              the value ``g1.theta(g2) = theta(g1, g2)`` 
              as a |Parity| object.
            * If ``g2`` is ``None`` (default), we return the value
              ``g1.theta() = theta(g1)``  as a |Cocode| object.
              Note that  ``g1.theta(g2) == g1.theta() & g2 .``

        The importance of the ``theta`` function comes from the
        fact that the multiplication of elements of the Parker loop
        is based on the cocycle. We embed the set of Golay code words  
        into  the set of positive Parker loop elements, which are 
        instances of class |PLoop|. 
 
        Let ``g1`` and ``g2`` be Golay code words of type |GCode|.
        Then  ``PLoop(g1)`` and ``PLoop(g2)`` are the corresponding
        positive Parker loop elements,  and  ``g1.theta(g2)`` is an 
        integer  modulo ``2`` of type |Parity|. We have:  
       
        ``PLoop(g1) * PLoop(g2) == (-1)**g1.theta(g2) * PLoop(g1 + g2) .``
        """
        if import_pending:
            complete_import()
        th = mat24.ploop_theta(self.value) 
        if g2 == None:
            complete_import()
            return Cocode(th)
        if isinstance(g2, GCode):
            return Parity(mat24.scalar_prod(g2.value, th))
        err = "Types %s is illegal for method theta()"
        raise TypeError(err % type(g2))

    def str(self):
        return "<GCode_%s>" % ihex(self.value & 0xfff, 3)
    __repr__ = str

    @classmethod
    def str_basis(cls):
       b = "Golay code basis g_i[0,...,23], i = 0,...,11"
       return str_basis(b, "g", GCODE_BASIS)

    @classmethod
    def show_basis(cls):
        print(cls.str_basis())







#######################################################################
# Class GcVector
#######################################################################

class GcVector:
    """Models a 24-bit vector in the underlying space of the Golay code.

    The :math:`2^{24}` bit vectors are numbered from
    ``0`` to ``0xffffff``. Bit ``i`` in that number refers to the
    coefficient of the  ``i``-th basis vector of the vector space,
    which may be :math:`0` or :math:`1`.

    :param value:

      This is the value of the bit vector to be returned. 

    :type value: see table below for legal types

    :return: A bit vector of 24 bit length.
    :rtype:  an instance of class |GcVector|

    :raise:
        * TypeError if ``type(value)`` is not in the table given above.
        * ValueError if ``value`` cannot be converted to an
          instance of class  |GcVector|.

    Depending on its type parameter **value** is  interpreted as follows:

    .. table:: Legal types for constructor of class ``GcVector``
      :widths: 20 80

      ===================== ================================================
      type                  Evaluates to
      ===================== ================================================
      ``int``               Here the bit vector with number ``value`` is
                            returned. 
                            ``0 <= value < 0x1000000`` must hold.

      ``list`` of ``int``   A list of integers  ``0 <= i < 24`` is 
                            interpreted as a list of bit positions
                            to be set in a bit vector. That bit vector 
                            is returned.

       class |GCode|        The bit vector corresponding to the Golay 
                            code word ``value`` is returned. 

       class |PLoop|        The Parker loop element ``value`` is converted
                            to a Golay code word (by dropping its sign)
                            and then to a bit vector. That bit vector is
                            returned.

       class |GcVector|     A deep copy of parameter ``value`` is returned.


      ``str``               Create random element depending on the string
                             | ``'r'``: Create arbitrary bit vector

      ===================== ================================================


    **Standard operations**
    
    Addition of bit vectors and scalar multiplication (with scalars of 
    type  ``int`` or |Parity|) are done as usual. The sum of a bit
    vector and a Golay code word (of type |GCode|) is a bit vector.
    The sum of a bit vector and a cocode word (of type |Cocode|) is a 
    cocode word of type |Cocode|.

    The bitwise and operators ``&``, ``|``, and ``~`` operate on two
    bit vectors as expected. If the other operand is a Golay code
    vector (of type |GCode|) then it is converted to a bit vector. 

    **Standard functions**

    ``len(v)`` returns the bit weight of the bit vector ``v``, i.e.
    the number of bits set in the vector ``v``. 
 
    """
    __slots__ = "value", "bit_list_", "bit_vector"
    parity_class = True
    sign = 1
    STR_VT = ["S%d","U%d","T%d","S%d+","U%d-"] 

    def __init__(self, value):
        if import_pending:
            complete_import()
        if isinstance(value, Integral):
            self.value = value & 0xffffff
        elif isinstance(value, GCode):
            self.value = mat24.gcode_to_vect(value.value)
        elif isinstance(value, GcVector):
            self.value =  value.value
        elif isinstance(value, PLoopIntersection):
            self.value = (mat24.gcode_to_vect(value.v1) & 
                mat24.gcode_to_vect(value.v2) & 0xffffff)
        elif isinstance(value, str):
            self.value = randint(0, 0xffffff)
            if 'e' in value and not 'o' in value:
                if mat24.bw24(self.value) & 1:
                    self.value ^= 1 << randint(0, 23)
            if 'o' in value and not 'e' in value:
                if not mat24.bw24(self.value) & 1:
                    self.value ^= 1 << randint(0, 23)
        else:
            self.value = as_vector24(value)
      

    def __add__(self, other):
        if isinstance(other, GcVector):
            return GcVector(self.value ^ other.value)
        elif isinstance(other, GCode):
            return GCode(self.value ^ mat24.gcode_to_vect(other.value))
        elif other == 0:
            return self
        else:
            if import_pending:
                complete_import()
            return Cocode(self).__add__(other)

    __radd__ = __add__
    __sub__ = __add__
    __rsub__ = __add__

    def __and__(self, other):
        if isinstance(other, GcVector):
            return GcVector(self.value & other.value)
        elif isinstance(other, GCode):
            return GCode(self.value & mat24.gcode_to_vect(other.value))
        else:
            if import_pending:
                complete_import()
            return Cocode(self).__and__(other)

    __rand__ = __and__

    def __or__(self, other):
        if isinstance(other, GcVector):
            return GcVector(self.value | other.value)
        elif isinstance(other, GCode):
            return GCode(self.value | mat24.gcode_to_vect(other.value))
        else:
            return NotImplemented

    __ror__ = __and__

    def __invert__(self):
        return GcVector(self.value ^ 0xffffff)

    def __eq__(self, other):
        return (isinstance(other, GcVector)  
            and ((self.value ^ other.value) & 0xffffff) == 0)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __pos__(self):
        return self

    __neg__ = __pos__

    def __mul__(self, other):
        if isinstance(other, Integral):
             return GcVector(self.value & -(other & 1))
        elif isinstance(other, AutPL):
            return GcVector(mat24.op_vect_perm(self.value, other.perm))
        else:
            return NotImplemented
      
    def __rmul__(self, other):
        if isinstance(other, Integral):
             return GcVector(self.value & -(other & 1))
        else:
            return NotImplemented
     
    def __truediv__(self, other):
        if isinstance(other, Integral):
            other = abs(other)
            if other == 1:
                return self
            w = mat24.bw24(self.value)
            if other == 2 and  w & 1 == 0:
                return Parity(w >> 1)
            elif other == 4 and  w & 3 == 0:
                return Parity(w >> 2)
            elif other in (2,4):
                err = "Don't know weight of Omega vector / %d"
                raise ValueError(err % other)
            else:
                raise TypeError(ERR_DIV4 % self.__class__)
        else:
            return NotImplemented

    def __mod__(self, other):
        return Parity(self.parity) % other


    def __len__(self):
        return mat24.bw24(self.value)

    @property
    def ord(self):
        return self.value & 0xffffff

    @property
    def parity(self):
        return mat24.bw24(self.value) & 1
        
    @property
    def cocode(self):
        """Return number of the cocode word corresponding to the vector.
       
        """
        return mat24.vect_to_cocode(self.value)

    @property
    def octad(self):
        """If the bit vector is an octad, return the number of that octad
        
        Complements of octads are also accepted as octads.

        :raise: 
            * ValueError if the bit vector is not an octad.
        """
        return mat24.vect_to_octad(self.value, 0)

    @property
    def gcode(self):
        """Return number of Golay code word corresponding to the bit vector 
        
        :raise: 
             * ValueError if the bit vector is not an Golay code word.
        """
        return mat24.vect_to_gcode(self.value)

    @property
    def ord(self):
        """Return the number of the bit vector.

        We have ``0 <= i < 0x1000000`` for the returned number ``i``.
        """
        return self.value & 0xffffff

        
    @property    
    def vector(self):
        """Return the number of the bit vector.

        We have ``0 <= i < 0x1000000`` for the returned number ``i``.

        Same as property ``ord``.
        """
        return  self.value & 0xffffff
        
    @property    
    def bits(self):
        """Return the Golay code word as a ``24``-bit vector.

        Returns a list with  ``24`` entries equal to
        ``0`` or ``1``.
        """
        try:
            return self.bit_vector
        except:
            self.bit_vector = ([(self.value >> i) & 1 for i in range(24)]) 
            return self.bit_vector

    @property    
    def bit_list(self):
        """Return the ordered list of positions of bits being set"""
        try:
            return self.bit_list_
        except AttributeError:
            self.bit_list_ = [i for i in range(24) if self.value & (1<<i)]
            return self.bit_list_

    def syndrome(self, i=None):
        """Return Golay code syndrome of the vector as a bit vector.

        :param i: ``None`` (default) or an integer ``0 <= i < 24`` 
                  used to select a syndrome.
                 
        :return: Syndrome of a Golay cocode element as a bit
                 vector of type |GcVector|.

        :raise:
           *  ValueError if argument ``i`` is ``None`` and the
              syndrome is not unique.

        ``v.syndrome(i)`` is equivalent to ``Cocode(v).syndrome(i)``.
        """
        if i is None: i = 24
        return GcVector(mat24.syndrome(self.value, i))

    def all_syndromes(self):
        """Return list of all Golay code syndromes of the vector.

        :return: List of all Golay code syndromes of the bit vector
                 as a list of vectors of type |GcVector|. Such a list
                 contains either six vectors of bit weight 4, or one
                 vector of bit weight less than 4.
        """
        return [GcVector(x) for x in mat24.all_syndromes(self.value)]

    def syndrome_list(self, i=None):
        """Return syndrome of cocode element as list of bit positions.

        :param i: ``None`` (default) or an integer ``0 <= i < 24`` 
                  used to select a syndrome.
                 
        :return: Syndrome of a Golay cocode element as a list
                 of at most four bit positions.

        :raise:
           *  ValueError if argument ``i`` is ``None`` and the
              syndrome is not unique.

        The syndrome of the Golay cocode element is calculated
        in the same way as in method **syndrome**. But here
        the result is returned as an ordered list of bit positions
        corresponding to the bits which are set in the syndrome.

        ``v.syndrome_list(i)`` is equivalent to 
        ``Cocode(v).syndrome_list(i)``.
        """
        if i is None: i = 24
        syn = mat24.syndrome(self.value, i)
        return [j for j in range(24) if (1 << j) & syn]


    def vtype(self, as_int = False):
        """Compute type of a bit vector ``v`` in :math:`\mathbb{F}_2^{24}`

        Here the type of a vector ``v`` is its orbit under the action
        of the Mathieu group :math:`M_{24}`. These orbits are denoted
        as in Figure 10.1 in :cite:`CS99`, Ch. 10.2.6. An orbit of
        weight <= 12 is

        * Special [S] if it contains or is contained in a octad (i.e
          in a Golay code word of weight 8),

        * Umbral [U] if it is not special and contains or is contained
          in a dodecad (i.e. in a Golay code word of weight 12),

        * Transversal [T] otherwise.

        A vector of weight > 12 is special, umbral, or transversal,
        if its complement has that property.

        A vector of weight 12  is extraspecial [S12+] if it contains
        three octads, and penumbral [U12-] if it has Hamming
        distance two from a dodecad.

        The type is returned as a string as in :cite:`CS99`. E.g. 'U7'
        means an umbral vector of weight 7, and 'S19' means a special
        vector of weight 19. 'T12', 'S12+', and 'U12-' means a
        transversal, extraspecial, and penumbral vector of weight 12,
        respectively.

        If the argument ``as_int`` of this funtion is ``True`` then
        the type is returned as an integer as in the C function
        ``mat24_vect_type`` in file ``mat24_functions.c``.
        """
        vt = mat24.vect_type(self.value)
        if as_int:
            return vt
        return self.STR_VT[vt >> 5] % (vt & 31)

    def str(self):
        return "<GcVector_%s>" % ihex(self.value & 0xffffff, 6)
    __repr__ = str

        