r"""We deal with octads and certain subsets of octads called suboctads.

An *octad* is a  Golay code word of weight ``8``. The Golay code 
has ``759`` octads. A signed octad is a Parker loop element 
corresponding to an octad or a complement of an octad. Signed 
octads are used for indexing some unit vectors of the 
representation :math:`\rho` of the monster :math:`\mathbb{M}`.

Function ``Octad()`` returns a signed octad as an instance of 
class |PLoop|. Arguments of function ``Octad()`` are interpreted
in the same way as arguments of the constructor for class |PLoop|,
but the function raises `ValueError` if the result is not a
(possibly complemented) signed octad.

We use some lexical order for numbering the ``759`` octads,
which we do not describe in detail.  
Due to the well-known properties of the Golay code we can create 
a random octad and display its number as follows:

.. code-block:: python
    
        from random import sample
        from mmgroup import Octad
        # Select 5 random integers between 0 and 23
        int_list = sample(range(24), 5)
        # Complete these 5 integers to a (unique) octad o
        o = Octad(int_list)
        # Display octad o and its number
        print("Octad", o.bit_list, "has number", o.octad)


A *suboctad* of an octad ``o`` is an element of the Golay cocode
:math:`\mathcal{C}^*` of even parity which can be represented
as a subset of ``o``. Each octad has ``64`` suboctads.
A pair (octad, suboctad) is used for indexing some basis
vectors in the representation :math:`\rho`. Class |SubOctad|
models such a pair. The constructor of  class |SubOctad|
takes two arguments ``octad`` and ``suboctad``. The first argument
evaluates to a  signed octad as in function ``Octad()``. The
second argument evaluates to a suboctad, see class |SubOctad|
for details.

The raison d'etre of a  suboctad is indexing a basis vector in
the representation  :math:`\rho`. For this purpose we need a pair 
of integers refering to the octad and the suboctad. For an instance 
``so`` of class |SubOctad| that pair is given by
``(so.octad, so.suboctad)``.

For an octad ``o`` and a suboctad ``s`` precisely one of the
pairs ``(o,s)`` and ``(~o,s)`` is actually used as an index for 
:math:`\rho`. The pair ``(~o,s)`` is used if the cocode element 
corresponding to ``s`` has minimum weight ``2``; the pair 
``(o,s)`` is used if that element has minimum weight ``0`` or 
``4``. See :cite:`Con85` or :cite:`Seysen20` for background. In 
this  implementation both, ``(o,s)`` and ``(~o,s)``, refer to the 
same basis vector of :math:`\rho`.
"""



from functools import reduce
from operator import __xor__
from numbers import Integral, Number
from random import randint


try:
    # Try importing the fast C function
    from mmgroup import mat24 
except (ImportError, ModuleNotFoundError):
    # Use the slow python function if the C function is not available
    from mmgroup.dev.mat24.mat24_ref import  Mat24
    mat24 = Mat24



from mmgroup.structures.auto_group import AbstractGroupWord
from mmgroup.structures.parity import Parity
from mmgroup.structures.parse_atoms import ihex

from mmgroup.structures.autpl import AutPL, AutPlGroup

from mmgroup.structures.gcode import GCode, GcVector
from mmgroup.structures.cocode import Cocode
from mmgroup.structures.ploop import PLoop


ERR_RAND = "Illegal string for constricting type %s element" 

ERR_DIV4 = "%s object may be divided by 2 or 4 only"
ERR_DIV2 = "%s object may be divided by 2 only"



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
# Function Octad
#######################################################################







def Octad(octad):
    """Return a (signed) octad as an element of the Parker loop. 

    :param octad:  the value of the octad to be returned. 
    :type octad: see table below for legal types

    :return:  a (signed) octad
    :rtype:   an instance of class |PLoop|

    :raise:
        * TypeError if ``type(octad)`` is not in the table given above.
        * ValueError if ``octad`` cannot be converted to a
          (possibly negated and complemented) octad.


    Depending on its type parameter **octad** is  interpreted as follows:

    .. table:: Legal types for parameter ``octad``
      :widths: 20 80

      ===================== ================================================
      type                  Evaluates to
      ===================== ================================================
      ``int``               Here the (positive) octad with  number 
                            ``octad`` is returned. There are 759 octads
                            in the Golay code. 
                            So ``0 <= octad < 759`` must hold.
                            

      ``list`` of ``int``   Such a list is converted to a Golay code word,
                            see class |GCode|, and the corresponding 
                            (positive) Parker loop element is returned.

      class |GCode|         The corresponding 
                            (positive) Parker loop element is returned. 

      class |PLoop|         A deep copy of parameter ``octad`` is returned.

      class |GcVector|      This is converted to a Golay code word,
                            see class |GCode|, and the corresponding 
                            (positive) Parker loop element is returned.

      class |SubOctad|      The *octad* part of the |SubOctad| ``octad``  
                            is  returned. 

      ``str``               Create random element depending on the string
                             | ``'r'``: Create arbitrary octad

      ===================== ================================================

    A complement of an octad is also accepted; then the corresponding 
    Parker loop element is retured. The function raises ValueError
    if  parameter ``octad`` does not evaluate to an octad or a
    complement of an octad.      
    """
    if isinstance(octad, Integral):
        if  not 0 <= octad < 759:
            raise ValueError("Bad octad number")
        return PLoop(mat24.octad_to_gcode(octad) & 0xfff)
    elif isinstance(octad, str):
        return PLoop(mat24.octad_to_gcode(randint(0, 758)))
    elif isinstance(octad, PLoop):
        _ = octad.octad
        return octad
    elif isinstance(octad, GCode):
        _ = octad.octad
        return PLoop(octad.ord)
    elif isinstance(octad, SubOctad):
        return  octad.octad_
    else:
        v =  PLoop(octad)
        _ = v.octad
        return v





#######################################################################
# Class SubOctad
#######################################################################


class SubOctad():
    """Models a pair (octad, suboctad)

    :param octad:

      The first component of the pair *(octad, suboctad)* to be 
      created.

    :type octad: same as in function
      :py:class:`~mmgroup.Octad`

    :param suboctad:
    
      The second component of the pair 
      *(octad, suboctad)* to be created.

    :type suboctad: see table below for legal types


    :return: 

     the pair *(octad, suboctad)* 

    :rtype: an instance of class |SubOctad|

    :raise:
        * TypeError if  argument *octad* or *suboctad* is not 
          of correct type.
        * ValueError  argument *octad* or *suboctad* does not
          evaluate to an octad or to a correct suboctad,
          respectively.

    A *suboctad* is an even cocode element that can be represented
    as a subset of the octad given by the argument *octad*. 

    The raison d'etre of class |SubOctad| is that pairs
    *(octad, suboctad)* are used for indexing vectors in the
    representation of the monster group. Here we want to 
    number the octads from ``0`` to ``758`` and the suboctads
    form ``0`` to ``63``, depending on the octad. Note that
    every octad has ``64``  suboctads.

    Depending on its type parameter **suboctad** is  interpreted as follows:

    .. table:: Legal types for parameter ``suboctad``
      :widths: 20 80

      ===================== ================================================
      type                  Evaluates to
      ===================== ================================================
      ``int``               Here the suboctad with the number given 
                            in the argument is 
                            taken.  That numbering depends on the octad 
                            given in   the argument ``octad``. 
                            ``0 <= suboctad < 64`` must hold.                           

      ``list`` of ``int``   Such a list is converted to a bit vector
                            as in class |GcVector|,
                            and the cocode element corresponding to that
                            bit vector is taken.

       class |GCode|        The intersection of the octad given as the 
                            first argument and the Golay code word given
                            as the second argument is taken. 
  
       class |GcVector|     This is converted to a cocode element,
                            see class |Cocode|, and that cocode element 
                            is taken.

       class |Cocode|       That cocode element is taken as the suboctad.

       class |SubOctad|     The ``suboctad`` part of the |SubOctad| ``value``  
                            is  taken. 

      ``str``               Create random element depending on the string
                             | ``'r'``: Create arbitrary suboctad
      ===================== ================================================

    **Remark**

    For an object ``so`` of type |Suboctad| use ``Octad(so)`` to obtain its
    octad as a Parker loop element and ``Cocode(so)`` obtain its suboctad 
    as a cocode element. You may obtain more information from the properties
    of the object  ``so``.
    """
    __slots__ = "octad_", "o_value", "sign_", "suboctad_"
    parity = 0
    ERR_MUL = "SubOctad can be multiplied with 1 or -1 only"

    def __init__(self, octad, suboctad):
        self.octad_ = Octad(octad)
        self.o_value = self.octad_.octad
        self.sign_ = (-1)**(self.octad_.value >> 12)
        gcode =  self.octad_.gcode
        if isinstance(suboctad, str):
            self.suboctad_ = randint(0, 63)
        elif isinstance(suboctad, Integral):
            self.suboctad_ = suboctad & 0x3f
        elif isinstance(suboctad, GCode):
            value = mat24.ploop_cap(gcode, suboctad.value)
            self.suboctad_ = mat24.cocode_to_suboctad(value, gcode)
        else:
            value = Cocode(suboctad).cocode 
            self.suboctad_ = mat24.cocode_to_suboctad(value, gcode)

    def __mul__(self, other):
        if isinstance(other, Integral):
            if abs(other) == 1:
                return SubOctad(value * self.octad_ , 
                      self.suboctad_)
            else:
                raise ValueError(self.ERR_MUL) 
        elif isinstance(other, AutPL):
            octad = self.octad_ * other
            cocode = Cocode(self.cocode) * other
            return SubOctad(octad, cocode)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Integral):
            return self.__mul__(other)
        else:
            return NotImplemented

    def __truediv__(self, other):
        if  isinstance(other, Integral):
            if abs(other) == 1:
                self.__mul__(other)
            elif abs(other) == 2:
                return Parity(mat24.suboctad_weight(self.suboctad_))
            else:
                raise TypeError(ERR_DIV2 % type(self))
        else:           
            return NotImplemented

    def __abs__(self):
         return SubOctad(abs(self.octad_ ),  self.suboctad_)

    def __pos__(self):
        return self

    def __neg__(self):
        return SubOctad(-self.octad_ ,  self.suboctad_)

    def __invert__(self):
        return SubOctad(~self.octad_ ,  self.suboctad_)


    @property
    def sign(self):
        """Return the sign of the octad.

        This is ``1`` for a positive and ``-1`` for a negative octad.
        """
        return self.sign_

    @property
    def octad(self):
        """Return number of octad.
        
        We have  ``0 <= o < 759`` for the returned number ``o``.
        """
        return self.o_value

    @property
    def gcode(self):
        """Return number of Golay code word corresponding to the octad.

        We have ``0 <= g < 0x1000`` for the returned number ``g``.
        """
        return self.octad_.gcode
        
    @property
    def suboctad(self):
        """Return the number of the suboctad.

        We have ``0 <= s < 64`` for the returned number ``s``.

        Let ``[b0, b1,..., b7``] be the bits set in the octad of
        an instance of this class in natural order. The following 
        table shows the suboctad numbers for some suboctads given
        as cocode elements. More suboctad numbers can be obtained 
        by combining suboctads and their corresponding numbers with
        ``XOR``.

        .. table:: Suboctad numbers of some cocode elements
         :widths: 16 16 16 16 18 18

         =========== =========== =========== =========== =========== =========== 
         ``[b0,b1]`` ``[b0,b2]`` ``[b0,b3]`` ``[b0,b4]`` ``[b0,b5]`` ``[b0,b6]``
         ``s  = 1``  ``s  = 2``  ``s  = 4``  ``s  = 8``  ``s = 16``  ``s = 32``
         =========== =========== =========== =========== =========== =========== 

        E.g. ``[b0, b5, b6, b7]`` is equivalent to ``[b1, b2, b3, b4]`` 
        modulo the Golay code and has number ``s = 1 ^ 2 ^ 4 ^ 8 = 15``.
        """
        return self.suboctad_
    

    @property
    def cocode(self):
        """Return number of cocode word corresponding to the suboctad.

        We have ``0 <= c < 0x1000`` for the returned number ``c``.
        """
        return mat24.suboctad_to_cocode(
                 self.suboctad_, self.octad_.value)

    def as_tuple(self):
        """Return data as a tuple 

        The function returns a tuple
        ``(sign, 'T', octad, suboctad)``
        where ``sign``, ``octad``, and ``suboctad`` are the
        integers given by the corresponding properties in
        class |SubOctad|.

        Such a tuple is used for  indexing vectors in the
        representation of the monster group. 
        """
        return (self.sign_, "T", self.octad_no(), self.suboctad_)

    def __eq__(self, other):
        return (isinstance(other, SubOctad) and 
               self.sign_ ==  other.sign_ and
               self.o_value ==  other.o_value and
               self.suboctad_ ==  other.suboctad_)

    def __ne__(self, other):
        return not self.__equ__(other)

    def str(self):
        return "<SubOctad_%d_%s>" % (o_value, ihex(self.suboctad_, 2))
    __repr__ = str    

