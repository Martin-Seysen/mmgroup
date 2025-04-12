"""We deal with the cocode of the Golay code

Documentation see module mmgroup.structure.gcode
"""
from functools import reduce
from operator import __xor__
from numbers import Integral, Number
from random import randint


from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.parity import Parity
from mmgroup.structures.parse_atoms import ihex

from mmgroup.structures.gcode import as_vector24, str_vector24
from mmgroup.structures.gcode import GCode, GcVector



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
    global import_pending, mat24, COCODE_BASIS, PLoop
    global AutPL, AutPlGroup, XLeech2
    from mmgroup import mat24, COCODE_BASIS 
    from mmgroup.structures.ploop import PLoop
    from mmgroup.structures.autpl import AutPL, AutPlGroup
    from mmgroup.structures.xleech2 import XLeech2
    import_pending = False




#######################################################################
# Class Cocode 
#######################################################################




class Cocode():
    """This class models an element of the cocode of the Golay code.

    The Golay code is a binary ``[24, 12, 8]`` code. So its cocode has 
    :math:`2^{12}` elements, each of it is a :math:`24` bit long word
    (modulo the Golay code).

    The :math:`2^{12}` Golay cocode elements are numbered from
    ``0`` to ``0xfff``. Bit ``i`` in that number refers to the
    coefficient of the  ``i``-th basis vector of the Golay cocode,
    which may be :math:`0` or :math:`1`. The chosen basis of the
    cocode is the reciprocal basis of the chosen basis of the
    Golay code.

    :param value:

      This is the value of the Golay cocode word to be returned. 

    :type value: see table below for legal types

    :return: A Golay cocode element
    :rtype:  an instance of class |Cocode|

    :raise:
        * TypeError if ``type(value)`` is not in the table given below.


    Depending on its type parameter **value** is  interpreted as follows:

    .. table:: Legal types for constructor of class ``Cocode``
      :widths: 20 80

      ==================== ================================================
      type                 Evaluates to
      ==================== ================================================
      ``int``              Here the cocode element with number ``value`` 
                           is returned. ``0 <= value <= 0x1000`` must hold.

      ``list`` of ``int``  A list of integers  ``0 <= i < 24`` is 
                           interpreted as a list of bit positions to be
                           set in a bit vector. The cocode word
                           corresponding to that bit vector is returned.

      class |Cocode|       A deep copy of parameter ``value`` is 
                           returned. 

      class |XLeech2|      The *cocode* part of the |XLeech2| object  
                           is converted to a cocode element and 
                           returned. 

      ``str``              Create random element depending on the string
                              | ``'r'``: Create arbitrary cocode element
                              | ``'e'``: Create an even cocode element
                              | ``'o'``: Create an odd cocode element

      ==================== ================================================


    **Standard operations**
    
    Addition of cocode elements and scalar multiplication (with scalars 
    of  type  ``int`` or |Parity|) are done as usual.

    ``g & c`` and ``c & g`` are legal operations for a Golay code
    element ``g`` of type  |GCode| and a cocode element ``c`` of type  
    |Cocode|. The result of this operation is not well defined as 
    a bit vector, but it  has a well-defined parity. Thus
    ``g & c`` returns a parity object of type |Parity|.

    **Standard functions**

    Let ``c`` be a cocode element of type |Cocode|.
    ``len(c)`` returns the bit weight of a shortest representative
    of the Golay code element ``c``. 

    **Special functions**
  
    The expression ``c % 2`` is a shorthand for  ``Parity(len(c))``.
    It returns the parity of ``c`` as a parity
    object of type  |Parity|. 

    """
    __slots__ = "value"
    parity_class = True

    def __init__(self, value):
        if import_pending:
            complete_import()
        if isinstance(value, Integral):
            self.value = value & 0xfff
        elif isinstance(value, Cocode):
            self.value = value.value
        elif isinstance(value, GcVector):
            self.value = mat24.vect_to_cocode(value.value)
        elif isinstance(value, XLeech2):
            v = value.value
            self.value = (mat24.ploop_theta(v >> 12) ^ v) & 0xfff
        elif isinstance(value, str):
            if len(value) == 1 and value in "reo":
                self.value = randint(0, 0xfff)
                if "o" in value: self.value |= 0x800
                if "e" in value: self.value &= 0x7ff
            else:
                raise ValueError(ERR_RAND % 'Cocode')
        else:
            vector = as_vector24(value)
            self.value = mat24.vect_to_cocode(vector)


    def __eq__(self, other):
        return (isinstance(other, Cocode)  
            and (self.value ^ other.value) & 0xfff == 0)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __pos__(self):
        return self
        
    __neg__ = __pos__

    def __and__(self, other):
        if isinstance(other, GCode):
            return Parity(mat24.scalar_prod(self.value, other.value))
        else:
            return NotImplemented


    def __add__(self, other):
        if isinstance(other, Cocode):
            return Cocode((self.value ^ other.value) & 0xfff)
        elif isinstance(other, GcVector):
            return Cocode(self.value ^ mat24.vect_to_cocode(other.value))
        elif other == 0:
            return Cocode(self)
        else:
            return NotImplemented


    __sub__ = __add__
    __radd__ = __add__
    __rsub__ = __radd__

    def __mul__(self, other):
        if isinstance(other, Integral):
            mask = -(other & 1)
            return Cocode(self.value & mask & 0xfff)
        elif isinstance(other, AutPL):
            return Cocode(mat24.op_cocode_perm(self.value, other.perm))
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Integral):
            mask = -(other & 1)
            return Cocode(self.value & mask & 0xfff)
        else:           
            return NotImplemented

    def __len__(self):
        return mat24.cocode_weight(self.value)

    def half_weight(self):
        err ="Don't know the halved weight of that cocode vector"
        raise ValueError(err)

    def __truediv__(self, other):
        if isinstance(other, Integral):
            if abs(other) == 1:
                return self
            elif abs(other) == 2:
                return self.half_weight()
            else:
                raise ValueError(ERR_DIV2 % type(self))
        else:
            return NotImplemented

    def __mod__(self, other):
        return Parity(self.parity) % other

    @property
    def ord(self):
        """Return the number of the cocode element.

        We have ``0 <= i < 0x1000`` for the returned number ``i``.
        """
        return self.value & 0xfff
        
    @property
    def cocode(self):
        """Return the number of the cocode element.

        We have ``0 <= i < 0x1000`` for the returned number ``i``.

        Same as property ``ord``.
        """
        return self.value & 0xfff

    @property    
    def parity(self):
        """Return parity of the cocode element as a |Parity| object.

        """
        return (self.value >> 11) & 1


    def syndrome(self, i=None):
        """Return syndrome of cocode element as a bit vector.

        :param i: ``None`` (default) or an integer ``0 <= i < 24`` 
                  used to select a syndrome.
                 
        :return: Syndrome of a Golay cocode element as a bit
                 vector of type |GcVector|.

        :raise:
           *  ValueError if argument ``i`` is ``None`` and the
              syndrome is not unique.

        Any Golay cocode element has either a unique shortest
        representative of bit weight ``<= 3`` or precisely six
        shortest representatives of bit weight ``4`` which form
        a partition of the underlying set of the code. Such a
        partition is called a *tetrad*. 

        In coding theory, a shortest representative of a cocode 
        element is called a *syndrome*.
             
        If the syndrome is unique, the function returns that 
        syndrome. Otherwise it returns the syndrome containing
        bit at position ``i``.
        """
        if i is None: i = 24
        return GcVector(mat24.cocode_syndrome(self.value, i))

    def all_syndromes(self):
        """Return list of all Golay code syndromes of the element.

                 
        :return: List of all Golay code syndromes of the cocode element
                 as a list of vectors of type |GcVector|. Such a list
                 contains either six vectors of bit weight 4, or one
                 vector of bit weight less than 4. 
        """
        return [GcVector(x) 
            for x in mat24.cocode_all_syndromes(self.value)]

    def syndrome_list(self, i=None):
        """Return syndrome of cocode element as list of bit positions.

        :param i: ``None`` (default) or an integer ``0 <= i < 24`` 
                  used to select a syndrome.
                 
        :return: Syndrome of a Golay cocode element as a list
                 of at most four bit positions.

        :raise:
            * ValueError if argument ``i`` is ``None`` and the
              syndrome is not unique.

        The syndrome of the Golay cocode element is calculated
        in the same way as in method **syndrome**. But here
        the result is returned as an ordered list of bit positions
        corresponding to the bits which are set in the syndrome.
        """
        if i is None: i = 24
        return mat24.cocode_to_bit_list(self.value, i)

    def syndromes_llist(self):
        """Return shortest syndromes of cocode element as list of lists.

        The function returns the list of all shortest representatives
        of the cocode element as a list ``ll`` of lists.  Each entry
        of ``ll`` corresponds to a shortest representative. Such an
        entry is an ordered list of bit positions corresponding to
        the bits which are set in the representative. Output ``ll``
        is sorted in natural order.
        """
        try:
            return [mat24.cocode_to_bit_list(self.value, 24)]
        except:
            bl = mat24.cocode_to_sextet(self.value)
            return [bl[i:i+4] for i in range(0, 24, 4)]

    def str(self):
        return "<Cocode_%s>" % ihex(self.value & 0xfff, 3)
    __repr__ = str    


    @classmethod
    def str_basis(cls):
       b = "Golay cocode basis c_i[0,...,23], i = 0,...,11"
       return str_basis(b, "c", COCODE_BASIS)

    @classmethod
    def show_basis(cls):
        print(cls.str_basis())









#######################################################################
# Subclasses of class Cocode
#######################################################################



class PLoopIntersection(Cocode):
    """Models the intersection of two Golay code words as a cocode word.
   
    Not for public use!
    """
    __slots__ = "value", "v1", "v2"
    parity = 0

    def __init__(self, v1, v2):
        if import_pending:
            complete_import()
        assert isinstance(v1, GCode) and isinstance(v2, GCode) 
        self.v1, self.v2 = v1.value, v2.value
        self.value = mat24.ploop_cap(self.v1, self.v2)

    def half_weight(self):
        return Parity(mat24.ploop_comm(self.v1, self.v2))
    
    def str(self):
        return "<Cocode_%s>" % ihex(self.value & 0xfff, 3)
    __repr__ = str    




 