
from numbers import Integral, Number
from random import randint


#######################################################################
# Class Parity 
#######################################################################


class Parity():
    """This class models the parity of a bit vector.

    As one might expect, the parity of a bit vector may be odd or 
    even, and it is an element of the field ``GF(2)`` of two elements.

    This means that arithmetic operations ``+``, ``-``, and ``*`` 
    work as expected on instances of class |Parity|, or on 
    instances of this class combined with integers. 

    The values ``1**p`` and ``(-1)**p`` are defined as usual for an
    object ``p`` of class |Parity|.
    Similarly, ``g**p`` is defined  if ``g`` is an element of a group 
    defined in this package and has order ``1`` or ``2``. More precisely, 
    ``g`` must be an instance of  class ``AbstractGroupWord`` here.

    Bitwise operations ``&``, ``|``, and ``^`` are illegal on  
    instances of class  |Parity|.

    ``Parity(x)`` is defined if ``x`` an integer, an instance of class
    |Parity|, or an instance of any class defined is this package which 
    has a natural  parity, e.g. an instance of class |GcVector| 
    or |Cocode|.

    ``int(p)`` evaluates to ``0`` or ``1`` for an instance of 
    class |Parity|.


    Implementation remarks
    
      If an object ``x`` has an attribute or a property ``parity`` then
      ``Parity(x)`` means ``Parity(x.parity)`` and ``x + p`` means 
      ``Parity(x) + p`` for any  instance ``p`` of class Parity.

      If an object ``x`` has an attribute or a property ``parity_class`` 
      then  ``x * p``  and  ``p * x``  mean  
      ``Parity(parity_class(x) * int(p))``. If that attribute or a 
      property has value ``True`` then ``x * p``  and  ``p * x``  mean  
      ``Parity(x * int(P))``.
    """
    __slots__ = "value"
    names = ["<even>", "<odd>"]
    ERR_MOD2 = "%s object may be reduce mod 2 only"
    ERR_PARITY = "%s object does not have a parity"

    def __init__(self, value):
        if isinstance(value, Integral):
            self.value = value & 1
        elif isinstance(value, Parity):
            self.value = value.value & 1
        elif isinstance(value, str):
            self.value = randint(0, 1)
        else:
            try:
                self.value = int(value.parity) & 1
            except:
                raise TypeError(self.ERR_PARITY % type(value)) 

    def __int__(self):
        return self.value & 1

    @property
    def ord(self):
        """Return ``1`` if the parity is odd and ``0`` if it is even"""
        return self.value

    def __eq__(self, other):
        return (isinstance(other, Parity)  
            and (self.value ^ other.value) & 1 == 0)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __pos__(self):
        return self
        
    __neg__ = __pos__


    def __add__(self, other):
        if isinstance(other, Integral):
            return Parity((self.value + other) & 1)
        elif isinstance(other, Parity):
            return Parity((self.value + other.value) & 1)
        else:
            try:
                parity = other.parity
            except:
                return NotImplemented
            else:
                return Parity((self.value + parity) & 1)

    __sub__ = __add__
    __radd__ = __add__
    __rsub__ = __radd__

    def __mul__(self, other):
        if isinstance(other, Integral):
            return Parity(self.value & other & 1)
        elif isinstance(other, Parity):
            return Parity(self.value & other.value & 1)
        else:
            try:
                parity_class = other.parity_class
            except AttributeError:
                return NotImplemented
            else:
                if parity_class is True:
                   return other * self.value
                else:
                   return parity_class(other) * self.value

    __rmul__ = __mul__


    def __truediv__(self, other):
        if isinstance(other, Integral):
            if abs(other) == 1:
                return self
        else:
            return NotImplemented

    def __mod__(self, other):
        if isinstance(other, Integral):
            if other == 2:
                return self
            else:
                raise ValueError(ERR_MOD2 % type(self))
        else:
            return NotImplemented


    def __rpow__(self, other):
        if isinstance(other, Integral):
            if other == 1:
                return 1
            elif other == -1:
                return (-1) ** self.value
            else:
                err = "Basis must be +-1 for type Parity exponent"
                raise ValueError(err)
        else:
            return NotImplemented
 
    def __str__(self):
        return self.names[self.value & 1]

    __repr__ = __str__


