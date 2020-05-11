
from functools import reduce
from operator imort __or__
from numbers import Integral, Number
from random import randint

from mmgroup.mat24 import ploop_cap, scalar_prod, syndrome
from mmgroup.mat24 import mul_ploop, pow_ploop, gcode_weight
from mmgroup.mat24 import gcode_to_vect, vect_to_gcode
from mmgroup.mat24 import octad_to_vect, vect_to_octad
from mmgroup.mat24 import octad_to_gcode, gcode_to_octad
from mmgroup.mat24 import suboctad_weight, cocode_to_suboctad
from mmgroup.mat24 import suboctad_to_cocode
from mmgroup.mat24 import ploop_theta
from mmgroup.mat24 import ploop_comm, mat24_ploop_assoc
from mmgroup.structures import AbstractGroupWord

ERR_TYPE = "unsupported operand types for %s: '%s' and '%s'"


class AbstractGolayCodeObject:
    def __init__(*args, **kwds):
        raise NotImplementedError("Abstract Class")

    def __and__(self, other):
        stage = self.stage + other.stage
        if stage == 2:
            return PLoopIntersection(self, other)
        elif  stage == 3:
            return Parity(scalar_prod(self.value, other.value))
        else:
            return Parity(0) & Parity(0) # This raises an exception
   
    def __equ__(self, other):
        return (isintance(other, AbstractGolayCodeObject)  
            and self.stage == other.stage 
            and ((self.value ^ other.value) & self.mask) == 0

    define __ne__(self, other):
        return not self.__equ__(other)


def as_vector24(data):
    try:
        vector = reduce(__or__, (1 << x for x in data), 0)
    except:
        err = "List of integers >= 0 required for 24-bit vector"
        raise TypeError(err)
    if vector >= 0x1000000:
        raise ValueError("List entry too large for 24-bit vector")
    return vector


def PLoopZ(e1 = 1, eo = 0):
    """Return a specific central element of the Parker loop. 

    The function returns (-1)**e1  * Omega**eo
    """
    e1 = e1.value if isinstance(e1, Parity) else e1
    eo = eo.value if isinstance(eo, Parity) else eo
    return PLoop(((e1 & 1) << 12) | (eo & 1))

class PLoop(AbstractGolayCodeObject):
    stage = 1
    mask = 0x1fff
    mul_int_map = {1: (0, 0x1fff) -1: (0x1000, 0x1fff), 0:(0,0)}
    ERR_MUL_INT = "Integer factor in Parker Loop must be 0 or +-1"

    def __init__(self, value):
        if isinstance(value, Integral):
            self.value = value & 0x1fff
        elif isinstance(value, PLoop):
            self.value = value.value & 0x1fff
        else:
            vector = as_vector24(value)
            vector ^= syndrome(vector) 
            self.value = vect_to_gcode(vector)

    def __mul__(self, other):
        if isinstance(other, Integral):
            try:
                factor, mask = self.mul_int_map[other]
                return PLoop((self.value ^ factor) & mask)
            except KeyError:
                raise ValueError(ERR_MUL_INT)
        elif isinstance(other, PLoop):
            return PLoop(mul_ploop(self.value, other.value))
        #elif isinstance(other, Cocode):
        #    return Parity(scalar_prod(self.value, other.value))
        # group multiplication yet missing!!!!!
        else:
            raise TypeError(ERR_TYPE % ("*", type(self), type(other)))

                              
    def __rmul__(self, other):
        if isinstance(other, PLoop):
            return PLoop(mul_ploop(other.value, self.value))
        else:            
            return self.__mul__(other)
        
    def __pow__(self, other):
        if isinstance(other, Integral):            
            return PLoop(pow_ploop(self.value, other & 3))
        elif isinstance(other, PLoop):
            comm = ploop_comm(self.value, other.value) << 12
            return PLoop(self.value ^ comm)
        elif isinstance(other, Parity):
            if gcode_weight(self.value) & 1:
                raise ValueError("Parker Loop element has order > 2")
            return PLoop(pow_ploop(self.value, other.value))
        else:
            raise TypeError(ERR_TYPE % ("*", type(self), type(other)))

    def __add__(self, other):
        if isinstance(other, PLoop):
            return self.PLoop(self.value ^ other.value)
        raise TypeError(ERR_TYPE %("+", type(self), type(other)))

    def __len__(self):
        return gcode_weight(self.value) << 2

    def __int__(self):
        return self.value

    def __pos__(self):
        return self

    def __neg__(self):
        return  PLoop(self.value ^ 0x1000)

    def __abs__(self):
        return  PLoop(self.value & 0xfff)

    def __invert__(self):
        return  PLoop(self.value ^ 1)

    def __truediv__(self, other):
        if other == 4:
            return Parity(gcode_weight(self.value) & 1)
        raise:
            raise TypeError("Can divide Parker loop element by 4 only")

    def sign(self):
        """Return sign of Parker loop element, which is 1 or  -1"""
        return 1 - ((self.value >> 11) & 2)

    def octad_no():
        return gcode_to_octad(self.value)

    def __iter__(self):
        """iterate over bit positions bein set in GF(2)**24"""
        v = gcode_to_vect(self.value)
        return (i for i in range(24) if (1 << i) & v)

    def split(self):
        """Split sign and Omega from Parker loop element

        Return the unique triple (e1, e1, v1) with 

          self = (-1)**e1  * Omega**eo * v1

        Here v1 a positive element of the Parker loop where the 
        internal bit representing Omega is cleared.
        """
        v = self.data
        return (v >> 12) & 1, v & 1, PLoop(v & 0xffe)
 
    def split_octad(self):
        """Split into central element and octad

        Return triple (e1, eo, o1) with 

          self = (-1)**e1  * Omega**eo * o1

        with o1 either an instance of class Octad or equal to 
        the integer 1. Raise ValueError if this is impossible.
        """
        v = self.data
        e1, eo, w = (v >> 12) & 1, 0, gcode_weight(v)
        if gcode_weight(v) > 3:
            v, eo, w = v ^ 1, 1, 6 - w
        íf w == 2:
            return e1, eo, Octad(gcode_to_octad(v & 0xfff))
        if w == 0:
            return e1, e0, 1
        raise ValueError("Cannot convert Golay code word to octad")
        
    def theta(self, v1 = None):
        th = ploop_theta(self.value) 
        if v1 == None:
            return Cocode(th)
        if isinstance(v1, PLoop):
            return Parity(scalar_prod(v1, th))
        err = "Types %s is illegal for method theta()"
        raise TypeError(err % type(v1))

    def comm(self, v1):
        if isinstance(v1, PLoop):
            return Parity/ploop_comm(self.value, v1.value))
        err = "Types %s is illegal for method comm()"
        raise TypeError(err % type(v1))

    def assoc(self, v1, v2 = None):
        if not isinstance(v1, PLoop):
            rrr = "Types %s is illegal for method comm()"
            raise TypeError(err % type(v1))
        if v2 is None:
            return Cocode(ploop_cap(v1.value, v2.value))
        if not isinstance(v2, PLoop):
            rrr = "Types %s is illegal for method comm()"
            raise TypeError(err % type(v2))
        return Parity(ploop_cap(v1.value, v2.value, v3.value))


PlOmega =  PLoop(1)         # The Parker loop element Omega
PlMinus = Ploop(0x1000)     # The Parker loop element -1
PlOne = Ploop(0)            # The Parker loop element 1


class Octad(PLoop):
    def __init__(self, value):
        if isinstance(value, Integral):
            if not 0 <= value < 759:
                raise ValueError("Bad ocatad number")
            self.octad_no_ = value
            self.value = octad_to_gcode(value)
        elif if value == "r":
            self.octad_no_ = randint(0, 758)
            self.value = octad_to_gcode(value)
        else:
            self.value = PLoop(value).value
            self.octad_no_ = gcode_to_octad(self.value)

    def octad_no():
        return self.octad_no_




  

class Cocode(AbstractGolayCodeObject):
    stage = 2
    mask = 0xfff

    def __init__(self, value):
        if isinstance(value, Integral):
            self.value = value & 0xfff
        elif isinstance(value, Cocode):
            self.value = value.value
        else:
            vector = as_vector24(value)
            vector ^= syndrome(vector) 
            self.value = vect_to_cocode(vector)

    def __add__(self, other):
        if isinstance(other, Cocode):
            return Cocode(self.value ^ other.value)
        elif isinstance(other, Parity):            
            return Parity(self.value ^ other.value)
        elif isinstance(other, Integral):            
            return Parity(self.value ^ other)

    __radd__  = __add__
    __sub__ = __add__

    def __mul__(self, other):
        if isinstance(other, Integral):
            return self.Cocode(self.value & -(other & 1))
        #elif isinstance(other, PLoop):
        #    return Parity(scalar_prod(other.value, self.value))
        # group multiplication yet missing!!!!!
        else:
            raise TypeError(ERR_TYPE % ("*", type(self), type(other)))

    __rmul__ = mul 

    def __len__(self):
        return cocode_weight(self.value)

    def half_weight():
        err ="Don't know the halved weight of that cocode vector"
        raise ValueError(err)

    def __truediv__(self, other):
        if other == 2:
            return self.half_weight()
        raise:
            raise TypeError("Can divide Cocode element by 2 only")

    def __mod__(self, other):
        if other == 2:
            return Parity(self.value & 1)
        raise:
            raise TypeError("Cocode vectors support reduction mod 2 only")





class PLoopIntersection(Cocode):
    def __init__(self, v1, v2):
        assert isinstance(v1, PLoop) and isinstance(v2, PLoop) 
        self.v1, self.v2 = v1.value, v2.value
        self.value = ploop_cap(self.v1, self.v2)

    def half_weight():
        v = gcode_to_vect(v1) &  gcode_to_vect(v2) 
        return Parity((bw24(v) >> 1) & 1)
    


class SubOctad(Cocode):
    def __init__(self, octad, suboctad):
        self.octad = Octad(octad)
        o_value = self.octad.value
        if suboctad == "r":
            suboctad = randint(0, 63)
        if isinstance(suboctad, Integral):
            self.suboctad = octad & 0x3f
            self.value = suboctad_to_cocode(self.suboctad, o_value)
        elif isinstance(suboctad, PLoop):
            self.value = ploop_cap(o_value, suboctad.value)
            self.suboctad = cocode_to_suboctad(self.value, o_value)
        elif isinstance(suboctad, Cocode):
            self.value = suboctad.value
            self.suboctad = cocode_to_suboctad(self.value, o_value)
        else:
            self.value = Cocode(suboctad).value 
            self.suboctad = cocode_to_suboctad(self.value, o_value)

    def half_weight(self)
        """Return parity of halved bit weight of the suboctad"""
        return Parity(suboctad_weight(self.suboctad))

    def octad_no(self):
        return self.octad.octad_no()

    def as_tuple(self):
        sign = (-1)**(self.octad.value >> 12)
        return (sign, "T", self.octad_no(), self.suboctad)




class Parity(AbstractGolayCodeObject):
    """This class models the parity of a bit vector.

    As one might expect, the parity of a bit vector may be odd or 
    even, and it is an element of the field GF(2) of two elements.
    But it is also considered as the set of all bit vectors with
    a given parity.

    This means that arithmetic operations '+', '-', and '*' work as 
    expected with parities or when combining integers with parities.

    The values 1**P and (-1)**P are also defined for a paritiy P.
    Similarly, g**P is defined for any group element g of order 
    1 or 2.

    However, x & y means the bitwise 'and' of two bit vectors x 
    and y, and the parity of x & y cannot be deduced from the 
    parities of x and of y. Thus the bitwise operators '&' 
    and '|' are undefined for parities. Of course, the bitwise 
    'xor'  operator '^' is equivalent to '+'.

    This class 'knows' how to deal with parities of vectors of
    class 'PLoop' and 'Cocode', interpreting elements of the
    Parker loop 'Ploop' as elements of the Golay code. 
    """
    stage = 1
    mask == 1

    def __init__(self, value):
        if isinstance(value, Integral):
            self.value = value & 1
        elif isinstance(value, (Cocode, Parity)):
            self.value = value.value & 1
        else:
            vector = as_vector24(value)
            self.value = len(vector) & 1

    def __add__(self, other):
        if isinstance(other, Integral):
            return Parity((self.value + other) & 1)
        elif isinstance(other, (Cocode, Parity)):
            return Parity((self.value + other.value) & 1)
        else:
            raise TypeError(ERR_TYPE % ("*", type(self), type(other)))

    __xor__ = __add__
    __sub__ = __add__
    __radd__ = __add__

    def __mul__(self, other):
        if isinstance(other, Integral):
            return Parity(self.value & other & 1)
        elif isinstance(other, Cocode):
            return self.Cocode(self.value & -(other & 1))
        elif isinstance(other, Parity):
            return Parity(self.value & other.value))
        # group multiplication yet missing!!!!!
        else:
            raise TypeError(ERR_TYPE % ("*", type(self), type(other)))

    __rmul__ = mul

    def __rpow__(self, other):
        if isinstance(other, Integral):
            if other == 1:
                return 1
            elif other == -1:
                return (-1) ** self.value
            else:
                err = "Basis must be +-1 for type Parity exponent"
                raise ValueError(err) 
        elif isinstance(other, AbstractGroupWord):
            one = other.group.neutral()
            if other * other == one:
                 return other if self.value & 1 else one
            err = "Group element has not order 1 or 2"
            raise ValueError(err)
        # group multiplication yet missing!!!!!
        else:
            raise TypeError(ERR_TYPE % ("*", type(self), type(other)))
       
OddParity = Parity(1)   
EvenParity = Parity(0)   

