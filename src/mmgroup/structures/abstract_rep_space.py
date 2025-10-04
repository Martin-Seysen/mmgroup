from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import re
import warnings
import copy
from random import randint, randrange, sample

from copy import deepcopy
from numbers import Integral

from mmgroup.structures.parse_atoms import AtomDict 
from mmgroup.structures.parse_atoms import eval_atom_expression, ihex 


######################################################################
# Utility functions for scalars
######################################################################




def mod_inv(a, p):
    """Return x with x * a == 1 (mod p)  for integers a, p.

    p >= 1 must hold.
    Raise ZerodivsionError if inversion is not possible
    """
    assert p >= 1
    a1, a2, x1, x2 = p, a % p, 0, 1
    # Invert a modulo p with the extended Euclidean algorithm
    while a2 > 0:
        # assert x1 * a % p == a1 % p and  x2 * a % p == a2 % p
        q, r = divmod(a1, a2)
        a1, a2, x1, x2 = a2, r, x2, x1 - q * x2
    if a1 != 1:
        raise ZeroDivisionError("Cannot invert %d mod %d" % (a, p))
    # assert x1 * a % p == 1
    return x1 % p


def mod_exp(a, e, p):
    """Return a**e (mod p) for integers  a, e, p  with p >= 1.

    Raise ZeroDivsionError if this is not possible
    """
    assert p >= 1
    return pow(a, e, p) if e >= 0 else pow(mod_inv(a, p), -e, p) 


def mod_pwr2(e, p):
    """Return 2**e (mod p) for integers e, p with p >= 1.

    Raise ZeroDivsionError if this is not possible. This is
    faster than  mod_exp(2, e, p)  for negative e.
    """
    assert p >= 1
    if e >= 0:
        return pow(2, e, p)
    elif p & 1:
        return pow((p + 1) >> 1, -e, p)
    else:
        raise ZeroDivisionError("Cannot invert 2 mod %d" % p)
 

def mod_rand_invertible(p):
    """Return a random invertible element modulo p for  p > 1"""
    assert p > 1
    while True:  
        a = randint(1, p - 1)
        # Check with the Euclidean algorithm if a is invertible
        a2, a1 = p, a
        while a1 > 0:
            a2, a1 = a1, a2 % a1
        # Now a is invertible mod p iff a1 == 1
        if a2 == 1:
            return a


def mod_rand_unit(p):
    """Return a random unit 1 or -1 modulo p for  p > 1"""
    assert p > 1
    return 1 + (p - 2) * randint(0, 1)
    



######################################################################
# class AbstractRepVector
######################################################################



class AbstractRepVector(object):
    """Model an vector in an abstract representation vector space 

    Here a vector space is an instance of class AbstractRepSpace
    """
    __slots__ = "space"
    def __init__(self, space, *args, **kwds):
        self.space = space

    def __eq__(self, other):    
        return( isinstance(other, AbstractRepVector) 
            and  self.space == other.space
            and  self.space.equal_vectors(self, other)
        )

    def __ne__(self, other): 
        return not self.__eq__(other)

    def copy(self):
        """Return deep copy of the vector"""
        return self.space.copy_vector(self)

    def __imul__(self, other):
        if  isinstance(other, Integral):
            return self.space.imul_scalar(self, other) 
        else:
            return self.space.imul_group_word(self, other)

    def __mul__(self, other):
        return self.copy().__imul__(other)

    def __rmul__(self, other):
        return self.space.imul_scalar(self.copy(), other) 

    def __itruediv__(self, other):
        if  isinstance(other, Integral):
            p = self.p
            try:
                return self.__imul__(mod_inv(other, p)) 
            except ZeroDivisionError:
                if p != 0:
                    raise
                res = self
                if other < 0:
                    res = self.__imul__(-1) 
                    other = -other
                bl = p.bit_length()
                if p != (1 << bl) >> 1:
                    raise
                return res.__ilshift__(other)
                    
        else:
            other = self.space.group(other)
            return self.space.imul_group_word(self, other**(-1))

    def __truediv__(self, other):
        return self.copy().__itruediv__(other)


    def __ilshift__(self, other):
        try: 
            factor = mod_pwr2(other, self.p)
            return self.space.imul_scalar(self, factor) 
        except ZeroDivisionError:
            if p != 0:
                 raise

    def __lshift__(self, other):
        return self.copy().__ilshift__(other)

    def __irshift__(self, other):
        return self.__ilshift__(-other)

    def __rshift__(self, other):
        return self.copy().__ilshift__(-other)

    def __neg__(self):
        return self.space.imul_scalar(self.copy(), -1) 

    def __pos__(self):
        return self  

    def __iadd__(self, other):
        if other in self.space:
            return self.space.iadd(self, other)
        elif isinstance(other, Integral) and other==0:
            return self
        raise TypeError("Vectors must be in the same vector space")

    def __add__(self, other):
        return self.copy().__iadd__(other)

    __radd__ = __add__

    def __isub__(self, other):
        return self.__iadd__(other.__neg__())

    def __sub__(self, other):
        return self.copy().__iadd__(other.__neg__())

    def __rsub__(self, other):
        return self.__neg__.__iadd__(other)

    def __getitem__(self, index):
        return self.space.vector_get_item(self, index) 

    def __setitem__(self, index, value):
        return self.space.vector_set_item(self, index, value) 
        
    def reduce(self):
        return self.space.reduce(self)

    def str(self):
        try:
            return self.space.str_vector(self)
        except NotImplementedError:
            return super(AbstractRepVector, str)()

    __repr__ = str

######################################################################
# class AbstractRepSpace
######################################################################


class AbstractRepSpace(object):
    group = None                    # group operating on the space
    vector_type = AbstractRepVector # type of an vector in the space 
    atom_parser = {}                # see method parse()
    FRAME = re.compile(r"^(\w*)$")  # see method parse()

    def __init__(self, p, *args, **kwds):
        """Creating instances is only possible for concrete spaces

         
        """
        self.p = p


    ### The following methods must be overwritten ####################

    def zero(self):
        """Return the zero vector of the vector space"""
        raise NotImplementedError("No unit vector in abstract space")

    def unit(self,  *args):
        """Return a unit vector of the vector space

        Constructing a unit vector without any arguments should
        return the zero vector.
        """
        raise NotImplementedError("No unit vector in abstract space")

    def iadd(self, v1, v2):
        """Return sum v1 + v2 of vectors v1 and v2.

        Vector v1  may be destroyed, but not v2.

        This method is called for vectors v1, v2 of the space
        'self'  only.
        """
        raise NotImplementedError("No vector addition in abstract space")


    def imul_scalar(self, v1, a):
        """Return product a*v1 of scalar a and vector v1.

        v1 may be destroyed.

        This method is called for elements v1 of the space
        'self' only.

        Here scalars are integers. In case self.p > 0 scalars
        are integers modulo p.

        Here we a bit sloppy, not claiming that p is prime.
        So 'self' may also be a free module over the integers
        modulo p for composite p.
        """
        raise NotImplementedError("No scalar multiplication in abstract space")


    def imul_group_word(self, v1, g):
        """Return product v1*g of vector v1 and group word g.

        v1 may be destroyed.

        This method is called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.
        """
        raise NotImplementedError("No group multiplication in abstract space")




    def vector_get_item(self, v1, index):
        """Return v1[index] for a vector v1 in the vector space"""
        raise NotImplementedError("Cannot get component in abstract space")

    def vector_set_item(self, v1, index, value) :
        """Put v1[index] = value for a vector v1 in the vector space"""
        err = "Item assignment not supported in space of type '%s'"
        raise NotImplementedError(err % type(self))

    def equal_vectors(self, v1, v2):
        """Return True iff vectors v1 and v2 are equal 

        This method is called for elements v1 and v2 of the space
        'self' only.
        """
        raise NotImplementedError("Cannot test equality in abstract space")


    ### The following methods should be overwritten ###################

    def copy_vector(self, v1):
        """Return deep copy of group element v1"""
        v2 = deepcopy(v1)
        v2.space = v1.space
        return v2


    def iadd_tuple(self, v1, t):
        """Compute v += self.unit(*t)

        In a concrete space, this may be optimized considerably.
        """
        raise NotImplementedError("Not supported")
        return self.iadd(v1, self.unit(*t))


    def reduce(self, v):
        """Reduce the vector v which is an element of the space

        We assume that the representation of a vector
        is not always given in a unique form that we call
        the reduced form. 

        This method tries to achieve this goal. Vectors
        may be reduced by any operator, even by '==',
        without notice.
        """
        return v

    def str_vector(self, v):
        """Convert vector v to a string

        """
        raise NotImplementedError


    def ishift_right(self, v1, shift):
        """Return product v1 / 2**shift .

        v1 may be destroyed.

        This method is called for elements v1 of the space
        'self' in case ``v1.p == 0`` only.

        The function should return ``V1 * 2**(-shift)`` is this is 
        defined.
        """
        raise NotImplementedError("Division by two not supported")


    ### The following methods need not be overwritten #################


    def parse(self, s):
        """Convert the string 's' to a vector
        
        Here 's' may ba any valid python expression. Then certain 
        identifiers are recognized e.g. as unit vectors of the 
        space, and operations given by the string 's' are performed 
        on the values of these identifiers. The function tries to 
        convert the result  of the evaluated python expression to 
        a  vector using the standard conversion functions 
        in the __call__  method of this group.

        The identifiers to be recognized are given by the dictionary 
        'self.atom_parser'. This dictionary is empty in the abastract 
        space but it may be set to a specific directory in a 
        specific space.

        A suitable atom parser is AtomDict(self.unit), with AtomDict 
        taken from module parse_atoms.py. This parser maps certain 
        identifiers to tuples and it calls method self.unit() with
        the entries of that tuple as arguments.

        If the string 's' matches the match object self.FRAME, then
        the string self.FRAME.match(s)[1] is parsed instead of 's'.
        This feature may be used to unpack the string 's', e.g. 
        for reversing the effect of a formatting method for 
        vectors. See the standard module re.py for match objects.

        Of course the user may overwrite this method with his own
        parser.
        """
        m = self.FRAME.match(s)
        if m: s = m[1]
        return eval_atom_expression(s, self.atom_parser)

    ### The following methods should not be overwritten ###############

    def set_preimage(self, space, morphism):
        """Set a natural 'morphism' from 'space' to 'self'.

        Here 'morphism' must be a natural homomorphism that maps 
        elements of the space 'space' to elements of the space 
        'self'. 'space' must be an instance of class
        AbstractRepSpace.

        If v is an element of 'space' then self(v) returns
        morphism(v), which is an element of the space self.
        So we may call class 'self', which represents a vector
        space, to map elements of other spaces to the space 'self'. 
        """
        raise NotImplementedError("Not supported")
        self.preimages[space] = morphism


    def __contains__(self, other):
        """Return True iff other is an element of the space"""
        try:
            if not isinstance(other, AbstractRepVector): 
                 return False
            return other.space == self
        except:
            return False

           

    def __call__(self,  *args):
        """Construct a vector of the rep space 'self'.

        The returned vector is the sum of all arguments 
        given to the function. Here each argument v  is 
        converted to a vector as follows:

        - An element v of this space is not changed.

        - If v is an element of another space V and there is a 
          natural hmomorphism phi from V to this space (set by
          method set_preimage) then v is converted to phi(v).

        - A tuple v is converted to self.unit(*tuple).

        - The integer 0 then is converted to the zero vector.

        - If v is a string, we convert v to self.parse(v),
          and perform one of the actions listed above.

        If none of these condition is met, we raise TypeError.
        """
        raise NotImplementedError
        result = self.zero()
        for v in args:
            if isinstance(v, str):
                v = self.parse(v)
            if isinstance(v, tuple):
                self.iadd_tuple(result, v)
            elif isinstance(v, AbstractRepVector):
                if v.space != self:
                    try:
                        morphism = self.preimages[v.space]
                    except KeyError:
                        err = "Cannot map vector to given vector space"
                        raise TypeError(err)
                    v = morphism(v)
                self.iadd(result, v)
            elif v == 0:
                pass
            else:
                err = "Cannot convert '%s' object to a rep vector"
                raise TypeError(err % str(type(v)))
        return self.reduce(result)


















