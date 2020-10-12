from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


"""Modelling an abstract group

Yet to be documented

"""
import sys
import re
import warnings
from copy import deepcopy
from random import sample, randint
from numbers import Integral, Number

from mmgroup.structures.parse_atoms import  eval_atom_expression        
from mmgroup.structures.parity import Parity



####################################################################
####################################################################
### Class AbstractGroup and helpers for that class
####################################################################
####################################################################

####################################################################
### Class AbstractGroupWord
####################################################################


class AbstractGroupWord(object):
    """Model an element of an abstract group.

    Users should not refer to this class directly. They should create 
    a group as an instance of subclass of class AbstractGroup and use 
    the methods of that group for creating elements.

    The standard group operations '*' '/' (=right multiplication with
    the inverse) and '**' are implemented here.

    g1 ** g2  means  g2**(-1) * g1 * g2  for group elements g1, g2.

    Here a group is an instance of (a subclass of) class AbstractGroup. 

    For each word a group 'g' should be passed as a keyword argument 
    'group = g'. If a class of type 'AbstractGroup' contains one 
    instance only, the corresponding subclass of this class may 
    contain a class attribute 'group' referring to that group. Then 
    the user  may contruct elements of that group using the 
    constructor of that subclass of this class.
    """
    __slots__ = "group"
    def __init__(self, *args, **kwds):
        try:
            self.group = kwds['group']
        except:
            assert isinstance(self.group, AbstractGroup)

    # There is no need to modify an methods below this line.
    # You should overwrite the corresonding methods in the
    # subclasses of class AbstractGroup insead.
 
    def __eq__(self, other):    
        return( isinstance(other, AbstractGroupWord) 
            and  self.group == other.group
            and  self.group._equal_words(self, other)
        )

    def __ne__(self, other): 
        return not self.__eq__(other)

    def copy(self):
        """Return a deep copy of the group element"""
        return self.group.copy_word(self)

    def __imul__(self, other):
        """Implementation of the group multiplication"""
        g = self.group
        return g._imul(self, g._to_group(other))


    def __mul__(self, other):
        """Implementation of the group multiplication"""
        g = self.group
        return g._imul(g.copy_word(self), g._to_group(other))

    def __rmul__(self, other):
        """Implementation of the reverse group multiplication"""
        g = self.group
        if isinstance(other, Parity):
            return other
        return g._imul(g._to_group(other), self)
 

    def __itruediv__(self, other):
        """Implementation of the group division

        Here self / other    means    self * other**(-1) .
        """
        g = self.group
        return g._imul(self, g._invert(g._to_group(other)))

    def __truediv__(self, other):
        """Implementation of the group division

        Here self / other    means    self * other**(-1) .
        """
        g = self.group
        return g._imul(g.copy_word(self), g._invert(g._to_group(other)))

    def __rtruediv__(self, other):
        """Implementation of the reverse group division

        Here self / other    means    self * other**(-1) .
        """
        g = self.group
        return g._imul(g.copy_word(g._to_group(other)), g._invert(self))
      
    def __pow__(self, exp):
        """Implementation of the power operation

        This is exponentiation for integer eponents and conjugation
        if the exponent is a group element.
        """
        g = self.group
        if isinstance(exp, Integral):
            if exp > 0:
                res, start = g.copy_word(self), self
            elif exp == 0:
                return g.neutral()
            else:
                start, exp = g._invert(self), -exp
                res = g.copy_word(start) 
            for i in range(int(exp).bit_length() - 2, -1, -1):
                res = g._imul(res, res)
                if exp & (1 << i):
                    res = g._imul(res, start) 
            return res      
        elif isinstance(exp, AbstractGroupWord):
            e = self.group._to_group(exp) 
            return g._imul(g._imul(g._invert(e), self), e)
        elif isinstance(other, Parity):
            one = self.group.neutral()
            if self * self == one:
                return self if other.value & 1 else one
            raise ValueError("Group element has not order 1 or 2")
        else:
            return NotImplemented


    def reduce(self, copy = False):
        """Reduce a group element

        If group elements are implemented as words, some functions
        may produce unreduced words. This function  reduces the
        group element in place.

        Note that all operators return reduced words. Functions return
        reduced words unless stated otherwise. However, reducing all
        words representing the same group element to the same word may 
        be beyond the capabilties of a program. 

        If ``copy`` is set then a reduced copy of the element is 
        returned, in case that the input element is not already 
        reduced.
        """
        return self.group.reduce(self, copy)

    def str(self):
        """Represent group element as a string"""
        try:
            return self.group.str_word(self)
        except NotImplementedError:
            return super(AbstractGroupWord, str)()
    __repr__ = str

    def as_tuples(self):
        """Convert group element to a list of tuples

        For a group element ``g`` the following should hold:

        ``g.group.word(*g.as_tuples()) == g`` .

        So passing the tuples in the list returned by this method
        as arguments to ``g.group`` or to ``g.group.word`` 
        reconstructs the element ``g``.

        This shows how to convert a group element to a list of tuples
        and vice versa.
        """
        return self.group.as_tuples(self)




####################################################################
### Class AbstractGroup
####################################################################


class AbstractGroup(object):
    """Model an abstract group"""
    word_type = AbstractGroupWord  # type of an element (=word) in the group
    atom_parser = {}               # see method parse()
    FRAME = re.compile(r"^(\w*)$") # see method parse()
    conversions = {}               # see method __call__
    ERR_CONV = "Cannot convert '%s' object to group word"

    def __init__(self, *data, **kwds):
        """Creating instances is only possible for concrete groups

         
        """
        self.preimages = {}  # dictionary of  preimages

    ### The following methods must be overwritten ####################

    def atom(self, *args):
        """Return an atomic element of this group. 

        Calling this function without any arguments should return
        the neutral element of this group.
        """
        raise NotImplementedError("No atoms defined in abtract group")  

    def _imul(self, g1, g2):
        """Return product g1 * g2 of group elements g1 and g2.

        g1 may be destroyed but not g2.

        This method is called for elements g1 and g2 of the group
        'self' only. It should return the reduced product.
        """
        raise NotImplementedError("No multiplication in abstract group")

    def _imul_nonreduced(self, g1, g2):
        """ Non-reduced multiplication g1 * g2

        The result of this mathod the product g1 * g2. If group 
        elements are represented by words then the function returns  
        the concatenation of the words g1 and g2 without reducing it.

        Here g1 and g2 must be elements of the same group.

        The default implementation does not distinguish between
        non-reduced and reduced multipLication. 

        g1 may be destroyed but not g2.

        This method is called for elements g1 and g2 of the group
        'self' only.

        This method is used in method word() of this class for
        constructing a word of the group without reducing it.
        """
        return self._imul(g1, g2) 


    def _invert(self, g1):
        """Return inverse g1**(-1) of group element g1.

        g1 must not be destroyed.

        This method is called for elements g1 of the group
        'self' only. It should return the reduced inverse.
        """
        raise NotImplementedError("No inversion in abstract group")
        
    ### The following methods should be overwritten ###################

    def copy_word(self, g1):
        """Return deep copy of group element ``g1``"""
        g_copy = deepcopy(g1)
        # Even a deep copy of an element is still in the same group
        g_copy.group = g1.group
        return g_copy

    def _equal_words(self, g1, g2):
        """Return True iff elements g1 and g2 are equal 

        This method is called for elements g1 and g2 of the group
        'self' only.
		
        In concrete group this method should be overwritten with
        a comparison of the relevant attributes of g1 and g2.

        Caution:
        Non-reduced words may be considered unequal even if they
        represent the same element. Use g1.reduce() or g1 * 1 to
        obtain the reduced form of g1. See method reduce() for
        details.
        """
        return g1 == g2


    def reduce(self, g, copy = False):
        """Reduce the word ``g`` which is an element of the group

        We assume that the representation of a group element
        is not always given in a unique form that we call
        the reduced form. 

        This method tries to achieve this goal. Group elements
        are reduced by any operator, except for the
        ``==`` and  ``!=`` operators. 

        For test purposes, is is useful to obtain a group 
        element in non-reduced form. Applications should
        create reduced group elements only.

        One way to obtain avoid reduction is to call method 
        ``word()`` of this class with elements separated by 
        commas. Then no reduction takes place across the factors
        separated by commas. 

        If argument ``copy`` is True, a reduced copy of ``g``
        should be returned if ``g`` is not reduced.
        """
        return g




    def as_tuples(self, g):
        """Convert group element ``g`` to a list of tuples.

        The returned tuple should represent a reduced word.

        The sequence:: 

            l = g.group.as_tuples(g) 
            g1 = g.group(*l)

        should compute a group element ``g1`` with ``g1 == g``.
        """
        raise NotImplementedError("Abstract method")


    def str_word(self, g):
        """Convert group atom ``g`` to a string

        For an element ``g`` of this group ``g.group.str_word(g)``
        should be equivalent   to ``g.str()``.
        """
        raise NotImplementedError

                 
    ### The following methods need not be overwritten #################

    def neutral(self):
        """Return neutral element of the group"""
        return self.atom()


    def _embed_number(self, n):
        """Try to embed the number n into the group.

        The function returns the number n as a group element or
        it raises TypeError if this is not possible.
        For a matrix group it may e.g. return n times the unit matrix.

        By default, we only map the number 1 to the neutral element.
        """
        if n == 1:
            return self.neutral()
        raise TypeError("Cannot convert a number to a group element")

   

    def parse(self, s):
        """Convert the string 's' to a group element
        
        Here 's' may ba any valid python expression. Then certain 
        identifiers are recognized e.g. as atomic elements of the 
        group, and operations given by the string 's' are performs 
        on the values of these identifiers. The function tries to 
        convert the result  of the evaluated python expression to 
        a  group element using the standard conversion functions 
        in the __call__  method of this group.

        The identifiers to be recognized are given by the dictionary 
        'self.atom_parser'. This dictionyry is empty in the abastract 
        group but it may be set to a specific directory in a 
        specific group.

        A suitable atom parser is AtomDict(self.atom), with AtomDict 
        taken from module parse_atoms.py. This parser maps certain 
        identifiers to tuples and it calls method self.atom() with
        the entries of that tuple as arguments.

        If the string 's' matches the match object self.FRAME, then
        the string self.FRAME.match(s)[1] is parsed instead of 's'.
        This feature may be used to unpack the string 's', e.g.
        for reversing the effect of a formatting method for group
        elements. See the standard module re.py for match objects.

        Of course the user may overwrite this method with his own
        parser.
        """
        m = self.FRAME.match(s)
        if m: s = m[1]
        return eval_atom_expression(s, self.atom_parser)

    def sample(self, tags, n_samples = None, *args):
        """Return a random product of atoms with given tags

        A random product of  ``n_samples`` atomic elements of 
        the group  is returned. 
        Therfore a list of  ``n_samples`` tags is sampled from
        the string or list ``tags``. For each tag this method
        calls method ``atom()`` of this class with that tag as  
        the first argument in order to generate an atomic element
        of the group. It returns the product of all these atomic 
        elements. ``n_samples`` defaults to the length of the 
        string or list ``tags``.

        If (optional) arguments ``*args`` follow the argument 
        ``n_samples`` then  these arguments are passed to method 
        ``atom()`` after the first argument ``tag``. Otherwise 
        ``tag`` is the only argument passed to method ``atom()``.

        The product of the random atoms is not reduced.

        Method ``sample()`` is most useful if calling method 
        ``atom()`` with a single parameter ``tag`` returns a 
        random atomic element with that tag. 
        """
        n = len(tags) if n_samples is None else n_samples
        y = (self.atom(tag, *args) for tag in sample(tags, n))
        return self.word(*y)
             
    def rand_word(self, tags, minlen = 1, maxlen = 0, *args):
        """Return a random product of atoms with given ``tags``

        This is a variant of method ``sample()``.

        A product of at least ``minlen`` and at most ``maxlen`` 
        atoms is returned. Here all atoms are generated 
        independently with a random tag taken from the string 
        or list ``tags``. The value ``minlen`` defaults to ``1``
        and  ``maxlen`` defaults to ``minlen``.

        The product of the random atoms is not reduced.

        Additional optional arguments ``*args`` may be passed to 
        each call of method  ``atom()`` in the same way as in 
        method ``sample()``.
        """
        length = randint(minlen, max(minlen, maxlen))
        used_tags = (sample(tags,1)[0] for i in range(length))
        y = (self.atom(tag, *args) for tag in used_tags)
        return self.word(*y)

    def rand_atom(self, tags, *args):
        """Equivalent to ``self.sample(tags, 1, *args)``.

        Thus the generated random word consists of one atom.
        """
        return self.sample(tags, 1, *args)



    ### The following methods should not be overwritten ###############


    def convert_via_tuple(self, g1):
        """Convert g1 to group element via method as_tuples()

        This is an auxiliary method for method set_preimage().
        
        If parameter 'morphism' in method set_preimage() is set
        to the value 'tuple' the elements of the group specified
        in that method are mapped to elements of group 'self'
        as follows:

        g1 is converted to a list of tuples via method as_tuples().
        Then the tuples in that list are passed to the constructor
        of this class parameters.
        
        This conversion method may be a group homomorphism in 
        some, but not in all cases. 
        """
        return self(*(g1.reduce(True).as_tuples()))



    def set_preimage(self, group, morphism):
        """Set a natural morphism from anotehr group to this group.

        Here ``morphism`` must be a natural homomorphism that maps 
        elements of the group ``group`` to elements of this group. 
        ``group`` must be an instance of class ``AbstractGroup``.

        If ``g`` is an element of ``group`` then calling the
        constructor of this class with argument ``g`` return
        ``morphism(g)``, which is an element of this group.
        So we may call the constructur of this group to
        map elements of other groups to this group. 

        Here ``morphism`` should be natural morphism from 
        another group ``group`` to group modelled by this 
        isntance, obtained e.g. by inclusion or by factoring of 
        groups. The diagram of all morphisms should commute.

        Parameter ``morphism`` may be set to the value ``tuple``.
        Here we mean the type ``tuple``, not the string.
        Then an element of the group ``group`` is converted to
        a list of tuples, and that list of tuples in converted
        to an element of this group. See method
        ``convert_via_tuple()`` for details. Such a conversion 
        may be useful in some, but not in all cases.

        Multiplication is done from left to right. If ``a`` is 
        an element of group ``g_A`` and ``b`` is an element of 
        ``g_B``, and both, ``g_A(b)`` and ``g_B(a)``, are defined,  
        then ``a * b`` evaluates to ``a * g_A(b)``.
         
        We make no attempt to compute the transitive closure 
        of the morphisms defined by this function.
        """
        if morphism == tuple:
            morphism = self.convert_via_tuple
        self.preimages[group] = morphism



    def __contains__(self, other):
        """Return True iff 'other' is an element of the group"""
        try:
            if not isinstance(other, self.word_type): 
                 return False
            if other.group != self: 
                 return False
            return True
        except:
            return False

    def _to_group(self, g):
        """Try to convert element g to an element of the group g

        This function is used in the implementation of the group
        operation.
        """
        if isinstance(g, AbstractGroupWord):
            og = g.group
            if og == self:
                return g
            try:
                morphism = self.preimages[og]
            except KeyError:
                raise TypeError(self.ERR_CONV % type(g))
            return morphism(g)
        elif isinstance(g, Number):
            return self._embed_number(g)
        else:
            raise TypeError(self.ERR_CONV % type(g))


    def _cast(self, g):
        """Convert the object ``g`` to an element of this group

        This function tries the conversions on ``g``
        as documented in the '__call__' method.

        In case of success it returns an element of the group
        'self'. Otherwise it raises TypeError.
        """
        g1 = g
        if isinstance(g, str):
            g = self.parse(g)
        if isinstance(g, AbstractGroupWord):
            og = g.group
            if og == self:
                return g
            for group, morphism in  self.preimages.items():
                if og == group:
                    return morphism(g)
            raise TypeError(self.ERR_CONV % type(g))
        elif isinstance(g, tuple):
            return self.atom(*g)
        elif g == 1:
            return self.neutral()
        elif type(g) in self.conversions:
            return self.conversions[type(g)](self, g)
        raise TypeError(self.ERR_CONV % type(g))
           


    def __call__(self, *args):
        """Convert args to group elements and return their product

        For each argument g the function tries the following 
        conversions on g:

        - If g is an element of this group then we return g.

        - If g is an element of another group G and there is a 
          natural hmomorphism phi from G to this group (set by
          method set_preimage) then we return phi(g).

        - If g is a tuple then we return self.atom(*tuple)

        - If g is equal to 1 then we return the neutral 
          element of the group.

        - If g is a string, we replace g by self.parse(g),
          and we try to convert the result to a group element
          using the conversion methods listed above.

        - If type(g) is contained in the dictionary
          'cls.conversions' of the class, then g is replaced
          cls.conversions[type(g)](self, g). Conversion should 
          be considered as fixed built-in functions of a class
          representing a group.

        If none of these condition is met, we raise TypeError.
        """
        if (len(args)):
            result = self._cast(args[0])
            for g in args[1:]:
                result = self._imul(result,  self._cast(g))
            return result
        return self.neutral()

    
    def word(self, *args):
        """Construct a word (i.e. an element) of the group ``G``.
 
        Here ``G`` is an instance of this class representing a group.

        The returned group element is the product of all arguments 
        given to the function. Here each argument ``g``  is converted 
        to a group element by calling ``G(g)``. So ``G.word(g)`` is
        equivalent to ``G(g)``.

        If elements of the group ``G`` are represented as words then
        ``G(g1, g2)`` computes ``G(g1) * G(g2)`` and tries to 
        represent the result as a *reduced* word in ``G``.
        
        On the other hand, ``G.word(g1, g2)`` simply concatenates
        the words ``G(g1)`` and ``G(g2)`` in order to represent the
        product  ``G(g1) * G(g2)``. In some cases this feature can
        be useful for testing.
        """
        if (len(args)):
            result = self._cast(args[0])
            for g in args[1:]:
                result = self._imul_nonreduced(result,  self._cast(g))
            return result
        return self.neutral()
    




