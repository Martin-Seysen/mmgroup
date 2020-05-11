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
    """
    def __init__(self, group, *args, **kwds):
        self.group = group

    # There is no need to modify an methods below this line.
    # You should overwrite the corresonding methods in the
    # subclasses of class AbstractGroup insead.
 
    def __eq__(self, other):    
        return( isinstance(other, AbstractGroupWord) 
            and  self.group == other.group
            and  self.group.equal_words(self, other)
        )

    def __ne__(self, other): 
        return not self.__eq__(other)

    def copy(self):
        """Return deep copy of element"""
        return self.group.copy_word(self)

    def __imul__(self, other):
        """Implementation of the group multiplication"""
        g = self.group
        return g.imul(self, g.to_group(other))


    def __mul__(self, other):
        """Implementation of the group multiplication"""
        g = self.group
        return g.imul(g.copy_word(self), g.to_group(other))

    def __rmul__(self, other):
        """Implementation of the reverse group multiplication"""
        g = self.group
        return g.imul(g.to_group(other), self)
 

    def __itruediv__(self, other):
        """Implementation of the group division

        Here self / other    means    self * other**(-1) .
        """
        g = self.group
        return g.imul(self, g.invert(g.to_group(other)))

    def __truediv__(self, other):
        """Implementation of the group division

        Here self / other    means    self * other**(-1) .
        """
        g = self.group
        return g.imul(g.copy_word(self), g.invert(g.to_group(other)))

    def __rtruediv__(self, other):
        """Implementation of the reverse group division

        Here self / other    means    self * other**(-1) .
        """
        g = self.group
        return g.imul(g.copy_word(g.to_group(other)), g.invert(self))
      
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
                start, exp = g.invert(self), -exp
                res = g.copy_word(start) 
            for i in range(int(exp).bit_length() - 2, -1, -1):
                res = g.imul(res, res)
                if exp & (1 << i):
                    res = g.imul(res, start) 
            return res      
        elif isinstance(exp, AbstractGroupWord):
            e = self.group.to_group(exp) 
            return g.imul(g.imul(g.invert(e), self), e)


    def __imatmul__(self, other):
        """Non-reduced multiplication g1 @ g2 of elements g1, g2

        g1 @ g2 represents the product g1 * g2. If group elements
        are represented by words then the function returns the 
        concatenation of the words g1 and g2 without reducing it.

        Here g1 and g2 must be elements of the same group.

        Non-reduced multiplication should be used for testing only.
        """
        g = self.group
        if not isinstance(other, AbstractGroupWord or other.group != g):
            raise TypeError("Can concatenete elements of same group only")
        return g.imul_nonreduced(self, other)

    def __matmul__(self, other):
        """Non-reduced multiplication g1 @ g2 of elements g1, g2

        g1 @ g2 represents the product g1 * g2. If group elements
        are represented by words then the function returns the 
        concatenation of the words g1 and g2 without reducing it.

        Here g1 and g2 must be elements of the same group.

        Non-reduced multiplication should be used for testing only.
        """
        g = self.group
        if not isinstance(other, AbstractGroupWord or other.group != g):
            raise TypeError("Can concatenete elements of same group only")
        return g.imul_nonreduced(g.copy_word(self), other)

    def inv(self):
        """Non-reduced inversion of an group element.

        The result represents the inverse of the given group element.
       
        If group elements are interpreted as words, the function
        returns just the reversed word of the inverses of the
        atoms of the input word. Here we make no attempt to 
        reduce the inverse word.

        Non-reduced inversion should be used for testing only.
        """
        return self.group.invert_nonreduced(self)

    def reduce(self, copy = False):
        """Reduce a group element

        If 'copy' is set, a reduced copy of the element is returned
        if the element is not reduced

        If group elements are implemented as word, some functions
        may produce unreduced words. This function  reduces the
        word 'self' in place.

        Note that all operators (except '@') return reduced words.
        """
        return self.group.reduce(self, copy)

    def str(self):
        try:
            return self.group.str_word(self)
        except NotImplementedError:
            return super(AbstractGroupWord, str)()
    __repr__ = str

    def as_tuples(self):
        """Convert group element to a list of tuples"""
        return self.group.as_tuples(self)


################################################### #################
### Class AbstractGroup
####################################################################


class AbstractGroup(object):
    """Model an abstract group"""
    word_type = AbstractGroupWord  # type of an element (=word) in the group
    atom_parser = {}               # see method parse()
    FRAME = re.compile(r"^(\w*)$") # see method parse()

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

    def imul(self, g1, g2):
        """Return product g1 * g2 of group elements g1 and g2.

        g1 may be destroyed but not g2.

        This method is called for elements g1 and g2 of the group
        'self' only. It should return the reduced product.
        """
        raise NotImplementedError("No multiplication in abstract group")


    def invert(self, g1):
        """Return inverse g1**(-1) of group element g1.

        g1 must not be destroyed.

        This method is called for elements g1 of the group
        'self' only. It should return the reduced inverse.
        """
        raise NotImplementedError("No inversion in abstract group")
        
    ### The following methods should be overwritten ###################

    def copy_word(self, g1):
        """Return deep copy of group element g1"""
        g_copy = deepcopy(g1)
        # Even a deep copy of an element is still in the same group
        g_copy.group = g1.group
        return g_copy

    def equal_words(self, g1, g2):
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
        """Reduce the word g which is an element of the group

        We assume that the representation of a group element
        is not always given in a unique form that we call
        the reduced form. 

        This method tries to achieve this goal. Group elements
        are reduced by any operator, except for the
        '==',  '!=', and '@' operators. 

        For test purposes, is is useful to obtain a group 
        element in non-reduced form. Applications should
        create reduced group elements only.

        One way to obtain avoid reduction is to call
        method word() with elements separated by commas.
        Then no reduction takes place across the factors
        separated by commas. 

        If argument 'copy' is True, a reduced copy of g
        must be returned if g is not reduced.

        The '@' operation multiplies elements without reducing
        the product. For a group element g, method g.inv() 
        returns a non-reduced inverse of g.
        """
        return g

    def imul_nonreduced(self, g1, g2):
        """ Non-reduced multiplication g1 @ g2

        g1 @ g2 represents the product g1 * g2. If group elements
        are represented by words then the function returns the 
        concatenation of the words g1 and g2 without reducing it.

        Here g1 and g2 must be elements of the same group.

        (g1 @ g2).reduce()  should be equal to g1 * g2.

        The default implementation does not distinguish between
        non-reduced and reduced multipLication. 

        g1 may be destroyed but not g2.

        This method is called for elements g1 and g2 of the group
        'self' only.
        """
        return self.imul(g1, g2) 


    def invert_nonreduced(self, g1):
        """Return non-reduced inverse g1**(-1) of group element g1.

        The result represents the inverse of the given group element.
       
        If group elements are interpreted as words, the function
        returns just the reversed word of the inverses of the
        atoms of the input word. Here we make no attempt to 
        reduce the inversed word.

        g1 must not be destroyed.

        This method is called for elements g1 of the group
        'self' only.
        """
        return self.invert(g1)


    def as_tuples(self, g):
        """Convert group element g to a list of tuples.

        The tuple should represent a reduce word.

        The sequence 

            l = self.as_tuples(g); g1 = self(*l)

        should compute an element g1 = g of the group 'self'.
        """
        raise NotImplementedError("Abstract method")


    def str_word(self, g):
        """Convert group atom g to a string

        """
        raise NotImplementedError

                 
    ### The following methods need not be overwritten #################

    def neutral(self):
        """Return neutral element of the group"""
        return self.atom()


    def embed_number(self, n):
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
        """Return a random product of atoms

        Assuming that the first parameter of self.atom() is a
        tag, the function samples 'n_samples' tags from the
        list 'tags' of tags and generates one atom with each of
        these tags. It returns the product of these tags. 

        'n_samples' defaults to the length of the list 'tags'.
        If (optional) arguments follow the argument 'tags'
        these arguments are passed to self.atom() after the
        first argument 'tag'. Otherwise 'tag' is the only 
        argument passed to method self.atom().

        The product of the random atoms is not reduced.

        It is useful to create a random atom with a given
        tag if that tag is the only argument passed to
        function self.atom().
        """
        n = len(tags) if n_samples is None else n_samples
        y = (self.atom(tag, *args) for tag in sample(tags, n))
        return self.word(*y)
             
    def rand_word(self, tags, minlen = 1, maxlen = 0, *args):
        """Return a random product of atoms from 'tags'

        A product of at least 'minlen' and at most 'maxlen' atoms
        is returned. Here all atoms are generated independently
        with a random tag taken from the string 'tags'.

        Optional arguments '*args' may be passed to the generator
        of each atom in the same was as in methos sample().
        """
        length = randint(minlen, max(minlen, maxlen))
        used_tags = (sample(tags,1)[0] for i in range(length))
        y = (self.atom(tag, *args) for tag in used_tags)
        return self.word(*y)

    def rand_atom(self, tags, *args):
        """Equivalent to self.sample(tags, 1, *args).

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
        """Set a natural 'morphism' from 'group' to 'self'.

        Here 'morphism' must be a natural homomorphism that maps 
        elements of the group 'group' to elements of the group 
        'self'. 'group' must be an instance of class
        AbstractGroup.

        If g is an element of 'group' then self(g) returns
        morphism(g), which is an element of the group self.
        So we may call class 'self', which represents a group, to
        map elements of other groups to the group 'self'. 

        Here 'morphism' should be natural morphism from 'group'
        to 'self', obtained e.g. by inclusion or by factoring of 
        groups. The diagram of all morphisms should commute.

        Parameter 'morphism' may be set to the value 'tuple'.
        Here we mean the type 'tuple', not the string.
        Then an element of the group 'group' is converted to
        a list of tuples, and that list of tuples in converted
        to an element of the group 'self'. See method
        convert_via_tuple() for details. Such a conversion may 
        be useful in some, but not in all cases.

        Multiplication is done from left to right. If a is an
        element of group g_A and b is an element of g_B, and
        both, g_A(b) and  g_B(a), are defined,  then a * b  
        evaluates to  a * g_A(b).
         
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

    def to_group(self, g):
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
                raise NotImplementedError
            return morphism(g)
        elif isinstance(g, Number):
            return self.embed_number(g)


    def cast(self, g):
        """Convert the object 'g' to an element of this group

        This function tries the conversions on 'g'
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
            raise TypeError("Cannot map group element to given group")
        elif isinstance(g, tuple):
            return self.atom(*g)
        elif g == 1:
            return self.neutral()
        print("""\nTrying to convert type %s object
%s to group element\n""" % (type(g1), g1)
        )
        raise TypeError("Cannot convert an argument to a group element")
           


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

        If none of these condition is met, we raise TypeError.
        """
        result = self.neutral()
        for g in args:
            result = self.imul(result,  self.cast(g))
        return result

    
    def word(self, *args):
        """Construct a word (i.e. an element) of the group.

        The returned group element is the product of all arguments 
        given to the function. Here each argument g  is converted 
        to a group element by calling self(g).

        See the calling method __call__() of this class for
        supported conversions.
        
        We deliberately do not reduce the product of the arguments.
        """
        result = self.neutral()
        for g in args:
            result = self.imul_nonreduced(result, self.cast(g))
        return result
    




