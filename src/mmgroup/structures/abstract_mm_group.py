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


from mmgroup.structures.abstract_group import AbstractGroup
from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.construct_mm import iter_mm
from mmgroup.structures.construct_mm import iter_strings_from_atoms
from mmgroup.structures.construct_mm import iter_tuples_from_atoms



####################################################################
####################################################################
### Class AbstractGroup and helpers for that class
####################################################################
####################################################################

####################################################################
### Class AbstractGroupWord
####################################################################


class AbstractMMGroupWord(AbstractGroupWord):
    """Model an element of an abstract group.

    This is an abstract base class. Subclasses of this class
    correspond to (preimages of) susbgroups of the monster generated
    by some or all of the standard generators of the monster group.
    Instances of such subclasses are elements of the corresponding 
    groups.

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
    __slots__ = []
    group_name = "AbstractMMGroup"
    def __init__(self, tag = None, atom = None):
        self.from_data(self.group, tag, atom)
  


    def reduce(self):
        """Reduce a group element

        If group elements are implemented as words, some functions
        may produce unreduced words. This function  reduces the
        group element in place.

        Note that all operators return reduced words. Functions return
        reduced words unless stated otherwise. However, reducing all
        words representing the same group element to the same word may 
        be beyond the capabilties of a program. 
        """
        return self.group.reduce(self)

    @property
    def mmdata(self):
        """Return the internal rpresentation of the element

        This method returns the internal representation of an element
        of the monster group as a numpy array of unsigend 32-bit
        integers.

        The internal representation is described in section 
        *The monster group* in the *API reference*.
        """
        raise NotImplementedError("Abstract class")


    def str(self):
        """Represent group element as a string"""
        try:
            return self.group.str_word(self)
        except NotImplementedError:
            return super(AbstractGroupWord, str)()
    __repr__ = str

    def raw_str(self):
        """Represent group element as a string

        In contrast to the standard method ``str``, this method
        returns a string containing the raw data of the group
        element only.
        """
        return self.group.raw_str_word(self)


    def as_tuples(self):
        """Convert group element to a list of tuples

        For a group element ``g`` the following should hold:

        ``self.__class__(self.as_tuples()) == self`` .

        This shows how to convert a group element to a list of tuples
        and vice versa.
        """
        return list(iter_tuples_from_atoms(self.mmdata))



def is_mmgroup_word(g):
    return isinstance(g,AbstractMMGroupWord) and g.group.is_mmgroup 

####################################################################
### Class AbstractGroup
####################################################################


class AbstractMMGroup(AbstractGroup):
    """Model an abstract group"""
    word_type = AbstractMMGroupWord  # type of an element (=word) in the group
    group_name = "AbstactMMGroup"
    is_mmgroup = True

    def __init__(self, *data, **kwds):
        """Creating instances is only possible for concrete groups

         
        """

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

        This method is called for elements g1 and g2 of the grou 'self'
        only. It should return the (possibly unreduced) reduced product.
        """
        raise NotImplementedError("No multiplication in abstract group")

    def _invert(self, g1):
        """Return inverse g1**(-1) of group element g1.

        g1 must not be destroyed.

        This method is called for elements g1 of the group 'self' only. 
        It should return the (possibly unreduced) inverse.
        """
        raise NotImplementedError("No inversion in abstract group")
        
    ### The following methods should be overwritten ###################

    def copy_word(self, g1):
        """Return deep copy of group element ``g1``"""
        raise NotImplementedError("No copy constructor in abstract group")

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
        raise NotImplementedError("No equality check in abstract group")


    def reduce(self, g):
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
        """
        return g


    ### The following methods need not be overwritten #################



    def raw_str_word(self, g):
        """Convert group atom ``g`` to a string

        For an element ``g`` of this group ``g.group.raw_str_word(g)``
        should be equivalent   to ``g.raw_str()``.
        """
        atoms = g.mmdata
        s = "*".join(iter_strings_from_atoms(atoms, abort_if_error=0))
        return s if s else "1"
                 

    def str_word(self, g):
        """Convert group atom ``g`` to a string

        For an element ``g`` of this group ``g.group.str_word(g)``
        should be equivalent   to ``g.str()``.
        """
        return "%s<%s>" % (self.group_name, self.raw_str_word(g))

    def neutral(self):
        """Return neutral element of the group"""
        assert self.word_type.group == self
        return self.word_type()


    def _embed_number(self, n):
        """Try to embed the number n into the group.

        The function returns the number n as a group element or
        it raises TypeError if this is not possible.
        For a matrix group it may e.g. return n times the unit matrix.

        By default, we only map the number 1 to the neutral element.
        """
        raise ValueError("Nixda!!!") # TODO!!!
        if n == 1:
            return self.neutral()
        raise TypeError("Cannot convert a number to a group element")

   




    ### The following methods should not be overwritten ###############



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


    def __call__(self, *args):
        """Convert args to group element

        """
        assert self.word_type.group == self
        return self.word_type(*args)
    



AbstractMMGroupWord.group = AbstractMMGroup

