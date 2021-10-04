"""An exprimental free group with relations


""" 
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import re
from random import randint

from mmgroup.structures.parse_atoms import TaggedAtom
from mmgroup.tests.groups.auto_group import AutoGroup
from mmgroup.tests.groups.auto_group import AutoGroupMulError

####################################################################
####################################################################
### An experimental free group with relations for testing
####################################################################
####################################################################


def f_const(c):
    """Return a function that returns list of tags c on any input"""
    def f(*args):
        return [TaggedAtom(b) for b in c]
    return f

def _failure(*args):
    raise AutoGroupMulError




class FreeGroup(AutoGroup): 

    def __init__(self, generators = "a", commute="", rewrite = {}, bad_commute=""):
        """Model a free group with n generators

        Generators are named a,b,c..., their inverses are A,B,C... .

        If commute is a string then all characters in that string
        commute.
        rewrite is optional dictionary with keys and values both
        being strings. Then the items of the dictionary are
        interpreted as rewriting rules. 

        bad_commute gives an experimental bad matching rule for 
        making all characters in bad_commute commuting. Here every
        commuting pair is matched, and the substitution function
        raises AutoGroupMulError if the pair should should not be
        echanged.
        """
        rules = {}        
       
        for arg in [generators, commute, bad_commute]:
            assert arg == "" or arg.islower()
        for g in generators:
            G = g.upper()
            rules[G + g] = rules[g + G] = f_const("")
           
        for c1 in commute + commute.upper():
            c1l = c1.lower()
            for c2 in commute + commute.upper():
                c2l = c2.lower()
                if c1 != c2:
                    if c1l > c2l:
                        rules[c1 + c2] = f_const(c2 + c1)
                    elif c1 in bad_commute and c2 in bad_commute: 
                        rules[c1 + c2] = _failure
 
        super(FreeGroup, self).__init__({}, rules, {})

    def equal_atoms(self, atom1, atom2):
        return atom1.tag == atom2.tag


    def reduce_atom(self, atom):
        return atom

    def iter_inv(self, atom):
        t = atom.tag   
        yield TaggedAtom(t.lower() if t.isupper() else t.upper())

             
    def __call__(self, w=""):
        """Return the word of the free group given by the string w"""
        assert isinstance(w, str)
        res = self.word_type(*[TaggedAtom(c) for c in w], group = self)
        return res.reduce()

    def str_word(self, word):
        x = "".join((str(atom.tag) for atom in word.iter_generators()))
        return x if x else "<1>"
		

