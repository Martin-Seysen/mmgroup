from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import sys
import warnings
from random import sample, randint
import re


from copy import deepcopy
from collections.abc import Sequence


from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.abstract_group import AbstractGroup





####################################################################
####################################################################
### Class AutoGroup and helpers for that class
####################################################################
####################################################################



####################################################################
### Exceptions provided for the user
####################################################################


"""Exception for a user-defined rewrite rule

The function implementing a rewrite rule may detect that it cannot
improve the input word given to it as an argument.

In this case it must throw the following excption AutoGroupMulError
to indicate that the application of that rule has failed.

A typical example is the case where a rule is applied to all 
elements of a subgroup G to simplify these elements. Then it
will try to simplify any element e of G and it will return the 
empty word if e is the neutral element of G. But if e cannot be
further simplfied, we must return an indication of failure to 
prevent infinite reapplication of the same rule.
""" 
class AutoGroupMulError(Exception):
    """Multiplication of group elements is not possible"""
    pass
    
####################################################################
### Class ReverseTaggedRewriteNode and helpers for that class
####################################################################



class ReverseTaggedRewriteNode(object):
    def __init__(self, parent = None, tag = None):
        self.children = {}
        self.action = None
        self.parent = parent
        if not parent is None:
            assert not tag in self.parent.children
            parent.children[tag] = self

    def child(self, tag):
        try:
            return self.children[tag]
        except KeyError:
            return ReverseTaggedRewriteNode(self, tag)

    def set_action(self, action):
        if not self.action is None:
            raise ValueError("Duplicate rewrite action")
        self.action = action



class TaggedRewriteSystem(object):
    def __init__(self, rules):
        self.word = None
        self.rules = ReverseTaggedRewriteNode()
        self.translation_rules = {}
        for tags, rule in rules.items():
            self.add_rule(tags, rule)

    def add_rule(self, tags, rule):
        if len(tags) == 1:
            assert not tags[0] in self.translation_rules
            self.translation_rules[tags[0]] = rule
        else:
            node = self.rules
            for tag in reversed(tags):
                node = node.child(tag)
            node.set_action(rule)

    def didactic_subst_word(self, group, word):
        """

        Didactic version: starts with shortest possible rule, and
        orders priority of rules by ascending length.
        """
        node = self.rules
        for i, w in enumerate(reversed(word)):
            try:
                node = node.children[w.tag]
            except KeyError:
                return 0, None
            if node.action:
                try:
                    return i + 1, node.action(group, word[-i - 1:])
                except AutoGroupMulError:
                    continue
        return 0, None
                     
    def subst_word(self, group, word):
        """

        Pragmatic version: starts with rule of length 2, then it tries 
        rules of length 3, 4, 5,... , and finally, it tries rule of
        length 1 if present.

        Dealing too much with rules of length 1 is a bit like self
        satifsction and causes too much overhead.
        """
        node = self.rules
        for i, w in enumerate(reversed(word)):
            try:
                node = node.children[w.tag]
            except KeyError:
                break
            if node.action:
                try:
                    return i + 1, node.action(group, word[-i - 1:])
                except AutoGroupMulError:
                    continue
        if len(word) and word[-1].tag in self.translation_rules:
            rule = self.translation_rules[word[-1].tag]
            try:
                return 1, rule(group, word[- 1:])
            except AutoGroupMulError:
                return 0, None
        return 0, None




####################################################################
####################################################################
### Class AutoGroup and helpers for that class
####################################################################
####################################################################






####################################################################
### Class AutoGroupWord
####################################################################




class AutoGroupWord(AbstractGroupWord):
    """TODO: Yet to be documented     

    """

    def __init__(self,  *generators, **kwds):
        """Create word a in group from the generators given as arguments.

        Users should not refer to this class directly.
        """
        self.group = kwds['group']
        assert isinstance(self.group, AutoGroup)
        self.seq =  []
        self.reduced = 0
        for generator in generators:
            g = self.group.reduce_generator(generator)
            if g:
                self.seq.append(g)

    def __getitem__(self, i):
        return AutoGroupWord(self.group, self.seq[i])

    def __len__(self):
        return len(self.seq)

    def iter_generators(self):
        for generator in self.seq:
            yield generator

    def append_word(self, word):
        self.seq.extend(word.seq)
        

    def iter_inv(self):
        for generator in reversed(self.seq):
            yield from self.group.iter_inv(generator) 

    def str(self):
        return self.group.str_word(self)

    __repr__ = str

    def to_list(self):
        """Deprecated!

        For compatiblity with older functions only!
        """
        raise NotImplementedError("Deprecated")
        return [x.to_tuple() for x in self.seq]


####################################################################
### Class AutoGroup
####################################################################


ABORT_MANY_REWRITES = 50000
ABORT_MANY_REWRITES_TXT = (
  "Reduction in class AutoGroup aborted after %s rewrites" % (
                         ABORT_MANY_REWRITES))



class AutoGroup(AbstractGroup):   

    FRAME = re.compile(r"^(\w*)$")
    STR_FORMAT = "%s"

 
    word_type = AutoGroupWord  # type of an element (=word) in the group 

    def __init__(self, creator, rules, inverter):
        """ TODO: Yet to be documented     


        """
        self.rewrite_system = TaggedRewriteSystem(rules)
        self.substituter = self.rewrite_system.subst_word
        self.creator = creator
        self.inverter = inverter
        super(AutoGroup, self).__init__()


    ################################################################
    # There is no need to overwrite any methods below this line
    ################################################################


    def generator(self, *data):
        """Convert tuple *data to a group word"""
        if len(data):
            creator = self.creator[data[0]](*data)
            return self.word_type(*creator, group=self)
        return self.word_type(group=self)

        

    def reduce_generator(self, generator):
        return generator.reduce() if generator else None



    def str_word(self, g1):
        if len(g1.seq):
            return self.STR_FORMAT  % "*".join(map(str, g1.seq))
        return self.STR_FORMAT % "1"


    def reduce(self, word):
        if word.reduced >= len(word.seq):
            return word
        w = word.seq
        gap = 0
        pos = word.reduced + 1
        min_gap = 64
        n_iterations = 0
        #print("start", w)
        while True:
            length, subst = self.substituter(self, w[:pos])
            #print("subst", length,  w[:pos], pos, subst, gap,  w[pos+gap:])
            if length:
                pos -= length
                gap += length - len(subst)
                if  gap < 0:
                    insert = max(min_gap, -gap)
                    w[pos:pos] = [None] * insert
                    gap += insert
                    min_gap *= 2
                subst = [self.reduce_generator(x) for x in subst]
                w[pos + gap : pos + gap + len(subst)] = subst 
            while pos + gap < len(w) and not w[pos + gap]:
                gap += 1
            if pos + gap >= len(w):
                w[pos:] = []
                word.seq = w
                word.reduced = pos
                #print("result", w)
                return word
            w[pos] = w[pos + gap]
            pos += 1
            #print("after:", w[:pos], pos, gap,  w[pos+gap:])
            n_iterations += 1 
            if n_iterations == ABORT_MANY_REWRITES:    
                warnings.warn(ABORT_MANY_REWRITES_TXT)
                w[pos:pos+gap] = []
                word.seq = w
                word.reduced = pos
                return word

    def copy_word(self, g1):
        """Return deep copy of group element g1"""
        result = self.neutral()
        result.seq = deepcopy(g1.seq)
        result.reduced = g1.reduced
        return result


    def _equal_words(self, g1, g2):
        """Return True iff elements g1 and g2 are equal 

        This method is called for elements g1 and g2 of the group
        'self' only.
        """
        g1 = self.reduce(g1.copy())
        g2 = self.reduce(g2.copy())
        tup1 = self.iter_to_tuples(g1)
        tup2 = self.iter_to_tuples(g2)
        equ = (x == y for x, y in zip(tup1, tup2))
        return len(g1) == len(g2) and all(equ)


    def _imul(self, g1, g2):
        """Return product g1*g2 of group elements g1 and g2.

        g1 may be destroyed but not g2.

        This method is called for elements g1 and g2 of the group
        'self' only.
        """
        assert g1.group == g2.group == self, (g1.group, g2.group)
        g1.seq += g2.seq[:]
        return self.reduce(g1)


    def iter_inv(self, generator):
        """yield inverse of the generator

        The function must yield one or more generators of type
        TaggedAtom, so that the product of these generators is the 
        inverse of the input 'generator'.
        """
        yield from self.inverter[generator.tag](self, generator)
		
    def _invert(self, g1):
        """Return product g1**-1 for group elements g1.

        g1 must not be destroyed.

        This method is called for elements g1 of the group
        'self' only.
        """
        result = self.neutral()
        assert g1.group == self
        result.seq = [
           inv  for g1_generator in reversed(g1.seq)
               for inv in self.iter_inv(g1_generator)
        ]
        return result.reduce()


    def iter_to_tuples(self, g1):
        assert g1.group == self
        return (t for generator in g1.seq for t in generator.as_tuples())


    def as_tuples(self, g1):
        return list(self.iter_to_tuples(g1))




       