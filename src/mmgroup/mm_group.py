from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import sys
import os
import collections
import re
import warnings
from numbers import Integral
import numpy as np
from random import randint, sample


from mmgroup.structures.parse_atoms import ihex, TaggedAtom
from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.abstract_group import AbstractGroup
from mmgroup.structures.parse_atoms import  AtomDict      
#from mmgroup.structures.parse_atoms import eval_atom_expression        
#from mmgroup.structures.parse_atoms import TaggedAtom

from mmgroup.mat24 import MAT24_ORDER, pow_ploop

from  mmgroup.mm import mm_group_mul_words


###########################################################################
# Word class for the group MM
###########################################################################


class MMGroupWord(AbstractGroupWord):
    MIN_LEN = 16
    def __init__(self, group, data, reduced = 0):
        self.group = group
        self.length = len(data)
        self._data = np.array(data, dtype = np.uint32) 
        self.reduced = reduced 
        self._extend(self.MIN_LEN)
                  
    def _extend(self, length):
        len_ = len(self._data)
        if length > len_:
            ap = np.zeros(max(length - len_, len_), dtype = np.uint32)
            self._data = np.append(self._data, ap)
             
    @property
    def data(self):
        return self._data[:self.length]
        
    def __len__(self):
        return self.length
        
    def __getitem__(self,i):
        if 0 <= i < self.length:
            return self._data[i] 
        if -self.length <= i < 0:
            return self._data[i + self.length]
        raise IndexError
        
    def is_reduced(self):
        return self.length == self.reduced

###########################################################################
# Atoms for the group M
###########################################################################


tag_dict = {
        "d": 0x10000000, 
        "p": 0x20000000, 
        "x": 0x30000000, 
        "y": 0x40000000, 
        "t": 0x50000000, 
        "l": 0x60000000, 
}




tags = " dpxytl"

def gen_d(tag, d = "r"):
    cocode = d
    if isinstance(d, str):
       cocode = randint(int('n' in d), 0xfff) 
       if "o" in d and not "e" in d:
           cocode |= 1 
       if "e" in d and not "o" in d:
           ccocode &= ~1
    return [0x10000000 + (cocode & 0xfff)]

def gen_p(tag, perm = "r"):
    if isinstance(perm, str):
        perm = randint(0, MAT24_ORDER-1)
    else:
        if not 0 <= perm < MAT24_ORDER:
            raise ValueError("Bad permutation number for Mathieu group")
    return [0x20000000 + perm]

def gen_xy(tag, r = "r"):
    pl = randint(0, 0x1fff) if isinstance(r, str) else r & 0x1fff
    return [ tag_dict[tag] + pl]

def gen_z(tag, r = "r"):
    pl = randint(0, 0x1fff) if isinstance(r, str) else r & 0x1fff
    pl = pow_ploop(pl, 3)
    return [tag_dict['y'] + pl, tag_dict['z'] + pl]

def gen_tl(tag, r = "r"):
    e = r
    if isinstance(r, str):
        e = randint(int('n' in r), 2) 
    return  [ tag_dict[tag] + e % 3] 


gen_tag_dict = {
        "d": gen_d, 
        "p": gen_p, 
        "x": gen_xy, 
        "y" :gen_xy, 
        "z" :gen_z, 
        "t": gen_tl, 
        "l": gen_tl, 
}


def gen_atom(*data):
    """Return list of integers representing element of monster group

    Yet to be better documented!
    """
    if len(data) == 0:
         return []
    try: 
        gen_function = gen_tag_dict[data[0]]
    except KeyError:
        err = "Illegal tag %s for MM group atom"
        raise ValueError(err % data[0])
    return gen_function(*data)


###########################################################################
# The class representing the group MM
###########################################################################


class MMGroup(AbstractGroup):
    word_type = MMGroupWord
    tags, formats = " dpxytl", [None, ihex, str, ihex, ihex, str, str]
    atom_parser = {}               # see method parse()
    FRAME = re.compile(r"^M?\<(\w*)\>$") # see method parse()
    STR_FORMAT = r"M<%s>"

    def __init__(self):
        """ TODO: Yet to be documented     


        """
        super(MMGroup, self).__init__()
        self.atom_parser = AtomDict(self.atom)

    def atom(self, *data):
        """Convert arguments *data to a group element"""
        return self.word_type(self, gen_atom(*data))


    def as_tuples(self, g):
        assert g.group == self
        # g = g.reduce(copy = True)
        data = g.data
        if len(data) and max(data) >= 0x70000000:
            raise ValueError("Illegal group element")
        return [(tags[a >> 28], a & 0xfffffff)  
                   for a in data if (a >> 28)]

    @classmethod
    def str_atom(cls, a):
        itag = (a >> 28) & 0xF
        if itag in [0, 8]: 
            return "1"
        if itag >= 8:
            return "(1/%s)" % cls.str_atom(a ^ 0x80000000)
        try:
            tag = cls.tags[itag]  
        except:
            return "(unknown)"
        fmt = cls.formats[itag]
        return tag + "_" + fmt(a & 0xfffffff)
        
   
    def str_word(self, g, fmt = None):
        s = "*".join(map(self.str_atom, g.data)) if len(g) else "1"
        return (fmt if fmt else self.STR_FORMAT) % s

    def reduce(self, g1, copy = False):
        l1 = g1.length
        if g1.reduced < l1:
            if copy:
                g1 = self.copy_word(g1)
            l_tail = l1 - g1.reduced
            g1._extend(l1 + l_tail + 1)
            g1._data[l1 : l1 + l_tail] = g1._data[g1.reduced : l1]
            tail =  g1._data[l1:]
            l1 = mm_group_mul_words(g1._data, g1.reduced, tail, l_tail, 1)
            g1.reduced = g1.length = l1
        return g1

    def imul_nonreduced(self, g1, g2):
        l1, l2 = g1.length, g2.length
        g1._extend(l1 + l2)
        g1._data[l1 : l1 + l2] = g2._data[:l2]
        g1.length = l1 + l2
        return g1
        
    def imul(self, g1, g2):
        l1, l2 = g1.length, g2.length
        g1._extend(l1 + 2 * l2 + 1)
        g1._data[l1 : l1 + l2] = g2._data[:l2]
        l1 += l2
        l_tail = l1 - g1.reduced
        g1._data[l1 : l1 + l_tail] = g1._data[g1.reduced : l1]
        tail = g1._data[l1:]
        l1 = mm_group_mul_words(g1._data, g1.reduced, tail, l_tail, 1)
        g1.reduced = g1.length = l1
        return g1

    def invert_nonreduced(self, g1):
        return self.word_type(self, np.flip(g1.data) ^ 0x80000000)

    def invert(self, g1):
        return self.reduce(self.invert_nonreduced(g1))

    def copy_word(self, g1):
        return self.word_type(self, g1.data, g1.reduced)

    def equal_words(self, g1, g2):
        d1, d2 = g1.reduce(True).data, g2.reduce(True).data
        return len(d1) == len(d2) and (d1 == d2).all()



