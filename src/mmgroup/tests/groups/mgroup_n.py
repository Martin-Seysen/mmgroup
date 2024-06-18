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


from mmgroup.structures.parse_atoms import  AtomDict      
#from mmgroup.structures.parse_atoms import eval_atom_expression        
from mmgroup.structures.parse_atoms import TaggedAtom
from mmgroup.structures.parse_atoms import ihex
from mmgroup.structures.construct_mm import iter_mm       
from mmgroup.structures.construct_mm import load_group_name     
from mmgroup.structures.construct_mm import add_to_embedded_classes     
from mmgroup.structures.abstract_mm_group import AbstractMMGroup


from mmgroup.tests.groups.auto_group import AutoGroupWord
from mmgroup.tests.groups.auto_group import AutoGroup
from mmgroup.tests.groups.auto_group import AutoGroupMulError


from mmgroup import mat24




###########################################################################
# Atoms for the group N 
###########################################################################


n_rules = 0
n_inv = 0 
n_symbols = 0




class AutPLoopAtom(TaggedAtom):
    __slots__ = 'tag', 'm24num', 'perm', 'rep', 'cocode'
    def __init__(self, tag, cocode, perm=0):
        self.tag = tag
        if isinstance(cocode, Integral):
            self.cocode = cocode & 0xfff
            if perm:
                if isinstance(perm, Integral):
                    self.m24num = perm
                else:
                    self.m24num = mat24.perm_to_m24num(perm)
                self.perm = mat24.m24num_to_perm(self.m24num)
                self.rep = mat24.perm_to_autpl(cocode, self.perm)
            else:
                self.m24num =  0
                self.perm = mat24.m24num_to_perm(0)
                self.rep = mat24.cocode_to_autpl(cocode)
        else:
            self.rep = cocode
            self.perm = mat24.autpl_to_perm(self.rep)
            self.cocode = mat24.autpl_to_cocode(self.rep) 
            self.m24num = mat24.perm_to_m24num(self.perm)

    def reduce(self):
        self.cocode  &= 0xfff
        return None if self.cocode == self.m24num == 0 else self


    def as_atoms(self):
        yield 0x10000000 + self.cocode
        yield 0x20000000 + self.m24num

    def as_tuples(self):
        yield "d", self.cocode
        yield "p", self.m24num

    def str(self):
        prod = []
        c = self.cocode & 0xfff
        if c:  
            prod.append("d_%s" % ihex(c)) 
        if self.m24num: 
            prod.append("p_%d" % self.m24num)
        s = "*".join(prod)
        return s if s else "1"
    __repr__ = str

tag_dict = {
        "x": 0x30000000, 
        "y": 0x40000000, 
        "t": 0x50000000, 
        "l": 0x60000000, 
}

class xy_Atom(TaggedAtom):
    __slots__ = 'tag', 'pl'
    format_tuple = (str, ihex)

    def __init__(self, tag, pl):
        self.tag = tag
        self.pl = pl & 0x1fff

    def reduce(self):
        self.pl  &= 0x1fff
        return None if self.pl == 0 else self

    def as_atoms(self):
        yield tag_dict[self.tag] + self.pl

    def as_tuples(self):
        yield self.tag, self.pl


class tl_Atom(TaggedAtom):
    __slots__ = 'tag', 'exp'
    format_tuple = (str, str)

    def __init__(self, tag, exp):
        self.tag = tag
        self.exp = exp % 3

    def reduce(self):
        self.exp  %= 3
        return None if self.exp == 0 else self

    def as_atoms(self):
        yield tag_dict[self.tag] + self.exp

    def as_tuples(self):
        yield self.tag, self.exp




###########################################################################
# Generating atoms of the group N from tuples
###########################################################################

# deprecated

def gen_d(tag, c = "r"):
    cocode = c
    if isinstance(c, str):
       cocode = randint(0,0xfff)
       if "o" in c and not "e" in c:
           cocode |= 0x800 
       if "e" in c and not "o" in c:
           cocode &= ~0x800
    yield AutPLoopAtom('p', cocode & 0xfff, 0)

def gen_p(tag, perm = "r"):
    if isinstance(perm, str):
        perm = randint(0, mat24.MAT24_ORDER-1)
    else:
        if not 0 <= perm < mat24.MAT24_ORDER:
            raise ValueError("Bad permutation number for Mathieu group")
    yield AutPLoopAtom('p', 0, perm)

def gen_xy(tag, r = "r"):
    pl = randint(0, 0x1fff) if isinstance(r, str) else r & 0x1fff
    yield xy_Atom(tag, pl & 0x1fff)


def gen_z(tag, r = "r"):
    pl = randint(0, 0x1fff) if isinstance(r, str) else r & 0x1fff
    pl = mat24.pow_ploop(pl, 3)
    yield xy_Atom('y', pl)
    yield xy_Atom('x', pl)



def gen_tl(tag, e = "r"):
    if isinstance(e, str):
        e = randint( (0 if "z" in e else 1), 2) 
    else:
        e = e % 3
    yield tl_Atom(tag, e) 




###########################################################################
# Rules deleting neutral elements
###########################################################################


def rule_x0(group, word):
    x = word[0].pl & 0x1fff
    if x == 0:
        return []
    else:
        raise AutoGroupMulError

def rule_p0(group, word):
    d, p = word[0].cocode & 0xfff, word[0].m24num
    if (d | p)  == 0:
        return []
    else:
        raise AutoGroupMulError


def rule_tl0(group, word):
    exp = word[0].exp & 3
    if exp == 0:
        return []
    else:
        raise AutoGroupMulError


###########################################################################
# Rules for the group N 
###########################################################################



def rule_pp(group, word):
    global n_rules
    n_rules += 1
    rep = mat24.mul_autpl(word[0].rep, word[1].rep)
    return  [AutPLoopAtom('p', rep)]


def rule_yp(group, word):
    global n_rules
    n_rules += 1
    y, p = word[0], word[1]
    pl = mat24.op_ploop_autpl(y.pl, p.rep)
    if p.cocode & 0x800:
        pl = mat24.pow_ploop(pl, 3)
        return [p, xy_Atom('y', pl), xy_Atom('x', pl)]
    else:
        return [p, xy_Atom('y', pl)]


def rule_xp(group, word):
    global n_rules
    n_rules += 1
    x, p = word[0], word[1]
    pl = mat24.op_ploop_autpl(x.pl, p.rep)
    return [p, xy_Atom('x', pl)] 





 
# Extra relations given by Kernels K0 and K1: 
# K0:  y_(-1) = x_(-Omega), y_Omega = x_(-1), y_(-Omega) = x_Omega
# K1:  y_(-1) = x_(-Omega)
# Anything apart from kernel K0 isno longer supported!!!
dict_kernel0_yx = {0:0, 0x1000:0x1800, 0x800:0x1000, 0x1800:0x800, }
#dict_kernel1_yx = {0:0, 0x1000:0x1800}

dict_kernel0_xy = {0:0, 0x1000:0x800, 0x800:0x1800, 0x1800:0x1000, }
#dict_kernel1_xy = {0:0, 0x1000:0x800}


#kernel_data = [
#   (0x1800, dict_kernel0_yx, dict_kernel0_xy),  # for kernel K_0
#   (0x1000, dict_kernel1_yx, dict_kernel1_xy),  # for kernel K_1
#]



def kernel_rule_yx(group, word):
    global n_rules
    n_rules += 1
    y, x = word[0].pl & 0x1fff,  word[1].pl & 0x1fff
    mask, table_yx, table_xy = 0x1800, dict_kernel0_yx, dict_kernel0_xy
    if y & ~mask and not x & ~mask:
        y ^= table_xy[x]
        return  [xy_Atom('y', y)]
    else:
        x ^= table_yx[y & mask]
        y_new = y & ~mask 
        if x == 0:
            return  [xy_Atom('y', y_new)]
        elif y_new == y:
            raise AutoGroupMulError
        else:
            return [ xy_Atom('y', y_new),  xy_Atom('x', x) ]


def kernel_rule_y(group, word):
    global n_rules
    n_rules += 1
    y = word[0].pl & 0x1fff
    mask, table  =  0x1800, dict_kernel0_yx
    if y & mask == y:
        return [ xy_Atom('x', table[y]) ]
    else:
        raise AutoGroupMulError

    mask = y & mask
    x = table[y & mask]
    y_new = y & ~mask & 0x1fff
    #print("kernel", hex(y), hex(y_new))
    if y == 0:
        return []
    if y_new == y:
        raise AutoGroupMulError
    return [ xy_Atom('y', y_new),  xy_Atom('x', x) ]


    
def rule_xx(group, word):
    global n_rules
    n_rules += 1
    x1, x2 = word[0], word[1]
    delta = AutPLoopAtom('p', mat24.ploop_cap(x1.pl, x2.pl), 0)
    x = xy_Atom(x1.tag, mat24.mul_ploop(x1.pl, x2.pl))
    return [delta, x]

rule_yy = rule_xx


def rule_xy(group, word):
    global n_rules
    n_rules += 1
    x, y = word[0], word[1]
    delta = AutPLoopAtom('p', mat24.ploop_cap(x.pl, y.pl), 0)
    comm =  mat24.ploop_comm(x.pl, y.pl) << 12
    x1 = xy_Atom(x.tag, x.pl ^ comm)
    y1 = xy_Atom(y.tag, y.pl ^ comm)
    return [delta, y1, x1]



def rule_tt(group, word):
    global n_rules
    n_rules += 1
    exp = (word[0].exp + word[1].exp) % 3
    return  [tl_Atom(word[0].tag, exp)] if exp else []


def rule_pt(group, word):
    global n_rules
    n_rules += 1
    p, t = word[0], word[1]
    exp = 3 - t.exp  if p.cocode & 0x800 else t.exp
    t1 = tl_Atom('t', exp % 3)
    return [t1, p]

def rule_xt(group, word):
    global n_rules
    n_rules += 1
    x, t = word[0], word[1]
    assert t.exp in [1,2]
    if t.exp == 2:
        pl = mat24.pow_ploop(x.pl, 3)
        x, y = xy_Atom('x', pl), xy_Atom('y', pl)
        return [t, y, x] 
    else:
        y = xy_Atom('y', x.pl)
        return [t, y] 
  

def rule_yt(group, word):
    global n_rules
    n_rules += 1
    y, t = word[0], word[1]
    assert t.exp in [1,2]
    if t.exp == 1:
        pl = mat24.pow_ploop(y.pl, 3)
        x, y = xy_Atom('x', pl), xy_Atom('y', pl)
        return [t, y, x] 
    else:
        x = xy_Atom('x', y.pl)
        return [t, x] 


rule_ll = rule_tt




###########################################################################
# Inversions of atoms in the group N
###########################################################################


def iter_inv_p(group, atom):
    global n_inv
    n_inv += 1
    rep = mat24.inv_autpl(atom.rep)
    yield  AutPLoopAtom(atom.tag, rep) 


def iter_inv_xy(group, atom):
    global n_inv
    n_inv += 1
    res = mat24.pow_ploop(atom.pl, 3)
    yield  xy_Atom(atom.tag, res) 


def iter_inv_tl(group, atom):
    global n_inv
    n_inv += 1
    yield tl_Atom(atom.tag, -atom.exp % 3) 


    
###########################################################################
# The group N 
###########################################################################


ERR_ATOM = "Illegal atom %s in constuctor of class MGroupNWord"

def iter_generators_from_atom(atom):
    atom = int(atom)
    if atom & -0x80000000 == 0:
        tag = atom >> 28
        data = atom & 0xfffffff
        if tag == 0:    return
        elif tag == 1:  yield AutPLoopAtom('p', data, 0)
        elif tag == 2:  yield AutPLoopAtom('p', 0, data)
        elif tag == 3:  yield xy_Atom('x', data)
        elif tag == 4:  yield xy_Atom('y', data)
        elif tag == 5:  yield tl_Atom('t', data)
        elif tag == 6:  yield tl_Atom('l', data)
        elif data == 0: return
        else:
            raise ValueError(ERR_ATOM % hex(atom))
    else:
        raise ValueError(ERR_ATOM % hex(atom))
            
         
    


parse_functions = {  # deprecated!!!!
            'd' : gen_d,
            'p' : gen_p,
            'x' : gen_xy,         
            'y' : gen_xy,
            'z' : gen_z,
            't' : gen_tl,
            'l' : gen_tl,
}
 
  

class MGroupNWord(AutoGroupWord, AbstractMMGroup):
    def __init__(self,  tag = None, atom = None, *args, **kwds):
        self.seq =  []
        self.reduced = 0
        atoms = iter_mm(self.group, tag, atom)
        for a in atoms:
            for generator in iter_generators_from_atom(a):
                self.seq.append(generator)
            
    def _iter_atoms(self):
        for gen in self.seq:
            yield from gen.as_atoms()

    @property
    def mmdata(self):
        return np.fromiter(self._iter_atoms(), dtype = np.uint32)


class MGroupN(AutoGroup):
    """Models the subgroup N_0 of the Monster MM or one of its covers.

    We may also suport a specific generator l of order 3 of MM \\ N_0 
    so that the monster generated by N_0 and l. The group N_0 or its
    covers are modeled correctly according to their definition in [2].
       
    The present implementation does not model any relations containing
    the generator l. The raison d'etre of generator l is to support 
    operations of MM on its 196884-dimensional represention.
    """
    word_type = MGroupNWord
    is_mmgroup = True
    tags = "pxytl"   
    __instance = None

    def __new__(cls):
        if MGroupN.__instance is None:
             MGroupN.__instance = MGroupNWord.__new__(cls)
        return MGroupN.__instance
 
    rules = {
        r"pp": rule_pp,
        r"xx": rule_xx,
        r"yy": rule_xx,
        r"tt": rule_tt,
        r"ll": rule_ll,
        r"xy": rule_xy,
        r"xp": rule_xp,
        r"yp": rule_yp,
        r"pt": rule_pt,
        r"xt": rule_xt,
        r"yt": rule_yt,
        r"x":  rule_x0,
        r"p":  rule_p0,
        r"t":  rule_tl0,
        r"l":  rule_tl0,
    }

    kernel_rules = {
        r"yx": kernel_rule_yx,
        r"y": kernel_rule_y,
    }
    rules.update(kernel_rules)

    inverter = {
            'p' : iter_inv_p,
            'x' : iter_inv_xy,         
            'y' : iter_inv_xy,
            't' : iter_inv_tl,
            'l' : iter_inv_tl,
    }

    group_name = "MN"
    #kernel = 0  # for compatiblity with older version
    



    def __init__(self):
        """Create an instance of subgroup N_0 of the Monster
		
        Here N_0 is the supgroup the monster defined in 
        [Conw85] We follow the conventions in [Seys19] for
        calculating in N_0. We also define an atom (t,e) in N_0
        which maps to the generator xi**e in [Seys19]; but we
        implement no relations for xi apart from xi**3 = 1.		
		        """
        rules = self.rules 
        super(MGroupN, self).__init__(parse_functions, 
            rules, self.inverter)

    def __call__(self, *args, **kwds):
         return MGroupNWord(*args, **kwds)

    def neutral(self):
         return MGroupNWord()
  



StdMGroupN = MGroupN()
MGroupNWord.group = StdMGroupN
load_group_name(StdMGroupN)
add_to_embedded_classes(MGroupNWord)




