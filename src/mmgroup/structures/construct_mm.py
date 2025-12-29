"""Construct an element of the Monster group

The main function ``iter_mm()`` in this module constructs an element
of the Monster group from data structures as specifies in section
**The Monster group** of the **API reference**. This function
yields the entries of a numpy array of type ``np.uint32`` containing 
internal represntation of the construted element.

Functions or data structures in this module starting with an
underscore (``_``) should not be imported by other modules.
"""

import collections
import re
import warnings
from numbers import Integral
import numpy as np
from random import randint
from functools import partial

try:
    import mmgroup
    from mmgroup import mat24
    from mmgroup.mat24 import ploop_theta, pow_ploop, MAT24_ORDER
except:
    w = "Extension mmgroup.mat24 not found, package not functional!"
    mmgroup._warn(w)


from mmgroup.generators import rand_get_seed, gen_leech2_type
from mmgroup.generators import gen_rng_modp
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import mm_group_invert_word
from mmgroup.generators import mm_group_n_reduce_element
from mmgroup.generators import mm_group_n_clear
from mmgroup.generators import mm_group_n_mul_atom
from mmgroup.clifford12 import xsp2co1_rand_word_G_x0
from mmgroup.clifford12 import xsp2co1_rand_word_N_0

from mmgroup.structures.abstract_group import AbstractGroup
from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.parse_atoms import ihex, TaggedAtom, AtomDict
from mmgroup.structures.parse_atoms import  eval_atom_expression  
from mmgroup.structures.autpl import autpl_from_obj, rand_perm_num
from mmgroup.structures.cocode import Cocode
from mmgroup.structures.ploop import PLoop
from mmgroup.structures.autpl import AutPL
from mmgroup.structures.random_mm import iter_random_mm



###########################################################################
# Standard atoms for the group M
###########################################################################


_tag_dict = {
        "d": 0x10000000, 
        "p": 0x20000000, 
        "x": 0x30000000, 
        "y": 0x40000000, 
        "t": 0x50000000, 
        "l": 0x60000000, 
}


_STD_Q_ELEMENTS = {
   "v+": 0x200, "v-":0x1000200, 
   "Omega":0x800000, "-Omega":0x1800000,
   "-":0x1000000, "+":0,
   "omega": 0x400, "-omega": 0x1000400, 
}





ERR_ATOM_VALUE = "Illegal value of atom in constructor of monster group with tag '%s'"
ERR_ATOM_TYPE = "Illegal type '%s' of atom in constructor of monster group with tag '%s'"
ERR_TAG_VALUE = "Illegal tag '%s' in constructor of monster group"
ERR_TAG_TYPE = "Illegal type '%s' of tag in constructor of monster group"


def std_q_element(tag, r):
    try:
        e = _STD_Q_ELEMENTS[r]
    except KeyError:
        raise ValueError(ERR_ATOM_VALUE % tag) 
    if tag == 'q':
        return e
    if tag in 'xyz' and e & 0xfff == 0:
        return e >> 12
    if tag == 'd' and e & 0x1fff000 == 0:
        return e
    raise ValueError(ERR_ATOM_VALUE % tag) 
  



def _iter_d(tag, d):
    if isinstance(d, Integral ):
        yield 0x10000000 + (d & 0xfff)
    elif isinstance(d, str):
        cocode = randint(d == 'n', 0xfff) 
        if d == "o":
            cocode |= 0x800 
        elif d == "e":
            cocode &= ~0x800
        if  len(d) == 1 and d in "rnoe":
            yield 0x10000000 + (cocode & 0xfff)
        else:
            yield 0x10000000 + std_q_element(tag, d)            
    else:
        try:
            cocode = Cocode(d).cocode
            yield  0x10000000 + (cocode & 0xfff)
        except:
            raise TypeError(ERR_ATOM_TYPE % (type(d), 'd'))

def _iter_p(tag, perm):
    if isinstance(perm, Integral):
        if not 0 <= perm < MAT24_ORDER:
            err = "Tag 'p': bad permutation number for Mathieu group M24"
            raise ValueError(err)
        yield 0x20000000 + perm
    elif isinstance(perm, str):
        try:
            perm_num = rand_perm_num(perm)
        except ValueError:
            if perm == 'n':
                perm_num = randint(1, MAT24_ORDER-1)
            else:
                raise ValueError(ERR_ATOM_VALUE % 'p')
        yield 0x20000000 + perm_num
    elif isinstance(perm, AutPL):
        yield 0x10000000 + perm._cocode
        yield 0x20000000 + perm._perm_num
    else:
        try:
            cocode, perm_num = autpl_from_obj(0, perm)
            yield  0x10000000 + (cocode & 0xfff)
            yield  0x20000000 + perm_num
        except:
            raise TypeError(ERR_ATOM_TYPE % (type(perm), tag))


def _gen_ploop_element(tag, r):
    if isinstance(r, Integral):
        return  r & 0x1fff
    elif isinstance(r, str):
        r_num = randint(r == 'n', 0x1fff)
        if  len(r) == 1 and r in "rn":
            return r_num
        return std_q_element(tag, r)            
    else:
        try:
            return PLoop(r).ord
        except:
            raise TypeError(ERR_ATOM_TYPE % (type(r), tag))

def _iter_xy(tag, r):
    yield _tag_dict[tag] + _gen_ploop_element(tag, r)

def _iter_z(tag, r ):
    pl = pow_ploop(_gen_ploop_element(tag, r), 3)
    yield _tag_dict['x'] + pl
    yield _tag_dict['y'] + pl

def _iter_tl(tag, r):
    if isinstance(r, Integral):
        e = r % 3
    elif isinstance(r, str):
        e = randint(int('n' in r), 2) 
        if  len(r) > 1 or not r in "rn":
            raise ValueError(ERR_ATOM_VALUE % tag)            
    else:
        raise TypeError(ERR_ATOM_TYPE % (type(r), tag))
    if e:
        yield _tag_dict[tag] + e 


###########################################################################
# Extended atoms for the group M
###########################################################################


ERR_Q_x0 = "Atom for tag 'q' must represent element of subgroup Q_x0 of the monster"




def _iter_q(tag, r):
    if isinstance(r, Integral):
        e = r & 0x1ffffff
    elif isinstance(r, str):
        if  len(r) == 1 and r in "rn":
            e = randint(r == 'n', 0x1ffffff) 
        elif r in _STD_Q_ELEMENTS:
            e = _STD_Q_ELEMENTS[r]
        else:
            raise ValueError(ERR_ATOM_VALUE % tag) 
    else:
        try:
            e = r.as_Q_x0_atom()
            assert isinstance(e, Integral)
        except:
            raise TypeError(ERR_Q_x0)
    d = (e >> 12) & 0x1fff
    delta = (e ^ ploop_theta(d)) & 0xfff
    yield _tag_dict['x'] + d
    yield _tag_dict['d'] + delta




ERR_LEECH2 = "Atom for tag 'c' must represent a type-4 vector in Leech lattice mod 2"


def _rand_type4_vector():
    c = 0
    while gen_leech2_type(c) != 4:
        c = randint(0, 0xffffff) 
    return c

def _iter_c(tag, r):
    if isinstance(r, str):
        if r == 'r':
            c = _rand_type4_vector() 
        else:
            raise ValueError(ERR_LEECH2)
    elif isinstance(r, Integral):
        c = r & 0xffffff
    else:
        try:
            c = r.as_Q_x0_atom()
            assert isinstance(c, Integral)
        except:
            raise TypeError(ERR_LEECH2)
    if gen_leech2_type(c) != 4:
        raise ValueError(ERR_LEECH2)
    a = np.zeros(6, dtype = np.uint32)
    len_a = gen_leech2_reduce_type4(c, a)
    mm_group_invert_word(a, len_a)
    yield from a[:len_a]
     



###########################################################################
# Atoms generated from arrays of 32-bit integers
###########################################################################


ERR_TAG_A = "Atom for tag 'a' must be a list of unsigend 32-bit integers"

def _iter_a(tag, a):
    try:
        for x in a:
            yield int(x) & 0xffffffff
    except:
        raise TypeError(ERR_TAG_A)
       
###########################################################################
# Generating random atoms
###########################################################################


def _iter_rand_subgroup(tag, s):
    yield from iter_random_mm(s)
    
    

         

###########################################################################
# Generating atoms
###########################################################################


def _iter_neutral(*args):
    return
    yield 0

_gen_tag_dict = {
        "d": _iter_d, 
        "p": _iter_p, 
        "x": _iter_xy, 
        "y" :_iter_xy, 
        "z" :_iter_z, 
        "t": _iter_tl, 
        "l": _iter_tl, 
        "q": _iter_q, 
        "c": _iter_c, 
        "a": _iter_a, 
        "r": _iter_rand_subgroup,
        "1": _iter_neutral
}


def _iter_atom(tag = None, number = None):
    """Return list of integers representing element of monster group

    This is the workhorse for method MMGroup.atom().
    See that method for documentation.
    """
    try: 
        gen_function = _gen_tag_dict[tag]
    except KeyError:
        if isinstance(tag, str):
            raise ValueError(ERR_TAG_VALUE % tag)
        else:
             raise TypeError(ERR_TAG_TYPE % type(tag))
    yield from gen_function(tag, number)




###########################################################################
# Converting the input of a construtor of a monster element to atoms
###########################################################################


_FRAME = re.compile(r"^([A-Za-z][A-Za-z_0-9]*)?\<(.+)\>$") 

_EMBEDDED_GROUPS = [
]

_EMBEDDED_CLASSES = (
   AutPL, Cocode,
)

_TYPES_PARSED_FROM_LIST = _EMBEDDED_CLASSES + (
   Integral, str, AbstractGroupWord,
)

_ATOM_PARSERS = {
}

def add_to_embedded_classes(cls):
    global _EMBEDDED_CLASSES
    if not cls in _EMBEDDED_CLASSES:
        _EMBEDDED_CLASSES = _EMBEDDED_CLASSES + (cls,) 


def _element_from_atom(group, tag, *data):
     if tag.isalpha():
         return group(tag, *data)
     else:
         raise ERR_TAG_VALUE(tag)


def _iter_parse_mm_string(group, s):
    m = _FRAME.match(s)
    group_name, string = (m[1], m[2]) if m else (None, s)
    if not group_name:
        group_name = group
    try:
        atom_dict = _ATOM_PARSERS[group_name]
    except KeyError:
        if isinstance(group_name, str):
            err = "Cannot parse string of shape '%s<..>' in monster group" 
            raise ValueError(err % group_name)
        else:
            atom_dict = _ATOM_PARSERS["M"] 
    element = eval_atom_expression(string, atom_dict)
    if element != 1:
        yield from element.mmdata



def load_group_name(group, group_name = None):
    global _EMBEDDED_GROUPS, _ATOM_PARSERS
    assert isinstance(group, AbstractGroup)
    _EMBEDDED_GROUPS.append(group)
    atom_dict = AtomDict(partial(_element_from_atom, group))
    _ATOM_PARSERS[group] = atom_dict
    if group_name:
        assert isinstance(group_name, str) and group_name.isidentifier()
        _ATOM_PARSERS[group_name] = atom_dict

import_groups_pending = True

def import_groups():
    global import_groups_pending
    global  MM, MMGroup, AbstractMMGroupWord
    from mmgroup import MM, MMGroup
    from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord
    #load_group_name(MM().group, "M")  # TODO: Fixme!!!!
    import_groups_pending = False




def iter_mm(group, tag = None, atom = None, in_G_x0 = False):
    if isinstance(tag, str) and len(tag) == 1:
        if atom is None and tag == 'r':
            atom = 'G_x0' if in_G_x0 else 'M' 
        yield from _iter_atom(tag, atom)
        return
    if tag is None:
        return
    if import_groups_pending:
        import_groups()
    if isinstance(tag, AbstractMMGroupWord):
        yield from tag.mmdata
    elif isinstance(tag, _EMBEDDED_CLASSES):
        yield from tag.mmdata
    elif isinstance(tag, str):
        if group is None:
            group = MMGroup
        yield from _iter_parse_mm_string(group, tag)
    elif isinstance(tag, list):
        for element in tag:
            if isinstance(element, tuple):
                yield from _iter_atom(*element)
            elif isinstance(element, _TYPES_PARSED_FROM_LIST):
                yield from iter_mm(element)
            else:
                raise TypeError(ERR_TAG_TYPE % type(element))
    elif isinstance(tag, Integral):
        if tag == 1:
            return
        elif tag == -1 and in_G_x0:
            yield 0x30001000
            return
        err = "Cannot convert integer to monster group element"
        raise ValueError(err)
    else:
        raise TypeError(ERR_TAG_TYPE % type(tag))



def print_mm_atoms(tag = None, data = None, group = None):
    print([hex(x) for x in iter_mm(group, tag, data)])



###########################################################################
# Converting atoms to tuples
###########################################################################


ERR_ILLEGAL_ATOM = "Illegal atom %s in monster group element"

def iter_tuples_from_N_x0_element(element):
    mm_group_n_reduce_element(element)
    if element[1]: yield('y', int(element[1]))
    if element[2]: yield('x', int(element[2]))
    if element[3]: yield('d', int(element[3]))
    if element[4]: yield('p', int(element[4]))
    mm_group_n_clear(element)


def iter_tuples_from_atoms(atoms):
    N_x0_element = np.zeros(5, dtype = np.uint32)
    for a in atoms:
        a &= 0xffffffff
        t = a & 0x7fffffff
        if t < 0x50000000:
            mm_group_n_mul_atom(N_x0_element, a)
        elif a & 0xfffffff == 0:
            continue
        elif t < 0x70000000:
            tag, v = t >> 28, a & 0xfffffff
            v = (-v if a & 0x80000000 else v) % 3
            if v:
                yield from iter_tuples_from_N_x0_element(N_x0_element)
                yield ("tl"[tag-5], int(v))
        else:
            raise ValueError(ERR_ILLEGAL_ATOM % hex(a))
    yield from iter_tuples_from_N_x0_element(N_x0_element)


def print_mm_tuples(tag = None, data = None, group = None):
    tuples = iter_tuples_from_atoms(iter_mm(group, tag, data))
    print([hex(x) for x in tuples])


###########################################################################
# Converting atoms to strings
###########################################################################


ERR_ILLEGAL_ATOM = "Illegal atom %s in monster group element"

def _iter_strings_from_N_x0_element(element):
    mm_group_n_reduce_element(element)
    if element[1]: yield "y_%s" % ihex(element[1])
    if element[2]: yield "x_%s" % ihex(element[2])
    if element[3]: yield "d_%s" % ihex(element[3])
    if element[4]: yield "p_%d" % element[4]
    mm_group_n_clear(element)


def iter_strings_from_atoms(atoms, abort_if_error = True):
    N_x0_element = np.zeros(5, dtype = np.uint32)
    for a in atoms:
        a &= 0xffffffff
        t = a & 0x7fffffff
        if t < 0x50000000:
            mm_group_n_mul_atom(N_x0_element, a)
        elif a & 0xfffffff == 0:
            continue
        elif t < 0x70000000:
            tag, v = t >> 28, a & 0xfffffff
            v = (-v if a & 0x80000000 else v) % 3
            if v:
                yield from _iter_strings_from_N_x0_element(N_x0_element)
                yield "%s_%d" % ("tl"[tag-5], v)
        else:
            if abort_if_error:
                raise ValueError(ERR_ILLEGAL_ATOM % hex(a))
            yield "<Bad atom %s>" % hex(a)
    yield from _iter_strings_from_N_x0_element(N_x0_element)


def print_mm_string(tag = None, data = None, group = None):
    strings = iter_strings_from_atoms(iter_mm(group, tag, data))
    s = "*".join(strings)
    print(s if len(s) else "1")




