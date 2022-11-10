import collections
import re
import warnings
from numbers import Integral
import numpy as np
from random import randint
from functools import partial

try:
    from mmgroup import mat24
    from mmgroup.mat24 import ploop_theta, pow_ploop, MAT24_ORDER
except:
    w = "Extension mmgroup.mat24 not found, package not functional!"
    warnings.warn(w, UserWarning)


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
from mmgroup import Cocode, PLoop, AutPL



###########################################################################
# Standard atoms for the group M
###########################################################################


tag_dict = {
        "d": 0x10000000, 
        "p": 0x20000000, 
        "x": 0x30000000, 
        "y": 0x40000000, 
        "t": 0x50000000, 
        "l": 0x60000000, 
}


STD_Q_ELEMENTS = {
   "v+": 0x200, "v-":0x1000200, 
   "Omega":0x800000, "-Omega":0x1800000,
   "-":0x1000000, "+":0,
   "omega": 0x400, "-omega": 0x1000400, 
}




tags = " dpxytl"


ERR_ATOM_VALUE = "Illegal value of atom in constructor of monster group with tag '%s'"
ERR_ATOM_TYPE = "Illegal type '%s' of atom in constructor of monster group with tag '%s'"
ERR_TAG_VALUE = "Illegal tag '%s' in constructor of monster group"
ERR_TAG_TYPE = "Illegal type '%s' of tag in constructor of monster group"


def std_q_element(tag, r):
    try:
        e = STD_Q_ELEMENTS[r]
    except KeyError:
        raise ValueError(ERR_ATOM_VALUE % tag) 
    if tag == 'q':
        return e
    if tag in 'xyz' and e & 0xfff == 0:
        return e >> 12
    if tag == 'd' and e & 0x1fff000 == 0:
        return e
    raise ValueError(ERR_ATOM_VALUE % tag) 
  



def iter_d(tag, d):
    if isinstance(d, Integral ):
        yield 0x10000000 + (d & 0xfff)
    elif isinstance(d, str):
        cocode = randint(d == 'n', 0xfff) 
        if d == "o":
            cocode |= 1 
        elif d == "e":
            cocode &= ~1
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

def iter_p(tag, perm):
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


def gen_ploop_element(tag, r):
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

def iter_xy(tag, r):
    yield tag_dict[tag] + gen_ploop_element(tag, r)

def iter_z(tag, r ):
    pl = pow_ploop(gen_ploop_element(tag, r), 3)
    yield tag_dict['x'] + pl
    yield tag_dict['y'] + pl

def iter_tl(tag, r):
    if isinstance(r, Integral):
        e = r % 3
    elif isinstance(r, str):
        e = randint(int('n' in r), 2) 
        if  len(r) > 1 or not r in "rn":
            raise ValueError(ERR_ATOM_VALUE % tag)            
    else:
        raise TypeError(ERR_ATOM_TYPE % (type(r), tag))
    if e:
        yield  tag_dict[tag] + e 


###########################################################################
# Extended atoms for the group M
###########################################################################


ERR_Q_x0 = "Atom for tag 'q' must represent element of subgroup Q_x0 of the monster"




def iter_q(tag, r):
    if isinstance(r, Integral):
        e = r & 0x1ffffff
    elif isinstance(r, str):
        if  len(r) == 1 and r in "rn":
            e = randint(r == 'n', 0x1ffffff) 
        elif r in STD_Q_ELEMENTS:
            e = STD_Q_ELEMENTS[r]
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
    yield tag_dict['x'] + d
    yield tag_dict['d'] + delta




ERR_LEECH2 = "Atom for tag 'c' must represent a type-4 vector in Leech lattice mod 2"


def rand_type4_vector():
    c = 0
    while gen_leech2_type(c) != 4:
        c = randint(0, 0xffffff) 
    return c

def iter_c(tag, r):
    if isinstance(r, str):
        if r == 'r':
            c = rand_type4_vector() 
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

def iter_a(tag, a):
    try:
        for x in a:
            yield x & 0xffffffff
    except:
        raise TypeError(ERR_TAG_A)
       
###########################################################################
# Generating random atoms
###########################################################################


ERR_RAND_INTERN = "Internal error in generating random element of monster"

def iter_rand_N_0(in_N_x0, even):
    r"""Return random element of subgroup :math:`N_{0}`

    The function returns a uniform distributed random element
    of the subgroup :math:`N_{x}` of the monster of structure
    :math:`2^{2+11+22}.(M_{24} \times \mbox{Sym}_3)`. The group 
    :math:`N_0` is generated by the generators with tags
    ``x, y, d, p, t``. The function uses the internal random 
    generator of the ``mmgroup`` package.

    If parameter ``in_N_x0`` is nonzero then we compute a random
    element of the subgroup :math:`N_{x0}` of index 3 in :math:`N_0` 
    generated by the generators with tags ``x, y, d, p``.

    If parameter ``even`` is nonzero then we compute a random
    element of the  subgroup :math:`N_{\mbox{even}}` of index 2
    in :math:`N_{x}`  generated by the generators with
    tags ``x, y, d, p, t``, where all generators with tag ``d``
    correspond to even Golay cocode words.

    If both, ``in_N_x0`` and ``even``, are nonzero then we compute
    a random element
    of :math:`N_{xyz0} = N_{\mbox{even}} \cap N_{x0}`.
    """
    buf = np.zeros(10, dtype = np.uint32)
    seed = rand_get_seed()
    length = xsp2co1_rand_word_N_0(buf, in_N_x0, even, seed) 
    if not 0 <= length <= 10:
        raise ValueError(ERR_RAND_INTERN)
    yield from buf[:length]
     

def iter_rand_G_x0():
    r"""Return random element of subgroup :math:`G_{x0}`

    The function returns a uniform distributed random element
    of the subgroup :math:`G_{x0}`. The function uses the
    internal random generator of the ``mmgroup`` package.
    """
    buf = np.zeros(10, dtype = np.uint32)
    seed = rand_get_seed()
    length = xsp2co1_rand_word_G_x0(buf, seed) 
    if not 0 <= length <= 10:
        raise ValueError(ERR_RAND_INTERN)
    yield from buf[:length]



def iter_rand_mm(quality):
    r"""Return a random element of the monster group

    The function returns a random element of the monster group.
    Here ``quality`` means a measure for the quality of the
    randimization process, where a higher value menas that
    the distribution of the elements is closer to uniform.

    If ``quality`` is an integer ``k`` then a product containing
    ``k`` powers of the triality element it generated. Here the
    default value creates an almost uniform distribution.

    In future versions the default value of parameter ``quality``
    may correspond to the generation of a truly uniform
    distribution.
    """
    a = np.zeros(10, dtype = np.uint32)
    seed = rand_get_seed()
    len_a = xsp2co1_rand_word_G_x0(a, seed) 
    if not 0 <= len_a <= 10:
        raise ValueError(ERR_RAND_INTERN)
    yield from a[:len_a]
    for k in range(quality):
        yield 0x50000001 + gen_rng_modp(2, seed)
        c = 0
        while gen_leech2_type(c) != 4:
            c = gen_rng_modp(0x1000000, seed) 
        len_a = gen_leech2_reduce_type4(c, a)
        mm_group_invert_word(a, len_a)
        if not 0 <= len_a <= 6:
            raise ValueError(ERR_RAND_INTERN)
        yield from a[:len_a]



RAND_FUNCTIONS = {
    "M"      : (iter_rand_mm, 18),  
    "G_x0"   : (iter_rand_G_x0,),
    "N_0"    : (iter_rand_N_0, 0, 0), 
    "N_x0"   : (iter_rand_N_0, 1, 0), 
    "N_0_e"  : (iter_rand_N_0, 0, 1), 
    "N_x0_e" : (iter_rand_N_0, 1, 1), 
}

def iter_rand_subgroup(tag, s):
    if isinstance(s, Integral):
        yield from iter_rand_mm(s)
        return
    try:
        ff = RAND_FUNCTIONS[s]
    except KeyError:
        if isinstance(s, str):
            err = "Bad group name '%s' for tag 'r'"
            raise ValueError(err % s)
        else:
            err = "Atom for tag 'r' must be a string"
            raise TypeError(err)
    yield from ff[0](*(ff[1:]))
    
    

         

###########################################################################
# Generating atoms
###########################################################################


def iter_neutral(*args):
    return
    yield 0

gen_tag_dict = {
        "d": iter_d, 
        "p": iter_p, 
        "x": iter_xy, 
        "y" :iter_xy, 
        "z" :iter_z, 
        "t": iter_tl, 
        "l": iter_tl, 
        "q": iter_q, 
        "c": iter_c, 
        "a": iter_a, 
        "r": iter_rand_subgroup,
        "1": iter_neutral
}


def iter_atom(tag = None, number = None):
    """Return list of integers representing element of monster group

    This is the workhorse for method MMGroup.atom().
    See that method for documentation.
    """
    try: 
        gen_function = gen_tag_dict[tag]
    except KeyError:
        if isinstance(tag, str):
            raise ValueError(ERR_TAG_VALUE % tag)
        else:
             raise TypeError(ERR_TAG_TYPE % type(tag))
    yield from gen_function(tag, number)




###########################################################################
# Converting the input of a construtor of a monster element to atoms
###########################################################################


FRAME = re.compile(r"^([A-Za-z][A-Za-z_0-9]*)?\<(.+)\>$") 

EMBEDDED_GROUPS = [
]

EMBEDDED_CLASSES = (
   AutPL, Cocode,
)

TYPES_PARSED_FROM_LIST = EMBEDDED_CLASSES + (
   Integral, str, AbstractGroupWord,
)

ATOM_PARSERS = {
}

def add_to_embedded_classes(cls):
    global EMBEDDED_CLASSES
    if not cls in EMBEDDED_CLASSES:
        EMBEDDED_CLASSES = EMBEDDED_CLASSES + (cls,) 


def element_from_atom(group, tag, *data):
     if tag.isalpha():
         return group(tag, *data)
     else:
         raise ERR_TAG_VALUE(tag)


def iter_parse_mm_string(group, s):
    m = FRAME.match(s)
    group_name, string = (m[1], m[2]) if m else (None, s)
    if not group_name:
        group_name = group
    try:
        atom_dict = ATOM_PARSERS[group_name]
    except KeyError:
        if isinstance(group_name, str):
            err = "Cannot parse string of shape '%s<..>' in monster group" 
            raise ValueError(err % group_name)
        else:
            atom_dict = ATOM_PARSERS["M"] 
    element = eval_atom_expression(string, atom_dict)
    if element != 1:
        yield from element.mmdata



def load_group_name(group, group_name = None):
    global EMBEDDED_GROUPS, ATOM_PARSERS
    assert isinstance(group, AbstractGroup)
    EMBEDDED_GROUPS.append(group)
    atom_dict = AtomDict(partial(element_from_atom, group))
    ATOM_PARSERS[group] = atom_dict
    if group_name:
        assert isinstance(group_name, str) and group_name.isidentifier()
        ATOM_PARSERS[group_name] = atom_dict

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
        yield from iter_atom(tag, atom)
        return
    if tag is None:
        return
    if import_groups_pending:
        import_groups()
    if isinstance(tag, AbstractMMGroupWord):
        yield from tag.mmdata
    elif isinstance(tag, EMBEDDED_CLASSES):
        yield from tag.mmdata
    elif isinstance(tag, str):
        if group is None:
            group = MMGroup
        yield from iter_parse_mm_string(group, tag)
    elif isinstance(tag, list):
        for element in tag:
            if isinstance(element, tuple):
                yield from iter_atom(*element)
            elif isinstance(element, TYPES_PARSED_FROM_LIST):
                yield from iter_mm(element)
            else:
                raise TypeError(ERR_TAG_TYPE % type(element))
    elif isinstance(tag, Integral):
        if tag == 1:
            return
        elif tag == -1 and in_G_x0:
            yield 0x31000000
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
    if element[1]: yield('y', element[1])
    if element[2]: yield('x', element[2])
    if element[3]: yield('d', element[3])
    if element[4]: yield('p', element[4])
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
                yield ("tl"[tag-5], v)
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

def iter_strings_from_N_x0_element(element):
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
                yield from iter_strings_from_N_x0_element(N_x0_element)
                yield "%s_%d" % ("tl"[tag-5], v)
        else:
            if abort_if_error:
                raise ValueError(ERR_ILLEGAL_ATOM % hex(a))
            yield "<Bad atom %s>" % hex(a)
    yield from iter_strings_from_N_x0_element(N_x0_element)


def print_mm_string(tag = None, data = None, group = None):
    strings = iter_strings_from_atoms(iter_mm(group, tag, data))
    s = "*".join(strings)
    print(s if len(s) else "1")



###########################################################################
# Converting tuples to strings
###########################################################################


TUPLE_FMT_DICT = dict(zip("dpxytl", [ihex, str, ihex, ihex, str, str]))

def iter_strings_from_tuples(tuples):
     for tag, value in tuples:
          try:
              fmt = TUPLE_FMT_DICT[tag]
          except KeyError:
              if isinstance(tag, str):
                  err = "Illagal tag in tuple of monster group atoms"
                  raise ValueError(err)
              else:
                  err = "Tag in tuple of monster group atoms must be a string"
                  raise TypeError(err)
          yield "%s_%s" % (tag, fmt(value))

