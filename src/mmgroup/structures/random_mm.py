"""Construct a random element of the Monster group

The main function ``iter_reandom_mm()`` in this module constructs an 
random element of the Monster group from data structures as specified 
in section **The Monster group** of the **API reference**. This function
yields the entries of a numpy array of type ``np.uint32`` containing 
internal represntation of the constructed element.
"""

import collections
import re
import warnings
from numbers import Integral
import numpy as np
from random import randint
from functools import partial

try:
    from mmgroup import mat24
    from mmgroup.mat24 import MAT24_ORDER
    from mmgroup.mat24 import m24num_rand_local
    from mmgroup.mat24 import m24num_rand_adjust_xy
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
from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2


"""Bit semantics

# Fix x_{-1}  and the set  \{x_{\pm \Omega}\}
# and the following partition of the set \tilde{\Omega}
Bit  0: MAT24_RAND_2 # fixes \{2, 3\}  
Bit  1: MAT24_RAND_o # fixes \{0, \ldots,7 \}
Bit  2: MAT24_RAND_t # fixes \{\{8i,\ldots,8i+7\} \mid i < 3 \}
Bit  3: MAT24_RAND_s # fixes \{\{4i,\ldots,4i+3\} \mid i < 6 \}
Bit  4: MAT24_RAND_l # fixes \{\{2i, 2i+1\} \mid  4 \leq i < 12 \}
Bit  5: MAT24_RAND_3 # fixes \{1, 2, 3\}  

# Fix certain involutions
Bit  6 : # fix x_{-1}
Bit  7 : # fix set \{x_{-1}, x_{\pm \Omega}\}

Bit  8 : # fix set \{x_{-1}, x_{\pm \Omega}\} and parity of cocode word
Bit  9 : # Subgroup AutPL
Bit 10 : # Subgroup Q_x0
"""

_NO       = 0x20000000
_NOT_FULL = 0x40000000


SUBGOUP_MAP = {
   'M':      0,   
   'G_x0':   0x40,  'G_1': 0x40,
   'N_0':    0x80,  'G_2': 0x80,
   'G_3':    0x8      | _NO,
   'G_5t':   0x4      | _NO,
   'G_5l':   0x10     | _NO,
   'G_10':   0x2      | _NO,
   'B':      0x1      | _NOT_FULL,
   '2E_6':   0x20     | _NO,
   'H+'  :   0x41,
   'N_0_e':  0x180,
   'N_x0':   0xc0,
   'N_x0_e': 0x1c0,
   'Q_x0':   0x440    | _NO,
   'AutPL':  0x2C0    | _NO,
   'quick':  0     # for future optimizations
}



def _parse_group_description(s_in):
    parts = [t.strip() for t in s_in.split('&')]
    count = 0
    flags = 0
    for s in parts:
        try:
            mask = SUBGOUP_MAP[s]
            flags |= mask
            count += bool(mask)
        except KeyError:
            err = "Unknown subgroup description '%s'"
            raise ValueError(err % s)
        if flags & _NO:
            err = "Subgroup '%s' of the Monster not supported"
            raise ValueError(err % s) 
    if (flags & _NOT_FULL) and (count > 1):
        err = "Subgroup '%s' not supported"
        raise ValueError(err % s_in) 
    return flags

def _random_tag_pi(flags):
    if flags & 0x400 == 0:
        u_rand = randint(0, MAT24_ORDER - 1)
        pi = m24num_rand_local(flags, u_rand)
        return 0x20000000 + pi

def _iter_tags_yxdp(flags):
    # tag y
    if flags & 0x600 == 0:
        u_rand = randint(0, 0x1fff)
        y = m24num_rand_adjust_xy(flags, u_rand)
        yield 0x40000000 + y
    # tag x
    if flags & 0x200 == 0:
        u_rand = randint(0, 0x1fff)
        x = m24num_rand_adjust_xy(flags, u_rand)
        yield 0x30000000 + x
    # tag d
    d = randint(0, 0xfff)
    if flags & 0x180 == 0x180:
        d &= 0x7ff
    yield 0x10000000 + d
    # tag p
    if flags & 0x400 == 0:
        u_rand = randint(0, MAT24_ORDER - 1)
        pi = m24num_rand_local(flags, u_rand)
        yield 0x20000000 + pi


BETA = 0x200  # Standard short vector BETA in Leech lattice mod 2

def _rand_Co_2_coset_No():
    r"""Return number of random coset of Co_2 / Co_2 \cap N_x0

    If this function returns the number c then 
    MM('c', c) is a representative of a random coset of
    Co_2 / (Co_2 \cap N_x0)
    """
    while True:
        # generate short vector v2 in Leech lattice mod 2
        ve = randint(300, 98579)  
        vs = mm_aux_index_extern_to_sparse(ve)
        v2 = mm_aux_index_sparse_to_leech2(vs) 
        # Check if v2 is orthogonal to BETA in the real Leech lattice
        v4 = v2 ^ BETA
        if gen_leech2_type(v4) == 4:
            # Return v2 + BETA if this is the case
            return v4 



def _iter_coset_G_x0(flags):
    if flags & 0x780:
        # Generator \xi is not in the subgroup
        return
    relevant_flags = flags & 0x3f
    if relevant_flags & 0x8:
        # Then we deal with in G_3' = G_3 \cap G_x0
        # Here we have |G_3' / G_3' \cap N_x0| = 3
        e = randint(0, 2)
        yield 0x60000000 + e # Append power of \xi
    elif relevant_flags & ~1 == 0:
        # Then we deal with H = 2.B \cap G_x0 or with H = G_x0
        # Generate suitable type-4 vector in the Leech lattice mod 2
        if relevant_flags:
            # Generate a random type-4 vector c that is orthogonal
            # to the standard type-2 vector in the Leech lattice
            c = _rand_Co_2_coset_No()
        else:
            # Generate a uniform random type-4 vector
            c = 0
            while gen_leech2_type(c) != 4:
                c = randint(0, 0xffffff)
        # Append a word in G_x0 that maps Omega to c
        a = np.zeros(6, dtype = np.uint32)
        len_a = gen_leech2_reduce_type4(c, a)
        mm_group_invert_word(a, len_a)
        yield from a[:len_a]
    else:
        # For the other cases we don't have a really good strategy
        for i in range(4):
            e = randint(1,2)
            yield 0x60000000 + e         # Append nonzero power of \xi
            yield _random_tag_pi(flags)  # Append tag pi generator

# Coset representatives of G_3 / G_3 \cap G_x0. 
# A representative (e,f) means \tau**e * \xi**f.        
G3_COSETS = [(0,0), (1,0), (2,0), (1,1), (1,2), (2,2), (2,2)] 

def _iter_rand_mm_(flags, n_rounds = 0):
    # A subgroup H of the monster is given by argument ``flags``
    # Generate an element of H_0 =  N_x0 \cap H 
    yield from _iter_tags_yxdp(flags)
    # Put H_1 = <<\xi> \cap H_1, H_0>. Append a coset representative
    # of  H_1 / H_0  in order to extend H_0 to H_1 
    yield from _iter_coset_G_x0(flags)
    # Put H_2 = <<\tau> \cap H_2, H_1>. Append a coset representative
    # of  H_2 / H_1 in order to extend H_1 to H_2 
    if flags & 0x640:
        # Then the triality element \tau is not in the group
        return
    elif flags & 0x188:
        # Then appending a fixed number of generators suffices
        # Append power of \tau
        yield 0x50000000 + randint(0,2)
        if flags & 8:
            # We are in G_3 and there are 3 * 7 cosets
            coset_data = G3_COSETS[randint(0, 6)]
            yield 0x50000000 + coset_data[0] # \tau
            yield 0x60000000 + coset_data[1] # \xi
    else:
        # Standard strategy for large subgroups of the Monster
        if not n_rounds or n_rounds <= 0:
            n_rounds = 5 if flags & 0x3f else 8
        for i in range(n_rounds):
            yield 0x50000000 + randint(1,3)  # \tau
            yield _random_tag_pi(flags)
            yield from _iter_coset_G_x0(flags)




def iter_random_mm(s):
    if isinstance(s, Integral):
        yield from _iter_rand_mm_(0, s)
    else:
        flags = _parse_group_description(s)
        yield from _iter_rand_mm_(flags)



##################################################################
# Deprecated stuff
##################################################################


ERR_RAND_INTERN = "Internal error in generating random element of monster"

def _iter_rand_N_0(in_N_x0, even):
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
     

def _iter_rand_G_x0():
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



def _iter_rand_mm(quality):
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



_RAND_FUNCTIONS = {
    "M"      : (_iter_rand_mm, 18),  
    "G_x0"   : (_iter_rand_G_x0,),
    "N_0"    : (_iter_rand_N_0, 0, 0), 
    "N_x0"   : (_iter_rand_N_0, 1, 0), 
    "N_0_e"  : (_iter_rand_N_0, 0, 1), 
    "N_x0_e" : (_iter_rand_N_0, 1, 1), 
}



def ___iter_random_mm(s):
    if isinstance(s, Integral):
        yield from _iter_rand_mm(s)
        return
    try:
        ff = _RAND_FUNCTIONS[s]
    except KeyError:
        if isinstance(s, str):
            err = "Bad group name '%s' for tag 'r'"
            raise ValueError(err % s)
        else:
            err = "Atom for tag 'r' must be a string"
            raise TypeError(err)
    yield from ff[0](*(ff[1:]))
    
