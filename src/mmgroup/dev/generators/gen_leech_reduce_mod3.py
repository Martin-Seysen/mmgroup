"""Map vector in Leech lattice mod 3 to a representative of its type 

This is yet a stub!

"""

import sys
import math
from collections import defaultdict
from numbers import Integral

import numpy as np

import_pending = True

def import_all():
    global import_pending 
    if import_pending:
        global mat24, gen_leech3_op_vector_word
        if __name__ == "__main__":
            sys.path.append(r"../../..")
        try:
            from mmgroup import mat24
            from mmgroup.generators import gen_leech3_op_vector_word
            import_pending = False
        except (ImportError, ModuleNotFoundError, AssertionError):
            from mmgroup.dev.mat24.mat24_ref import  Mat24
            mat24 = Mat24    
            import_pending = False
        if __name__ == "__main__":
            sys.path.pop()


###############################################################
#  Compute intermedate results for debugging or testing
###############################################################


def str_vector3(a):
    """Convert vector in Leech lattice mod 3 to string"""
    supp = (a ^ (a >> 24)) & 0xffffff
    neg = a & supp      
    return "Leech3_%06xh_%06xh" % (supp, neg)

def weight_vector3(a):
    """Return weight of support of vector in Leech lattice mod 3"""
    supp = (a ^ (a >> 24)) & 0xffffff
    return mat24.bw24(supp)

def reduce_vector3(a):
    """Change vector in Leech lattice mod 3 to a standard form"""
    ovfl = (a & (a >> 24)) & 0xffffff
    return (a ^ ovfl ^ (ovfl << 24)) & 0xffffffffffff 

def mul_vector3(a, g):
    """Multiply vector in Leech lattice mod 3 by a group element

    """ 
    a1 = gen_leech3_op_vector_word(a, g, len(g))
    return reduce_vector3(a1)



###############################################################
# Find a good tetrad for reduction
###############################################################


def augment_bitvector(v, source, weight):
    """Add bits taken from ``source`` to bit vector ``v``

    The function removes bit from bit vector ``source`` and adds them
    to bit vector ``v`` until ``source`` is zero or ``v`` has bit
    weight ``weight``. It returns the augmented bit vector ``v``.
    """
    source &= 0xffffff & ~v
    w =  mat24.bw24(v)
    while w < weight and source:
        b0 = source & -source
        v |= b0
        source &= ~b0
        w += 1
    return v

def neg_points_hexadecad(gc, neg, sub):
    """Deal with subset of 'negative points' in a hexadecad

    Let ``gc`` be a bit vector that is a hexadecad (i.e. has weight 16)
    and that is in the Golay gode.  Let ``neg`` be a subset of bits 
    of ``gc`` that are considered as 'negative'. We assume that we may 
    change the signs of the bits in ``gc`` by flipping the signs of all 
    bits of a Golay code word.
    Here we ignore the signs of the bits outside of ``gc``. Then it can
    be shown that we can reduce the bit weight  of ``neg`` to 2 by a
    sequence of such sign flips (ignoring the bits outside ``gc``). 

    Let ``sub`` be a subset of ``gc`` such that the sigs of the bits
    in  ``sub`` can also be ignored. If sub is nonzero, we may further
    reduce bit weight of ``neg`` to 1 by sign flips as above, provided
    that we also ignore the signs of the bits in ``sub``.
    
    The function returns a modified subset ``neg`` of ``gc & ~sub``
    (obtained from input ``neg``  by sign flips as described above)
    such that the bit weight out the output ``neg`` is minimal.  
    """ 
    sub &= gc
    o = ~gc & 0xffffff
    neg = mat24.syndrome(neg, mat24.lsbit24(o)) & gc & ~sub
    w = mat24.bw24(neg)
    if w < 2 or (w == 2 and sub == 0):
        return neg
    pool = augment_bitvector(neg, sub, 3)
    return mat24.intersect_octad_tetrad(o, pool) & gc & ~pool &~sub

def cohexad(dodecad, duad, included):
    """Compute hexad from dodecad an duad.

    If ``dodecad`` is a dodecad and ``duad`` is a duad (i.e. a set
    of weight 2) disjoint to the dodecad then the complement of the
    dodecad contains precisely to (disjoint) hexads that complete
    the duad to a octad. The function returns one of these hexads.
    If  ``include`` is a subset of the complement of the dodecad of
    weight 1 then then hexad containing that subset is returned. 
    """
    hexad = mat24.cocode_as_subdodecad(
        mat24.vect_to_cocode(duad & 0xffffff),
        mat24.vect_to_gcode(dodecad & 0xffffff),
        0)
    if included & hexad == 0:
        hexad ^= dodecad
    return hexad

def bitvector_syndromes(v):
    """Compute Golay code syndromes of a vector and their weights

    Given a bit vector vector ``v``, the function returns a triple
    ``(w_v, w_syn, syn_list)``. Here ``w_v`` is the bit weight of
    vector ``v``, and ``w_syn`` is the common bit weight of all
    Golay code syndromes computed. ``syn_list`` is a list of pairs
    ``(w_add, syn)``. Any such pair denotes a Golay code syndrome
    ``syn``, such that ``syn & v`` has bit weight ``w_add``.

    Note that ``gc = v ^ syn`` is a Golay code word. If ``gc`` is
    a docecad (i.e. has bit weight 12), and there is another
    syndrome ``syn1`` such that ``v ^ syn1`` is not a dodecad,   
    then we make sure that ``v ^ syn_list[0][1]`` is not a docecd.  
    """
    w_v = mat24.bw24(v)
    syn_l = [(mat24.bw24(x & v), x) for x in mat24.all_syndromes(v)]
    w_syn = mat24.bw24(syn_l[0][1])
    w_add_0 = syn_l[0][0]
    w_gc = w_v + w_syn - 2 * w_add_0
    if w_gc == 12:                 # We don't like dodecads too much
        for i, (w_add, _) in enumerate(syn_l):
            if w_add != w_add_0 : 
                syn_l[i], syn_l[0] = syn_l[0], syn_l[i]
                break
    return w_v, w_syn, syn_l 

def reduce_even(v, w_syn, syn_list):
    """Do vector in Leech lattice mod3 with Hamming distance 0 or 4

    Let ``w_syn, syn_list`` be as obtained from applying function
    ``bitvector_syndromes`` to the bit vector ``v``. If ``w_syn`` is
    0 or 4 then the function returns a tetrad suitable for reducing
    vector ``v`` as described in function ``find_tetrad_leech_mod3``. 
    Otherwise the function returns 0 meaning that not tetrad has
    been found.
    """ 
    if w_syn == 4:
        w, tet = 0, 0
        for w0, syn in syn_list:
            if w0 == 4:
                return syn 
            if w < 4:
                tet |= syn & v
                w += w0
        return tet if w == 4 else 0
    # _, syn = syn_list[0]
    if w_syn == 0:
        return augment_bitvector(0, v, 4)
    # one might make it better here!
    return 0


def find_tetrad_16_large(gc, sub):
    """Do vector in Leech lattice mod3 with support close to a hexadecad

    Let ``gc`` be a bit vector that is a hexadecad (i.e. has weight 16)
    and that is in the Golay gode. The 16 bits of ``gc`` have a natural
    structure as an affine 4-dimensional svace ofer ``GF(2)`` Let 
    ``sub`` be a subset of ``gc`` of weight at most 4. If ``sub`` can
    be completed to an affine plane in ``gc`` then the function
    returns the four points of that affine plane as a tetrad. Otherwise
    the function returns a tetrad corrsponding to an affine plane in
    ``gc`` such that each plane in ``gc`` parallel to that tetrad
    intersects with ``sub`` in 0 or 2 points.  
    """ 
    #print("plane", hex(gc), hex(sub))
    if mat24.bw24(sub) == 3:
        tet = mat24.intersect_octad_tetrad(gc, sub) & gc
        return tet
    p0 = sub & -sub
    remain = sub & ~p0
    plane = mat24.intersect_octad_tetrad(gc, remain) & gc
    if sub == plane:
        return plane
    p1 = plane & ~remain
    assert mat24.bw24(p1) == 1
    line = p0 | p1 
    #print("plane", hex(gc), hex(sub), hex(p0), hex(p1), hex(remain))
    while remain:
        point = remain & -remain
        remain &= ~point
        tet = mat24.intersect_octad_tetrad(gc, point | line) & gc
        #print(hex(point), hex(tet))
        if mat24.bw24(tet & sub) == 2:
            return tet
    return 0

    
def find_tetrad_leech_mod3(a):
    r"""Find tetrad suitable for reduction in Leech lattice mod 3

    Let ``a`` be a vector in the Leech lattice mod 3 given in
    Leech lattice mod 3 encoding. Let the ``v`` be the support 
    of ``a``, i.e. the bit vector containing the the nonzero
    entries of ``a``.

    The function returns a tetrad ``t`` that may be used by
    function ``reduce_tetrad_leech_mod3`` to find an
    element ``g`` of the group :math:`\mbox{Co}_1` such that
    the support of ``a * g`` has bit weight less than the weight
    of ``v``. The function returns 0 if no such tetrad can be found.

    This procedure reduces the bit weight of ``v`` by

    - at least 12 if ``v`` has bit weight 20, 22, 23, or 24

    - at least 6 if ``v`` has bit weight 19 or more

    - at least 3 except in the following cases:

      - ``v`` has bit weight 5 or less than 4

      - ``v`` is an umbral hexad, heptad, nonad or undecad

      - ``v`` is a transversal octad

      - ``v`` is a special nonad and the unique octad contained
              in ``v`` contains an odd number of entries -1. 

    For characterizing the bit vector ``v`` we use the termminology
    in [CS99], Ch. 10.2.6.
    """
    if import_pending:
        import_all()
    v = (a ^ (a >> 24)) & 0xffffff # non-zero entries of vector
    neg = a & v                    # 'negative' entries of v
    w_v = mat24.bw24(v)
    if w_v >= 20:
        outside = (0xffffff & ~v)
        outside_bit = mat24.lsbit24(outside | 0x800000) 
        neg = mat24.syndrome(neg, outside_bit)
        tet = augment_bitvector(outside, neg, 4)
        return augment_bitvector(tet, 0xffffff, 4)
    if (w_v <= 5):
        return v if w_v == 4 else 0
    w_v, w_syn, syn_list = bitvector_syndromes(v)
    w_add, syn = syn_list[0]
    gc = v ^ syn                    # gc is a Golay code word
    w_gc = w_v + w_syn - 2 * w_add  # bit weight of gc
    add = syn & v                   # v = (gc | add) & ~sub, 
    sub = syn & ~v                  # syn = gc ^ v = add | sub
    w_sub = w_syn - w_add           # bit weight of sub
    if w_gc == 16 and w_sub < 4:
        # The case w_add = 4 (i.e. w = 20) has been done above
        neg = neg_points_hexadecad(gc, neg, sub)
        #print("~", hex(neg & gc), hex(gc), hex(sub), hex(add))
        if w_sub >= 3:
            #print("neg", w_sub, w_add, hex(gc), hex(sub), hex(add), hex(neg))
            sub = augment_bitvector(sub, neg, 4)
            return find_tetrad_16_large(gc, sub) 
        #print("gc", hex(gc), hex(sub), hex(add))
        tet = gc & mat24.intersect_octad_tetrad(gc, sub | add)
        assert tet
        return tet
    if w_gc == 8:
        if w_add == 1 and (w_sub | (mat24.bw24(gc & neg) & 1)):
            return 0
        tet = gc & mat24.intersect_octad_tetrad(gc, sub | add)
        if tet:
            return tet
    if w_sub & 1 == 0 and w_syn & 3 == 0:
        tet = reduce_even(v, w_syn, syn_list)
        return tet
    if w_gc == 12:
        # The cases w_sub + w_add == 4 and w_sub == w_add == 0 (mod 2)
        # have been done above.
        if w_add == 3:
            return augment_bitvector(add, ~gc, 4)
        if w_sub == 2 and w_add == 0:
            hexad = cohexad(~gc, sub, add)
            return augment_bitvector(0, hexad, 4)            
        # Here (w_sub, w_add) is (0, 1), (1, 0), (1, 1), (1, 2),
        # (2, 1), or (3, 0).
        if w_add == 2:
            hexad = cohexad(gc, add, sub)
            sub =  augment_bitvector(sub, hexad, 2)
            return hexad & ~sub
        if w_add == 1:
            sub = augment_bitvector(sub, gc, 2)
            hexad = cohexad(~gc, sub, add)
            if w_sub < 2:
                return augment_bitvector(add | sub, hexad, 4)
            else:
                hexad = (hexad ^ ~gc) & 0xffffff
                return augment_bitvector(sub, hexad, 4)
        # Here (w_sub, w_add) is (1, 0) or (3, 0). We cannot reduce
        # the case (3, 0). There may be a trick for reducing the case
        # (1, 0), i.e. an umbral undecad. Here we essentially have to
        # 'move' duads inside a tetrad. This is not implemented here. 
        return 0 
    raise ValueError("This should not happen")





###############################################################
# Reduce an element of the Leech lattice mod 3
###############################################################




def reduce_tetrad_leech_mod3(a, tetrad):
    r"""Reduce weight of a vector in the Leech lattice mod 3

    Let ``a`` be a vector in the Leech lattice mod 3 given in
    Leech lattice mod 3 encoding. Let the ``t`` be  a tetrad
    given as a bit vector. Let the weight of ``a`` be the number
    of nonzero bits in ``a``.

    The function tries to reduce the weight of ``a`` by applying
    the following procedure. It first computes a permmutation in the
    Mathieu group :math:`M_{24}` that  maps ``t`` to the standard
    tetrad (0,1,2,3). Then it applies a power :math:`\xi^e, e = \pm1`
    of the generator :math:`\xi`  of :math:`\mbox{Co}_1`.

    In case of success the function returns an element ``g``
    of :math:`\mbox{Co}_1` such that ``w * g`` has a smaller
    weight than ``g``. Here ``g`` is returned as a numpy array of
    generators as defined in file ``mmgroup_generators.h``.
    
    For each tetrad ``t'`` of the standard sextet the function
    computes a score for exponenets :math:`e = \pm 1`  as follows:

    The score for :math:`e = 1` is +1 if ``t'`` has weight 4 and
    an odd number of entries -1. The score for :math:`e = -1` is +1
    if ``t'`` has weight 4 and an even number of entries -1.
    The score for any :math:`e` is -1 if if ``t'`` has weight 1.
    In all other cases the score is 0.

    The function succeeds with exponent :math:`e` if the sum of the
    scores of all tetrads is strictly positive for that exponent.
    Otherwise it fails and raises an exception.

    We remark that in case of success the weight of (the support of)
    ``a`` is reduced by a multiple of 3.
    """
    g = np.array([0x20000000, 0x60000000], dtype = np.uint32)
    length, tetrad_bits =  mat24.vect_to_bit_list(tetrad)
    res, pi = mat24.perm_from_map(tetrad_bits[:length], range(4))
    assert res > 0
    g[0] += mat24.perm_to_m24num(pi)
    a = mul_vector3(a, g[:1]) 
    supp = (a ^ (a >> 24)) & 0xffffff
    neg = a & 0xffffff
    SCORE = 0x2000000100010110    
    n_plus = n_minus = n_bad = 0
    for i in range(0, 24, 4):
        tet = (supp >> i) & 15
        score_tet = (SCORE >> (4 * tet)) & 15
        if score_tet == 2:
            neg1 = (neg >> i) & tet
            wn = (0x6996 >> neg1) & 1
            n_plus += wn
            n_minus += 1 - wn
        else:
            n_bad += score_tet
    if n_plus >= n_minus:
        n, e = n_plus, 1
    else:
        n, e = n_minus, 2
    assert n > n_bad
    g[1] += e
    a = mul_vector3(a, g[1:]) 
    return g, a
 

###############################################################
# Special reduction for an umbral undecad
###############################################################
   


def neg_dodecad(v, syn, neg):
    """Adjust signs of a weight-11 vector in Leech lattice mod 3

    Let ``v`` be the support of a vector ``w`` in the Leech lattice
    mod 3, i.e. the bit vector containing the  nonzero entries 
    of ``w``. Here ``v``must be a umbral undecad, and the (disjoint)
    union of ``v`` and the singleton ``syn`` must be a umbral 
    dodecad. Let ``neg`` be the bit vector containing the bits
    where ``w`` has an entry -1.
    
    We may change the signs of a vector in ``w`` by flipping the
    signs of all entries corresponding to the bits of a Golay code
    word. The function returns a suitable Golay code wort ``gc`` (as
    a bit vector) such that flipping the signs in ``w`` corresponding
    to the bits in ``gc`` makes all entries of ``w`` nonegative. 
    """
    neg &= v
    neg1 = mat24.cocode_as_subdodecad(
        mat24.vect_to_cocode(neg),
        mat24.vect_to_gcode(0xffffff & ~(v | syn)),
        mat24.lsbit24(syn)
    )
    return neg | neg1
    


def reduce_umbral_undecad(a):
    r"""Reduce vector of weight 11 in the Leech lattice mod 3

    Let ``a`` be a vector in the Leech lattice mod 3 such that the
    support `´v`` of ``a`` (i.e. the bit vector containing the nonzero
    entries of ``a``) is an umbral undecad. The function returns an
    element ``g`` of :math:`\mbox{Co}_1` such that the support
    of ``w * g`` is a non-umbral undecad. Here ``g`` is returned as a
    numpy array of generators as defined in 
    file ``mmgroup_generators.h``.

    Element ``g`` first flips suitable bit of vector ``a``; then
    permutes the bits of ``a`` in a suitable way; and, finally, 
    it applies generator :math:`\xi`.  
    """
    g = np.zeros(0, dtype = np.uint32)
    v = (a ^ (a >> 24)) & 0xffffff # non-zero entries of vector
    w_v = mat24.bw24(v)
    syn = mat24.syndrome(v, 0)
    if syn & (syn - 1) or syn & v or mat24.bw24(v) != 11:
        return g, a
    g = np.array([0x20000000, 0x40000000, 0x60000001], dtype = np.uint32)
    tet = augment_bitvector(syn, v, 4)
    length, tetrad_bits =  mat24.vect_to_bit_list(tet)
    res, pi = mat24.perm_from_map(tetrad_bits[:length], range(4))
    assert res > 0
    g[0] += mat24.perm_to_m24num(pi)
    a = mul_vector3(a, g[:1]) 
    v = (a ^ (a >> 24)) & 0xffffff
    syn = mat24.syndrome(v, 0) 
    #print("u11", hex(mul_vector3(tet, g[:1])), hex(v))
    y = neg_dodecad(v, syn, a & v)
    n_swapped = 0
    gc_opp = 0xffffff ^ v ^ syn
    for i in range(4, 24, 4):
        t = v & (15 << i)
        if t and n_swapped < 3:
            duad = syn ^ (t & -t)
            y ^= cohexad(gc_opp, duad, 0) ^ duad
            n_swapped += 1
    g[1] ^= mat24.vect_to_gcode(y)
    a = mul_vector3(a, g[1:2]) 
    #print("u11 a", hex(v), hex(a & v))
    # return g[:2], a
    a = mul_vector3(a, g[2:])
    v = (a ^ (a >> 24)) & 0xffffff; 
    syn = mat24.syndrome(v, 0)
    #print("u11 b", hex(v), hex(syn))
    return g, a



###############################################################
# The final permutation
###############################################################

X = 24
POSITIONS = np.array([
    [X]*8,
    [0] + [X]*7,
    [2,3] + [X]*6,
    [1,2,3] + [X]*5, 
    [X]*8,
    [0, 4, 5, 6, 7] + [X]*3,
    [0, 4, 8, 12, 16, 20] + [X]*2, 
    [2, 3, 0, 1, 4, 5, 8] + [X],
    [0, 1, 2, 3, 4, 5, 8] + [X],
    [0, 1, 2, 3, 4, 5, 8] + [X],
], dtype = np.uint8)


UMBRAL9 = [0, 4, 8, 13,14,15, 17,18,19]


def make_sign_octads():
    d = {}
    OCTAD_MASK = 0x1ff
    SIGN_OCTADS = np.zeros(9, dtype = np.uint32)
    GOOD_OCTADS = dict( [(1 + (1 << i), i) for i in range(1,8)])
    done = 0
    for i in range(759):
        o = mat24.octad_to_vect(i)
        masked = o & OCTAD_MASK
        if masked in GOOD_OCTADS:
            index = GOOD_OCTADS[masked]
            SIGN_OCTADS[index] = o
            done |= 1 << index
    assert done == 0xfe, hex(done)
    SIGN_OCTADS[0] = SIGN_OCTADS[3]
    SIGN_OCTADS[8] = 0xff00 
    return SIGN_OCTADS

_SIGN_OCTADS = None

def sign_octads():
    global _SIGN_OCTADS
    if _SIGN_OCTADS is None:
        _SIGN_OCTADS = make_sign_octads()
    return _SIGN_OCTADS  


def reduce_sign(a):
    SIGN_OCTADS = sign_octads()
    v = (a ^ (a >> 24)) & 0xffffff
    neg = (a >> 24) & 0xffffff
    sign_vector = 0
    if v == 0x111111:
        for i in range(0, 24, 4):
            if neg & (1 << i):
                sign_vector ^= 0xeeeeee ^ (15 << i)
    else:
        for i in range(1, 9):
            if (neg >> i) & 1:
                sign_vector ^= SIGN_OCTADS[i]
        if (neg ^ sign_vector) & 1 and v & 8 == 0:
            sign_vector ^= SIGN_OCTADS[0]
    return mat24.vect_to_gcode(sign_vector)



def reduce_final_perm(a):
    g = np.array([0x20000000, 0x40000000], dtype = np.uint32)
    v = (a ^ (a >> 24)) & 0xffffff # non-zero entries of vector
    w_v = mat24.bw24(v)
    if w_v == 0:
        return g[:0], a
    if w_v <= 6:
        _, v_bits =  mat24.vect_to_bit_list(v)
        res, pi = mat24.perm_from_map(v_bits[:w_v], POSITIONS[w_v, :w_v])
        assert res > 0
        g[0] += mat24.perm_to_m24num(pi)
        a = mul_vector3(a, g[:1])
        g[1] ^= reduce_sign(a)
        a = mul_vector3(a, g[1:])
        return g, a
    syn = mat24.syndrome(v, 0)
    sub = syn & ~v
    add = syn & v
    w_core, core_bits =  mat24.vect_to_bit_list(v & ~add)
    core_bits = core_bits[:w_core]  
    w_sub, sub_bits =  mat24.vect_to_bit_list(sub)
    sub_bits = sub_bits[:w_sub]  
    if w_v == w_core == 9:
        # Then a is a (signed) umbral nonad
        pi = mat24.perm_from_dodecads(sub_bits + core_bits, UMBRAL9) 
        g[0] += mat24.perm_to_m24num(pi)
        a = mul_vector3(a, g[:1]) 
        v = (a ^ (a >> 24)) & 0xffffff
        assert v == 0xEEE000, hex(v)
        neg = (a >> 24) & 0xffffff
        neg1 = neg_dodecad(0xEEE111, 1, neg)
        gc = mat24.vect_to_gcode(neg | neg1)
        g[1] ^= gc
        a = mul_vector3(a, g[1:]) 
        assert a == 0xeee000
        return g, a
    w_add, add_bits =  mat24.vect_to_bit_list(add)
    add_bits = add_bits[:w_add]
    if w_add == 1 and  w_v + w_sub == 9:
        base_bits = (sub_bits + core_bits)[:6] + add_bits
        res, pi = mat24.perm_from_map(base_bits, POSITIONS[w_v, :7])
        assert res > 0  
        g[0] += mat24.perm_to_m24num(pi)
        a = mul_vector3(a, g[:1])
        g[1] ^= reduce_sign(a)
        a = mul_vector3(a, g[1:])
        return g, a
    return g[:0], a
        
     
###############################################################
# The final reduction function
###############################################################

def reduce_weight_leech_mod3(a):
    g = np.zeros(0, dtype = np.uint32)
    #w_start = weight_vector3(a)
    #print("w_start =", w_start)
    for i in range(6):
        tetrad = find_tetrad_leech_mod3(a)
        #print('tet', hex(tetrad))
        if mat24.bw24(tetrad) != 4:
            assert tetrad == 0
            if weight_vector3(a) >= 12: 
                print("too heavy!!!", weight_vector3(a))
                raise ValueError("too heavy")
            #print("done")
            return g, reduce_vector3(a)
        g1, a = reduce_tetrad_leech_mod3(a, tetrad)
        #print(str_vector3(a), ", weight =", weight_vector3(a))
        g = np.concatenate([g, g1]) 
    raise ValueError("WTF")


def reduce_leech_mod3(a):
    global set_a_reduced
    if import_pending:
        import_all()
    g, a = reduce_weight_leech_mod3(a)
    if mat24.bw24((a ^ (a >> 24)) & 0xffffff) == 11:
        g1, a = reduce_umbral_undecad(a)
        g2, a = reduce_weight_leech_mod3(a)
        g = np.concatenate([g, g1, g2])
    g3, a = reduce_final_perm(a)
    g = np.concatenate([g, g3])
    a = reduce_vector3(a)
    return g, a
     



###############################################################
# Sum up shortest vector in Leech lattice mod 3.
###############################################################

def binom(n, k):
    """Binomial coefficent"""
    return math.factorial(n) // math.factorial(k) // math.factorial(n - k)


def norms_24_mod_3():
    """Number of vectors of given norm in Leech lattice mod 3

    The function returns an array ``a`` such that ``a[n % 3]`` is
    the number of vectors of norm ``n`` in the Leech lattice mod 3. 
    """
    a = [0, 0, 0]
    for k in range(25):
        a[k % 3] += binom(24, k) << k
    assert a[1] == a[2]
    assert sum(a) == 3**24
    return a


# N_LEECH[t] is the number of vectors of type t in the Leech lattice.
# The names of these types are taken from [CS99], Ch. 10.3.3. 
# The numbers of the vectors of a given type are taken from the ATLAS.
C = 65520
N_LEECH = {
   '0' : 1,  '2' : 3 * C,  '3': 256*C,  '4': 6075 * C,  '5': 70656 * C, 
   '6_22' : 6900 * C,  '6_23' : 518400 * C,  '7' : 2861568 * C,
   '8_42' : 12295800 * C, '9_33' : 32972800 * C 
}

# Possible types of shortest vectors in the Leech lattice mod 3
TYPES_MOD3 = list(N_LEECH.keys())

N_SHARED = {
   '6_22' : 3,  '7' : 2,  '8_42' : 9,  '9_33' : 36
}


def n_shared(vtype):
    """Return No of shortest vectors of a type in Leech lattice mod 3

    Given a type ``vtype`` of a vector in the Leech lattice mod 3,
    the function returns the number of shortest possible vectors in
    the Leech lattice that map to a vector of that type in the
    Leech lattice mod 3. 

    Here ``vtype`` must be a string the the list ``TYPES_MOD3``.
    """
    return N_SHARED[vtype] if vtype in  N_SHARED else 1

def type_value(vtype):
    """Map the types in the list TYPES_MOD3 to their numeric values"""
    return int(vtype[0])


def check_types():
    """Check vector types in Leech lattice mod 3
    
    Basically, the function checks that vector types in Leech lattice
    mod 3 sum up as expected.
    """
    norms_mod3 = norms_24_mod_3()
    n_shortest_mod3 = [0, 0, 0]
    for x in TYPES_MOD3:
         t = type_value(x) % 3
         n_shortest_mod3[t] += N_LEECH[x] / n_shared(x)
    for i in range(3):
        assert n_shortest_mod3[i] == norms_mod3[i]



###############################################################
# Analysing test results
###############################################################



set_reduced = set()  # Reduced vectors in Leech lattice mod 3
maxlen = 0           # Max. length of reducing group element


def bit_vector_name(w_v, w_sub, w_add):
    MAMES = { (7,2,1): 'U', (8,1,1): 'T', (6,3,1):'U',
       (9,3,0):'U', (9,0,1):'S' }        
    if w_v == 0:
        return "(0)"
    if w_v < 6:
        t = 'S'
    else:
        t = MAMES[(w_v, w_sub, w_add)]
    return "%s_%d" % (t, w_v)  


def leech_name(w_v, w_sub, w_add):
    A_NAMES = ['0', '4', '2', '6_22',  0, '8_42', '6_23', '7', '5']
    if w_v < len(A_NAMES):
         return str(A_NAMES[w_v])
    DICT_NAMES = {(9,3,0) : '9_33', (9,0,1) : '3'}
    return str(DICT_NAMES[(w_v, w_sub, w_add)])
     
def group_name(w_v, w_sub, w_add):
    A_NAMES = ['2.Co_1', '2^{11}.M_{23}', 'Co_2', 'PSU_6(2):S_3',
               None,  '2^{1+8}.A_9', 'M_{24}', 'HS:2', 'McL']
    if w_v < len(A_NAMES):
         return str(A_NAMES[w_v])
    DICT_NAMES = {(9,3,0):'3^6:2.M_12', (9,0,1):'Co_3'}
    return str(DICT_NAMES[(w_v, w_sub, w_add)])
     
def analyse_final(g, a):
    global set_reduced, maxlen
    maxlen = max(len(g), maxlen)
    set_reduced.add(a)

def display_final_analysis():
    data = sorted(list(set_reduced))
    print("""Reduced vectors in Leech lattice mod 3:
Representative : (w, w_sub, w_add), set name; Leech type, Aut; N"""
    )
    for a in data:
        v = (a ^ (a >> 24)) & 0xffffff;
        neg = "[-]" if (a >> 24) & 0xffffff else "   "
        w_v, w_syn, syn_list = bitvector_syndromes(v)
        w_add, _ = syn_list[0]
        w_sub = w_syn - w_add
        v_name = bit_vector_name(w_v, w_sub, w_add)
        type_name = leech_name(w_v, w_sub, w_add)
        aut_name = group_name(w_v, w_sub, w_add)
        N = n_shared(type_name)
        print("0x%06x%s: (%d, %d, %d), %-4s;  type %-5s, %-15s ; %2d" %
           (v, neg, w_v, w_sub, w_add, v_name, type_name, aut_name, N))
    print("\nMax length of reduction word: %d" %
        maxlen)


###############################################################
# Tests
###############################################################


from random import randint, sample



MIN_TESTS = 10000


def reduce_leech_mod3_testdata(ntests = 1):
    ntests = max(int(ntests), MIN_TESTS)
    for i in range(9):
        yield (1 << i) - 1
    A_TESTS = [min(i * ntests // 384, 2 * binom(24, i) << i) for i in range(16)]
    A_TESTS += [16 * ntests // 384] * 9
    #print(A_TESTS)
    for i, n_len_tests in enumerate(A_TESTS):
        for j in range(n_len_tests):
            la = sample(range(24), i)
            rnd = randint(0, 0xffffff)
            v = sum(1 << (24 * ((rnd >> i) & 1) + x) for i, x in enumerate(la))
            yield v
    for i in range(ntests // 3):
        yield randint(0, (1 << 48) - 1)

def test_reduce_leech_mod3(ntests = 1000, f_test = None, verbose = 0):
    import_all()
    check_types()
    M = "Testing function reduce_leech_mod3() with about %d tests"
    print(M % max(int(ntests), MIN_TESTS))
    for i, a in enumerate(reduce_leech_mod3_testdata(ntests)):
        g, a1 = reduce_leech_mod3(a)
        a_g = mul_vector3(a, g)
        assert a_g == a1, [str_vector3(a_g), str_vector3(a1)]
        analyse_final(g, a1)
        if f_test:
            g_test, a1_test = f_test(a)
            assert reduce_vector3(a1_test) == a1
            assert  g_test == g
        if verbose > 1 and (i + 1) % 10000 == 0:
            print('\r%d       ' % (i+1), end = '')
            sys.stdout.flush()
    if verbose > 1:
        print("\r"+ 40 * ' ', end = "\r")
        sys.stdout.flush()
    if verbose:
        display_final_analysis()
    if len(set_reduced) == 10:
        print("Test passed")
    if len(set_reduced) > 10:
        ERR = "Could not reduce all test cases"
        raise ValueError(ERR)
    if len(set_reduced) < 10:
        MSG = "Test did not fail, but too few cases have been tested"
        print(MSG)



###############################################################
# Crated tables for generating C code
###############################################################




class MockupTables:
    directives = {}
    tables = {
        "Leech3SignOctads" : [0] * 10,
        "Leech3BitPositions" : POSITIONS,
    }
    def __init__(self, *args, **kwds):
        pass

     
class Tables(MockupTables):
    def __init__(self, *args, **kwds):
        import_all()
        self.tables = {}       
        self.tables.update(MockupTables.tables)
        self.tables["Leech3SignOctads"] = sign_octads()

###############################################################
# Main function
###############################################################


if __name__ == "__main__":
    test_reduce_leech_mod3(1.e4, None, verbose = 2)



       
