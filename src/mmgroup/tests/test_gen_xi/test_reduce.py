


from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, choices, shuffle, sample
from numbers import Integral
from collections import defaultdict

import numpy as np
import pytest

from mmgroup import XLeech2, Xsp2_Co1, PLoop, GCode, AutPL, Cocode, GcVector

from mmgroup import mat24
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_subtype
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_start_type4
from mmgroup.generators import gen_leech2_start_type24
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_reduce_type2_ortho
from mmgroup.generators import gen_leech2_reduce_type4

from mmgroup.tests.test_gen_xi.test_leech2_type import TYPE_DATA
from mmgroup.tests.test_gen_xi.test_leech2_type import rand_n_elem
from mmgroup.tests.test_gen_xi.test_leech2_type import xs_vector

# Standard vector in the Leech lattice mod 2 in Leech lattice encoding
# The standard fram \Omega
OMEGA = 0x800000
# The standard type-2 vector \beta
BETA = 0x200 

assert Cocode(BETA).syndrome() == GcVector(0xC)

#####################################################################
# Auxiliary functions
#####################################################################

def xi_reduce_odd_type4(v, verbose = 0):
    r"""Compute power of :math:`\xi` that reduces a vector ``v``

    Let ``v`` be a vector in the Leech lattice mod 2 in Leech
    lattice encoding. We assume that ``v`` is of subtype 0x43.

    We compute an exponent ``e`` such that :math:`\xi^e` maps 
    ``v`` to a vector of subtype 0x42 or 0x44.

    The function returns ``e`` if ``v`` is mapped to type 0x42
    and ``0x100 + e`` if ``v`` is mapped to type 0x44.  A negative
    return value indicates that no such exponent ``e`` exists.
    """
    assert v & 0x800         # Error if cocode part ov v is even
    coc = (v ^ mat24.ploop_theta(v >> 12)) & 0xfff
    # Obtain cocode as table of bit fields of 5 bits
    tab = mat24.syndrome_table(coc & 0x7ff)
    # Check if the syndrome bits are in 3 different MOG columns.
    # We first XOR bit field i with bit field (i-1)(mod 3)
    # and then zero the lowest two bits of each bit field. 
    tab ^= ((tab >> 5) & 0x3ff) ^ ((tab & 0x1f) << 10);
    tab &= 0x739c;
    # Now all three bit fields are nonzero iff the syndrome bits
    # are in three differnt columns. Next add 32 - 4 to each bit
    # field in order to produce a carry if the field is nonzero.
    tab += 0x739c;
    # Next we isolate the three carry bits 
    tab &= 0x8420;
    # Return -1 if all carry bits are set, i.e all syndrome bits
    # are in different columns.
    if (tab == 0x8420):
         return -1
    # Let scalar be the scalar product of the Golay part of v
    # with the standard tetrad \omega
    scalar = (v >> 22) & 1;
    # Exponent for element \xi of G_x0 is 2 - scalar
    exp = 2 - scalar;
    if verbose:
        w = gen_leech2_op_atom(v, 0x60000000 + exp) 
        print("Reducing c = %s, subtype %s, t=%s, e=%d, to v = %s, subtype %s"
            %  (hex(mat24.cocode_syndrome(coc,0)), hex(gen_leech2_subtype(v)),
            hex(tab), exp,  hex(mat24.gcode_to_vect(w >> 12)),
            hex(gen_leech2_subtype(w))
            ))
    # Return exponent for \xi in the lower 4 bits of the retrun value;
    # Return 0 in bit 8 if all syndrome bits of v are in the same
    # MOG column and 1 in bit 8 otherwise.
    return ((tab != 0) << 8) + exp

def xi_reduce_odd_type2(v, verbose = 0):
    r"""Compute power of :math:`\xi` that reduces a vector ``v``

    Let ``v`` be a vector in the Leech lattice mod 2 in Leech
    lattice encoding. We assume that ``v`` is of subtype 0x21.

    We type to compute an exponent ``e`` such that :math:`\xi^e` 
    maps ``v`` to a vector of subtype 0x22.

    The function returns ``e`` if such an eponent exists.  A negative
    return value indicates that no such exponent ``e`` exists.
    """
    assert v & 0x800         # Error if cocode part ov v is even
    # Let scalar be the scalar product of the Golay part of v
    # with the standard tetrad \omega
    scalar = (v >> 22) & 1;
    # Exponent for element \xi of G_x0 is 2 - scalar
    exp = 2 - scalar;
    return exp


def xi_reduce_octad(v, verbose = 0):
    r"""Compute power of :math:`\xi` that reduces a vector ``v``

    Let ``v`` be a vector in the Leech lattice mod 2 in Leech
    lattice encoding. We assume that ``v`` is of subtype 
    0x22, 0x42, or 0x44.

    We try to compute an exponent ``e`` such that :math:`\xi^e`  
    maps ``v`` to a vector of subtype 0x40 or 0x20.

    The function returns ``e`` if such an eponent exists.  A negative
    return value indicates that no such exponent ``e`` exists.

    Assume :math:`v_2 = \lambda_d + \lambda_\delta + \lambda_\epsilon`, 
    where :math:`d \in \mathcal{C}, \delta, \epsilon \in \mathcal{C}^*`, 
    with :math:`d, \delta` grey, even,  :math:`\epsilon` coloured. 
    The function returns
  
    :math:`e=0` if :math:`d=0 \pmod{\Omega}`,
     
    :math:`e=1` if :math:`\delta=\gamma(d) \pmod{\omega} `,

    :math:`e=2` if :math:`\delta=0  \pmod{\omega} `. 
  
    In all other cases there is no suitable exponent :math:`e`.
    """
    #print("Std12", hex(v))
    if v & 0x7ff800 == 0:
        return 0
    if v & 0x7f080f == 0: 
        return 1
    parity = -((0x6996 >> (v & 0xf)) & 1)
    v ^= ((v >> 12) ^ parity) & 0xf
    if v & 0x7f080f == 0: 
        return 2
    return -1 # no exponent found


def xi_reduce_dodecad(v, verbose = 0):
    r"""Compute power of :math:`\xi` that reduces a vector ``v``

    Let ``v`` be a vector in the Leech lattice mod 2 in Leech
    lattice encoding. We assume that ``v`` is of subtype 0x46.

    We compute an exponent ``e`` such that :math:`\xi^e` maps 
    ``v`` to a vector of subtype 0x44.

    The function returns ``e`` if such an eponent exists.  A negative
    return value indicates that no such exponent ``e`` exists.
    """
    # Let ``vect`` be the Golay code part of v as a bit vector.
    vect = mat24.gcode_to_vect(v >> 12)
    # Set bit 4*i of s if all bits 4*i, 4*i+1, 4*i+2, 4*i+3 of 
    # ``vect`` are equal, otherwise clear bit 4*i, for 0 <= i < 6.
    s1 = vect | (vect >> 2); s1 = s1 | (s1 >> 1);
    s0 = vect & (vect >> 2); s0 = s0 & (s0 >> 1); 
    s = (s0 | ~s1) & 0x111111;
    # If the Golay code part of v is a docecad then either no or two
    # bits in s are set. Fail if no bit in s is set.
    if (s == 0):
        return -1
    # Here two bits of s (in two different MOG columns) are set. 
    # Set all bits in a MOG column if one bit is set in that column.
    # Thus the bits being set in s form a grey even octad.
    s *= 15;
    # Let 'coc' be the cocode part of v
    coc = v ^ mat24.ploop_theta(v >> 12)
    # Compute scalar product of octad s and ``coc`` in ``scalar``
    tab = mat24.syndrome_table((mat24.recip_basis[0] ^ coc) & 0x7ff)
    scalar = s ^ (s >> (tab & 31)) ^  (s >> ((tab >> 5) & 31))
    scalar ^=  (s >> ((tab >> 10) & 31))
    scalar &= 1
    # The requested exponent is now equal to ``2 - scalar``.
    exp = 2 - scalar
    return exp

def apply_perm(v, src, dest, n, log_list, verbose = 0):
    r"""Apply permutation to vector in Leech lattice mod 2.
  
    The function computes a permutation :math:`\pi` that maps
    the entries of the array ``src`` of length ``n`` to
    the entries of the array ``dest`` (of the same length) in
    the given order. 

    Let :math:`v_2` be the vector in the Leech lattice mod  2 given 
    by parameter ``v2``. The function returns :math:`v_2 x_\pi`.
    Parameter ``v2`` and the return value are given in Leech
    lattice encoding.
  
    Parameter ``p_res`` points to an integer where the function 
    stores the element :math:`x_\pi` as a generator of the
    monster group as as described  in file ``mmgroup_generators.h``.
    That generator is stored with tag  ``MMGROUP_ATOM_TAG_IP`` so
    that we can compute the inverse of :math:`\pi` very 
    efficiently. 

    We compute the inverse of the lowest permutation (in lexical
    order) that maps ``dest[:n]`` to ``src[:n]``.
    """
    res, p =  mat24.perm_from_map(dest[:n], src[:n])
    assert res > 0, (res, dest[:n], src[:n])
    p_inv =  mat24.inv_perm(p)
    p_num =  mat24.perm_to_m24num(p)
    log_list.append(0xA0000000 + p_num)
    xd = (v >> 12) & 0xfff
    xdelta = (v ^ mat24.ploop_theta(xd)) & 0xfff
    m =  mat24.perm_to_matrix(p_inv)
    xd = mat24.op_gcode_matrix(xd, m)
    xdelta = mat24.op_cocode_perm(xdelta, p_inv)
    v_out = (xd << 12) ^ xdelta ^ mat24.ploop_theta(xd)
    if verbose:
        print("Apply permutation (mapping v to gcode %s):\n%s" %
          (hex(mat24.gcode_to_vect(v_out >> 12)), p_inv)
        )
    return v_out


#####################################################################
# Find and apply permutations in Mat24
#####################################################################

OCTAD = [0,1,2,3,4, 8,9]

def find_octad_permutation(v, result, verbose = 0):
    coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff
    w = mat24.gcode_weight(v >> 12) 
    vect = mat24.gcode_to_vect((v ^ ((w & 4) << 21)) >> 12)
    src = mat24.vect_to_list(vect, 5)
    syn =  mat24.cocode_syndrome(coc, src[0]) & ~vect
    if syn:
        v5 = (1 << src[0]) | (1 << src[1]) | (1 << src[2])
        v5 |= syn 
        special = mat24.syndrome(v5, 24);
        src = src[:3]
        src.append(mat24.lsbit24(special & vect))
        src.append(mat24.lsbit24(vect & ~(special | v5)))
        src.append(mat24.lsbit24(syn))  
        syn &= ~(1 << src[-1])
        src.append(mat24.lsbit24(syn))  
    return apply_perm(v, src, OCTAD, len(src), result, verbose)


#####################################################################
# Function reduce_type_2
#####################################################################

def reduce_type2(v, verbose = 1):
    r"""Map (orthgonal) short vector in Leech lattice to standard vector

    This is a python implementation of the C function
    ``gen_leech2_reduce_type2`` in file ``gen_leech.c``.
   
    Let ``v \in \Lambda / 2 \Lambda`` of type 2 be given by 
    parameter ``v`` in Leech lattice encoding. 

    Let ``beta`` be the short vector in the Leech  lattice propotional
    to  ``e_2 - e_3``, where ``e_i`` is  the ``i``-th basis vector
    of ``\{0,1\}^{24}``.
   
    Then the function constructs a ``g \in G_{x0}`` 
    that maps ``v`` to ``beta``.
 
    The element ``g`` is returned as a word in the generators
    of ``G_{x0}`` of length ``n \leq 6``. Each atom of the 
    word ``g`` is encoded as  defined in the header 
    file ``mmgroup_generators.h``. 

    The function stores ``g`` as a word of generators in the
    array ``pg_out`` and returns the length  ``n``  of that
    word. It returns a negative number in case of failure, 
    e.g. if ``v`` is not of type 2.
    """
    vtype = gen_leech2_subtype(v)
    if (vtype >> 4) != 2:
        raise ValueError("Vector is not short")
    result = []
    for _i in range(4):
        if verbose:
            coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff
            vt = gen_leech2_subtype(v)
            coc_anchor = 0
            if vt in [0x22]:
                w = mat24.gcode_weight(v >> 12) 
                vect = mat24.gcode_to_vect((v ^ ((w & 4) << 21)) >> 12)
                coc_anchor = mat24.lsbit24(vect)
            coc_syn = Cocode(coc).syndrome_list(coc_anchor)
            gcode = mat24.gcode_to_vect(v >> 12)
            print("Round %d, v = %s, subtype %s, gcode %s, cocode %s" % 
                (_i, hex(v & 0xffffff), hex(vt), hex(gcode),  coc_syn)
            )
        assert vtype == gen_leech2_subtype(v)
        if vtype == 0x21:
            exp = xi_reduce_odd_type2(v)
            vtype = 0x22
        elif vtype == 0x22:
            exp = xi_reduce_octad(v)            
            if exp < 0:
                w = mat24.gcode_weight(v >> 12) 
                vect = mat24.gcode_to_vect((v ^ ((w & 4) << 21)) >> 12)
                src = mat24.vect_to_list(vect, 4) 
                dest = [0,1,2,3]
                v = apply_perm(v, src, dest, 4, result, verbose)
                exp = xi_reduce_octad(v) 
                assert exp >= 0
            vtype = 0x20
        elif vtype == 0x20:
            exp =  0
            # map v to stadard cocode word [2,3]
            if v & 0x7fffff != 0x200:
                 syn = (mat24.cocode_syndrome(v, 0))
                 src = mat24.vect_to_list(syn, 2)
                 v  = apply_perm(v, src, [2,3], 2, result, verbose)
            # correct v2 if v2 is the cocode word [2,3] + Omega
            if v & 0x800000:
                atom = 0xC0000200  
                   # operation y_d such that d has odd scalar
                   # product with cocode word [2,3]
                v = gen_leech2_op_atom(v, atom)
                result.append(atom)
            assert v & 0xffffff == 0x200
            return np.array(result, dtype = np.uint32)
        else:
            raise ValueError("WTF")
        if exp: 
            exp = 0xE0000003 - exp
            v = gen_leech2_op_atom(v, exp)
            result.append(exp)
    raise ValueError("WTF1")

#####################################################################
# Test function reduce_type_2
#####################################################################


#####################################################################
# Test data for function reduce_type_2()


def type2_testdata(ntests):
    """The function yields type-2 vectors in the Leech lattice mod 2

    These type-2 vectors are given in **Leech lattice encoding**.
    """
    testdata = [
        0x200,
        0x1000200,
        0x800200,
        0x1800200,
    ]
    for t in testdata:
         yield t
    for ploop, cocode, subtype in TYPE_DATA:
        if subtype in [0x20, 0x22, 0x21]:
            v = xs_vector(ploop, cocode)
            yield v
            for i in range(ntests):
                g =  rand_n_elem()
                yield g.xsp_conjugate(v)

    for i in range(ntests):
        yield XLeech2('r', 2).ord

#####################################################################
# Testing function reduce_type_2


@pytest.mark.gen_xi
def test_reduce_type_2(ntests = 500, verbose = 0):
    """Test function ``reduce_type2`` """
    for n, v in enumerate(type2_testdata(ntests)):
        if verbose:
            print(" \nTest %d, v = %s, subtype  =%s" % 
                (n+1, hex(v), hex(gen_leech2_subtype(v)))
            )
        op = reduce_type2(v, verbose)
        w = gen_leech2_op_word(v, op, len(op))
        assert w & 0xffffff == BETA, hex(w)
        a = np.zeros(6, dtype = np.uint32)
        l = gen_leech2_reduce_type2(v, a)
        if l < 0:
            err = "Error %s in function gen_leech2_reduce_type2"
            raise ValueError(err % hex(l & 0xffffffff))
        a = a[:l]
        assert list(a) == list(op), ((a), (op))  

#####################################################################
# Function reduce_type_2_ortho
#####################################################################

def reduce_type2_ortho(v, verbose = 0):
    r"""Map (orthgonal) short vector in Leech lattice to standard vector

    This is a python implementation of the C function
    ``gen_leech2_reduce_type2_ortho`` in file ``gen_leech.c``.
   
    Let ``v \in \Lambda / 2 \Lambda`` of type 2 be given by 
    parameter ``v`` in Leech lattice encoding. 

    In the real Leech lattice, (the origin of) the vector ``v`` must
    be orthogonal to the standard short vector ``beta``. Here ``beta``
    is the short vector in the Leech  lattice  propotional
    to  ``e_2 - e_3``, where ``e_i`` is  the ``i``-th basis vector
    of ``\{0,1\}^{24}``.
   
    Let ``beta'`` be the short vector in the Leech lattice propotional
    to  ``e_2 + e_3``.  Then the function constructs a ``g \in G_{x0}`` 
    that maps ``v`` to ``beta'`` and fixes ``beta``.
 
    The element ``g`` is returned as a word in the generators
    of ``G_{x0}`` of length ``n \leq 6``. Each atom of the 
    word ``g`` is encoded as  defined in the header 
    file ``mmgroup_generators.h``. 

    The function stores ``g`` as a word of generators in the
    array ``pg_out`` and returns the length  ``n``  of that
    word. It returns a negative number in case of failure, 
    e.g. if ``v`` is not of type 2,  or not orthogonal 
    to ``beta'`` in the real Leech lattice.
    """
    vtype = gen_leech2_subtype(v)
    if (vtype >> 4) != 2:
        raise ValueError("Vector is not short")
    if gen_leech2_type(v ^ 0x200) != 4:
        raise ValueError("Vector not orthogonal to standard vector")
    result = []
    for _i in range(4):
        if verbose:
            coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff
            vt = gen_leech2_subtype(v)
            coc_anchor = 0
            if vt in [0x22]:
                w = mat24.gcode_weight(v >> 12) 
                vect = mat24.gcode_to_vect((v ^ ((w & 4) << 21)) >> 12)
                coc_anchor = mat24.lsbit24(vect)
            coc_syn = Cocode(coc).syndrome_list(coc_anchor)
            gcode = mat24.gcode_to_vect(v >> 12)
            print("Round %d, v = %s, subtype %s, gcode %s, cocode %s" % 
                (_i, hex(v & 0xffffff), hex(vt), hex(gcode),  coc_syn)
            )
        assert vtype == gen_leech2_subtype(v)
        if vtype == 0x21:
            exp = xi_reduce_odd_type2(v)
            vtype = 0x22
        elif vtype == 0x22:
            exp = xi_reduce_octad(v)            
            if exp < 0:
                w = mat24.gcode_weight(v >> 12) 
                vect = mat24.gcode_to_vect((v ^ ((w & 4) << 21)) >> 12)
                if vect & 0x0c:
                    vect &= ~0x0c
                    src = mat24.vect_to_list(vect, 2) + [2,3]
                    dest = [0,1,2,3]
                else:
                    src = [2,3] + mat24.vect_to_list(vect, 3)
                    v5 = (1 << src[2]) | (1 << src[3]) | (1 << src[4])
                    v5 |= 0x0c
                    special = mat24.syndrome(v5, 24);
                    src.append(mat24.lsbit24(special & vect))  
                    dest = [2,3,4,5,6,7]
                v = apply_perm(v, src, dest, len(src), result, verbose)
                exp = xi_reduce_octad(v) 
                assert exp >= 0
            vtype = 0x20
        elif vtype == 0x20:
            if ((v & 0xffffff) == 0x800200):
                return np.array(result, dtype = np.uint32)
            syn = (mat24.cocode_syndrome(v, 0)) & ~0xc
            if syn and syn != 3:
                src =  mat24.vect_to_list(syn, 2) + [2, 3]
                v = apply_perm(v, src, [0, 1, 2, 3], 4, result, verbose)
            exp = 2 - ((v >> 23) & 1);
        else:
            raise ValueError("WTF")
        if exp: 
            exp = 0xE0000003 - exp
            v = gen_leech2_op_atom(v, exp)
            result.append(exp)
    raise ValueError("WTF1")


#####################################################################
# Test function reduce_type_2_ortho
#####################################################################




def type2_ortho_testdata(ntests):
    """Yields special type-2 vectors in the Leech lattice mod 2

    These type-2 vectors are orthogonal to the vector BETA in the 
    real Leech lattice and given in **Leech lattice encoding**.
    """
    for v in type2_testdata(2 * ntests):
        if gen_leech2_type(v ^ BETA) == 4:
            yield v



@pytest.mark.gen_xi
def test_reduce_type_2_ortho(ntests = 500, verbose = 0):
    """Test function ``reduce_type2_ortho`` """
    for n, v in enumerate(type2_ortho_testdata(ntests)):
        if verbose:
            print(" \nTest %d, v = %s, subtype  =%s" % 
                (n+1, hex(v), hex(gen_leech2_subtype(v)))
            )
        op = reduce_type2_ortho(v, verbose)
        w = gen_leech2_op_word(v, op, len(op)) 
        assert w & 0xffffff == 0x800200, hex(w)
        b = gen_leech2_op_word(BETA, op, len(op)) 
        assert b == BETA
        a = np.zeros(6, dtype = np.uint32)
        l = gen_leech2_reduce_type2_ortho(v, a)
        if l < 0:
            err = "Error %s in function gen_leech2_reduce_type2_ortho"
            raise ValueError(err % hex(l & 0xffffffff))
        a = a[:l]
        assert list(a) == list(op), ((a), (op))  
    
    

#####################################################################
# Function reduce_type_4_std
#####################################################################

LSTD = [0,1,2,3,4,5]

def reduce_type4_std(v, verbose = 0):
    r"""Map type-4 vector in Leech lattice to standard vector

    This is (almost) a python implementation of the C function
    ``gen_leech2_reduce_type4`` in file ``gen_leech.c``.
   
    Let ``v \in \Lambda / 2 \Lambda`` of type 4 be given by 
    parameter ``v`` in Leech lattice encoding. 

    Let ``Omega`` be the type- vector in the Leech  lattice 
    corresponding to the standard coordinate frame in the Leech
    lattice.
   
    Then the function constructs a ``g \in G_{x0}`` 
    that maps ``v`` to ``Omega``.
 
    The element ``g`` is returned as a word in the generators
    of ``G_{x0}`` of length ``n \leq 6``. Each atom of the 
    word ``g`` is encoded as  defined in the header 
    file ``mmgroup_generators.h``. 

    The function stores ``g`` as a word of generators in the
    array ``pg_out`` and returns the length  ``n``  of that
    word. It returns a negative number in case of failure, 
    e.g. if ``v`` is not of type 4.

    We remark that the C function ``gen_leech2_reduce_type4`` 
    treats certain type-4  vectors ``v`` in a special way
    as indicated in function ``reduce_type4``.
    """
    if verbose:
        print("Transforming  type-4 vector %s to Omega" %
            hex(v & 0x1ffffff) )
    vtype = gen_leech2_subtype(v)
    result = []
    for _i in range(5):
        coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff
        if verbose:
            vt = gen_leech2_subtype(v)
            coc_anchor = 0
            if vt in [0x42, 0x44]:
                w = mat24.gcode_weight(v >> 12) 
                vect = mat24.gcode_to_vect((v ^ ((w & 4) << 21)) >> 12)
                coc_anchor = mat24.lsbit24(vect)
            coc_syn = Cocode(coc).syndrome_list(coc_anchor)
            gcode = mat24.gcode_to_vect(v >> 12)
            print("Round %d, v = %s, subtype %s, gcode %s, cocode %s" % 
                (_i, hex(v & 0xffffff), hex(vt), hex(gcode),  coc_syn)
            )
        assert vtype == gen_leech2_subtype(v)
        if vtype == 0x48:
            if verbose:
                 res = list(map(hex, result))
                 print("Transformation is\n", res)
            return np.array(result, dtype = np.uint32)
        elif vtype == 0x40:
            if v & 0x7ffbff:
                syn = mat24.cocode_syndrome(coc, 0)
                src = mat24.vect_to_list(syn, 4)
                v = apply_perm(v, src, LSTD, 4, result, verbose)
                #print("after type 40", hex(v),  Cocode(v).syndrome(0))
            exp = 2 - ((v >> 23) & 1)
            vtype = 0x48
        elif vtype in [0x42, 0x44]:
            exp = xi_reduce_octad(v)            
            if exp < 0:
                v = find_octad_permutation(v, result, verbose)
                exp = xi_reduce_octad(v) 
                assert exp >= 0
            vtype = 0x40
        elif vtype == 0x46:
            exp = xi_reduce_dodecad(v, verbose)
            if exp < 0:
                vect = mat24.gcode_to_vect(v >> 12)
                src = mat24.vect_to_list(vect, 4)
                v = apply_perm(v, src, LSTD, len(src), result, verbose)
                exp = xi_reduce_dodecad(v, verbose) 
                assert exp >= 0
            vtype = 0x44
        elif vtype == 0x43:
            exp = xi_reduce_odd_type4(v, verbose)
            if exp < 0:
                vect = mat24.gcode_to_vect(v >> 12)
                syn = mat24.cocode_syndrome(coc, 24)
                src = mat24.vect_to_list(syn, 3)
                #print("coc list", src)
                v = apply_perm(v, src, LSTD[1:], len(src), result, verbose)
                exp = xi_reduce_odd_type4(v, verbose)
                assert exp > 0
            vtype = 0x42 + ((exp & 0x100) >> 7)
            exp &= 3
        else:
            raise ValueError("WTF")
        if exp: 
            exp = 0xE0000003 - exp
            v_old = v
            v = gen_leech2_op_atom(v, exp)
            assert v & 0xfe000000 == 0, (hex(v_old), hex(exp), hex(v))
            result.append(exp)
    raise ValueError("WTF1")


#####################################################################
# Test function reduce_type_4_std
#####################################################################


#####################################################################
# Test data for function reduce_type_4_std)


def type4_testdata(ntests):
    """The function yields type-4 vectors in the Leech lattice mod 2

    These type-4 vectors are given in **Leech lattice encoding**.
    """
    testdata = [
        0x800000,
        0x1800000,
    ]
    for t in testdata:
         yield t
    for ploop, cocode, subtype in TYPE_DATA:
        if subtype in [0x40, 0x42, 0x44, 0x46, 0x43]:
            v = xs_vector(ploop, cocode)
            yield v
            for i in range(ntests):
                g =  rand_n_elem()
                yield g.xsp_conjugate(v)

    for i in range(ntests):
        yield XLeech2('r', 4).ord

    for v in type2_ortho_testdata(ntests):
        yield v ^ BETA

#####################################################################
# Testing function  data for function reduce_type_4_std()

@pytest.mark.gen_xi
def test_reduce_type_4_std(ntests = 500, verbose = 0):
    """Test function ``reduce_type4_std`` """
    for n, v in enumerate(type4_testdata(ntests)):
        if verbose:
            print(" \nTest %d, v = %s, subtype  =%s" % 
                (n+1, hex(v), hex(gen_leech2_subtype(v)))
            )
        op = reduce_type4_std(v, verbose)
        w = gen_leech2_op_word(v, op, len(op)) 
        assert w & 0xffffff == 0x800000, hex(w)


#####################################################################
# Test function reduce_type_4
#####################################################################




def reduce_type4(v, verbose = 0):
    r"""Map type-4 vector in Leech lattice to standard vector

    This is  a python implementation of the C function
    ``gen_leech2_reduce_type4`` in file ``gen_leech.c``.
   
    Let ``v \in \Lambda / 2 \Lambda`` of type 4 be given by 
    parameter ``v`` in Leech lattice encoding. 

    Let ``Omega`` be the type- vector in the Leech  lattice 
    corresponding to the standard coordinate frame in the Leech
    lattice. Let ``beta`` be the short vector in the Leech  
    lattice propotional  to  ``e_2 - e_3``, where ``e_i`` is  
    the ``i``-th basis vector  of ``\{0,1\}^{24}``.
   
    Then the function constructs a ``g \in G_{x0}`` 
    that maps ``v`` to ``Omega``. If ``v + beta`` is of type 2
    and orthogonal to ``beta`` in the real Leech lattice then the 
    returned element ``g`` also fixes ``beta``.
 
    The element ``g`` is returned as a word in the generators
    of ``G_{x0}`` of length ``n \leq 6``. Each atom of the 
    word ``g`` is encoded as  defined in the header 
    file ``mmgroup_generators.h``. 

    The function stores ``g`` as a word of generators in the
    array ``pg_out`` and returns the length  ``n``  of that
    word. It returns a negative number in case of failure, 
    e.g. if ``v`` is not of type 4.
    """
    vtype = gen_leech2_subtype(v)
    if (vtype >> 4) != 4:
        err = "Leech lattice vector must be of type 4"
        raise ValueError(err)
    vtype_beta = gen_leech2_subtype(v ^ BETA)
    if (vtype_beta >> 4) == 2:
        return reduce_type2_ortho(v ^ BETA, verbose)
    return reduce_type4_std(v, verbose)
     

#####################################################################
# Testing function  data for function reduce_type_4_std()


@pytest.mark.gen_xi
def test_reduce_type_4(ntests = 500, verbose = 0):
    """Test function ``reduce_type4`` """
    for n, v in enumerate(type4_testdata(ntests)):
        if verbose:
            print(" \nTest %d, v = %s, subtype  =%s" % 
                (n+1, hex(v), hex(gen_leech2_subtype(v)))
            )
        op = reduce_type4(v, verbose)
        w = gen_leech2_op_word(v, op, len(op)) 
        assert w & 0xffffff == 0x800000, hex(w)
        if gen_leech2_type(v ^ BETA) == 2:
            b = gen_leech2_op_word(BETA, op, len(op)) 
            assert b & 0xffffff == BETA, hex(b)
        a = np.zeros(6, dtype = np.uint32)
        l = gen_leech2_reduce_type4(v, a)
        if l < 0:
            err = "Error %s in function gen_leech2_reduce_type4"
            raise ValueError(err % hex(l & 0xffffffff))
        a = a[:l]
        assert list(a) == list(op), ((a), (op))  
            

