import sys
from collections import OrderedDict

sys.path.append(r"C:\Data\projects\MonsterGit\src")

from random import randint, shuffle, sample
from functools import reduce
from operator import __or__, __xor__
from collections import defaultdict
from collections.abc import Iterable
import re

import numpy as np
import pytest

from mmgroup import Xsp2_Co1, PLoop, GCode, AutPL, Cocode, SubOctad
from mmgroup import Octad, Parity
from mmgroup import mat24
from mmgroup.mat24 import MAT24_ORDER, gcode_weight
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import gen_leech2_op_word




def suboctad_type(octad, w, coc):
    """Return suboctad type.

    Let ``octad`` be an octad, i.e. a bit vector of length 8. Let
    w = 1 (mod 2) if ``octad`` denotes an octad and ``w = 0`` 
    if ``octad`` denotes a complemented octad. Let ``coc`` be
    an even cocode vector in cocode representation.

    The function returns 0 in bit 1 of the return value if the 
    cocode word ``coc`` can be written as a subset of the octad,
    and 1 in bit 1 otherwise.

    the function returns ``1 + w + bit_weight(coc)/2`` in bit 0
    of the return value.

    Then Leech lattice vector ``x_octad + x_coc`` is of subtype 
    0x22 if the return value is zero. Otherwise it is of
    subtype 0x44 (or 0x46) if bit 1 of the return value is 
    (or 1).   
    """
    lsb = mat24.lsbit24(octad)            # ls bit of octad
    syn = mat24.cocode_syndrome(coc, lsb) # cocode syndrome
    wsub = octad & syn == syn             # wsub = 1 if coc is suboctad
    cw = mat24.cocode_weight(coc) >> 1    # cw = cocode_weight(v) / 2
    return 2 * (1 - wsub) + ((w ^ cw ^ 1) & 1)


def leech2_start_type4(v):
    """Return subtype of a Leech lattice frame ``v`` used for reduction

    The function returns the subtype of a vector ``v`` of type 4 in
    the Leech lattice modulo 2. Parameter ``v2`` must be in Leech 
    lattice encoding. The function returns the subtype of ``v`` that
    will be used for reduction in function ``gen_leech2_reduce_type4``.

    In that function we take care of the special case that ``v + v0``
    is of type 2 for a specific short vector ``v0``. 
    
    A simpler (but slower) implementation of thhis function is:

    If ``v ^ v0`` is of type 2 the return the subtype of ``v ^ v0``.
    Otherwise return the subtype of ``v``.

    The function returns 0 if ``v`` is equal to  ``Omega`` and
    a negative value if ``v`` has not type 4.

    This is a refernece implementation for function 
    ``gen_leech2_start_type4()`` in file ``gen_leech.c``.
    """
    if v & 0x7ff800 == 0:
        # Then v or v + Omega is an even cocode element.
        # Return 0 if v == Omega and -1 if v ==  0.
        if v & 0x7fffff == 0:
            return 0 if v & 0x800000 else -1 
        # Let w be the cocode weight. Return -2 if w == 2.
        if mat24.cocode_weight(v) != 4:
            return -2
        # Here v has type 4. Let v23 be the standard type-2 vector.
        # Return 0x20 if v ^ v23 has type 2 and 0x40 otherwise.
        return 0x20 if mat24.cocode_weight(v ^ 0x200) == 2 else 0x40  
    if mat24.scalar_prod(v >> 12, v):
        # Then v has type 3 and we return -3
        return -3
    if v & 0x800:
        # Then the cocode word 'coc' of v is odd.
        coc = (v ^ mat24.ploop_theta(v >> 12)) & 0xfff
        syn = mat24.cocode_syndrome(coc)
        # If 'coc' has weight 1 then v is of type 2 and we return -2.
        if (syn & (syn - 1)) == 0:
            return -2
        # Here v is of type 4. 
        # Return 0x21 if v ^ v23 is of type 2 and 0x43 otherwise.
        if (syn & 0xc) == 0xc and (v & 0x200000) == 0:
            return 0x21
        return 0x43
    w = gcode_weight(v >> 12) # Weight of Golay code part divided by 4
    if w == 3:
        # Then the Golay code part of v is a docecad and we return 0x46.
        return 0x46
    # Here the Golay code part of v is a (possibly complemented) octad.
    # Add Omega to v if Golay part is a complemented octad.
    v ^= (w & 4) << 21
    # Put w = 1 if that part is an octad and w = 0 otherwise.   
    w = (w >> 1) & 1

    # Let 'octad' be the octad in the Golay code word in vector rep.
    octad = mat24.gcode_to_vect(v >> 12)
    coc = v ^ mat24.ploop_theta(v >> 12)  # cocode element of v
    # Return -2 if v is of type 2.
    sub = suboctad_type(octad, w, coc)
    if sub == 0:
        return -2
    # Return 0x22 if v ^ v23 is shsort
    if suboctad_type(octad, w, coc ^ 0x200) == 0:
        return 0x22
    # Otherwise return the subtype of v
    return 0x44 if sub & 2 else 0x42



#######################################################################
# Reduce a type-4 vector
#######################################################################

# Deprecated!!!!

def apply_perm(v, p, log_list):
    p_inv =  mat24.inv_perm(p)
    p_inv_num =  mat24.perm_to_m24num(p_inv)
    log_list.append(0xA0000000 + p_inv_num)
    xd = (v >> 12) & 0xfff
    xdelta = (v ^ mat24.ploop_theta(xd)) & 0xfff
    m =  mat24.perm_to_matrix(p)
    xd = mat24.op_gcode_matrix(xd, m)
    xdelta = mat24.op_cocode_perm(xdelta, p)
    return (xd << 12) ^ xdelta ^ mat24.ploop_theta(xd)


def apply_perm_map(v, src, dest, length, log_list):
    p =  mat24.perm_from_map(src, dest, length)
    return apply_perm(v, p, log_list)


def apply_perm_heptad(v, src, dest, log_list):
    apply_perm_map(v, src, dest, 7, log_list)


def apply_xi(v, exp, log_list):
    exp = 0xE0000003 - exp
    v = gen_leech2_op_atom(v, exp)
    assert v & 0xfe000000 == 0
    log_list.append(exp)
    return v


#######################################################################
# Find suitble power of xi for reducing dodecads and octads
#######################################################################




ERR_STD12_SPACE = "xi**e does not map Leech vector to standard-12 subspace"

def xi_grey_octad_to_tetrad(v):
    """Find exponent of ``xi`` mapping Leech vector ``v`` a subspace 

     

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
    



def xi_dodecad_to_octad(v):
    ok = 0; scalar = 0;
    for i in range(0, 24, 4):
        w = (v >> i) & 15
        good = (0x8001 >> i) & 1
        ok |= good
        scalar ^= w & -good
    scalar = (0x6996 >> scalar) & 1
    return  (2 - scalar) | (ok - 1)
    


STD_OCTAD = GCode([0,1,2,3,4,5,6,7]).ord
print("STD_OCTAD =", hex(STD_OCTAD))

RED_STD_HEPTAD = [0,1,2,3,4,5,8]
RED_IMG_CASE_44 = [0,1,2,3,8,9,4]
RED_STD_DODECAD = [0,1,2,3,8,9,14,15,18,19]


#######################################################################
# Reducing a type-4 vector
#######################################################################



def reduce_type4(v):
    result = []
    src = [None] * 24; dest = [None] * 8
    # Deal with vectors 0 and Omega first
    if (v & 0x7fffff) == 0:
        return ((v >> 23) & 1) - 1;
    # Abort if vector has type 3
    if (mat24.scalar_prod(v >> 12, v)):
        return -3;
    # next deal with an odd cocode vector
    if (v & 0x800):
        coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff
        # Obtain cocode as table of 3 bit fileds of length 5
        tab = mat24.syndrome_table(coc & 0x7ff)
        # Error if cocode part has weight 1
        if (tab >= 0x300):
            return -6;
        # Permute if not all syndome entries are in the same MOG column 
        if ((sy ^  (syn >> 10)) & 0x18):
            syn = (1 << (tab & 31)) + (1 << ((tab >> 5) & 31));
            syn += (1 << ((tab >> 5) & 31));
            assert bw24(syn) == 3
            src[0] = dest[0] = 2;
            src[1] = dest[1] = 3;
            w = 3;
            if (syn & 0xc == 0xc):
                src[2] = lsb24(syn & ~0xc);
                dest[2] = 1;   
            elif (syn & 0xc == 0x0):
                src[2:5] = mat24.vect_to_list(syn, 3);
                dest[2:5] = [5,6,7]
                w = 5;
            else:
                src[:3] = mat24.vect_to_list(syn, 3);
                dest[:3] = [1,2,3]
            v = apply_perm_map(v, src, dest, w, result)
        exp = 2 - ((v >> 22) & 1);
        v = apply_xi(v, exp, result);
        assert gen_leech2_type(v) == 0x42
    
    # Now v is of type 2 or 4 and has an even cocode part.
    # Next we deal the case that the Golay code part of v 
    # is a dodecad.         
    gcode = (v >> 12) & 0xfff
    if mat24.gcode_weight(gcode) & 1:
        if xi_dodecad_to_octad(v) < 0:
            vect = mat24.gcode_to_vect(gcode & 0x7ff);
            src[:4] = mat24.vect_to_list(vect,  4);
            dest[:4] = [0,1,2,3]
            v = apply_perm_map(v, src, dest, 4, result)
        exp = xi_dodecad_to_octad(v)
        assert exp > 0
        v = apply_xi(v, exp, result);
        assert gen_leech2_type(v) == 0x44
        gcode = (v >> 12) & 0xfff

    # Now v is of type 2 or 4 and has an even cocode part.
    # The Golay Code part of v has weight divisible by 8.
    # Next we deal the case that the Golay code part of v 
    # is a (possibly complemented) octad.   
    if (gcode & 0x7ff):
        # Permutation vector if required
        if (xi_grey_octad_to_tetrad(v) < 0):
            # Obtain octad bit vector
            octad = mat24.gcode_to_vect(gcode);
            # correct a complemented octad vector
            octad ^= -((mat24.gcode_weight(gcode) >> 2) & 1);
            octad &= 0xffffff;
            # obtain cocode vector. If cocode vector has weight 4
            # then set bit at position 2
            coc = (v ^  mat24.ploop_theta(gcode)) & 0xfff;
            syn = mat24.cocode_syndrome(coc, 2);
            # remove octad bits from cocode vector
            baby_syn = syn & ~octad;
            src[0] = dest[0] = 2;
            src[1] = dest[1] = 3;
            if baby_syn == 0:
                # Then bits [2,3] are inside the octad    
                octad &= ~0xc; # remove bits [2,3] from octad
                src[2:4] = mat24.vect_to_list(octad, 2); 
                w = 4;
            elif False:  #WTF
                # Then bits [2,3] are outside the octad and
                # Cocode(v) ^ Cocode([2,3]) is inside the octad
                src[2:5] = mat24.vect_to_list(octad, 3);  
                coc = vect_to_cocode(0xc + (1 << src[2]) 
                      + (1 << src[3]) + (1 << src[4]));
                src[5] =  mat24.syndrome_table(coc & 0x7ff) & 0x31;
                dest[2:6] = [4,5,6,7]
            else:
                # Ingnore bits [2,3] and map octad to grey octad.
                # If any syndrome is outside the octad then map
                # it into one tetrad.
                syn = mat24.cocode_syndrome(coc, mat24.lsb24(octad));
                syn &= ~octad;
                src[0:4] = mat24.vect_to_list(octad, 4);
                w = 4;
                if syn:
                    src[4:6] = mat24.vect_to_list(syn, 2);
                    coc = vect_to_cocode(syn + (1 << src[0]) 
                        + (1 << src[1]) + (1 << src[2]));
                    src[3] =  mat24.syndrome_table(coc & 0x7ff) & 0x31;
                    dest[0:6] = [0,1,2,3,4,5]
                    w = 6 
            v = apply_perm_map(v, src, dest, w, result)
        exp = xi_grey_octad_to_tetrad(v)
        assert exp > 0
        v = apply_xi(v, exp, result);
        assert gen_leech2_type(v) == 0x40

    # Now v is of type 2 or 4 and has an even cocode part.
    # The Golay Code part of v is 0 or Omega.
    # Permute v if its cocode is not grey and even
    if mat24.cocode_weight(v) != 4:
        return -4
    if v & 0x7ffbf0:
        syn = mat24.cocode_syndrome(coc, 2);
        # YET TO BE DONE!!!!!!
            
    # YET TO BE DONE!!!!!!    


    #########################################
    for _i in range(5):
        gc = mat24.gcode_to_vect(v >> 12)
        coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff
        vtype = gen_leech2_type(v)
        #print("round", _i, ", type", hex(vtype), hex(gc), hex(Cocode(coc).syndrome(0).ord))
        if vtype == 0x48:
            return np.array(result, dtype = np.uint32)
        elif vtype == 0x40:
            if v & 0x7ffbff:
                sextet = mat24.cocode_to_sextet(coc)
                src =  sextet[:6] + sextet[8:9]
                v = apply_perm_heptad(v, src, RED_STD_HEPTAD, result)
                #print("after type 40", hex(v),  Cocode(v).syndrome(0))
            exp = 2 - ((v >> 23) & 1)
            #print(hex(v),  Cocode(v).syndrome(0))
            #assert v & 0xffffff == 0x800000, hex(v)
        elif vtype == 0x42:
            exp = map_to_standard12_xi_exp(v)
            if exp < 0:
                w, src = mat24.vect_to_bit_list(gc)
                if (w == 16):
                    src = src[16:22] + src[0:1]
                else:
                    src[6] = src[8]
                v = apply_perm_heptad(v, src, RED_STD_HEPTAD, result)
                exp = map_to_standard12_xi_exp(v)
                assert exp >= 0
        elif vtype == 0x44:
            exp = map_to_standard12_xi_exp(v)
            if exp < 0:
                w, oct = mat24.vect_to_bit_list(gc)
                if (w == 16):
                    gc ^= 0xffffff
                    oct = oct[16:]  # + oct[:16]
                #print("\noct =", oct[:9])
                syn =  mat24.cocode_syndrome(coc, oct[0]) & ~gc
                w, syn_list = mat24.vect_to_bit_list(syn)
                #print("syn_list =", syn_list[:2])
                assert w == 2
                a = sum([1 << i for i in oct[:5]])
                a1 = a ^ (1 << syn_list[0]) ^ (1 << syn_list[1])
                #print("a1 =", hex(a1))
                special = mat24.syndrome(a1, 24)
                #print(hex(special))
                a2 = a & ~special
                w, a2_list = mat24.vect_to_bit_list(a2)
                #print("a2_list", a2_list[:4])
                assert w == 4
                src = a2_list[:4] + syn_list[:2]
                src.append(mat24.lsbit24(a & special))
                #print(src)
                v = apply_perm_heptad(v, src, RED_IMG_CASE_44, result)
                exp = map_to_standard12_xi_exp(v)
                assert exp >= 0
        elif vtype == 0x46:
            w, dodecad = mat24.vect_to_bit_list(gc)
            perm = mat24.perm_from_dodecads(dodecad, RED_STD_DODECAD) 
            v = apply_perm(v, perm, result) 
            coc = mat24.ploop_theta(v >> 12) ^ v
            scalar = mat24.scalar_prod(STD_OCTAD, coc)
            exp = 2 - scalar
            # map to a value v with  gen_leech2_type(v) == 0x44 
        elif vtype == 0x43:
            cv = mat24.cocode_syndrome(coc, 24)
            w, cv_list = mat24.vect_to_bit_list(cv)
            assert w == 3, (w, cv_list)
            oct = sum([1 << i for i in cv_list[:5]])
            syn =  mat24.syndrome(oct, 24)
            oct ^= syn
            w, octad_list =  mat24.vect_to_bit_list(oct)
            assert w == 8
            w, rem_list = mat24.vect_to_bit_list(oct & ~cv)
            assert w == 5
            src = cv_list[:3] + rem_list[:3] + octad_list[8:9]
            v = apply_perm_heptad(v, src, RED_STD_HEPTAD, result)
            #print("after type 40", hex(v),  Cocode(v).syndrome(0))
            exp = 2 - ((v >> 22) & 1)
        else:
            raise ValueError("WTF")
        if exp: 
            exp = 0xE0000003 - exp
            v = gen_leech2_op_atom(v, exp)
            assert v & 0xfe000000 == 0
            #print(hex(v), "->",  "li_%d" % e,  "->", hex(res))
            result.append(exp)
    raise ValueError("WTF1")



#######################################################################
# Reduce a type-4 vector with a C function
#######################################################################


def reduce_type4_c(v):
    result = np.zeros(10, dtype = np.uint32)
    l = gen_leech2_reduce_type4(v, result)
    assert l >= 0, hex(l + 0x10000)
    return result[:l]



def conj_leech_word(v, g):
    return gen_leech2_op_word(v, g, len(g))

#######################################################################
# Create test data
#######################################################################


def reduce_type4_testdata(very_large = False):
    testdata = [
        (0x800, 0),
        (0x800, [0,1,2,3]),
        (0, [0,1,2,3]),
    ]
    for gc, cocode in testdata:
        v = leech_vector(gc, cocode )
    for i in range(10):
        coc = list(range(24))
        shuffle(coc)
        yield randint(0,1)*0x800000 + Cocode(coc[:4]).ord
    for i in range(100):
        so = SubOctad(randint(0, 759-1), randint(0, 64)) 
        # print(hex(so.gcode), hex(so.cocode), Cocode(so).syndrome_list(0))
        v = leech_vector(so.gcode ^ 0x800, so.cocode)
        assert gen_leech2_type(v) == 0x42
        yield v
    for i in range(100):
        so = SubOctad(randint(0, 759-1), randint(0, 64))
        oct = Octad(so.octad).bit_list
        more = set()
        while len(more) < 2:
            rnd = randint(0, 23)
            if rnd not in oct:
                more.add(rnd)
        #print("more", more)
        coc = Cocode(so) + Cocode(list(more)) 
        if randint(0,1):
            so = ~so 
        #print("\n",hex(so.gcode), oct, coc.syndrome_list(oct[0]))
        v = leech_vector(so.gcode, coc)
        assert gen_leech2_type(v) == 0x44,  gen_leech2_type(v) 
        yield v
    for i in range(100):
        p_found = False
        while not p_found:
            gc = GCode("r")
            p_found = len(gc) == 12
        odd = 0
        while not odd:
            coc = Cocode("e")
            odd = int(gc & coc)
        v = leech_vector(gc, coc)
        if gen_leech2_type(v) == 0x36:
            v ^= 0x800000
        assert gen_leech2_type(v) == 0x46,  gen_leech2_type(v) 
        yield v
    for i in range(100):
        l = list(range(24))
        shuffle(l)
        v = leech_vector(randint(0, 0xfff), l[:3])
        if gen_leech2_type(v) != 0x43:
             v ^= 0x800000
        assert gen_leech2_type(v) == 0x43,  gen_leech2_type(v) 
        yield v
    if very_large:
        for v in range(0x1000000):
            if gen_leech2_type(v) & 0xf0 == 0x40:
                yield v


#######################################################################
# The main test
#######################################################################


def _test(verbose = 0, c_function = 1):
    for i, v in enumerate(reduce_type4_testdata(True)):
        #print("")
        if verbose:
            print( hex(gen_leech2_type(v)), hex(v))
        result = reduce_type4(v)
        assert conj_leech_word(v, result) & 0xffffff == 0x800000
        if verbose:
            print([hex(x) for x in result])
        if c_function:
            result_c = reduce_type4_c(v)
            if verbose:
                print([hex(x) for x in result_c])
            assert (result == result_c).all(), (result, result_c)
        if i & 0xffff == 0xffff:
             print(".", end = "", flush = True)
        #print(hex(v), [hex(x) for x in result])
           
#_test()


