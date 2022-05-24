


from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint 
from functools import reduce

import numpy as np
import pytest

from mmgroup import XLeech2 
from mmgroup import mat24
from mmgroup.generators import gen_leech2_subtype
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_reduce_n
from mmgroup.generators import gen_leech2_reduce_n_rep

from mmgroup.dev.generators.gen_leech_reduce_n import MAP_VECTOR
from mmgroup.dev.generators.gen_leech_reduce_n import TABLE_OCTAD
from mmgroup.dev.generators.gen_leech_reduce_n import TABLE_DODECAD
from mmgroup.dev.generators.gen_leech_reduce_n import STD_DODECAD
from mmgroup.dev.generators.gen_leech_reduce_n import CPL_DODECAD

from mmgroup.tests.test_gen_xi.test_leech2_type import rand_n_elem
from mmgroup.tests.test_gen_xi.test_reduce import apply_perm
from mmgroup.tests.test_gen_xi.test_reduce import find_octad_permutation


#####################################################################
# Auxiliary functions
#####################################################################

OCTAD_PLUS = [0,1,2,3,4, 8,9]


def find_octad_permutation_odd(v, result, verbose = 0):
    """ Find a suitable permutation for an octad.

    Similar to function ``find_octad_permutation`` in module
    ``mmgroup.dev.generators.gen_leech_reduce_n``.
    Here ``v, o, c`` are as in that function; but the scalar
    product of ``o`` and ``c`` must be 1. Apart from that
    operation is as in function ``find_octad_permutation``.

    We compute a permutation that maps octad ``o`` to the standard
    octad (0,1,2,3,4,5,6,7). If the cocode part ``c`` of ``v`` is 
    not a suboctad of octad ``o`` then we map (one shortest 
    representative of) ``c`` into the set (0,1,2,3,...7,8). 
    """
    coc = (v ^ mat24.ploop_theta(v >> 12)) & 0xfff
    w = mat24.gcode_weight(v >> 12) 
    vect = mat24.gcode_to_vect((v ^ ((w & 4) << 21)) >> 12)
    src = mat24.vect_to_list(vect, 5)
    if mat24.cocode_weight(coc) == 4:
        sextet = mat24.cocode_to_sextet(coc) 
        for i in range(0, 24, 4):
            syn = (1 << sextet[i]) | (1 << sextet[i+1])
            syn |= (1 << sextet[i+2]) | (1 << sextet[i+3])
            special = syn & vect
            if special & (special - 1): 
                break
    else:
        syn = mat24.cocode_syndrome(coc, 24)
    src.append(mat24.lsbit24(syn & ~vect))  
    return apply_perm(v, src, OCTAD_PLUS, 6, result, verbose)


def op_y_x(v, table, out):
    v = gen_leech2_op_atom(v, out[0])
    c = (v ^ table[0]) & 0x7ff
    y = 0
    for i in range(11):
        y ^=  table[i + 1] & -((c >> i) & 1)
    out.append(0xC0000000 + y)
    v = gen_leech2_op_atom(v, out[1])
    if (((v >> 23) ^ table[12]) & 1):
        y =  table[13]
        v = gen_leech2_op_atom(v, 0xC0000000 + y)
        out[1] ^= y         
    c = (v ^ table[0]) & 0x7ff
    x = table[14] if (v & 0x1000000) else 0
    out.append(0x90000000 + x) #!!!







def map_vector(v):
    return MAP_VECTOR[gen_leech2_subtype(v)] if v else 0


#####################################################################
# Python implementation of function to be tested
#####################################################################


def ref_gen_leech2_reduce_n(v, verbose = 0):
    vtype = gen_leech2_subtype(v);
    subtype = vtype & 0xf;

    out = []

    if subtype & 1:
        coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff;
        syn = mat24.cocode_syndrome(coc)
        src = mat24.vect_to_list(syn, mat24.bw24(syn))
        assert len(src)  in [1,3]
        lst = [1,2,3] if subtype == 3 else [0]
        apply_perm(v, src, lst, len(src), out)
        v = gen_leech2_op_atom(v, out[0])
        out.append(0xC0000000 + ((v >> 12) & 0x7ff))
        v = gen_leech2_op_atom(v, out[1])
        out.append(0xB0000000 + ((v >> 13) & 0x800))
    elif subtype == 6:
        gv = (v >> 12) & 0xfff 
        vect = mat24.gcode_to_vect(gv)
        src = mat24.vect_to_list(vect, mat24.bw24(vect))
        assert len(src) == 12
        dest = STD_DODECAD
        if (vtype == 0x36):
            coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff;
            w =  mat24.bw24(mat24.cocode_as_subdodecad(coc, gv))
            if w & 2:
               dest = CPL_DODECAD 
        pi = mat24.perm_from_dodecads(dest, src)
        out.append(0xA0000000 + mat24.perm_to_m24num(pi))
        op_y_x(v, TABLE_DODECAD, out)
    elif subtype in [2,4]:
        if vtype == 0x34:
            find_octad_permutation_odd(v, out)
        else:
            find_octad_permutation(v, out)
        op_y_x(v, TABLE_OCTAD, out)
    elif subtype in [0,8]:
        if ((v & 0x7ff) == 0):
            out.append(0xA0000000)
            out.append(0xC0000000)
            x = 0x800 if v & 0x1800000 == 0x1800000 else 0
            out.append(0x90000000 + x)
        else:
            syn = mat24.cocode_syndrome(v & 0x7ff, 0)
            src = mat24.vect_to_list(syn, mat24.bw24(syn))
            j = mat24.bw24(syn) & 2;
            lst, y0  = ([2,3], 0x200) if j else ([0,1,2,3], 0x400)
            apply_perm(v, src, lst, len(lst), out);
            v = gen_leech2_op_atom(v, out[0])
            y = y0 if v & 0x800000 else 0
            out.append(0xC0000000 + y)
            v = gen_leech2_op_atom(v, out[1])
            x = y0 if v & 0x1000000 else 0
            out.append(0xB0000000 + x)
    else:
        raise ValueError("Bad subtype " + hex(vtype))
    assert len(out) == 3
    return vtype, np.array(out, dtype = np.uint32)






#####################################################################
# Test data for testing function gen_leech2_reduce_n
#####################################################################


def reduce_n_testdata(ntests = 100):
    """The function yields test vectors in the Leech lattice mod 2

    These test vectors are given in **Leech lattice encoding**;
    they are used for testing function ``gen_leech2_reduce_n``.
    """
    for c in [0,0x200]:
        for d in range(0, 0x2000000, 0x800000):
             yield c + d
    testdata = [
    ]
    for t in testdata:
         yield t
    for i in range(ntests):
        v = randint(0, 0x7ff)
        o = randint(0, 3)
        yield (o << 23) + v
    for v in MAP_VECTOR.values():
        for i in range(ntests):
            g =  rand_n_elem()
            yield g.xsp_conjugate(v)
    for type_ in (2,3,4):      
        for i in range(ntests):
            yield XLeech2('r', type_).ord

#####################################################################
#####################################################################
# Testing function gen_leech2_reduce_n
#####################################################################
#####################################################################




#####################################################################
# Test transversal in function gen_leech2_reduce_n
#####################################################################



@pytest.mark.gen_xi
def test_std_subtypes():
    for subtype, v in [(0,0)] + list(MAP_VECTOR.items()):
        subtype1 =  gen_leech2_subtype(v)
        if subtype1 != subtype:
            print("Subtype obtained: %s, expected: %s" % 
                (hex(subtype1), hex(subtype)))
            print("gcode", gcode)
            print("cocode", cocode)
            err = "Error in transversal of N_x0"
            raise ValueError(err)
    for subtype, v in list(MAP_VECTOR.items()):
        v_c = gen_leech2_reduce_n_rep(subtype)
        assert v == v_c


#####################################################################
# Test function gen_leech2_reduce_n
#####################################################################

@pytest.mark.gen_xi
def test_reduce_n(ntests = 100, verbose = 0):
    """Test function ``gen_leech2_reduce_n`` """
    for n, v in enumerate(reduce_n_testdata(ntests)):
        ref_subtype = gen_leech2_subtype(v)
        if verbose:
            print(" \nTest %d, v = %s, subtype  = %s" % 
                (n+1, hex(v), hex(ref_subtype))
            )
        subtype, op = ref_gen_leech2_reduce_n(v)
        ok = True 
        err = "Unknown error"

        ok1 =  ref_subtype == subtype 
        if ok and not ok1:
            err  = "Function gen_leech2_reduce_n returns wrong subtype"
        ok &= ok1

        v_mapped = gen_leech2_op_word(v, op, len(op))
        ok2 = v_mapped == map_vector(v)
        if ok and not ok2:
            err  = "Wrong mapping in function gen_leech2_reduce_n"
        ok &= ok2

        ok3 = len(op) == 3
        tags = [x >> 28 for x in op]
        ok3 &= tags[:2] == [0xa, 0xc] and tags[2] in [9, 0xb]
        if ok and not ok3:
            err  = "Error in computed operation on vector"
        ok &= ok3

        if  (op[2] & 0xf0000000) == 0xB and v != 0x18000000:
            ok4 = op[2] & 0x800 == 0
        else:
            ok4 = True
        if ok and not ok4:
            err  = "Computed operation on vector must be even"
        ok &= ok4

        op_c = np.zeros(3, dtype = np.uint32)
        subtype_c = gen_leech2_reduce_n(v, op_c) 
        #op_c[1] &= ~2 # Produce an error for verifying the test
        ok_C_function = subtype_c == subtype
        ok_C_function &= (op_c == op).all()
        if ok and not ok_C_function:
            err  = "Error in C function"
        ok &= ok_C_function
    

        if verbose or not ok:
            if not ok:
                print(" \nTest %d, v = %s, subtype  = %s" % 
                    (n+1, hex(v), hex(ref_subtype))
                )
            print("transformation", op)
            print("map expected:", hex(map_vector(v)), ", obtained",
                hex(v_mapped) )
            if not ok_C_function:
                 print("Subtype from C function is %s, expected %s" %
                     (hex(subtype_c), hex(subtype))
                 )
                 print("Operation from C function:", op_c)
                 print("Expected:                 ", op)
            if not ok:
                raise ValueError(err)







