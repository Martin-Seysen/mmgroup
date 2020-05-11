from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from functools import reduce
from operator import __xor__

import pytest

from random import randint, shuffle

from mmgroup import mat24
mat24fast = mat24


from mmgroup.bitfunctions import bit_mat_mul, v2
from mmgroup.dev.mat24.mat24theta import theta_to_basis_vector


from mmgroup.dev.mat24.mat24_ref import Mat24


#####################################################################
# Test computation of cocycle theta on basis vectors
#####################################################################

@pytest.mark.mat24
@pytest.mark.parametrize("gc", [Mat24, mat24fast])
def test_theta_grey(gc):
    """ yet to be documented!!!! """
    for v1 in range(0, 0x1000, 0x400):
        for v2 in range(0x10):
            v0 = v1 + v2
            v = gc.gcode_to_vect(v0)
            cc0 = gc.ploop_theta(v0) 
            cc = gc.cocode_to_vect(cc0)
            ccref = theta_to_basis_vector(v)
            assert cc == ccref 


@pytest.mark.mat24
@pytest.mark.parametrize("gc", [Mat24, mat24fast])
def test_theta_colored(gc):
    """ yet to be documented!!!! """
    for v0 in range(0,0x400,0x10):
        v = gc.gcode_to_vect(v0)
        cc0 = gc.ploop_theta(v0)   & 0xc0f
        cc = gc.cocode_to_vect(cc0)
        ccref = ((v | (v >> 1) | (v >> 2)) & 0x111111)
        assert cc == ccref , map(hex, [v,cc,ccref])


#####################################################################
# Test functions related to the multiplication in the Parker loop
#####################################################################


@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref", [(mat24fast, Mat24)])
def test_Parker_loop(gc, ref):
    if gc is None: gc = Mat24
    for i in range(300):
        # test power map
        v1 = randint(0,0x1fff)
        c1 = gc.gcode_to_vect(v1)
        coc = gc.ploop_cocycle(v1,v1) 
        assert coc == (gc.bw24(c1) >> 2) & 1, map(hex, [i,v1, c1, coc])
        sq = gc.mul_ploop(v1,v1)
        assert sq == coc << 12, map(hex, [i, v1, c1, sq])
  
        # test commutator
        v2 = randint(0,0x1fff)
        c2 = gc.gcode_to_vect(v2)
        coc = gc.ploop_cocycle(v1,v2) 
        assert coc == gc.bw24(gc.ploop_theta(v1) & v2 & 0xfff) & 1
        comm = coc ^ gc.ploop_cocycle(v2,v1)
        assert comm == (gc.bw24(c1 & c2) >> 1) & 1, map(hex, [i,v1,c1,v2,c2,comm])
        assert comm == gc.ploop_comm(v1, v2) 
        pr = gc.mul_ploop(v1,v2)
        assert pr == v1 ^ v2 ^ (coc << 12), map(hex,[i,v1,c1,v2,c2,coc])  

        # test associator
        v3 = randint(0,0x1fff)
        c3 = gc.gcode_to_vect(v3)
        pr1 = gc.mul_ploop(pr,v3)
        pr2 = gc.mul_ploop(v1, gc.mul_ploop(v2,v3))
        computed_assoc = (pr1 ^ pr2) >> 12
        ref_assoc = gc.bw24(c1 & c2 & c3) & 1
        assert computed_assoc == ref_assoc
        other_assoc = gc.ploop_assoc(v1, v2, v3)
        assert other_assoc == ref_assoc

        if ref: 
            assert coc == ref.ploop_cocycle(v1,v2) 
            assert pr == ref.mul_ploop(v1,v2)
            assert gc.ploop_theta(v1) == ref.ploop_theta(v1)

        # test intersection
        cap = gc.ploop_cap(v1, v2)
        ref_cap = gc.vect_to_cocode(c1 & c2)
        assert cap == ref_cap

        # test exponentiation
        v4 = randint(0,0x1fff)
        v4_inv = gc.pow_ploop(v4,3)
        assert gc.mul_ploop(v4, v4_inv) == 0, map(hex, [v4, v4_inv, gc.mul_ploop(v4, v4_inv)])
        assert gc.pow_ploop(v4,4) == 0
        c4 = gc.gcode_to_vect(v4)
        assert gc.pow_ploop(v4,2) == (gc.bw24(c4) & 4) << 10
        assert gc.pow_ploop(v4,5) == v4
  
    print( "Parker Loop test passed" )



#####################################################################
# Test function mat24.ploop_solve
#####################################################################


def rand_ploop_solve_array(length):
    c = randint(0, 0xfff)
    a = [(1 << i) ^ ((c >> i) << 12) for i in range(12)]
    for i in range(150):
        j1, j2 = randint(0,11), randint(0,11)
        if j1 != j2:
            a[j1] ^= a[j2]
    shuffle(a)
    for i in range(12, length):
        c = randint(0, 0xfff) 
        v = reduce(__xor__, 
             ((a[i] & -((c>>i) & 1)) for i in range(12)))
        a.append(v)
    return a[:length]



def rand_ploop_solve_bad_array(length):
    a = [0x1000]
    while len(a) < length:
        a.append(randint(0, 0x1fff))
    c = randint(0, (1 << len(a)) - 1)
    for i in range(1, len(a)):
        a[0] ^= a[i] & -((c >> i) & 1)
    shuffle(a)
    return a



def ploop_solve_testcases():
    for length in list(range(14)) + [80, 161]:
        yield rand_ploop_solve_array(length)

def ploop_solve_bad_testcases():
    for length in list(range(3)) + [10, 11, 12, 13, 80]:
        yield rand_ploop_solve_bad_array(length)



@pytest.mark.mat24
@pytest.mark.parametrize("gc", [mat24fast])
def test_ploop_solve(gc):
    for m in ploop_solve_testcases():
        c = gc.ploop_solve(m)
        for v in m:
            assert gc.scalar_prod(c & 0xfff, v) == (v >> 12) & 1
    for m in ploop_solve_bad_testcases():
        with pytest.raises(ValueError):
            gc.ploop_solve(m)

     


