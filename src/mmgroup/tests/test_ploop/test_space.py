from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint, shuffle, sample

from mmgroup import mat24
from mmgroup import GCode, Cocode, PLoop, PLoopZ
from mmgroup import Octad, SubOctad, GcVector 
from mmgroup import Parity
from mmgroup.mat24 import MAT24_ORDER

from mmgroup.tests.spaces.spaces import MMTestSpace 
from mmgroup.mm_space import characteristics
from mmgroup import MM0, MMV

space = MMV
group = MM0

#####################################################################
# Test tags A[i, i] and D
#####################################################################


@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [10])
def test_AD(n_cases):
    for p in characteristics():
        sp = space(p)
        for _ in range(n_cases):
            i0 = randint(0, 23)
            scalar = randint(1, p-1) + p * randint(-5, 5)
            c0 = Cocode([i0])
            v0 = GcVector([i0])
            ref = sp('A', i0, i0)
            for j0 in [i0, c0, v0]:
                for j1 in [i0, c0, v0]:
                    assert sp('A', j0, j1) ==  ref
                    mmv = sp()
                    mmv['A', j0, j1] =  scalar
                    assert mmv == scalar * ref
                    assert mmv['A', j0, j1] == scalar % p
                d1 = sp('D', j0)
                assert sp('D', j0) ==  ref
                mmv = sp()
                mmv['D', j0] =  scalar
                assert mmv == scalar * ref
                assert mmv['D', j0] == scalar % p



#####################################################################
# Test tags A, B, and C 
#####################################################################


@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [20])
def test_ABC(n_cases):
    for p in characteristics():
        sp = space(p)
        for _ in range(n_cases):
            i1 = i0 = randint(0, 23)
            while i1 == i0:
                i1 = randint(0, 23)
            scalar = randint(1, p-1) + p * randint(-5, 5)
            for tag in "ABC":
                c0 = Cocode([i0])
                c1 = Cocode([i1])
                v0 = GcVector([i0])
                v1 = GcVector([i1])
                ref = sp(tag, i0, i1)
                for j0 in [i0, c0, v0]:
                    for j1 in [i1, c1, v1]:
                        assert sp(tag, j0, j1) ==  ref
                        mmv = sp()
                        mmv[tag, j0, j1] =  scalar
                        assert mmv == scalar * ref
                        assert mmv[tag, j0, j1] == scalar % p
                    

#####################################################################
# Test tag T
#####################################################################


@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [20])
def test_T(n_cases):
    for p in characteristics():
        sp = space(p)
        for _ in range(n_cases):            
            o = randint(0, 758)
            so = randint(0, 63)
            coc = Cocode(SubOctad(o, so))
            v = GcVector (coc.syndrome(randint(0,23)))
            sgn = randint(0, 1)
            pl = PLoopZ(sgn) * Octad(o)
            gc = GCode(pl)
            scalar = randint(1, p-1) + p * randint(-5, 5)
            assert gc == GCode(pl)
            ref = (-1)**sgn * sp('T', o, so)
            for sign, d in ((sgn,o), (sgn,gc), (0, pl)):
                for par2 in (so, coc, v):
                    assert (-1)**sign * sp('T', d, par2) == ref
                    mmv = sp()
                    mmv['T', d, par2] = (-1)**sign * scalar
                    assert mmv == scalar * ref
                    signed_scalar = (-1)**sign * scalar % p
                    assert mmv['T', d, par2] == signed_scalar


#####################################################################
# Test tags X, Y, and Z
#####################################################################


@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [20])
def test_XYZ(n_cases):
    for p in characteristics():
        sp = space(p)
        for _ in range(n_cases):
            d = randint(0, 0x1fff)
            i = randint(0, 23)
            c = Cocode([i])
            v = GcVector([i])
            pl = PLoop(d)
            gc = GCode(d & 0xfff)
            d0 = d & 0x7ff
            scalar = randint(1, p-1) + p * randint(-5, 5)
            for tag in "XYZ":
                s2 = s1 = d >> 12
                if tag == "Y": s1 ^= d >> 11
                ref = (-1)**s1 * sp(tag, d0, i)
                for sign, dd in ((s1, d0), (s2, gc), (0, pl)):
                    for j in [i, c, v]:
                        assert (-1)**sign * sp(tag, dd, j) == ref
                        mmv = sp()
                        mmv[tag, dd, j] = (-1)**sign * scalar
                        assert mmv == scalar * ref
                        signed_scalar = (-1)**sign * scalar % p
                        assert mmv[tag, dd, j] == signed_scalar



 