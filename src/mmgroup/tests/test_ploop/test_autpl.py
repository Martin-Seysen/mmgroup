from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint, shuffle, sample

from mmgroup import mat24
from mmgroup import GCode, Cocode, PLoop, PLoopZ
from mmgroup import Octad, GcVector 
from mmgroup import AutPL
from mmgroup import Parity




#####################################################################
# Generate random elements of AutPL
#####################################################################



def rand_u7():
    """Generate a random umbral heptad"""
    o = Octad(sample(range(24), 5))
    return sample(o.bit_list[:], 6) + sample(((~o).bit_list)[:], 1)

def rand_autpl_u7():
    """Generate random Permutation in Mat24 from umbral heptads"""
    return AutPL(0, zip(rand_u7(), rand_u7()), Cocode("r"))



#####################################################################
# Test generation of permutations
#####################################################################



def perm_mapping(length):
    """Create a mapping defining a unique permutation p in Mat24.

    Return a triple (src, dest, perm_num) with src, dest being
    lists of the given 'length' satisfying:

        p[src[i]] = dest[i] for 0 <= i < length .    

    for a unique permutation p in Mat24. p has the internal number 
    perm_num. 7 <= length < 24 must hold. The function does not 
    use any internal numbering or random process of the mmgroup 
    package for creating the permutation p.
    """
    p  = rand_autpl_u7()  # p is a random permutation in AutPL
    perm = p.perm
    assert length >= 7
    while True:
        src = sample(range(24), length)
        if length >= 9 or len(Cocode(src)) >= 2:
            # The the permutation is determined by its effect on src
            return  src, [perm[i] for i in src], p.perm_num
        



def create_perm_testvectors():
    for i in range(40):
        for length in range(7, 25):
            yield perm_mapping(length)





@pytest.mark.ploop
def test_create_group_permutation():
    for i, (src, dest, num) in enumerate(create_perm_testvectors()):
        p1 = AutPL(0, zip(src, dest))
        p1.check()
        p1c = p1.copy().check()
        assert p1c == p1
        assert p1.perm_num == num
        assert p1 == AutPL(0, dict(zip(src, dest))).check()
        assert p1 == AutPL(0, num).check()
        assert p1 == AutPL(p1).check()
        assert p1.cocode == 0
        perm = p1.perm
        assert isinstance(perm, list)
        assert min(perm) == 0 and max(perm) == 23
        for x, y in zip(src, dest):
            assert perm[x] == y
        assert p1 == AutPL(0, perm)
        coc = Cocode(randint(1,0xfff))
        p2 = AutPL(coc, 0) * AutPL(0, perm)
        p2.check()
        assert p2 == (AutPL(coc) * AutPL(0, perm)).check()
        assert p2 == AutPL(coc, perm).check()
        assert p2 != p1
        assert p2.perm_num == p1.perm_num
        assert p2.cocode == coc.cocode
        p2c = AutPL(0, perm) * AutPL(coc, 0).check()
        assert p2.perm_num == p1.perm_num
        assert Cocode(p2c.cocode) * p2 == coc
        assert p2 == AutPL(Cocode(p2.cocode), 0) * AutPL(0, p2.perm)


#####################################################################
# Test Group multiplication and inverse
#####################################################################

@pytest.mark.ploop
def test_group_op(n_cases = 100):
    print("")
    for i in range(n_cases):
        p1 = AutPL('r', 'r')
        p2 = rand_autpl_u7() 
        p2 *= AutPL(Cocode(randint(0, 0xfff)))
        p12 = (p1 * p2).check()
        assert p12.rep == mat24.mul_autpl(p1.rep, p2.rep)
        if i < 1:
            print(p1, "*", p2, "=", p12)
        p3 =  AutPL('r', 'r')
        assert (p1 * p2) * p3 == p1 * (p2 * p3)
        p1i = (p1 ** -1).check()
        assert p1 * p1i == p1i * p1 == AutPL()
        

#####################################################################
# Test Group operation on Parker loop
#####################################################################


def autpl_testwords(n_cases = 50):
    for i in range(n_cases):
        if i % 4 == 0:
            yield rand_autpl_u7()
        else:
            yield AutPL('r', 'r') 


def ploop_testvectors(n_cases = 50):
    for i in range(n_cases):
        yield PLoop(randint(0, 0x1fff))



@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [100])
def test_group_ploop(n_cases):
    print("")
    for i, (p, g) in enumerate(zip(ploop_testvectors(n_cases),
            autpl_testwords(n_cases))):
        pg = p * g
        if i < 1:
            print(p, "*", g, "=", pg)
        pg_ref = PLoop(mat24.op_ploop_autpl(p.ord, g.rep)) 
        assert pg == pg_ref
        assert p/4 == pg/4 == parity(len(p)/4)


#####################################################################
# Test Group operation on Cocode and on vector space
#####################################################################

def cocode_testvectors(n_cases = 50):
    for i in range(n_cases):
        yield Cocode(randint(0, 0xfff))

def get_testvectors(n_cases = 50):
    for i in range(n_cases):
        yield GcVector(randint(0, 0xffffff))


@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [100])
def test_group_ploop(n_cases):
    print("")
    for i, (c, v, g) in enumerate(zip(cocode_testvectors(n_cases),
           get_testvectors(n_cases), autpl_testwords(n_cases))):
        cg = c * g
        vg = v * g
        if i < 1:
            print(c, "*", g, "=", cg)
            print(v, "*", g, "=", vg)
        cg_ref = Cocode(mat24.op_cocode_perm(c.ord, g.perm)) 
        assert cg == cg_ref
        assert len(cg) == len(c)
        vg_ref = GcVector(mat24.op_vect_perm(v.ord, g.perm)) 
        assert vg == vg_ref
        assert len(vg) == len(v)
        if len(v) & 1 == 0:
            assert v/2 == vg/2 == Parity(len(v) >> 1)
        if len(v) & 3 == 0:
            assert v/4 == vg/4 == Parity(len(v) >> 2)
        assert v % 2 == vg % 2 == Parity(len(v))




#####################################################################
# Test construction of group words with permutations and cocode
#####################################################################


@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [100])
def test_group_from_perm(n_cases):
    for i in range(n_cases):
        h1 = rand_u7()
        h2 = rand_u7()
        autp = AutPL(0, zip(h1, h2))
        assert autp == AutPL(0, dict(zip(h1, h2)))
        assert autp.perm == mat24.perm_from_heptads(h1, h2)
        assert autp.cocode == 0
        perm_num = autp.perm_num
        assert perm_num == mat24.perm_to_m24num(autp.perm)
        assert autp == AutPL(0, mat24.perm_from_heptads(h1, h2))
        assert autp == AutPL(autp)
        assert autp == AutPL(0, autp.perm_num)
        assert autp == AutPL(0, zip(h1, h2))
        coc_num = randint(1, 0xfff)
        coc = Cocode(coc_num)
        assert coc != Cocode(0)
        assert coc.cocode == coc_num
        im_coc = coc * autp
        assert type(im_coc) == type(coc)
        assert AutPL(im_coc) ==  AutPL(coc)**autp
        aut_cp = AutPL(coc) * autp
        assert aut_cp == AutPL(coc_num, perm_num) 
        if coc_num and perm_num:
            assert aut_cp.as_tuples() == [('d', coc_num), ('p', perm_num)]        
        assert autp * AutPL(im_coc) ==  aut_cp
        assert type(aut_cp) == type(autp)
        assert aut_cp.perm_num == perm_num 
        assert aut_cp.cocode == coc_num 
        assert Parity(aut_cp) == Parity(coc) 
        assert aut_cp.parity == coc.parity  == Parity(coc).ord
        assert autp == AutPL() * autp == autp *  AutPL()
        with pytest.raises(TypeError):
            autp * coc
        with pytest.raises(TypeError):
            autp * Parity(randint(0,9))
        with pytest.raises(TypeError):
            autp * randint(2,9)
        with pytest.raises(TypeError):
            randint(2,9) * autp
        with pytest.raises(TypeError):
            autp * PLoop(randint(0, 0x1fff)) 



#####################################################################
# Test generation of random elements
#####################################################################

        
@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [50])
def test_rand_group_word(n_cases):
    for i in range(n_cases):
        aut_cp = AutPL('r', 'r') 
        perm_num = aut_cp.perm_num
        coc_num = aut_cp.cocode 
        assert aut_cp == AutPL(coc_num, perm_num) 
        assert 0 <= perm_num < mat24.MAT24_ORDER
        assert 0 <= coc_num < 0x1000
        aut_cpe = AutPL('e', 'r') 
        assert aut_cpe.parity == 0
        aut_cpe = AutPL('o', 'r')
        assert aut_cpe.parity == 1



        