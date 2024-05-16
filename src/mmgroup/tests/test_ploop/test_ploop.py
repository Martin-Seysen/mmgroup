from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint, shuffle, sample

from mmgroup import mat24
from mmgroup import GCode, Cocode, PLoop, PLoopZ
from mmgroup import Octad, SubOctad, GcVector 
from mmgroup import AutPL
from mmgroup import Parity


#####################################################################
# Test class Parity
#####################################################################


@pytest.mark.ploop
def test_parity():
    odd, even = Parity(3), Parity(-4)
    assert odd != even
    assert Parity(-311) == odd
    assert int(odd) == 1
    assert int(even) == 0
    assert (-1)**odd == -1
    assert (-1)**even == 1
    assert 1**odd == 1**even == 1
    assert odd + 1 == odd - 1 == 1 + odd == odd + odd == 7 - odd == even
    assert odd * odd == -odd == odd
    assert odd * even == even * odd == 2 * odd == even * 3 == even
    print("\nPossible parities are:", odd,",", even)
    with pytest.raises(ValueError):
        2**even
    with pytest.raises(TypeError):
        even & even
    cc1, cc0 = Cocode(0x823),  Cocode(0x7ef)
    assert cc1 + even == even + cc1 == cc0 + odd == odd + cc0 == odd
    assert cc1 + odd == odd + cc1 == cc0 + even == even + cc0 == even
    assert odd * cc1 == cc1 * odd == cc1 
    ccn = Cocode(0)
    assert even * cc1 == cc1 * even == even * cc0 == cc0 * even == ccn
    pl = PLoop(0x1234)
    gc = GCode(pl)
    assert even * pl == even * gc == GCode(0)
    assert odd * pl == odd * gc == gc
    aut =  AutPL('e', 'r')
    assert odd * aut == odd
    assert even * aut == even

#####################################################################
# Test Parker loop
#####################################################################


def ploop_creation_testvectors(n_cases = 100):
    i = 0
    while i < n_cases:
        v = randint(0, 0xffffff)
        try:
            syn = mat24.syndrome(v)
            sign = randint(0, 1)
            i += 1
            yield sign, v, syn
        except ValueError:
            pass




@pytest.mark.ploop
def test_create_Parker_loop():
    print("\nTesting creation of Parker loop elements")
    for sign, v, syn in ploop_creation_testvectors():
        # Test power map, inveriosn, sign, theta, and conversion to GCode
        vlist = [i for i in range(24) if (v >> i) & 1]
        shuffle(vlist)
        p1 = (-1)**sign * PLoop(vlist)
        ploop = mat24.vect_to_gcode(v ^ syn) + (sign << 12)
        p2 = PLoop(ploop)
        assert p1 == p2
        g1 = GCode(ploop)
        p3 = (-1)**sign * PLoop(g1)
        assert p1 == p3
        assert p1 == PLoop(p3)
        assert abs(p1) == PLoop(GCode(vlist))
        assert abs(p1) == PLoop(GcVector(vlist))
    print( "Parker Loop creation test passed" )



@pytest.mark.ploop
def test_Parker_loop():
    print("\nTesting operation on Parker loop")
    for i in range(200):
        # Test power map, inveriosn, sign, theta, and conversion to GCode
        n1 = randint(0,0x1fff)
        p1 = PLoop(n1)
        assert p1.ord == n1 == p1.gcode + 0x1000 * (1 - p1.sign) / 2
        if (i < 2):
            print("Testing Parker loop element", p1, ", GCode = ", GCode(p1))
        assert len(p1) == mat24.gcode_weight(n1) << 2
        assert p1.theta() == Cocode(mat24.ploop_theta(n1))
        assert p1/4 == p1.theta(p1) == Parity(mat24.ploop_cocycle(n1,n1))
        assert p1**2 == PLoopZ(p1/4) == (-PLoopZ())**(p1/4)
        assert p1**(-5) == PLoopZ(p1/4) * p1 == (-1)**(p1/4) * p1 == 1/p1
        assert (1/p1).ord ^ n1 == (mat24.gcode_weight(n1) & 1) << 12
        assert -1 / p1 == -p1**3 == -(1/p1)
        assert p1 * (-1/p1) == PLoopZ(1) == -PLoopZ()
        assert abs(p1).ord == GCode(p1).ord == p1.ord & 0xfff
        assert p1 != GCode(p1)
        assert +p1 == 1 * p1 == p1 * 1 == p1
        assert p1 != -p1 and p1 != ~p1
        assert (-p1).ord == p1.ord ^ 0x1000
        assert (~p1).ord == p1.ord ^ 0x800
        s, o, p1_pos = p1.split()
        assert s == (p1.ord >> 12) & 1
        assert o == (p1.ord >> 11) & 1
        assert p1.sign == (-1)**s
        assert p1_pos.ord == p1.ord & 0x7ff
        assert PLoopZ(s, o) * p1_pos == p1
        assert PLoopZ(0, o) * p1_pos == abs(p1) 
        assert PLoopZ(s+1, o) * p1_pos == -p1 
        assert -p1 == -1 * p1 == p1 * -1 == p1 / -1 == -1 / p1**-5
        assert PLoopZ(s, 1+o) * p1_pos == ~p1
        assert PLoopZ(1+s, 1+o) * p1_pos == ~-p1 == -~p1 == ~p1 / -1
        assert 2*p1 == GCode(0) ==  -2 *  GCode(p1)
        assert -13*p1 == GCode(p1) == 7 * GCode(p1)
        if len(p1) & 7 == 0:
            assert p1**Parity(1) == p1
            assert p1**Parity(0) == PLoop(0)
        else:
            with pytest.raises(ValueError):
                p1**Parity(randint(0,1))
        assert p1.bit_list == mat24.gcode_to_bit_list(p1.value & 0xfff)
        assert p1.bit_list == GcVector(p1).bit_list
        assert p1.bit_list == GCode(p1).bit_list
        assert PLoop(GcVector(p1)+0) == PLoop(GCode(p1)+0) == abs(p1)
        assert p1 + 0 == 0 + p1 == GCode(p1) + 0 == 0 + GCode(p1)
        assert Cocode(GcVector(p1)) == Cocode(0)
        assert p1/2 == Parity(0)

        # Test Parker loop multiplication and commutator
        n2 = randint(0,0x1fff)
        p2 = PLoop(n2)
        coc = Cocode(mat24.ploop_cap(n1, n2))
        if (i < 1):
            print("Intersection with", p2, "is", coc)
        p2inv = p2**-1
        assert p1 * p2 == PLoop(mat24.mul_ploop(p1.ord, p2.ord))
        assert p1 / p2 == p1 * p2**(-1)
        assert p1 + p2 == p1 - p2 == GCode(p1 * p2) == GCode(n1 ^ n2)
        assert (p1 * p2) / (p2 * p1) == PLoopZ((p1 & p2)/2)
        assert p1 & p2 == coc
        assert p1.theta() == Cocode(mat24.ploop_theta(p1.ord))
        assert p1.theta(p2) == Parity(mat24.ploop_cocycle(p1.ord, p2.ord))
        assert (p1 & p2)/2 == p1.theta(p2) + p2.theta(p1) 
        assert p1 & p2 == p1.theta() + p2.theta() + (p1 + p2).theta()
        assert int((p1 & p2)/2) == mat24.ploop_comm(p1.ord, p2.ord) 
        assert GcVector(p1 & p2) ==  GcVector(p1) &  GcVector(p2)      
        assert ~GcVector(p1 & p2) ==  ~GcVector(p1) |  ~GcVector(p2) 
        assert Cocode(GcVector(p1 & p2)) == p1 & p2       

        # Test associator
        n3 = randint(0,0x1fff)
        p3 = PLoop(n3)
        assert p1 * p2 * p3 / (p1 * (p2 * p3)) == PLoopZ(p1 & p2 & p3)
        assert int(p1 & p2 & p3) == mat24.ploop_assoc(
                p1.ord, p2.ord, p3.ord)
        i = randint(-1000, 1000)
        par = Parity(i)
        s3 = ((p3 & p1) & p2) + par
        assert s3 ==  ((p3 & p1) & p2) + par 
        assert s3 ==  i + (p1 & (p2 & p3))
        assert s3  == par + (p1 & (p2 & p3))
  
        # Test some operations leading to a TypeError
        with pytest.raises(TypeError):
            p1 & p2 & p3 & p1
        with pytest.raises(TypeError):
            coc & coc
        with pytest.raises(TypeError):
            GCode(p1) * GCode(p2)
        with pytest.raises(TypeError):
            1 / GCode(p2)
        with pytest.raises(ValueError):
           coc / 4
        with pytest.raises(TypeError):
           p1 * coc 
        with pytest.raises(TypeError):
           coc * p1
        types = [GcVector, GCode, Cocode, PLoop]
        for type_ in types:
            with pytest.raises(TypeError):
                int(type_(0))

    print( "Parker Loop test passed" )




@pytest.mark.ploop
def test_cocode():
    print("")
    for i in range(200):
        # Test power map, inveriosn, sign, theta, and conversion to GCode
        n1 = randint(0,0x1fff)
        p1 = PLoop(n1)
        ccvector = randint(0, 0xffffff)
        coc = mat24.vect_to_cocode(ccvector)
        cclist =  [i for i in range(24) if (ccvector >> i) & 1]
        cc1 = Cocode(cclist)
        cc2 = Cocode(coc)
        if i < 1:
            print("\nTesting", GcVector(ccvector), ", cocode =", cc1)
        assert cc1 == cc2
        u = Parity(mat24.scalar_prod(p1.value, cc1.value))
        assert p1 & cc1 == u == cc1 & p1 == u * 1 == u + 0
        par = Parity(randint(0,1))
        assert cc1 + par == par + cc1 == cc1.value//0x800 + par
        assert cc1 % 2 == Parity(cc1)
        assert len(cc1) == mat24.cocode_weight(cc1.value)
        if len(cc1) < 4:
            syndrome = mat24.cocode_syndrome(cc1.value)
            assert cc1.syndrome().value ==  syndrome
            syn_from_list = sum(
                1 << i for i in GcVector(ccvector).syndrome_list())
            assert syn_from_list == syndrome
        i = randint(0,23)
        assert cc1.syndrome(i).value ==  mat24.cocode_syndrome(cc1.value, i)
        syndrome_list = cc1.syndrome(i).bit_list
        assert len(cc1) == len(syndrome_list)
        assert syndrome_list == mat24.cocode_to_bit_list(cc1.value, i)
        syn1 = [GcVector(x) for x in mat24.all_syndromes(ccvector)]
        syn2 = GcVector(ccvector).all_syndromes()
        assert cc1.all_syndromes() == syn1 == syn2

@pytest.mark.ploop
def test_octads():
    print("")
    OMEGA = ~PLoop(0)
    for i in range(200):
        no = randint(0, 758)
        o = Octad(no)
        assert o == PLoop(mat24.octad_to_gcode(no))
        assert o == Octad(sample(o.bit_list, 5))
        assert o.octad == no
        o_signbit, o_cpl = randint(0,1), randint(0,1)
        signed_o = PLoopZ(o_signbit, o_cpl) * o
        assert signed_o.octad == no
        assert signed_o.sign == (-1)**o_signbit
        assert o.gcode == mat24.octad_to_gcode(no)
        assert signed_o.gcode ==  (o * PLoopZ(0, o_cpl)).gcode
        assert signed_o.split_octad() == (o_signbit, o_cpl, o)
        nsub =  randint(0, 63)
        sub = (-1)**o_signbit * SubOctad(no, nsub)
        sub1 = (-1)**o_signbit * SubOctad(o, nsub)
        sub_weight = mat24.suboctad_weight(nsub)
        assert sub == sub1
        assert o == (-1)**o_signbit * Octad(sub) * OMEGA**sub_weight 
        assert sub.octad_number() == no
        ploop, cocode = sub.isplit()
        sign, gcode = (-1)** (ploop >> 12), ploop & 0xfff
        assert gcode == o.gcode ^ 0x800 * sub_weight  
        #assert sub.suboctad == nsub 
        assert sign == signed_o.sign == (-1)**o_signbit
        coc = mat24.suboctad_to_cocode(nsub, no)
        assert cocode == coc 
        assert Cocode(sub) == Cocode(coc)
        assert sub == SubOctad(sub.sign * o, Cocode(coc))
        assert sub ==  sub.sign  * SubOctad(no, Cocode(coc))

        o2 =  Octad(randint(0, 758))
        sub2 = SubOctad(o, o2)
        assert Cocode(sub2) == Cocode(o & o2)
        assert len(Cocode(sub2))//2 & 1 == int((o & o2)/2)

        t_sign, t_tag, t_i0, t_i1 = sub.vector_tuple()
        assert t_tag == 'T'
        assert t_sign == sign, (hex(sub.value), t_sign, o_signbit)
        assert t_i0 == no
        assert t_i1 == nsub



