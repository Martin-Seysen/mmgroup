from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint

from mmgroup import mat24
mat24fast = mat24


from mmgroup.bitfunctions import bit_mat_mul, v2
from mmgroup.dev.mat24.mat24theta import theta_to_basis_vector


from mmgroup.dev.mat24.mat24_ref import Mat24



@pytest.mark.mat24
def test_basis_is_same():
    for i in range(24):
        assert Mat24.basis[i] == mat24fast.basis[i]
        assert Mat24.recip_basis[i] == mat24fast.recip_basis[i]


@pytest.mark.mat24
@pytest.mark.parametrize("gc", [Mat24, mat24fast])
def test_basis(gc):
    def test_vects():
        for i in range(24): yield gc.basis[i]
        for i in range(80): yield randint(0,0xffffff)
    for v in test_vects():
        ve = gc.vect_to_vintern(v)
        assert ve == bit_mat_mul(v, gc.recip_basis),  map(hex,[i,v,ve])
        ved = gc.vintern_to_vect(ve)
        assert ved == bit_mat_mul(ve, gc.basis),  map(hex,[i,v,ve,ved])
        assert ved == v , map(hex,[i,v,ve,ved])


def lsbit24_testcases():
    for i in range(32):
        i2 = 1 << i
        yield i2
        yield i2-1
        yield i2+1
    yield 0xffffffff
    for i in range(100):
        yield randint(0,0xffffff)


@pytest.mark.mat24
@pytest.mark.parametrize("gc", [Mat24, mat24fast])
def test_lsbit24(gc):
    for x in lsbit24_testcases():
        try:
            v_ref = min(v2(x),24)
        except:
            v_ref = 24
        v = gc.lsbit24(x)              
        assert v == v_ref, (hex(x), v, v_ref)
    print( "test of lsbit determination passed" )


def spread_testcases():
    data = [ (7, 0, 0), (0x311111, 0x101ffe, 0x2e), (0xffffff,)*3 ]
    for d in data:
        yield d
    for k in range(20):
        for i in range(1,4):
            mask = 0xfffffff
            for j in range(i): mask &= randint(0, 0xffffff)
            x = randint(0, 0xffffff)
            l = [(x >> i) & 1 for i in range(24) if (mask >> i) & 1] 
            y = sum(z << i for i, z in enumerate(l))
            yield (mask, x, y)
 

@pytest.mark.mat24
@pytest.mark.parametrize("gc", [Mat24, mat24fast])
def test_spread(gc):
    for mask, x, y in  spread_testcases():
        res = gc.extract_b24(x, mask)
        assert res == y, map(hex, [mask, x, y, res])
        x &= mask
        res = gc.spread_b24(y, mask)
        assert res == x, map(hex, [mask, y, x, res])
    print( "test of 24 bit spread/extract passed" )




def syndrome_testcases():
    """yields test cases (codewordNo, cocode, t) for Mat24.syndrome

    codewordNo is 12-bit random number indicating a Golay code word
    cocode is a list of 0...4 bit position where to modify the code word
    t is an entry of cocode if len(cocode) = 4 and a random bit position
    (or 24) otherwise.
    """
    testcases = [
         (1, [1], 3), (0x400,[0],0), (0xd06,[0],0)
    ]
    for t in testcases: 
        yield t
    ll = [0,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4]
    for i in range(2000):
        codewordNo = randint(0, 0xfff)
        l = ll[randint(0,len(ll)-1)]
        cocode = []
        while len(cocode) < l:
            i = randint(0,23)
            if not i in cocode:   cocode.append(i)
        t = cocode[0] if len(cocode) == 4 else min(24, randint(0,29))
        yield codewordNo, cocode, t



@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref", [(mat24fast, Mat24)])
def test_syndrome(gc, ref):
    for i in range(24):
        for j in range(i+1,24):
            for k in range(j+1,24):
                v = (1<<i) + (1<<j) + (1<<k)
                if gc.syndrome(v ^ 0xffffff) != v:
                    print( i,j,k, binbv(v), binbv( LexGolayCode.syndrome(v)) )
                    raise ValueError( "decoding error" )

    for codewordNo, cocode, t in syndrome_testcases():
        codeword = gc.gcode_to_vect(codewordNo)
        decoded_codewordNo = gc.vect_to_gcode(codeword) 
        weight =  gc.bw24(codeword)
        assert weight  == 4 * gc.gcode_weight(codewordNo), (
             hex(codeword), gc.gcode_weight(codewordNo) )
        assert decoded_codewordNo == codewordNo
        if codewordNo: 
            assert gc.bw24(codeword) in [8,12,16,24], map(hex,[codewordNo,codeword])
        synd = sum([1 << i for i in cocode])
        vector = codeword ^ synd
        encoded = gc.vect_to_vintern(vector)
        decoded = gc.vintern_to_vect(encoded)
        assert vector == decoded

        syndrome = gc.syndrome(vector,t)
        assert syndrome  == synd , map(hex, [codeword, synd,t, syndrome])
        coc = gc.vect_to_cocode(vector)
        assert coc == gc.vect_to_cocode(synd)    
        sample_v =  gc.cocode_to_vect(gc.vect_to_cocode(vector))
        assert coc == gc.vect_to_cocode(sample_v) 
        if ref:
            assert  codeword == ref.gcode_to_vect(codewordNo)
            assert  decoded_codewordNo == ref.vect_to_gcode(codeword)
            assert  syndrome == ref.syndrome(vector,t) 
            assert  coc == ref.vect_to_cocode(vector)
            assert  weight  == 4 * ref.gcode_weight(codewordNo)
        if len(cocode) == 4:
            with pytest.raises(ValueError):
                gc.syndrome(vector,24)
    print( "Golay code syndome decoding test passed" )


def cocode_weight_testcases(gc):
   yield 0
   for i in range(24):
       yield gc.vect_to_cocode(1 << i)
   for i in range(500):
       yield randint(0, 0xfff)



@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref", [(mat24fast, Mat24)])
def test_cocode_weight(gc, ref):
    for c1 in cocode_weight_testcases(gc):
        w = gc.bw24(gc.cocode_syndrome(c1, 0))
        assert ref.cocode_weight(c1) ==  w
        assert gc.cocode_weight(c1) ==  w

def octad_testcases():
    testcases = [
        1,  0x800, 0x903
    ]
    for c in testcases: yield c
    for i in range(1000):
        yield randint(0,0xfff)   


@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref", [(mat24fast, Mat24)])
def test_octads(gc, ref):    
    for c in octad_testcases():
        try:
            oct =  gc.gcode_to_octad(c)
        except:
            oct = None
        v = gc.gcode_to_vect(c)
        if oct is None:
            assert gc.bw24(v) in [0,12,16,24], (map(hex,[c, v]))
            continue
        assert 0 <= oct < 759, (hex(c), hex(v), oct)
        assert gc.bw24(v) == 8
        assert gc.vect_to_octad(v) == oct
        assert gc.gcode_to_octad(c) == oct
        assert gc.octad_to_vect(oct) == v, (hex(c), oct, hex(v))
        assert gc.octad_to_gcode(oct) == c, (oct, hex(c))
        bitlist, dummy = gc.vect_to_bit_list(v)
        assert v == sum([1 << i for i in bitlist])
        if ref:
            assert  oct == ref.gcode_to_octad(c)
            assert ref.vect_to_octad(v) == oct
            assert ref.gcode_to_octad(c) == oct
            assert ref.octad_to_vect(oct) == v, (hex(c), oct, hex(v))
            assert ref.octad_to_gcode(oct) == c, (oct, hex(c))

    print( "Golay code octad test passed" )




@pytest.mark.mat24
@pytest.mark.parametrize("gc", [Mat24, mat24fast])
def test_theta_blackwhite(gc):
    """ yet to be documented!!!! """
    for v1 in range(0, 0x1000, 0x800):
        for v2 in range(0x20):
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
    for v0 in range(0,0x800,0x20):
        v = gc.gcode_to_vect(v0)
        cc0 = gc.ploop_theta(v0)   & 0x81f
        cc = gc.cocode_to_vect(cc0)
        ccref = ((v | (v >> 1) | (v >> 2)) & 0x111111)
        assert cc == ccref , map(hex, [v,cc,ccref])




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

        test_theta_blackwhite(gc)
        test_theta_colored(gc)

        # test exponentiation
        v4 = randint(0,0x1fff)
        v4_inv = gc.pow_ploop(v4,3)
        assert gc.mul_ploop(v4, v4_inv) == 0, map(hex, [v4, v4_inv, gc.mul_ploop(v4, v4_inv)])
        assert gc.pow_ploop(v4,4) == 0
        c4 = gc.gcode_to_vect(v4)
        assert gc.pow_ploop(v4,2) == (gc.bw24(c4) & 4) << 10
        assert gc.pow_ploop(v4,5) == v4
  
    print( "Parker Loop test passed" )




def Mat24_testcases():
    """yields pairs (v,k) for testing the implementation of Mat24

    v is a Golay code vector. k is an integer representing an element 
    of the Mathieu group Mat24. You may use the Mat24 method 
    m24num_to_perm to convert k to a permutation.

    """
    testdata = [
       (0,  0),
       (0xdb1235, 115873693 ),
       (0xffffff, 244823040-1 ),
       (0,  244823040-1),
    ]
    for t in testdata: yield t
    for i in range(300):
        yield randint(0, 0xffffff), randint(0,244823040-1)


@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref, small", [(mat24fast, Mat24, True)])
def test_Mat24(gc, ref, small):
    """test implementation of the Mathieu group Mat24"""
    for i, (v, k) in enumerate( Mat24_testcases()):
        p =  gc.m24num_to_perm(k)
        if ref and i < 20:
            p_ref = ref.m24num_to_perm(k)
            assert p == p_ref, (p, p_ref, k)
        assert sum(p) == 12*23
        v1 = gc.op_vect_perm(v, p)
        assert gc.vect_type(v) == gc.vect_type(v1), (
            map(hex,[v,k,v1]), i, gc.vect_type(v),gc.vect_type(v1))
        k1 = gc.perm_to_m24num(p)
        if ref and i < 20:
            k1_ref = ref.perm_to_m24num(p)
            assert k1 == k1_ref, (i, k1, k1_ref)
        assert k == k1,  (map(hex, [v, k, k1]),i)
        assert gc.perm_check(p) == 0
        #if (i & 0x3ff) == 0x3ff: print( '.', end ="" )

        vc = randint(0,0xfff)
        m = gc.perm_to_matrix(p)
        vc1 = gc.vect_to_gcode(
               gc.op_vect_perm(gc.gcode_to_vect(vc), p)   )
        vc2 = gc.op_gcode_matrix(vc, m)
        vc2a = gc.op_gcode_perm(vc, p)
        assert vc1 == vc2, (map(hex,[vc,vc1,vc2]),p,map(hex,m))
        assert vc1 == vc2a

        vc3 = gc.op_cocode_perm(vc, p)
        vc4 = gc.op_vect_perm(gc.cocode_to_vect(vc), p)
        vc4a = gc.vect_to_cocode(vc4)
        assert vc3 == vc4a

        p1 = gc.matrix_to_perm(m)
        assert p == p1

        if small and i > 100:
            break
        if ref and i < 100:
            assert p ==  ref.m24num_to_perm(k), (hex(k),p,ref.m24num_to_perm(k))
            assert k1 == ref.perm_to_m24num(p)
            assert vc2 == ref.op_gcode_matrix(vc, m)
    print( "test of the Mathieu group operation passed" )
    





def Mat24_group_testcases():
    """yields pairs (v,k1, k2) for testing the group operation of Mat24

    v is a Golay code vector. k1, k2 are integers representing an element 
    of the Mathieu group Mat24. You may use the Mat24 method 
    Mat24_int_to_perm to convert k1, k2 to a permutation.

    """
    testdata = [
    ]
    for t in testdata: yield t
    mmax = 244823040-1
    for i in range(40):
        yield randint(0, 0xffffff), randint(0,mmax), randint(0,mmax)



@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref, small", [(mat24fast, Mat24, True)])
def test_group_Mat24(gc, ref, small):
    """test implementation of the Mathieu group operations in Mat24"""
    if gc is None: gc = Mat24
    for i, (v, k1, k2) in enumerate( Mat24_group_testcases()):
        p1 =  gc.m24num_to_perm(k1)
        p2 =  gc.m24num_to_perm(k2)
        p12 = gc.mul_perm(p1, p2)
        v1 =  gc.op_vect_perm(v, p12)
        v2 =  gc.op_vect_perm(gc.op_vect_perm(v, p1),p2)
        assert v1 == v2
        p2i = gc.inv_perm(p2)
        p1a = gc.mul_perm(p12, p2i)
        assert p1 == p1a
    print( "test of the Mathieu group implementation passed" )


 

AUT_PL_SMALL = 20

def Aut_Pl_testcases(gc, small = False):
    """yield tuples (c1, p1, c2, p2, v1, v2) for testing Automorphisms of Pl

    (c1, p1) and (c2, p2) are pairs of a cocode element and an integer
    representation to form an automorphism of the Parker loop.
    v1 and v2 are elements of the Parker loop given as 13-bit integers   
    """
    Id = 0   # number of identity permutation 
    testcases = [
        (0, Id, 0, Id, 0x111, 0x222),
        (0x55, 567234, 0x8a5, None, 0x111, 0x222),
        (0x356, 567234, 1, None, 0x876, 0x9ab),
        (0, 1, 0, Id, 0x111, 0x222),
    ]
    for t in testcases:
        yield t
    mmax = gc.MAT24_ORDER-1
    def rand_aut():
        return randint(0, 0xfff), randint(0, mmax)
    def rand_v1v2():
        return randint(0, 0x1fff), randint(0, 0x1fff)
    for i in range(AUT_PL_SMALL if small else 200):
        yield  rand_aut() + rand_aut() +  rand_v1v2()
  

@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref, small", [(mat24fast, Mat24, True)])
def test_Aut_Pl(gc, ref, small):
    """test Automorphism operations of the parker loop"""
    def to_aut(c, p, check = True):
        p = gc.m24num_to_perm(0 if p is None else p)
        m = gc.perm_to_autpl(c, p)
        if check and ref:
             mref = ref.perm_to_autpl(c, p)
             assert m == mref, (map(hex,m), map(hex,mref))
        return m

    def autpl_op_tuple(tuple_, a, g):
        return tuple([g.op_ploop_autpl(v, a) for v in tuple_])

    for i, (c1, p1, c2, p2, v1, v2) in enumerate(Aut_Pl_testcases(gc, small)):
        p = gc.m24num_to_perm(p1)
        m = gc.perm_to_matrix(p)
        check = i < AUT_PL_SMALL
        a1, a2 = to_aut(c1, p1, check),  to_aut(c2, p2, check) 
   
        # test composition and decomposition
        assert c1 == gc.autpl_to_cocode(a1)
        assert p == gc.autpl_to_perm(a1)
        a1_0 = gc.perm_to_autpl(0, p)
        assert a1 == gc.mul_autpl(gc.cocode_to_autpl(c1), a1_0)

        # v1, v2 are Parker loop vectors
        # a1, a2 are Parker loop automorphisms
        # test (v1 * a1) * (v2 * a1) =  (v1 * v2) * a1 
        v3 = gc.mul_ploop(v1, v2)
        if check and ref:
            assert  v3 == ref.mul_ploop(v1, v2)
        v11, v12, v13 = autpl_op_tuple((v1, v2, v3), a1, gc)
        if check and ref:
            assert (v11,v12,v13) == autpl_op_tuple([v1,v2,v3], a1, ref)
        v13a = gc.mul_ploop(v11, v12)
        assert v13 == v13a, (i, map(hex,[v13, v13a]))

        # test (v * a1) * a2 = v * (a1 * a2) for v = v1, v2, v3
        v21, v22, v23 = autpl_op_tuple((v11, v12, v13), a2, gc)
        v23b = gc.mul_ploop(v21, v22)
        assert  v23 == v23b
        a12 = gc.mul_autpl(a1, a2)
        if ref and check:
            assert a12 == ref.mul_autpl(a1, a2)
        v21a, v22a, v23a = autpl_op_tuple((v1, v2, v3), a12, gc)
        assert (v21, v22, v23) == (v21a, v22a, v23a)

        # test inversion: check that (a1 * a2) * a2**(-1) == a2
        a2i = gc.inv_autpl(a2)
        a1n = gc.mul_autpl(a12, a2i)
        assert a1 == a1n, (map(hex,a1), map(hex,a1n))
        # test inversion function perm_to_autpl_inv()
        perm2 = gc.m24num_to_perm(0 if p2 is None else p2)
        perm2i, a2ia = gc.perm_to_iautpl(c2, perm2)
        assert a2i == a2ia 
        assert perm2i ==  gc.inv_perm(perm2)     
        
    print( "Test of Automorphisms of Parker Loop passed" )


def suboctad_testcases(gc):
    """yield testcases (octad, suboctad, status) for suboctad functions

    octad is an octad, suboctad is a suboctad, both in vector
    representation. status is 0 if test should not fail.
    status = 2  if octad is not a good octad
    status = 1  is suboctad is not a suboctad of the octad
    """
    testdata = [
        (0x82b114, 0xa50, 1),
        (0xa9c50, 0x29850, 0), 
        (0xffff, 0, 2), (0,1,2), (0x00ff00, 0x001000,1)
    ]
    for d in testdata:
        yield d
    for i in range(100):
        oct = gc.gcode_to_vect(randint(0,0xfff)) 
        sub = gc.cocode_syndrome(randint(0, 0xffffff), gc.lsbit24(oct))
        if gc.bw24(oct) != 8: 
            status = 2
        elif sub & oct != sub or gc.bw24(sub) & 1:
            status = 1
        else:
            status = 0
        yield  oct, sub, status
    for i in range(200):
        oct = gc.octad_to_vect(randint(0,758)) 
        sub = randint(0, 0xffffff) & oct
        status = gc.bw24(sub) & 1
        while (i > 50 and status):
            sub = randint(0, 0xffffff) & oct
            status = gc.bw24(sub) & 1
        yield oct, sub, status
                


@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref", [(mat24fast, Mat24)])
def test_suboctads(gc, ref):
    for octad, suboctad, status in suboctad_testcases(gc):
        v = gc.vect_to_gcode(octad)
        c = gc.vect_to_cocode(suboctad)
        last_u_sub = 0
        if status == 0:
            u_sub = gc.cocode_to_suboctad(c,v)
            c1 = gc.suboctad_to_cocode(u_sub,v)
            syn1 = gc.cocode_syndrome(c1, gc.lsbit24(octad))
            assert c == c1, map(hex,[octad, suboctad, u_sub, syn1, status])
            w = gc.suboctad_weight(u_sub)
            assert w == (gc.bw24(suboctad) >> 1) & 1
            scal_ref = (gc.suboctad_weight(u_sub ^ last_u_sub)
              ^ gc.suboctad_weight(u_sub) ^ gc.suboctad_weight(last_u_sub))
            assert gc.suboctad_scalar_prod(u_sub, last_u_sub) == scal_ref
            if ref:
                assert u_sub == ref.cocode_to_suboctad(c,v)
                assert c == ref.suboctad_to_cocode(u_sub,v)
                assert w == ref.suboctad_weight(u_sub)
                assert ref.suboctad_scalar_prod(u_sub, last_u_sub) == scal_ref
            last_sub = u_sub
        else:
            try:
                gc.cocode_to_suboctad(c,v)
            except:
                #print("gc.cocode_to_suboctad failed as expected")
                pass
            else:
                syn1 = gc.cocode_syndrome(suboctad, gc.lsbit24(octad))
                print ("Bad suboctad", map(hex,[octad, suboctad, syn1, status]))
                raise ValueError("cocode_to_suboctad() should fail bud didn't.")
            if status == 2:
                try:
                    gc.suboctad_to_cocode(randint(0,63),v)
                except:
                    pass
                else:
                    print ("Bad octad", map(hex,[octad, suboctad, v, c, status]))
                    raise ValueError("suboctad_to_cocode() should fail bud didn't.")
    print("Golay code suboctad test passed")




