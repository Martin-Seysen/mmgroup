from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint

from mmgroup import mat24
mat24fast = mat24


from mmgroup.bitfunctions import bit_mat_mul, v2
from mmgroup.dev.mat24.mat24theta import theta_to_basis_vector


from mmgroup.dev.mat24.mat24_ref import Mat24

#####################################################################
# Test basis of Golay code and cocode against its reference
#####################################################################

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



#####################################################################
# Test function lsbit()
#####################################################################


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


#####################################################################
# Test functions extract_b24() and spread_b24()
#####################################################################


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


#####################################################################
# Test functions related to Golay code conversions
#####################################################################

def gcode_testcases():
    """yields Golay code words as test cases
    """
    yield 0
    for i in range(500):
        yield randint(0, 0xfff)


@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref", [(mat24fast, Mat24)])
def test_gcode(gc, ref):
    for codewordNo in gcode_testcases():
        codeword = gc.gcode_to_vect(codewordNo)
        decoded_codewordNo = gc.vect_to_gcode(codeword) 
        weight =  gc.bw24(codeword)
        assert weight  == 4 * gc.gcode_weight(codewordNo), (
             hex(codeword), gc.gcode_weight(codewordNo) )
        assert decoded_codewordNo == codewordNo
        lg, gcode_bits = gc.vect_to_bit_list(codeword)
        gcode_bits = gcode_bits[:lg]
        assert gcode_bits == gc.gcode_to_bit_list(codewordNo)
        if codewordNo: 
            assert gc.bw24(codeword) in [8,12,16,24], map(hex,[codewordNo,codeword])
        if ref:
            assert gcode_bits == ref.gcode_to_bit_list(codewordNo)
            assert  codeword == ref.gcode_to_vect(codewordNo)
            assert  decoded_codewordNo == ref.vect_to_gcode(codeword)
            assert  weight  == 4 * ref.gcode_weight(codewordNo)
    print("Golay code test passed")


#####################################################################
# Test functions related to Golay code syndrome
#####################################################################


def syndrome_testcases(ntests = 2000):
    """yields test cases (codewordNo, cocode, t) for Mat24.syndrome

    codewordNo is 12-bit random number indicating a Golay code word
    cocode is a list of 0...4 bit position where to modify the code word
    t is an entry of cocode if len(cocode) = 4 and a random bit position
    (or 24) otherwise.
    """
    testcases = [
          (1, [1], 3), (0x400,[0],0), (0xd06,[0],0)
        # (0x800, [0x800], 3), (0x200,[0],0), (0x683,[0],0)
    ]
    for t in testcases: 
        yield t
    ll = [0,1,1,2,2,2,3,3,3,3,3,4,4,4,4,4]
    for i in range(ntests):
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
        ## The following stuff is already done in test_gcode()
        #decoded_codewordNo = gc.vect_to_gcode(codeword) 
        #weight =  gc.bw24(codeword)
        #assert weight  == 4 * gc.gcode_weight(codewordNo), (
        #     hex(codeword), gc.gcode_weight(codewordNo) )
        #assert decoded_codewordNo == codewordNo
        
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
        assert gc.cocode_to_bit_list(coc, t) == sorted(cocode)
        if ref:
            assert  ref.cocode_to_bit_list(coc, t) == sorted(cocode)
            assert  syndrome == ref.syndrome(vector,t), (hex(syndrome), t) 
            assert  coc == ref.vect_to_cocode(vector)
            ## The following stuff is already done in test_gcode()
            #assert  codeword == ref.gcode_to_vect(codewordNo)
            #assert  decoded_codewordNo == ref.vect_to_gcode(codeword)
            #assert  weight  == 4 * ref.gcode_weight(codewordNo)
        if len(cocode) == 4:
            with pytest.raises(ValueError):
                gc.syndrome(vector,24)
            t = gc.cocode_to_sextet(coc)
            assert t == ref.cocode_to_sextet(coc)
            assert sum([1 << i for i in t]) == 0xffffff
            for i in range(0, 24, 4):
                v1 = sum([1 << j for j in t[i:i+4]])
                assert gc.vect_to_cocode(v1) == coc
    print( "Golay code syndome decoding test passed" )


#####################################################################
# Test functions related to lists Golay code syndromes
#####################################################################


@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref", [(mat24fast, Mat24)])
def test_syndrome(gc, ref):
    for codewordNo, cocode, t in syndrome_testcases(500):
        codeword = gc.gcode_to_vect(codewordNo)
        synd = sum([1 << i for i in cocode])
        syn_list = gc.all_syndromes(codeword ^ synd)
        syn_list_c = gc.cocode_all_syndromes(gc.vect_to_cocode(synd))
        #print(hex(codeword), cocode, [hex(x) for x in syn_list])
        assert syn_list_c == syn_list
        syn_len = 6 if len(cocode) == 4 else 1
        assert len(syn_list) == syn_len
        for v in syn_list:
            # Next instruction fails if synd ^ v not in Golay code
            gc.vect_to_gcode(synd ^ v)


#####################################################################
# Test functions related to Golay cocode weight
#####################################################################


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


#####################################################################
# Test functions related to octad numbering
#####################################################################


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
        assert gc.octad_to_gcode(oct) == c, (oct, hex(c))
        assert gc.octad_to_vect(oct) == v, (hex(c), oct, hex(v))
        ll, bitlist = gc.vect_to_bit_list(v)
        bitlist = bitlist[:ll]
        assert v == sum([1 << i for i in bitlist])
        if ref:
            assert  oct == ref.gcode_to_octad(c)
            assert ref.vect_to_octad(v) == oct
            assert ref.gcode_to_octad(c) == oct
            assert ref.octad_to_vect(oct) == v, (hex(c), oct, hex(v))
            assert ref.octad_to_gcode(oct) == c, (oct, hex(c))

    print( "Golay code octad test passed" )






#####################################################################
# Test functions related to suboctads 
#####################################################################


def suboctad_testcases(gc):
    """yield testcases (octad, suboctad, status) for suboctad functions

    octad is an octad, suboctad is a suboctad, both in vector
    representation. status is 0 if test should not fail.
    status = 2  if octad is not a good octad
    status = 1  is suboctad is not a suboctad of the octad
    """
    testdata = [
      #   (0x00ffff, 0x110000, 0),
      #   (0x00ffff, 0x001111, 1),
      #   (0xeed2dd, 0x481001, 1),
    ]
    for d in testdata:
        yield d
    for i in range(10000):
        oct = gc.gcode_to_vect(randint(1,0xfff)) 
        short_oct = oct if gc.bw24(oct) in [8, 24] else oct ^ 0xffffff
        cocodev = randint(0, 0xffffff)
        sub = gc.cocode_syndrome(cocodev,  gc.lsbit24(short_oct))
        if gc.bw24(short_oct) == 8: 
            status = sub & short_oct != sub or gc.bw24(sub) & 1
        else:
            status = 2
        yield  oct, sub, status
    for i in range(200):
        oct = gc.octad_to_vect(randint(0,758)) 
        sub = randint(0, 0xffffff) & oct
        status = gc.bw24(sub) & 1
        while (i > 50 and status):
            sub = randint(0, 0xffffff) & oct
            status = gc.bw24(sub) & 1
        sign = randint(-1, 0)
        yield (oct ^ sign) & 0xffffff, sub, status
                


@pytest.mark.mat24
@pytest.mark.parametrize("gc, ref", [(mat24fast, Mat24)])
def test_suboctads(gc, ref):
    for octad, suboctad, status in suboctad_testcases(gc):
        #print("Testing suboctad", hex(octad), hex(suboctad), status)
        v = gc.vect_to_gcode(octad)
        c = gc.vect_to_cocode(suboctad)
        last_u_sub = 0
        if status == 0:
            o = gc.vect_to_octad(octad, 0)
            assert 0 <= o < 759 
            u_sub = gc.cocode_to_suboctad(c,v)
            assert 0 <= u_sub < 1 << 16
            u_sub &= 0x3f
            c1 = gc.suboctad_to_cocode(u_sub, o)
            syn1 = gc.cocode_syndrome(c1, gc.lsbit24(octad))
            assert c == c1, map(hex,[octad, suboctad, u_sub, syn1, status])
            w = gc.suboctad_weight(u_sub)
            assert w == (gc.bw24(suboctad) >> 1) & 1
            scal_ref = (gc.suboctad_weight(u_sub ^ last_u_sub)
              ^ gc.suboctad_weight(u_sub) ^ gc.suboctad_weight(last_u_sub))
            assert gc.suboctad_scalar_prod(u_sub, last_u_sub) == scal_ref
            if ref:
                ref_u_sub = ref.cocode_to_suboctad(c,v)
                assert u_sub == ref_u_sub & 0x3f, (hex(u_sub), hex(ref_u_sub))
                assert c == ref.suboctad_to_cocode(u_sub, o)
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
                syn1 = gc.syndrome(suboctad, gc.lsbit24(octad))
                print ("Bad suboctad", list(map(hex,[octad, suboctad, syn1, status])))
                raise ValueError("cocode_to_suboctad() should fail bud didn't.")
            if status == 2:
                continue # no reasonabale action do do in current setup
                """
                try:
                    gc.suboctad_to_cocode(randint(0,63), )
                except:
                    pass
                else:
                    print ("Bad octad", list(map(hex,[octad, suboctad, v, c, status])))
                    raise ValueError("suboctad_to_cocode() should fail bud didn't.")
                """
    print("Golay code suboctad test passed")





