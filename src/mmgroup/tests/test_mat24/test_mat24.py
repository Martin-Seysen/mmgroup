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
# Test operation of the group Mat24 on Golay code and cocode
#####################################################################



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
        gc.perm_check(p)
        #if (i & 0x3ff) == 0x3ff: print( '.', end ="" )

        vc = randint(0,0xfff)
        m = gc.perm_to_matrix(p)
        vc1 = gc.vect_to_gcode(
               gc.op_vect_perm(gc.gcode_to_vect(vc), p)   )
        vc2 = gc.op_gcode_matrix(vc, m)
        vc2a = gc.op_gcode_perm(vc, p)
        assert vc1 == vc2, (list(map(hex,[vc,vc1,vc2])),p,list(map(hex,m)))
        assert vc1 == vc2a

        vc3 = gc.op_cocode_perm(vc, p)
        vc4 = gc.op_vect_perm(gc.cocode_to_vect(vc), p)
        vc4a = gc.vect_to_cocode(vc4)
        assert vc3 == vc4a, (hex(vc3), hex(vc4), hex(vc4a), p)

        p1 = gc.matrix_to_perm(m)
        assert p == p1

        if small and i > 100:
            break
        if ref and i < 100:
            assert p ==  ref.m24num_to_perm(k), (hex(k),p,ref.m24num_to_perm(k))
            assert k1 == ref.perm_to_m24num(p)
            assert vc2 == ref.op_gcode_matrix(vc, m)
            p1 = [p[1], p[0]] + p[2:]
            with  pytest.raises(ValueError):
                gc.perm_check(p1)
    print( "test of the Mathieu group operation passed" )
    



#####################################################################
# Test group operation of the group Mat24 
#####################################################################


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


 
#####################################################################
# Test operation of the group Aut_PLoop on the Parker loop 
#####################################################################

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
@pytest.mark.parametrize("gc, ref, small", [(mat24fast, Mat24, False)])
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



