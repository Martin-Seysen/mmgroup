
from random import randint
import numpy as np

import pytest


from mmgroup import MM, MMSpace, Cocode, GCode, AutPL
#from general import get_case, leech_type, span

from mmgroup import mat24
from mmgroup.mat24 import MAT24_ORDER, gcode_weight

from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_op_word

from mmgroup.clifford12 import leech_matrix_2A_axis_type

#########################################################################
## Use: C:\Data\projects\scripts\Leech_lattice\reduce.py!!!!!!!!!!!
#########################################################################

#a, v = get_case("2A")

V = MMSpace(15)



def syndrome(v, point):
    cc = mat24.theta(v >> 12) ^ v
    return mat24.cocode_to_vector(cc, 0)


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

def apply_perm_heptad(v, src, dest, log_list):
    p =  mat24.perm_from_heptads(src, dest)
    return apply_perm(v, p, log_list)

def apply_perm_map(v, src, dest, log_list):
    res, p =  mat24.perm_from_map(src, dest)
    assert res > 0
    return apply_perm(v, p, log_list)



def apply_xi(v, e, log_list):
    e1 =  e % 3 + 0x60000000
    res = gen_leech2_op_atom(v, e1)
    assert res & 0xfe000000 == 0
    #print(hex(v), "->",  "li_%d" % e,  "->", hex(res))
    log_list.append(e1)
    return res

RED_STD_HEPTAD = [0,1,2,3,4,5,8]


def map_to_standard12_xi_exp(v):
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


def reduce_type2(v, sign = 1):
    result = []
    for _i in range(5):
        gc = mat24.gcode_to_vect(v >> 12)
        coc = (v ^  mat24.ploop_theta(v >> 12)) & 0xfff
        vtype = gen_leech2_type(v)
        if vtype == 0x21:
            exp = 2 - ((v >> 22) & 1)
        elif vtype == 0x22:
            exp = map_to_standard12_xi_exp(v)
            if exp < 0:
                w, src = mat24.vect_to_bit_list(gc)
                if (w == 16):
                    src = src[16:22] + src[0:1]
                elif w == 8:
                    src[6] = src[8]
                else:
                    raise ValueError("WTF1")
                v = apply_perm_heptad(v, src, RED_STD_HEPTAD, result)
                exp = map_to_standard12_xi_exp(v)
                assert exp >= 0
        elif vtype == 0x20:
            exp =  0
            # map v to stadard cocode word [2,3]
            syn = (mat24.cocode_syndrome(v,0))
            src = [i for i in range(24) if syn & (1 << i)]
            #print("cc", hex(v), syn, src)
            if src != [2,3]:
                v = apply_perm_map(v, src, [2,3], result)
            # correct v if v is the cocode word [2,3] + Omega
            if v & 0x800000:
                atom = 0xC0000200  
                   # operation y_d such that d has odd scalar
                   # product with cocode word [2,3]
                v = gen_leech2_op_atom(v, atom)
                result.append(atom)
            assert v & 0xffffff == 0x200
            if sign and v & 0x1000000:
                atom = 0xB0000200  
                   # operation y_d such that d has odd scalar
                   # product with cocode word [2,3]
                v = gen_leech2_op_atom(v, atom)
                result.append(atom)
                assert v  == 0x200
            return result
        else:
            raise ValueError("WTF")
        if exp: 
            exp = 0xE0000003 - exp
            v = gen_leech2_op_atom(v, exp)
            result.append(exp)
    raise ValueError("WTF1")
  


#########################################################################
## Use: C:\Data\projects\scripts\Leech_lattice\reduce.py!!!!!!!!!!!
#########################################################################
      

def reduce_type2_fast(v2, sign = 1):
    res = np.zeros(10, dtype = np.uint32)
    length = gen_leech2_reduce_type2(v2, sign, res)
    assert length >= 0, (hex(v2), hex(length))
    return list(res[:length])



v_start_tuple = ("I", 3, 2)
v_start = V(v_start_tuple)
v_opp = v_start * MM(('x', 0x200))


def make_testcases():
    V = v_start.space
    yield v_start.copy()
    yield V(  ("I", 11, 9) )
    for i in range(30):
        yield v_start * MM.rand_N_0(in_N_x0 = True)
    for i in range(200):
        yield v_start * MM.rand_G_x0()


@pytest.mark.involution
def test_reduce_type2(verbose = 0):
    for i, w in enumerate(make_testcases()):
        if verbose: print("\nTest", i+1)
        if verbose: print("Reducing", w)
        w1 = leech_matrix_2A_axis_type(15, w.data) & 0xffffff
        w_op = reduce_type2(w1)
        w_op_fast = reduce_type2_fast(w1)
        ok = w_op == w_op_fast
        if verbose or not ok: print("Op:  ", [hex(x) for x in w_op])
        if not ok:
            print("Fast:  ", [hex(x) for x in w_op_fast] )
            print("Op obtained from:", hex(w1))
            raise ValueError("Fast function reduce_type2 failed")
        assert len(w_op) <= 6, len(w_op)
        mm_op = MM.from_data(w_op)
        w_op = np.array(w_op, dtype = np.uint32)
        v1 = gen_leech2_op_word(w1, w_op, len(w_op)) & 0xffffff
        if verbose: print("Reduced v", hex(v1))
        assert v1 == 0x200
        vv = w * mm_op
        assert vv in [v_start, v_opp], (vv, v_start, v_opp)
        if verbose: print(vv)
    
    


 