
from random import randint, shuffle
import numpy as np

import pytest


from mmgroup import MM0, MMSpace, Cocode, GCode, AutPL, MMV
#from general import get_case, leech_type, span

from mmgroup import mat24
from mmgroup.mat24 import MAT24_ORDER, gcode_weight

from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_reduce_type2_ortho
from mmgroup.generators import gen_leech2_op_word

from mmgroup.clifford12 import leech_matrix_2A_axis_type



V = MMV(15)

#########################################################################
## Auxiliary functions
########################################################################



def apply_perm(v, p, log_list):
    r"""Apply permutation to vector in Leech lattice mod 2.
  
    Let ``pi`` be the permutation given by the array ``p`` 
    as a permutation on the set  ``\{0,\ldots,23\}``. Let 
    ``v_2`` be the vector in the Leech lattice mod  2 given 
    by parameter ``v2``. The function returns ``v_2 * pi``
    Parameter ``v2`` and the return value are given in Leech
    lattice encoding.
  
    The function computes  the permutation ``x_pi`` to the list
    ``log_list``. Here ``x_pi`` is encoded as a generator of the
    monster group as as described  in file ``mmgroup_generators.h``.
    That generator is stored with tag  ``MMGROUP_ATOM_TAG_IP`` so
    that we can compute the inverse of ``pi`` very efficiently. 
    """
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
    r"""A variant of function ``apply_perm``.
    
    Here permutation ``pi``  must be given as a mapping from
    an umbral heptad referred by ``p_src`` to an umbral heptad 
    referred by ``p_dest``. Here ``p_src`` and ``p_dest`` define
    the  permutation ``pi`` as in function 
    ``mat24_perm_from_heptads`` in file ``mat24_functions.c c``.
    """
    p =  mat24.perm_from_heptads(src, dest)
    return apply_perm(v, p, log_list)


def apply_perm_map(v, src, dest, log_list):
    r"""Another variant of function ``apply_perm``.
  
    Here permutation ``pi``  must be given as a mapping from any
    subset of ``\{0,1\}^{24}`` of length ``n`` referred by ``p_src`` 
    to another subset referred by ``p_dest``. Here ``p_src``,
    ``p_dest``,  and ``n``  specify  a  permutation ``pi`` as in 
    function   ``mat24_perm_from_heptads`` in 
    file ``mat24_functions.c c``.

    Caution: Permutation ``pi`` might not be specified uniquely 
    by  the given mapping!
    """
    res, p =  mat24.perm_from_map(src, dest)
    assert res > 0
    return apply_perm(v, p, log_list)




RED_STD_HEPTAD = [0,1,2,3,4,5,8]


def map_to_standard12_xi_exp(v):
    r"""Find exponent of ``xi`` mapping Leech vector ``v`` a subspace 

    Let ``v_2`` be the vector in the Leech lattice mod 2 given 
    by parameter ``v2`` in Leech lattice encoding. Let  
    ``\Omega, \omega, \gamma(.)`` be as in [Sey20].
  
    The function  tries to find an exponent ``e`` such that 
    ``v_2 \xi^e`` is  equal to an element 
    ``\lambda_\delta \pmod{\lambda_\Omega}`` of
    ``\Lambda / 2 \Lambda``, 
    ``\delta \in \mathcal{C}^*``, ``\delta`` even. The 
    function returns ``e`` if such an  ``e`` exists and -1
    otherwise.
  
    Assume ``v_2 = \lambda_d + \lambda_\delta + \lambda_\epsilon, 
    d \in \mathcal{C}, \delta, \epsilon \in \mathcal{C}^*``, 
    ``d, \delta`` grey, even,  ``\epsilon`` coloured. 
    The function returns
  
    ``e=0`` if ``d=0 \pmod{\Omega}``,
     
    ``e=1`` if ``\delta=\gamma(d) \pmod{\omega} ``,

    ``e=2`` if ``\delta=0  \pmod{\omega} ``. 
  
    In all other cases there is no suitable exponent ``e``.
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



###########################################################################
## Reduce a type-2 vector in the Leech lattice mod 2 
###########################################################################


def reduce_type2(v2, sign = 1):
    r"""Map short vector in Leech lattice to standard frame

    This is a python implementation of the C function
    ``gen_leech2_reduce_type2`` in file ``gen_leech.c``.
   
    Let ``v_2 \in \Lambda / 2 \Lambda`` of type 2 be given by 
    parameter ``v2`` in Leech lattice encoding. Then the function 
    constructs a ``g \in G_{x0}`` that maps ``v_2`` to the 
    standard vector  ``v_0``. Here  ``v_0`` is the short 
    vector the Leech lattice propotional  to  ``e_2 - e_3``, 
    where ``e_i`` is the ``i``-th basis vector 
    of ``\{0,1\}^{24}``.
  
    The element ``g`` is returned as a word in the generators
    of ``G_{x0}`` of length ``n \leq 6``. Each atom of the 
    word ``g`` is encoded as  defined in the header 
    file ``mmgroup_generators.h``. 

    The function stores ``g`` as a word of generators in the
    array ``pg_out`` and returns the length  ``n``  of that
    word. It returns a negative number in case of failure, 
    e.g. if ``v_2`` is not of type 2.
  
    If ``sign`` is not zero then ``v_2`` is interpreted as an
    element of the extraspecial group ``2^{1+24}`` an the 
    operation ``$g`` maps ``v_2`` to the positive vector ``v_2``. 
    """
    result = []
    for _i in range(5):
        gc = mat24.gcode_to_vect(v2 >> 12)
        coc = (v2 ^  mat24.ploop_theta(v2 >> 12)) & 0xfff
        vtype = gen_leech2_type(v2)
        if vtype == 0x21:
            exp = 2 - ((v2 >> 22) & 1)
        elif vtype == 0x22:
            exp = map_to_standard12_xi_exp(v2)
            if exp < 0:
                w, src = mat24.vect_to_bit_list(gc)
                if (w == 16):
                    src = src[16:22] + src[0:1]
                elif w == 8:
                    src[6] = src[8]
                else:
                    raise ValueError("WTF1")
                v2 = apply_perm_heptad(v2, src, RED_STD_HEPTAD, result)
                exp = map_to_standard12_xi_exp(v2)
                assert exp >= 0
        elif vtype == 0x20:
            exp =  0
            # map v2 to stadard cocode word [2,3]
            syn = (mat24.cocode_syndrome(v2, 0))
            src = [i for i in range(24) if syn & (1 << i)]
            #print("cc", hex(v2), syn, src)
            if src != [2,3]:
                v2 = apply_perm_map(v2, src, [2,3], result)
            # correct v2 if v2 is the cocode word [2,3] + Omega
            if v2 & 0x800000:
                atom = 0xC0000200  
                   # operation y_d such that d has odd scalar
                   # product with cocode word [2,3]
                v2 = gen_leech2_op_atom(v2, atom)
                result.append(atom)
            assert v2 & 0xffffff == 0x200
            if sign and v2 & 0x1000000:
                atom = 0xB0000200  
                   # operation y_d such that d has odd scalar
                   # product with cocode word [2,3]
                v2 = gen_leech2_op_atom(v2, atom)
                result.append(atom)
                assert v2  == 0x200
            return result
        else:
            raise ValueError("WTF")
        if exp: 
            exp = 0xE0000003 - exp
            v2 = gen_leech2_op_atom(v2, exp)
            result.append(exp)
    raise ValueError("WTF1")
  




#########################################################################
## Wrapper for the C function ``gen_leech2_reduce_type2``
#########################################################################
      

def reduce_type2_fast(v2, sign = 1):
    r"""Wrapper for the C function ``gen_leech2_reduce_type2``"""
    res = np.zeros(10, dtype = np.uint32)
    length = gen_leech2_reduce_type2(v2, sign, res)
    assert length >= 0, (hex(v2), hex(length))
    return list(res[:length])




###########################################################################
## Reduce (orthogonal) type-2 vector in the Leech lattice mod 2 
###########################################################################


RED_23_ORTHO =  [2, 3, 4, 5, 6, 7]
RED_23_23014 =  [2, 3, 0, 1, 4]


def reduce_type2_ortho(v2):
    r"""Map (orthgonal) short vector in Leech lattice to standard vector

    This is a python implementation of the C function
    ``gen_leech2_reduce_type2_ortho`` in file ``gen_leech.c``.
   
    Let ``v_2 \in \Lambda / 2 \Lambda`` of type 2 be given by 
    parameter ``v2`` in Leech lattice encoding. 

    In the real Leech lattice, (the origin of) the vector ``v_2`` must
    be orthogonal to the standard short vector ``v_0``. Here ``v_0``
    is the short vector in the Leech  lattice  propotional
    to  ``e_2 - e_3``, where ``e_i`` is  the ``i``-th basis vector
    of ``\{0,1\}^{24}``.
   
    Let ``v_1`` be the short vector in the Leech lattice propotional
    to  ``e_2 + e_3``.  Then the function constructs 
    a ``g \in G_{x0}`` that maps ``v_2`` to ``v_1``.
 
    The element ``g`` is returned as a word in the generators
    of ``G_{x0}`` of length ``n \leq 6``. Each atom of the 
    word ``g`` is encoded as  defined in the header 
    file ``mmgroup_generators.h``. 

    The function stores ``g`` as a word of generators in the
    array ``pg_out`` and returns the length  ``n``  of that
    word. It returns a negative number in case of failure, 
    e.g. if ``v_2`` is not of type 2,  or not orthogonal 
    to ``v_1`` in the real Leech lattice.
    """
    if (gen_leech2_type(v2) >> 4) != 2:
        raise ValueError("Vector is not short")
    if (gen_leech2_type(v2 ^ 0x200) >> 4) != 4:
        raise ValueError("Vector not orthogonal to standard vector")
    result = []
    for _i in range(4):
        gc = mat24.gcode_to_vect(v2 >> 12)
        vtype = gen_leech2_type(v2)
        #print(_i, hex(v2), hex(gc), hex(coc), hex(vtype))
        if vtype == 0x21:
            exp = 2 - ((v2 >> 22) & 1)
        elif vtype == 0x22:
            exp = map_to_standard12_xi_exp(v2)
            if exp < 0:
                if (mat24.gcode_weight(v2 >> 12) & 4):
                    gc ^= 0xffffff
                assert mat24.bw24(gc) == 8
                if gc & 0xc:
                    gc &= ~0xc
                    _, oct = mat24.vect_to_bit_list(gc)
                    oct = [2, 3] + oct[:3]
                    v2 = apply_perm_map(v2, oct, RED_23_23014, result)
                else:
                    a = 0
                    _, oct = mat24.vect_to_bit_list(gc)
                    for j in range(5): a += 1 << oct[j];
                    a &= ~mat24.syndrome(a ^ 0xc, 0)
                    _, oct = mat24.vect_to_bit_list(a); 
                    oct = [2, 3] + oct[:4]
                    v2 = apply_perm_map(v2, oct, RED_23_ORTHO, result)
                exp = map_to_standard12_xi_exp(v2)
                assert exp >= 0
        elif vtype == 0x20:
            if ((v2 & 0xffffff) == 0x800200):
                return result
            syn = (mat24.cocode_syndrome(v2, 0)) & ~0xc
            if syn and syn != 3:
                c = [2,3] + [i for i in range(24) if syn & (1 << i)][:2]
                v2 = apply_perm_map(v2, c, RED_23_23014[:4], result)
            exp = 2 - ((v2 >> 23) & 1);
        else:
            raise ValueError("WTF")
        if exp: 
            exp = 0xE0000003 - exp
            v2 = gen_leech2_op_atom(v2, exp)
            result.append(exp)
    raise ValueError("WTF1")
  
#########################################################################
## Wrapper for the C function ``gen_leech2_reduce_type2_ortho``
#########################################################################
      

def reduce_type2_ortho_fast(v2):
    r"""Wrapper for the C function ``gen_leech2_reduce_type2_ortho``"""
    res = np.zeros(1000, dtype = np.uint32)
    length = gen_leech2_reduce_type2_ortho(v2, res)
    assert length >= 0, (hex(v2), hex(length))
    return list(res[:length])




#########################################################################
## Testing the C function ``gen_leech2_reduce_type2``
#########################################################################



v_start_tuple = ("I", 3, 2)
v_start = V(*v_start_tuple)
v_opp = v_start * MM0('x', 0x200)


def make_testcases():
    #V = v_start.space
    yield v_start.copy()
    yield V("I", 11, 9)
    for i in range(30):
        yield v_start * MM0("r", "N_x0")
    for i in range(200):
        yield v_start * MM0("r", "G_x0")


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
        mm_op = MM0('a', w_op)
        w_op = np.array(w_op, dtype = np.uint32)
        v1 = gen_leech2_op_word(w1, w_op, len(w_op)) & 0xffffff
        if verbose: print("Reduced v", hex(v1))
        assert v1 == 0x200
        vv = w * mm_op
        assert vv in [v_start, v_opp], (vv, v_start, v_opp)
        if verbose: print(vv)
    
    

#########################################################################
## Testing the C function ``gen_leech2_reduce_type2_ortho``
#########################################################################

v_start_0 =  0x200 
assert Cocode(v_start_0).syndrome_list() == [2,3]
v_ortho_start = 0x800200


def rand_pi_mat22():
    r"""Generate certain 'random' element of AutPL

    The function generates a random element ``e`` of the automorphism
    group ``AutPL`` of the Parker loop, such that the permutation
    in the Mathieu ``M24`` group corresponding to ``e`` preserves the 
    set  ``\{2, 3\}``.
    """
    pi = AutPL('r', 'r')
    perm = pi.perm
    l_src = perm[2:4]
    shuffle(l_src)
    _, a = mat24.perm_from_map(l_src, [2, 3])
    return  pi  *  AutPL('r', a)  

def rand_Co2(quality = 5):
    r"""Generate certain 'random' element in a subgroup of ``G_x0``

    Let ``v_0`` be the element of the subgroup  ``Q_x0`` of ``G_x0``
    corresponding to the Golay cocode element ``[2,3]``.

    The function generates a random element of the centralizer
    of ``v_0`` in ``G_x0``.
    """
    a = MM0()
    for i in range(quality):
         pi = rand_pi_mat22()
         x1 = randint(0, 0xfff) & ~0x200
         y1 = randint(0, 0xfff) & ~0x200
         e = randint(0, 2)
         a *= MM0([('p', pi), ('x', x1), ('y', y1), ('l', e)]) 
    return a 


def mul_v2_mm(v, g):
    """Transfrom element of Leech lattice mod 2

    Let ``v`` be an element of the leech lattice mod 2, stored as
    an integer in **Leech lattice encoding**. The the function 
    transforms ``v`` with the element ``g`` of ``G_x0`` and returns
    the result.

    More precisely, ``v`` is intepreted as an element of the 
    extraspecial 2-group ``Q_x0``, and the function conjugates
    ``v`` with ``g``.
    """
    data = g.mmdata
    return gen_leech2_op_word(v, data, len(data))

def v_ortho_g(g):
    return mul_v2_mm(v_ortho_start, g)

def rand_v_ortho(quality = 5):
    """Return random short vector orthogonal to ``v_start_0``.

    Here ``v_start_0`` is the short vector in the Leech lattice
    mod 2 corresponding to the  Golay cocode element ``[2,3]``.
    The function returns a random short vector ``v`` in the
    Leech lattice mod 2 that is orthogonal to ``v_start_0`` in
    the real Leech lattice. The short vector ``v`` is returned 
    in **Leech lattice encoding** with a random sign bit.
    """
    return v_ortho_g(rand_Co2(quality))





def make_testcases_ortho():
    """Yield test cases for function ``reduce_type2_ortho_fast``

    Such a testcase is a short vector in the Leech lattice mod 2
    as returned by function ``rand_v_ortho`` in this module.
    """
    data0 = [0x1157a1a]
    for d in data0:
        yield d

    data = [
        [],
        [('y', 0x500)],
        [('p', rand_pi_mat22())],
    ]
    for d in data:
        yield v_ortho_g(MM0(d))
    for q in range(1,5):
        for i in range(200):
            yield rand_v_ortho(q)


@pytest.mark.involution
def test_reduce_type2_ortho(verbose = 0):
    """Test the C function ``reduce_type2_ortho_fast``"""
    for i, w in enumerate(make_testcases_ortho()):
        if verbose: print("\nTest", i+1)
        if verbose: print("Reducing", hex(w))
        w_op = reduce_type2_ortho(w)
        if verbose:
            print([hex(x) for x in w_op])
        
        w_op_fast = reduce_type2_ortho_fast(w)
        ok = w_op == w_op_fast
        if verbose or not ok: 
            print("Op:  ", [hex(x) for x in w_op])
            print("Fast:", [hex(x) for x in w_op_fast] )
        if not ok:
            print("Op obtained from:", hex(w))
            raise ValueError("Fast function reduce_type2 failed")
        
        assert len(w_op) <= 6, len(w_op)
        mm_op = MM0('a', w_op)

        v1 = mul_v2_mm(w, mm_op)
        if verbose: print("Reduced v", hex(v1))
        assert v1 & 0xffffff == 0x800200,  (hex(w), hex(v1))
    
        img = mul_v2_mm(v_start_0, mm_op)
        ok_vstart = ((v_start_0 ^ img) & 0xffffff) == 0
        assert ok_vstart, (hex(v_start_0), hex(img))


