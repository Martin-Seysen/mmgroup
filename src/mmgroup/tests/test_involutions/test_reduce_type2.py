
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



V = MMSpace(15)

#########################################################################
## Auxiliary functions
########################################################################



def apply_perm(v, p, log_list):
    """Apply permutation to vector in Leech lattice mod 2.
  
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
    """A variant of function ``apply_perm``.
    
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


def reduce_type2(v, sign = 1):
    """Map short vector in Leech lattice to standard frame

    This is a python implementation of the C function
    ``gen_leech2_reduce_type2`` in file ``gen_leech.c``.
   
    Let ``v_2 \in \Lambda / 2 \Lambda`` of type 2 be given by 
    parameter ``v2`` in Leech lattice encoding. Then the function 
    constructs a ``g \in G_{x0}`` that maps ``v_2`` to the 
    standard vector  ``v_0`` which corresponds to the Golay cocode
    word  ``e_2 + e_3``, where ``e_i`` is the ``i``-th basis
    vector of ``\{0,1\}^{24}``.
  
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
    operation f$g`` maps ``v_2`` to the positive vector ``v_2``. 
    """
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
## Wrapper for the C function ``gen_leech2_reduce_type2``
#########################################################################
      

def reduce_type2_fast(v2, sign = 1):
    r"""Wrapper for the C function ``gen_leech2_reduce_type2``"""
    res = np.zeros(10, dtype = np.uint32)
    length = gen_leech2_reduce_type2(v2, sign, res)
    assert length >= 0, (hex(v2), hex(length))
    return list(res[:length])


#########################################################################
## Testing the C function ``gen_leech2_reduce_type2``
#########################################################################



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
    
    


 