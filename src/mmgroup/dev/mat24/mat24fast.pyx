# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from libc.stdint cimport uint32_t, uint16_t, uint8_t, int32_t;




include "mat24_functions.pxd"






MAT24_ORDER =  244823040 


###########################################################################
# Some general bit operations
###########################################################################


def lsbit24(uint32_t v1):
    return mat24_lsbit24(v1) 


def bw24(uint32_t v1):
    return mat24_bw24(v1) 


###########################################################################
# Basis of Golay cocode and code
###########################################################################

cdef int _i 
basis =  [MAT24_BASIS[_i] for _i in range(24)]
recip_basis =  [MAT24_RECIP_BASIS[_i] for _i in range(24)]



###########################################################################
# Conversion between bit vectors of GF(2)**24
###########################################################################


def vect_to_bit_list(uint32_t v1):
    cdef uint8_t a[24]
    cdef uint32_t w
    w = mat24_vect_to_bit_list(v1, &a[0])
    return w, list(a[:24])



def vect_to_list(uint32_t v1, uint32_t u_len):
    cdef uint8_t a[24]
    cdef uint32_t w
    w = mat24_vect_to_list(v1, u_len, &a[0])
    return list(a[:w])



def extract_b24(uint32_t v1, uint32_t u_mask):
    return mat24_extract_b24(v1, u_mask)


def spread_b24(uint32_t v1, uint32_t u_mask):
    return mat24_spread_b24(v1, u_mask)



###########################################################################
# Conversion between representations of GF(2)**24, Golay code, etc.
###########################################################################



def vect_to_vintern(uint32_t v1): 
    return mat24_vect_to_vintern(v1)



def vintern_to_vect(uint32_t v1):
    return  mat24_vintern_to_vect(v1)


def vect_to_cocode(uint32_t v1):
    return mat24_vect_to_cocode(v1)


def gcode_to_vect(uint32_t v1):
    return mat24_gcode_to_vect(v1)


def cocode_to_vect(uint32_t c1):
    return mat24_cocode_to_vect(c1)


def vect_to_gcode(uint32_t v1):
    cdef uint32_t cn
    cn = mat24_vect_to_gcode(v1)
    if cn & 0xff000000:
        err = "Bit vector is not a Golay code word"
        raise ValueError, err
    return cn


def gcode_to_octad(uint32_t v1, uint32_t u_strict = 1):
    cdef uint32_t cn
    cn = mat24_gcode_to_octad(v1, u_strict)
    if cn & 0xff000000:
        err = "Golay code word is not an octad"
        raise ValueError, err
    return cn


def vect_to_octad(uint32_t v1, uint32_t u_strict = 1):
    cdef uint32_t v
    v = mat24_vect_to_octad(v1, u_strict)
    if v & 0xff000000:
        err = "Bit vector is not an octad"
        raise ValueError, err
    return v




def octad_to_gcode(uint32_t u_octad):
    cdef uint32_t cn
    cn = mat24_octad_to_gcode(u_octad)
    if cn & 0xff000000:
        err = "Illegal octad number"
        raise ValueError, err
    return cn


def octad_to_vect(uint32_t u_octad):
    cdef uint32_t cn
    cn = mat24_octad_to_vect(u_octad)
    if cn & 0xff000000:
        err = "Illegal octad number"
        raise ValueError, err
    return cn




###########################################################################
# Golay code syndromes and weights
###########################################################################


def syndrome_table(c1):
    return MAT24_SYNDROME_TABLE[c1 & 0x7ff]

def cocode_syndrome(uint32_t c1,  uint32_t u_tetrad = 24):
    cdef uint32_t res
    res = mat24_cocode_syndrome(c1, u_tetrad)
    if res == 0xffffffff:
        raise ValueError("Golay code syndrome is not unique")
    return res


def syndrome(uint32_t v1,  uint32_t u_tetrad = 24):
    cdef uint32_t res
    res =  mat24_syndrome(v1, u_tetrad)
    if res == 0xffffffff:
        raise ValueError("Golay code syndrome is not unique")
    return res



def vect_type(v1):
    cdef uint32_t res
    res = Mat24Sub_vtype(v1)
    return ( res >> 5, (res >> 2) & 7, res  & 3 )




def gcode_weight(uint32_t v1):
    return mat24_gcode_weight(v1)


def gcode_to_bit_list(uint32_t  v1):
    cdef uint8_t data[24]
    cdef uint32_t length
    length = mat24_gcode_to_bit_list(v1, &data[0])
    return list(data[:length])



def cocode_weight(uint32_t  c1):
    return mat24_cocode_weight(c1)


def cocode_to_bit_list(uint32_t  c1,  uint32_t u_tetrad = 24):
    cdef uint8_t data[4]
    cdef uint32_t length
    length = mat24_cocode_to_bit_list(c1, u_tetrad, &data[0])
    if length == 0xffffffff:
        raise ValueError("Golay code syndrome is not unique")
    return list(data[:length])
    
    
def cocode_to_sextet(uint32_t  c1):
    cdef uint8_t data[24]
    cdef uint32_t result
    result = mat24_cocode_to_sextet(c1, &data[0])
    if result == 0xffffffff:
        raise ValueError("Golay cocode word is not a sextet")
    return list(data[:24])


############################################################################
# Scalar product of Golay code and cocode
############################################################################



def scalar_prod(uint32_t v1, uint32_t c1):
    return mat24_scalar_prod(v1, c1)

############################################################################
# Conversion from and to suboctads
############################################################################


def suboctad_to_cocode(uint32_t u_sub, uint32_t v1):
    cdef uint32_t res
    res =  mat24_suboctad_to_cocode(u_sub, v1)
    if res & 0xfffff000:
        raise ValueError("Attempt to compute suboctad of a non-octad")
    return res


def cocode_to_suboctad(uint32_t c1, uint32_t v1):
    cdef uint32_t res
    res =  mat24_cocode_to_suboctad(c1, v1)
    if res & 0xfffff000:
        raise ValueError("Octad/suboctad mismatch")
    return res


def suboctad_weight(uint32_t u_sub):
    return mat24_suboctad_weight(u_sub)


def suboctad_scalar_prod(uint32_t u_sub1, uint32_t u_sub2):
    return mat24_suboctad_scalar_prod(u_sub1, u_sub2)

###########################################################################
# Represent a cocode element as a subset of a docecad
###########################################################################

def cocode_as_subdodecad(uint32_t c1, uint32_t v1, uint32_t u_single = 24):
    cdef uint32_t res
    res =  mat24_cocode_as_subdodecad(c1, v1, u_single)
    if res & 0xff000000:
        err = "Cocode element is not a subset of the docecad"
        raise ValueError(err)
    return res


###########################################################################
# Parker Loop
###########################################################################


def ploop_theta(uint32_t v1):
    return mat24_ploop_theta(v1)


def ploop_cocycle(uint32_t v1, uint32_t v2):
    return mat24_ploop_cocycle(v1, v2)



def mul_ploop(uint32_t v1, uint32_t v2):
    return mat24_mul_ploop(v1, v2) 


def pow_ploop(uint32_t v1, uint32_t u_exp):
    return mat24_pow_ploop(v1, u_exp) 



def ploop_comm(uint32_t v1, uint32_t v2):
    return mat24_ploop_comm(v1, v2)


def ploop_cap(uint32_t v1, uint32_t v2):
    return mat24_ploop_cap(v1, v2)


def ploop_assoc(uint32_t v1, uint32_t v2, uint32_t v3):
    return mat24_ploop_assoc(v1, v2, v3)


cdef enum:
     ploop_solve_size = 64  

def ploop_solve(a):
    cdef uint32_t pa[ploop_solve_size + 13]
    cdef unsigned int i
    cdef unsigned int res = 0
    cdef unsigned int length = len(a)
    cdef unsigned int start = 0
    cdef unsigned int d
    cdef unsigned int result_length = 0
    while start < length:
        d = min(length - start, ploop_solve_size)
        for i in range(d):
            pa[result_length + i] = a[start + i]
        res = mat24_ploop_solve(pa, result_length + d)
        result_length = res >> 16
        start += d
    if res & 0x1000:
        err = "Cannot correct signs of all Parker loop elements"
        raise ValueError(err)
    return res 

###########################################################################
# Mathieu group M24: conversion of representations
###########################################################################



def perm_complete_heptad(p_io):
    cdef uint8_t p1[24]
    cdef uint8_t *p_p1 = p1
    cdef int i
    for i in range(24): p1[i] = 24 if p_io[i] is None else int(p_io[i]) 
    cdef uint32_t res = mat24_perm_complete_heptad(p_p1)
    if res != 0:
        raise ValueError("Vector is not an umbral heptad")
    for i in range(24): p_io[i]  = p1[i]
    return res



def perm_check(p1):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef int i
    cdef uint32_t acc = 0
    cdef uint32_t entry
    if len(p1) != 24: 
        err = "Permutation in group Mat24 must have length 24"
        raise ValueError(err)
    for i in range(24): 
        entry = int(p1[i])
        p1a[i] = entry
        acc |= entry 
    if acc  & -32 or mat24_perm_check(p_p1a):
        err =   "Permutation is not in group Mat24" 
        raise ValueError(err)


def perm_complete_octad(p_io):
    cdef uint8_t p1[8]
    cdef uint8_t *p_p1 = p1
    cdef int i
    p1[5] = 24
    for i in range(5): p1[i] = p_io[i] 
    if p_io[5] is not None: p1[5] = p_io[5]
    cdef uint32_t res = mat24_perm_complete_octad(p_p1)
    if res:
         err = "Cannot complete a vector to an octad"
         raise ValueError(err)      
    for i in range(8): p_io[i]  = p1[i]


def perm_from_heptads(h1, h2):
    cdef uint8_t h1a[7]
    cdef uint8_t *p_h1a = h1a
    cdef uint8_t h2a[7]
    cdef uint8_t *p_h2a = h2a 
    cdef uint8_t p1[24]
    cdef uint8_t *p_p1 = p1
    cdef int i
    for i in range(7):
        h1a[i] = int(h1[i]) 
        h2a[i] = int(h2[i]) 
    cdef uint32_t res = mat24_perm_from_heptads(p_h1a, p_h2a, p_p1)
    if res:
        err = "Cannot construct permutation in Mat24 from heptads"
        raise ValueError(err)
    return [p1[i] for i in range(24)] 


def perm_from_map(h1, h2):
    cdef uint8_t h1a[24]
    cdef uint8_t *p_h1a = h1a
    cdef uint8_t h2a[24]
    cdef uint8_t *p_h2a = h2a 
    cdef uint8_t p1[24]
    cdef uint8_t *p_p1 = p1
    cdef uint32_t length = len(h1)
    if len(h2) != length:
        err = "Arrays in function perm_from_map() have different lengths"
        raise ValueError(err)
    if length > 24:
        err = "Array in function perm_from_map() is too long"
        raise ValueError(err)
    for i in range(length):
        h1a[i] = int(h1[i]) 
        h2a[i] = int(h2[i]) 
    cdef uint32_t res = mat24_perm_from_map(p_h1a, p_h2a, length, p_p1)
    if res < 0:
        err = "Illegal entry in array in function perm_from_map()"
        raise ValueError(err)
    return res, [p1[i] for i in range(24)] 
        


def m24num_to_perm(uint32_t u_m24):
    cdef uint8_t p1[24]
    cdef uint8_t *p_p1 = p1
    cdef int i
    if mat24_m24num_to_perm(u_m24, p_p1):
        err = "Illegal number for a permutation in Mat24"
        raise ValueError, err
    return [p1[i] for i in range(24)]




def perm_to_m24num(p1):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef int i
    for i in range(24): p1a[i] = int(p1[i]) 
    return  mat24_perm_to_m24num(p_p1a)


def perm_to_matrix(p1):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef uint32_t m1[12]
    cdef uint32_t *p_m1 = m1
    cdef int i
    for i in range(24): p1a[i] = int(p1[i])
    mat24_perm_to_matrix(p_p1a, p_m1)
    return [m1[i] for i in range(12)]


def matrix_to_perm(m1):
    cdef uint32_t m1a[12]
    cdef uint32_t *p_m1a = m1a
    cdef uint8_t p1[24]
    cdef uint8_t *p_p1 = p1
    cdef int i
    for i in range(12): m1a[i] = int(m1[i])
    mat24_matrix_to_perm(p_m1a, p_p1)
    return [p1[i] for i in range(24)]

def matrix_from_mod_omega(m1):
    cdef uint32_t m1a[12]
    cdef uint32_t *p_m1a = m1a
    cdef int i
    for i in range(12): m1a[i] = int(m1[i])
    mat24_matrix_from_mod_omega(m1a)
    for i in range(12): m1[i] = int(m1a[i])


###########################################################################
# Mathieu group M24: Mapping a dodecad
###########################################################################



def perm_from_dodecads(d1, d2):
    cdef uint8_t d1a[9]
    cdef uint8_t *p_d1a = d1a
    cdef uint8_t d2a[9]
    cdef uint8_t *p_d2a = d2a 
    cdef uint8_t p1[24]
    cdef uint8_t *p_p1 = p1
    cdef int i
    for i in range(9):
        d1a[i] = int(d1[i]) 
        d2a[i] = int(d2[i]) 
    cdef uint32_t res = mat24_perm_from_dodecads(p_d1a, p_d2a, p_p1)
    if res:
        err = "Cannot construct permutation in Mat24 from dodecads"
        raise ValueError(err)
    return [p1[i] for i in range(24)] 




###########################################################################
# Mathieu group M24: operation of group elements
###########################################################################


def op_vect_perm(uint32_t v1, p1):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef int i
    for i in range(24): p1a[i] = int(p1[i])
    return mat24_op_vect_perm(v1, p_p1a)



def op_gcode_matrix(uint32_t v1, m1):
    cdef uint32_t m1a[12]
    cdef uint32_t *p_m1a = m1a
    cdef int i
    for i in range(12): m1a[i] = int(m1[i])
    return mat24_op_gcode_matrix(v1,  p_m1a)


def op_gcode_perm(uint32_t v1, p1):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef int i
    for i in range(24): p1a[i] = int(p1[i])
    return mat24_op_gcode_perm(v1, p_p1a)




def op_cocode_perm(uint32_t c1, p1):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef int i
    for i in range(24): p1a[i] = int(p1[i])
    return mat24_op_cocode_perm(c1, p_p1a)



def mul_perm(p1, p2):
    cdef uint8_t pp1[24]
    cdef uint8_t *p_pp1 = pp1
    cdef uint8_t pp2[24]
    cdef uint8_t *p_pp2 = pp2
    cdef uint8_t pp3[24]
    cdef uint8_t *p_pp3 = pp3
    cdef int i
    for i in range(24): 
        pp1[i] = int(p1[i])
        pp2[i] = int(p2[i])
    mat24_mul_perm(p_pp1, p_pp2, p_pp3)
    return [pp3[i] for i in range(24)]


def inv_perm(p1):
    cdef uint8_t pp1[24]
    cdef uint8_t *p_pp1 = pp1
    cdef uint8_t pp2[24]
    cdef uint8_t *p_pp2 = pp2
    cdef int i
    for i in range(24): 
        pp1[i] = int(p1[i])
    mat24_inv_perm(p_pp1, p_pp2)
    return [pp2[i] for i in range(24)]


###########################################################################
# Automorphisms of the Parker Loop
###########################################################################
 


def autpl_set_qform(m_io):
    cdef uint32_t m_ioa[12]
    cdef uint32_t *p_m_ioa = m_ioa
    cdef int i
    for i in range(12): m_ioa[i] = int(m_io[i])
    mat24_autpl_set_qform(p_m_ioa)
    for i in range(12): m_io[i] = int(m_ioa[i])
    return m_io



def perm_to_autpl(c1, p1):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef uint32_t m1[12]
    cdef uint32_t *p_m1 = m1
    cdef int i
    for i in range(24): p1a[i] = int(p1[i])
    mat24_perm_to_autpl(int(c1), p_p1a, p_m1)
    return [m1[i] for i in range(12)]


def cocode_to_autpl(c1):
    cdef uint32_t m1[12]
    cdef uint32_t *p_m1 = m1
    mat24_cocode_to_autpl(int(c1), p_m1)
    return [m1[i] for i in range(12)]



def autpl_to_perm(m1):
    cdef uint32_t m1a[12]
    cdef uint32_t *p_m1a = m1a
    cdef uint8_t p1[24]
    cdef uint8_t *p_p1 = p1
    cdef int i
    for i in range(12): m1a[i] = int(m1[i])
    mat24_autpl_to_perm(p_m1a, p_p1)
    return [p1[i] for i in range(24)]


def autpl_to_cocode(m1):
    cdef uint32_t m1a[12]
    cdef uint32_t *p_m1a = m1a
    cdef int i
    for i in range(12): m1a[i] = int(m1[i])
    return mat24_autpl_to_cocode(p_m1a)



def op_ploop_autpl(v1, m1):
    cdef uint32_t m1a[12]
    cdef uint32_t *p_m1a = m1a
    cdef int i
    for i in range(12): m1a[i] = int(m1[i])
    return mat24_op_ploop_autpl(v1, p_m1a)


def mul_autpl(m1, m2):
    cdef uint32_t m1a[12]
    cdef uint32_t *p_m1a = m1a
    cdef uint32_t m2a[12]
    cdef uint32_t *p_m2a = m2a
    cdef uint32_t m3a[12]
    cdef uint32_t *p_m3a = m3a
    cdef int i
    for i in range(12): 
        m1a[i] = int(m1[i])
        m2a[i] = int(m2[i])
    mat24_mul_autpl(m1a, m2a, m3a)
    return [m3a[i] for i in range(12)]


def inv_autpl(m1):
    cdef uint32_t m1a[12]
    cdef uint32_t *p_m1a = m1a
    cdef uint32_t m2a[12]
    cdef uint32_t *p_m2a = m2a
    cdef int i
    for i in range(12): 
        m1a[i] = int(m1[i])
    mat24_inv_autpl(m1a, m2a)
    return [m2a[i] for i in range(12)]



def perm_to_iautpl(c1, p1):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef uint32_t m1[12]
    cdef uint32_t *p_m1 = m1
    cdef uint8_t pi[24]
    cdef uint8_t *p_pi = pi
    cdef int i
    for i in range(24): p1a[i] = int(p1[i])
    mat24_perm_to_iautpl(int(c1), p_p1a, p_pi, p_m1)
    return [pi[i] for i in range(24)], [m1[i] for i in range(12)]






###########################################################################
# Auxiliary functions for the Monster group
#
# These functions are not available in class mmgroup.dev.mat24_ref.Mat24
###########################################################################


def perm_to_net(p1):
    cdef uint8_t p[24]
    cdef uint8_t *p_p = p
    cdef uint32_t res[9]
    cdef uint32_t *p_res =res
    cdef int i
    for i in range(24):
        p[i] = p1[i]
    mat24_perm_to_net(p_p, p_res)
    return [res[i] for i in range(9)]





def op_all_autpl(m1):
    cdef uint32_t pm[12]
    cdef uint32_t *p_pm = pm
    cdef uint16_t res[2048]
    cdef uint16_t *p_res = res
    cdef int i
    for i in range(12):
        pm[i] = m1[i]
    mat24_op_all_autpl(p_pm, p_res)
    return [res[i] for i in range(2048)]



def op_all_cocode(c1):
    cdef uint8_t res[2048]
    cdef uint8_t *p_res = res
    cdef uint32_t c1i = c1
    cdef int i
    mat24_op_all_cocode(c1i, p_res)
    return [res[i] for i in range(2048)]




###########################################################################
# Functions from file mat24_random.ske
###########################################################################



cdef extern from "mat24_functions.h":
    uint32_t  MAT24_RAND_2, MAT24_RAND_o, MAT24_RAND_t
    uint32_t  MAT24_RAND_s, MAT24_RAND_l, MAT24_RAND_3

MAT24_RAND = {
    '2' : MAT24_RAND_2,
    'o' : MAT24_RAND_o,
    't' : MAT24_RAND_t,
    's' : MAT24_RAND_s,
    'l' : MAT24_RAND_l,
    '3' : MAT24_RAND_3,
}



def complete_rand_mode(uint32_t u_mode):
    cdef uint32_t res = mat24_complete_rand_mode(u_mode)
    return res


def perm_in_local(p1):
    cdef uint8_t p1a[24]
    cdef int i
    for i in range(24): p1a[i] = int(p1[i])
    cdef uint32_t res = mat24_perm_in_local(p1a)
    return res

ERR_MAT24_RANDOM = "Generation of random permutation in Mat24 has failed"

def perm_rand_local(uint32_t u_mode, uint32_t u_rand):
    cdef uint8_t p1a[24]
    cdef uint8_t *p_p1a = p1a
    cdef int32_t res = mat24_perm_rand_local(u_mode, u_rand, p_p1a)
    if res:
        raise ValueError(ERR_MAT24_RANDOM)
    cdef int i
    return [p1a[i] for i in range(24)]

def m24num_rand_local(uint32_t u_mode, uint32_t u_rand):
    cdef int32_t res = mat24_m24num_rand_local(u_mode, u_rand)
    if res < 0:
        raise ValueError(ERR_MAT24_RANDOM)
    return res

def m24num_rand_adjust_xy(uint32_t u_mode, uint32_t v):
    cdef int32_t res = mat24_m24num_rand_adjust_xy(u_mode, v)
    if res < 0:
        raise ValueError(ERR_MAT24_RANDOM)
    return res




########################################################################
########################################################################
# The following stuff is experimental and not documented officially 
########################################################################
########################################################################



cpdef uint32_t Mat24Sub_vtype(uint32_t x):
    """returns type of bit vector in GF(2)**24. 

    This is either (w,d,y), 0 <= d <= 3,  or (w,4,t)

    Here w is the bit weight of vector v. In the first case, there
    is a unique syndrome s(x) of x of weight d, 0 <= d <= 3, and
    the weight of x & s(x) is equal to y in case w <= 12 and to
    d - y in case w >= 12.
    In the second case the syndrome may be any of a tetrad of 
    six vectors s_i of length 4 with disjoint support. In case 
    w = 12, the value  t is the number of vectors s_i completely 
    contained in x. For other values of w the value t is zero.
  
    It is well known that all vectors of the same type are
    in one orbit of the automorphism group of the Golay code.

    The tuple is returned in the form  (w << 5 ) + (d << 2) + y .
    """
    cdef uint32_t w, s, d, y,  rest
    w = mat24_bw24(x)
    s = mat24_syndrome(x,0)
    d = mat24_bw24(s)
    if d < 4: 
        y =  mat24_bw24(x & s) 
        if w >= 12: y = d - y
    elif w != 12: 
        y = 0
    else: 
        y, rest = 0 + (s == x & s),  x & ~s & 0xffffff 
        while rest:
            s =  mat24_syndrome(x, mat24_lsbit24(rest))
            y +=  s == x & s  
            rest &= ~s
    return ( (w << 5 ) + (d << 2) + y )



def Mat24Sub_count_vtypes():
     cdef uint32_t a[8000]
     cdef uint32_t n
     for n in range(8000):
         a[n] = 0 
     for n in range(0x1000000):
         a[Mat24Sub_vtype(n)] += 1     
     d =  {}
     for n in range(8000):
        if a[n] > 0:
            d[ (n>>5, (n>>2) & 7, n & 3) ] = a[n]
     return d

     
       
