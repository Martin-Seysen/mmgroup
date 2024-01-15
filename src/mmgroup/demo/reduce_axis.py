
from random import choice


from mmgroup.demo import Mm, Leech2, MmV15
from mmgroup.demo.reduce_sub import mat15_norm
from mmgroup.demo.reduce_sub import mat15_rank_3
from mmgroup.demo.reduce_sub import mat15_apply
from mmgroup.demo.reduce_sub import map_type4_to_Omega
from mmgroup.demo.reduce_sub import map_type2_to_standard
from mmgroup.demo.reduce_sub import map_feasible_type2_to_standard
from mmgroup.demo.reduce_sub import find_triality_element_for_axis
from mmgroup.demo.reduce_sub import vect15_S
from mmgroup.demo.reduce_sub import leech2_span
from mmgroup.demo.reduce_sub import leech2_rad
from mmgroup.demo.reduce_sub import find_in_Nx0



def get_axis_type(v):
    """Yet to be documented
    """
    assert isinstance(v, MmV15) and v.p == 15
    norm = mat15_norm(v)
    AXIS_TYPES_KNOWN_FROM_NORM = {
        2:'8B', 3:'4C', 5:'6C', 10:'12C', 13:'4B'
    }
    if norm in AXIS_TYPES_KNOWN_FROM_NORM:
        return AXIS_TYPES_KNOWN_FROM_NORM[norm]
    elif norm == 4:
        r3, l2  = mat15_rank_3(v, 2)
        if r3 == 23:
            if l2.type == 2:
                M_v_l2 = mat15_apply(v, l2)
                if M_v_l2 == 4:
                    return '2A'
                elif M_v_l2 == 7:
                    return '6A'
                else:
                    raise ValueError("Vector is not an axis")
            else:
                raise ValueError("Vector is not an axis")
        elif r3 == 2:
            return '10A'
        else:
            raise ValueError("Vector is not an axis")
    elif norm == 8:
        r3, _,  = mat15_rank_3(v, 0)
        if r3 == 8:
            return '2B'
        elif r3 == 24:
            return '10B'
        else:
            raise ValueError("Vector is not an axis")
    elif norm == 14:
        r3, _, = mat15_rank_3(v, 0)
        if r3 == 8:
            return '6F'
        elif r3 == 23:
            return '4A'
        else:
            raise ValueError("Vector is not an axis")
    else:
        raise ValueError("Vector is not an axis")



def axis_leech2_vectors(v, axis_type):
    """Yet to be documented
    """
    assert isinstance(v, MmV15) and v.p == 15
    if axis_type == '2A':
        return []    
    if axis_type == '2B':
        return leech2_span(vect15_S(v, 4))
    if axis_type == '4A':
        _, l2  = mat15_rank_3(v, 0)
        return [l2]
    if axis_type in ['4B','4C']:
        return leech2_rad(vect15_S(v, 1))
    if axis_type == '6A':
        _, l2  = mat15_rank_3(v, 2)
        return [x + l2 for x in vect15_S(v, 5)]
    if axis_type == '6C':
        return leech2_span(vect15_S(v, 3))
    if axis_type in ['6F', '12C']:
        return leech2_rad(vect15_S(v, 7))
    if axis_type == '8B':
        S1 = (vect15_S(v, 1))
        l2 = choice(S1)
        return [l2 + x for x in S1]
    if axis_type == '10A':
        S1 = (vect15_S(v, 1))
        S3 = (vect15_S(v, 3))
        l2 = choice(S3)
        return [l2 + x for x in S1]
    if axis_type == '10B':
        return leech2_rad(vect15_S(v, 4))
       


TARGET_AXES_TYPES = {
   '2B' : ['2A'],
   '4A' : ['2A'],
   '4B' : ['2B'],
   '4C' : ['2B'],
   '6A' : ['4A'],
   '6C' : ['4A'],
   '6F' : ['4C'],
   '8B' : ['4A'],
   '10A' : ['6A'],
   '10B' : ['4B', '4C'],
   '12C' : ['4B', '6A'],
}


G_MINUS = Mm('x', 0x200)
assert MmV15('v-') * G_MINUS ==  MmV15('v+')

def reduce_axis(v):
    r"""Yet to be documented
    """
    v1 = v.copy()
    g = Mm(1)
    while True:
        axis_type = get_axis_type(v1)
        if axis_type == '2A':
            break
        leech2_vectors = axis_leech2_vectors(v1, axis_type)
        type4_vectors = [l2 for l2 in leech2_vectors if l2.type == 4]
        l2 = choice(type4_vectors)
        g_Gx0 = map_type4_to_Omega(l2)
        v1 = v1 * g_Gx0
        g = g * g_Gx0
        assert v * g == v1
        target_axes_types = TARGET_AXES_TYPES[axis_type]
        g_tau = find_triality_element_for_axis(v1, target_axes_types)
        v1 = v1 * g_tau
        g = g * g_tau
        assert v * g == v1
    
    _, l2 = mat15_rank_3(v1, 2)
    g_Gx0 = map_type2_to_standard(l2)
    g = g * g_Gx0
    v1 = v1 * g_Gx0
    if v1 != MmV15('v+'):
        v1 = v1 * G_MINUS
        g = g * G_MINUS
    assert v1 == MmV15('v+')
    return g       

 

BETA = Leech2('beta')
OMEGA = Leech2('Omega')

def reduce_baby_axis(v):
    r"""Yet to be documented
    """
    v1 = v.copy()
    g = Mm(1)
    while True:
        axis_type = get_axis_type(v1)
        if axis_type == '2A':
            break
        leech2_vectors = axis_leech2_vectors(v1, axis_type)
        feasible_type2_vectors = [
            l2 + BETA for l2 in leech2_vectors if
            l2.type == 4 and (l2 + BETA).type == 2
        ]
        l2 = choice(feasible_type2_vectors)
        g_Gx0 = map_feasible_type2_to_standard(l2)
        v1 = v1 * g_Gx0
        g = g * g_Gx0
        assert v * g == v1
        target_axes_types = TARGET_AXES_TYPES[axis_type]
        g_tau = find_triality_element_for_axis(v1, target_axes_types)
        v1 = v1 * g_tau
        g = g * g_tau
        assert v * g == v1
        
    _, l2 = mat15_rank_3(v1, 2)
    if l2 == BETA:
        assert v1 == MmV15('v-')
        return g

    g_Gx0 = map_feasible_type2_to_standard(l2)
    g = g * g_Gx0
    v1 = v1 * g_Gx0
    for e in 1, 2:
        v2 = v1 * Mm('t', e)
        if v2 == MmV15('v-'):
            return g * Mm('t', e)

    raise ValueError("Cannnot reduce axis for baby monster")



def reduce_G_x0(v):
    v1 = v.copy()
    r3, l2 = mat15_rank_3(v1, 0)
    assert r3 == 23
    g1 = map_type4_to_Omega(l2) 
    v1 = v1 * g1
    g2 = find_in_Nx0(v1)
    return g1 * g2
    

def reduce_monster_element(g):
    r"""Yet to be documented
    """
    v_plus = MmV15('v+') * g
    h = reduce_axis(v_plus)
    v_minus = MmV15('v-') * g * h
    h = h * reduce_baby_axis(v_minus)
    from mmgroup import MM
    assert (MM(g*h)).in_G_x0()

    v_1 = MmV15('v1') * g * h
    h = h * reduce_G_x0(v_1)
    return h
    

