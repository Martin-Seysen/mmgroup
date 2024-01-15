"""Yet to be documented"""

from mmgroup.demo import Mm, Leech2, MmV15
from mmgroup.demo.reduce_sub import mat15_rank_3
from mmgroup.demo.reduce_sub import map_type4_to_Omega
from mmgroup.demo.reduce_sub import find_in_Nx0
from mmgroup.demo.reduce_axis import reduce_axis
from mmgroup.demo.reduce_axis import reduce_baby_axis




def reduce_G_x0(v):
    r"""Yet to be documented
    """
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
    v_1 = MmV15('v1') * g * h
    h = h * reduce_G_x0(v_1)
    return h
    


def check_reduce_monster_element():
    r"""Yet to be documented
    """
    g = Mm('r', 14)
    assert g.count_triality_elements == 14
    h = reduce_monster_element(g)
    assert MmV15('v1') * g * h == MmV15('v1')
    assert h.count_triality_elements <= 7
    

