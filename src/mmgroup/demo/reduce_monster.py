r"""This module demonstrates the reduction algorithm for the Monster

The main function **reduce_monster_element** in this module reduces
an element of the Monster group. Given an word :math:`g` of
generators of the Monster of an arbitrary length, the function
returns an word :math:`h` of generators of the Monster of length
at most 7, such that :math:`g \cdot h` is the neutral element
of the Monster. Here the length of a word is is the number of
triality elements :math:`\tau^{\pm 1}` contained in that word.
"""

from mmgroup.demo import Mm, Leech2, MmV15
from mmgroup.demo.reduce_sub import mat15_rank_3
from mmgroup.demo.reduce_sub import map_type4_to_Omega
from mmgroup.demo.reduce_sub import find_in_Nx0
from mmgroup.demo.reduce_axis import reduce_axis
from mmgroup.demo.reduce_feasible import reduce_feasible_axis




def reduce_G_x0(v):
    r"""Reduce an element g of the subgroup G_x0 of the Monster

    :param v: Vector in representation of the Monster, see below
    :type v: class MmV15
    :return: Element g of G_x0 as defined below
    :rtype:  class Mm

    Let v be a vector in the representation of the Monster
    satisfying the conditition

    v * g = MmV15('v1')

    for an unknown element g of the group G_x0. The function computes
    g if such an element g exists. Otherwise it raises an error.
    """
    v2 = v.copy()                           # Local copy of v

    # Compute kernel of the matrix corresponding to part 300_x of v;
    # that kernel must contain a unique Leech lattice vector l2 of type 4
    r3, l2 = mat15_rank_3(v2, 0)
    assert r3 == 23                         # Rank of 300_x must be 23

    # Compute element g1 of G_x0 that transforms l2 to \lambda_\Omega
    g1 = map_type4_to_Omega(l2) 
    v2 = v2 * g1                            # Transform v2 with g1

    # Now v2 * g2 = MmV15('v1') for a (uinque) element g2 of the group N_x0
    g2 = find_in_Nx0(v2)                    # Finding g2 is easy
    assert v * g1 * g2 ==  MmV15('v1')      # This is what we expect
    return g1 * g2
    

def reduce_monster_element(g):
    r"""Reduce an element g of the Monster

    :param g: Element of the Monster of arbitrary length
    :type g:  class Mm
    :return: Element h of length at most 7 with g * h == 1
    :rtype:  class Mm
    """
    # Compute h such that g * h is in the subgroup H^+ of the Monster
    v_plus = MmV15('v+') * g
    h = reduce_axis(v_plus)                   # Now g * h is in H^+
    # Compute h such that g * h is in the subgroup H of G_x0
    v_minus = MmV15('v-') * g * h
    h = h * reduce_feasible_axis(v_minus)     # Now g * h is in H
    # Reduce the element g * h in the group G_x0
    v_1 = MmV15('v1') * g * h
    h = h * reduce_G_x0(v_1)                  # Now g * h is 1
    return h
    


def check_reduce_monster_element():
    r"""Generate a random element g of the Monster and test reduction of g
    """
    g = Mm('r', 14)                           # Random Monster element g
    assert g.count_triality_elements == 14    # Check length of g
    h = reduce_monster_element(g)             # Reduce g; the result is h
    assert MmV15('v1') * g * h == MmV15('v1') # Check that g * h is 1
    assert h.count_triality_elements <= 7     # Length of h must be <= 7
    

