r"""Module **mmgroup.demo.reduce_axis** demonstrates the reduction of an axis


Here the *axes* are certain vectors in the representation
:math:`\rho_{15}` the Monster. The axes are in a one-to-one
correspondence with the left cosets of the subgroup :math:`H^+` of
structure :math:`2.B` of the Monster. Thus mapping an arbitrary
axis **v** to the standard axis :math:`v^+`  corresponding to the
subgroup :math:`H^+` is equivalent to mapping an element of the
Monster to an element of :math:`H^+` by right multiplication.
For background we refer to :cite:`Seysen22`, Section 7.

The process of finding an element :math:`g` that maps axis **v**
to :math:`v^+` is called *reduction* of the axis  **v**. Function
*reduce_axis* in this module reduces an axis  **v**.

We name an orbit of an axis under the action of the subgroup
:math:`G_{x0}` by a string as in :cite:`Seysen22`, Section 8.2.
In the sequel an *orbit* of an axis means an
orbit under the action of :math:`G_{x0}`. Function
**axis_orbit** in this module computes the orbit of an axis.

For reducing an axis we first multiply an axis with a suitable
element of :math:`G_{x0}` in order to obtain a 'nice' axis
in the same orbit. Roughly speaking, an axis is 'nice' if
multiplication with a suitable power of the triality element
:math:`\tau` of the Monster maps that axis into a 'simpler'
orbit.

This way we may repeatedly multiply an axis first with an
element of :math:`G_{x0}` and then with a power of :math:`\tau`,
leading to a 'simpler' orbit in each step of the reduction
process. Possible sequences of orbits obtained during such a
reduction process are shown in Figure 2 in
:cite:`Seysen22`, Section 8.3.

For each axis :math:`{\bf v}` there is a set :math:`U_4({\bf v})`
of vectors in the Leech lattice mod 2 such that for any
:math:`g \in G_{x0}` the axis :math:`{\bf v} g` is 'nice', if
there is an :math:`l_2 \in U_4({\bf v})` with
:math:`l_2 g = \lambda_\Omega`. For background we refer to
:cite:`Seysen22`, Section 8.3. Function **compute_U** in this
module is used to compute the set :math:`U_4({\bf v})`. More
precisely, calling **compute_U(v)** returns a list of vectors in
the Leech lattice mod 2; and :math:`U_4({\bf v})` is the set of
type-4 vectors in that list. Function **map_type4_to_Omega**
in module **mmgroup.demo.reduce_sub** computes :math:`g` from
:math:`l_2`.

At the end of that process we obtain an axis in the orbit
'2A'. That orbit also contains the axis :math:`v^+`. Mapping
an arbitrary axis in orbit '2A' to the standard axis
:math:`v^+` is easy.
"""

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


def axis_orbit(v):
    """Compute the orbit of an axis

    :param v: The axis to be checked
    :type v: class MmV15
    :return: Name of the orbit of axis v
    :rtype: str
    """
    ORBITS_KNOWN_FROM_NORM = {
        2:'8B', 3:'4C', 5:'6C', 10:'12C', 13:'4B'
    }
    assert isinstance(v, MmV15) and v.p == 15

    # Compute norm of the 24 times 24 matrix 300_x contained in vector v (mod 15)
    norm = mat15_norm(v)

    if norm in ORBITS_KNOWN_FROM_NORM:
        # If there is only one orbit with that norm, return that orbit
        return ORBITS_KNOWN_FROM_NORM[norm]
    elif norm == 4:
        # If 300_x has norm 4 then compute rank of matrix 
        # M = 300_x - 1 * M_1, where M_1 is the unit matrix
        r3, l2  = mat15_rank_3(v, 2)
        if r3 == 23:
            # If M has rank 23 then we check the kernel of M
            # contains a vector l2 of type 2 in Leech lattice.
            if l2.type == 2:
                # compute y = transposed(l2) * M * l2 (mod 15)
                y = mat15_apply(v, l2)
                if y == 4:
                    return '2A'  # orbit is '2A' if y == 4
                elif y == 7:
                    return '6A'  # orbit is '6A' if y == 7
                else:
                    raise ValueError("Vector is not an axis")
            else:
               raise ValueError("Vector is not an axis")
        elif r3 == 2:
            return '10A'         # orbit is '10A' if M has rank 2
        else:
            raise ValueError("Vector is not an axis")
    elif norm == 8:
        # If 300_x has norm 8 then compute rank of matrix 300_x
        r3, _,  = mat15_rank_3(v, 0)
        if r3 == 8:
            return '2B'          # orbit is '2A' if rank is 8
        elif r3 == 24:
            return '10B'         # orbit is '10B' if rank is 8
        else:
            raise ValueError("Vector is not an axis")
    elif norm == 14:
        # If 300_x has norm 14 then compute rank of matrix 300_x
        r3, _, = mat15_rank_3(v, 0)
        if r3 == 8:
            return '6F'          # orbit is '6F' if rank is 8
        elif r3 == 23:
            return '4A'          # orbit is '4A' if rank is 23 
        else:
            raise ValueError("Vector is not an axis")
    else:
        raise ValueError("Vector is not an axis")



def compute_U(v):
    """Compute a set of Leech lattice vectors from an axis v

     For a given axis v the function computes the set U(v) of
     vectors in the Leech lattice mod 2, as defined above. 

    :param v: The axis v
    :type v: class MmV15
    :return: The set U(v) of vectors in the Leech lattice (mod 2)
    :rtype: list[Leech2]
    """
    assert isinstance(v, MmV15) and v.p == 15
    orbit = axis_orbit(v)      # Compute the orbit of axis v
    if orbit == '2A':
        return []    
    if orbit == '2B':
        return leech2_span(vect15_S(v, 4))
    if orbit == '4A':
        _, l2  = mat15_rank_3(v, 0)
        return [l2]
    if orbit in ['4B','4C']:
        return leech2_rad(vect15_S(v, 1))
    if orbit == '6A':
        _, l2  = mat15_rank_3(v, 2)
        return [x + l2 for x in vect15_S(v, 5)]
    if orbit == '6C':
        return leech2_span(vect15_S(v, 3))
    if orbit in ['6F', '12C']:
        return leech2_rad(vect15_S(v, 7))
    if orbit == '8B':
        S1 = (vect15_S(v, 1))
        l2 = choice(S1)   # l2 is a random element of S1
        return [l2 + x for x in S1]
    if orbit == '10A':
        S1 = (vect15_S(v, 1))
        S3 = (vect15_S(v, 3))
        l2 = S3[0]        # S3 is a singleton here
        return [l2 + x for x in S1]
    if orbit == '10B':
        return leech2_rad(vect15_S(v, 4))
    raise ValueError('Unknown orbit')
       


TARGET_ORBITS = {
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



def reduce_axis(v):
    r"""Return element of Monster reducing an axis v to the standard axis v^+

    :param v: The axis to be reduced
    :type v: class MmV15
    :return: Element g of the Monster with v * g = v^+
    :rtype: class Mm

    Here v^+ is the standard axis. 
    """ 
    v1 = v.copy()                  # local copy of axis v
    g = Mm(1)                      # the neutral element of the Monster

    # We accumulate in g the element of the Monster that transforms v

    # Map axis to a 'simpler' orbit; see documentation of the module
    while True:
        orbit = axis_orbit(v1)     # orbit of current axis v
        if orbit == '2A':
            break                  # done if we ar in orbit '2A'

        # Compute the set U_4(v) and select a random element of that set
        leech2_vectors = compute_U(v1)
        type4_vectors = [l2 for l2 in leech2_vectors if l2.type == 4]
        l2 = choice(type4_vectors) # a random element of U_4(v)

        # Find a Monster element g1 that maps v1 to a 'nice' axis
        # and map v1 to that 'nice' axis
        g1 = map_type4_to_Omega(l2) 
        v1 = v1 * g1               # Transfrom v1 with g1
        g = g * g1
        assert v * g == v1         # Correctness condition for loop

        # Find a Monster element g_tau that maps v1 to a 'simpler' axis,
        # and map v1 to that 'simpler' axis
        target_orbits = TARGET_ORBITS[orbit] # possible 'simpler' orbits
        g_tau = find_triality_element_for_axis(v1, target_orbits)
        v1 = v1 * g_tau            # Transform v1 with g_tau
        g = g * g_tau
        assert v * g == v1         # Correctness condition for loop
    
    # Now axis v has been transformed to an axis v1 in orbit '2A'
    # Map the short Leech lattice vector l2 corresponding to axis v1
    # to the standard short vector \lambda_\beta.
    _, l2 = mat15_rank_3(v1, 2)     
    g_Gx0 = map_type2_to_standard(l2)
    g = g * g_Gx0
    v1 = v1 * g_Gx0

    assert v1 in [MmV15('v+'), MmV15('v-')] # Here v1 must be v^+ or v^-
    if v1 != MmV15('v+'):         
        G_MINUS = Mm('negate_beta') # Correct v1 if it is equal to v^-
        v1 = v1 * G_MINUS
        g = g * G_MINUS
    assert v1 == MmV15('v+')        # This is what we expect
    return g       

 


def reduce_baby_axis(v):
    r"""Yet to be documented
    """
    BETA = Leech2('beta')
    OMEGA = Leech2('Omega')
    v1 = v.copy()
    g = Mm(1)
    while True:
        orbit = axis_orbit(v1)
        if orbit == '2A':
            break
        leech2_vectors = compute_U(v1)
        feasible_type2_vectors = [
            l2 + BETA for l2 in leech2_vectors if
            l2.type == 4 and (l2 + BETA).type == 2
        ]
        l2 = choice(feasible_type2_vectors)
        g_Gx0 = map_feasible_type2_to_standard(l2)
        v1 = v1 * g_Gx0
        g = g * g_Gx0
        assert v * g == v1
        target_orbits = TARGET_ORBITS[orbit]
        g_tau = find_triality_element_for_axis(v1, target_orbits)
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



    

