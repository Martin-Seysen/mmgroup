r"""Module **mmgroup.demo.reduce_feasible**: reduction of a feasible axis

Yet to be documented!!

Here the *feasible axes* are certain axes in the representation
:math:`\rho_{15}` the Monster, see :cite:`Seysen22`, Section 9.1
for details. They are in a one-to-one correspondence with the
left cosets of the subgroup :math:`H` of the group
:math:`H^+` defined in :cite:`Seysen22`.

Thus mapping an arbitrary feasible axis **v** to the standard
feasible axis :math:`v^-`  corresponding to the subgroup :math:`H`
is equivalent to mapping an element of :math:`H^+` to an element
of :math:`H` by right multiplication.

The process of finding an element :math:`g` that maps a feasible
axis **v** to :math:`v^-` is called *reduction* of the feasible
axis  **v**. Function **reduce_feasible_axis** in this module
reduces an feasible axis  **v**.

The algorithm in function **reduce_feasible_axis** for reducing
a feasible axis is discussed in :cite:`Seysen22`, Section 9.
This algorithm is quite similar to the algorithm in function
**reduce_axis** for reducing a general axis. Analogies
between these two algorithms are summarized in Table 2 in
that section.

More details yet to be documented!!


"""



from random import choice


from mmgroup.demo import Mm, Leech2, MmV15
from mmgroup.demo.reduce_axis import axis_orbit
from mmgroup.demo.reduce_axis import compute_U
from mmgroup.demo.reduce_axis import TARGET_ORBITS
from mmgroup.demo.reduce_sub import map_feasible_type2_to_standard
from mmgroup.demo.reduce_sub import find_triality_element_for_axis
from mmgroup.demo.reduce_sub import mat15_rank_3



def reduce_feasible_axis(v):
    r"""Return element of Monster reducing a feasible axis v

    Here reducing a fesible axis means reduction to the 
    standard feasible axis v^-.

    :param v: The feasible axis to be reduced
    :type v: class MmV15
    :return: Element g of the Monster with v * g = v^-
    :rtype: class Mm
    """
    BETA = Leech2('beta')    # Vector \lambda_\beta in leech lattice mod 2
    OMEGA = Leech2('Omega')  # Vector \lambda_\Omega in leech lattice mod 2
    v1 = v.copy()            # local copy of the feasible axis v
    g = Mm(1)                # the neutral element of the Monster

    # In g we will accumulate the element of the Monster that transforms v

    # Map axis to a 'simpler' orbit
    while True:
        orbit = axis_orbit(v1)
        if orbit == '2A':
            break

        # Compute the set U_f(v) and select a random element of that set
        U = compute_U(v1)
        U_f = [
            l2 + BETA for l2 in U if
            l2.type == 4 and (l2 + BETA).type == 2
        ]
        l2 = choice(U_f)     # a random element of U_f(v)


        # Find a Monster element g1 that maps v1 to a 'nice' axis
        # and map v1 to that 'nice' (feasible) axis
        g1 = map_feasible_type2_to_standard(l2)
        v1 = v1 * g1
        g = g * g1
        assert v * g == v1   # Correctness condition for loop

        # Find a Monster element g_tau that maps v1 to a 'simpler' axis,
        # and map v1 to that 'simpler' (feasible) axis
        target_orbits = TARGET_ORBITS[orbit]
        g_tau = find_triality_element_for_axis(v1, target_orbits)
        v1 = v1 * g_tau
        g = g * g_tau
        assert v * g == v1   # Correctness condition for loop
        
    # Now v has been transformed to an axis v1 in orbit '2A0' or '2A1'.
    # Compute the short Leech lattice vector l2 = \lambda(v2).
    _, l2 = mat15_rank_3(v1, 2)

    if l2 == BETA:
        # If l2 is \lambda_beta then v1 is in orbit '2A1' and we are done
        assert v1 == MmV15('v-')
        return g

    # Map v1 to an axis v1 with \lambda(v1) = \lambda\beta + \lambda\Omega
    g1 = map_feasible_type2_to_standard(l2)
    g = g * g1
    v1 = v1 * g1

    # Now a power of the triality element maps v1 to the
    # standard feasible axis v^-
    for e in [1, -1]:
        v2 = v1 * Mm('t', e)
        if v2 == MmV15('v-'):
            g = g * Mm('t', e)
            assert v * g == MmV15('v-')   # This is what we expect
            return g

    raise ValueError("Vector is not a feasible axis")



    

