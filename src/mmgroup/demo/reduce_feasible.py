r"""Module **mmgroup.demo.reduce_feasible**: reduction of a feasible axis


Here the *feasible axes* are certain axes in the representation
:math:`\rho_{15}` of the Monster, see :cite:`Seysen22`, Section 9.1
for details. They are in a one-to-one correspondence with the
left cosets of the subgroup :math:`H = G_{x0} \cap H^+` in the
group :math:`H^+` defined in :cite:`Seysen22`.

Thus mapping an arbitrary feasible axis **v** to the standard
feasible axis :math:`v^-`  corresponding to the subgroup :math:`H`
is equivalent to mapping an element of :math:`H^+` to an element
of :math:`H` by right multiplication.

The process of finding an element :math:`g` in :math:`H^+` that maps
a feasible axis **v** to :math:`v^-` is called *reduction* of the
feasible axis  **v**. Function **reduce_feasible_axis** in this
module reduces a feasible axis  **v**.

The algorithm in function **reduce_feasible_axis** for reducing
a feasible axis is discussed in :cite:`Seysen22`, Section 9.
This algorithm is quite similar to the algorithm in function
**reduce_axis** for reducing a general axis. Analogies
between these two algorithms are summarized in Table 2 in
that section.

An orbit of :math:`H^+` under the action of :math:`H` will be called
a :math:`H`-orbit. As in function **reduce_axis**, we repeatedly
multiply a feasible axis first with an element of :math:`H` and then
with a power of :math:`\tau`. This way we will map the axis into a
'simpler' :math:`H`-orbit in each step of the reduction process.

Details of that reduction process are discussed in :cite:`Seysen22`,
Section 9. Possible sequences of :math:`H`-orbits obtained during
such a reduction process are shown in Figure 3 in that section.

An orbit of the Monster under the action of :math:`G_{x0}` will be
called a :math:`G`-orbit. Note that disjoint :math:`H`-orbits
lie in disjoint :math:`G`-orbits, with the following exceptions.
:math:`H`-orbits **'2A0'** and **'2A1'** are in :math:`G`-orbit
**'2A'**; and :math:`H`-orbits **'2B0'** and **'2B1'** are in
:math:`G`-orbit **'2B'**.

The reduction process mentioned above terminates if the
feasible axis is in the :math:`H`-orbit **'2A0'** or **'2A1'**.
Orbit **'2A1'** is the singleton :math:`\{v^-\}`, so that we
are done in that case. Transforming a feasible axis from orbit
**'2A0'** to orbit **'2A1'** is easy.

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

    Here reducing a feasible axis means reduction to the
    standard feasible axis v^-.

    :param v: The feasible axis to be reduced
    :type v: class MmV15
    :return: Element g of subgroup H^+ of the Monster with v * g = v^-
    :rtype: class Mm
    """
    BETA = Leech2('beta')    # Vector \lambda_\beta in leech lattice mod 2
    OMEGA = Leech2('Omega')  # Vector \lambda_\Omega in leech lattice mod 2
    v1 = v.copy()            # Local copy of the feasible axis v
    g = Mm(1)                # Neutral element of the Monster

    # In g we will accumulate the element of H^+ that transforms v

    # Map axis to a 'simpler' orbit
    while True:
        orbit = axis_orbit(v1)
        if orbit == '2A':   # Done if we are in orbit '2A1' or '2A0'
            break

        # Compute the set U_f(v), as defined in [Sey22], Section 9.2;
        # and select a random element l2 of that set
        U = compute_U(v1)
        U_f = [
            l2 + BETA for l2 in U if
            l2.type == 4 and (l2 + BETA).type == 2
        ]
        l2 = choice(U_f)     # A random element of U_f(v)


        # Find a Monster element g1 that maps v1 to a 'nice' axis
        # and map v1 to that 'nice' (feasible) axis
        g1 = map_feasible_type2_to_standard(l2)
        v1 = v1 * g1         # Transform v1 with g1
        g = g * g1
        assert v * g == v1   # Correctness condition for loop

        # Find a Monster element g_tau that maps v1 to a 'simpler' axis,
        # and map v1 to that 'simpler' (feasible) axis
        target_orbits = TARGET_ORBITS[orbit]
        g_tau = find_triality_element_for_axis(v1, target_orbits)
        v1 = v1 * g_tau      # Transform v1 with g_tau
        g = g * g_tau
        assert v * g == v1   # Correctness condition for loop
        
    # Now v has been transformed to an axis v1 in orbit '2A0' or '2A1'.
    # Compute the short Leech lattice vector l2 = \lambda(ax(v2)).
    _, l2 = mat15_rank_3(v1, 2)

    if l2 == BETA:
        # If l2 is \lambda_beta then v1 is in orbit '2A1' and we are done
        assert v * g == MmV15('v-')       # This is what we expect
        return g

    # Map v1 to an axis with \lambda(ax(v1)) = \lambda\beta + \lambda\Omega
    g1 = map_feasible_type2_to_standard(l2)
    g = g * g1
    v1 = v1 * g1             # Transform v1 with g1

    # Now a power of the triality element maps v1 to the
    # standard feasible axis v^-
    for e in [1, -1]:
        v2 = v1 * Mm('t', e)
        if v2 == MmV15('v-'):
            g = g * Mm('t', e)
            assert v * g == MmV15('v-')   # This is what we expect
            return g

    raise ValueError("Vector is not a feasible axis")



    

