import numpy as np

from mmgroup import MMGroup, MMGroupWord
from mmgroup.mm_order import check_mm_order
from mmgroup.structures.xsp2_co1 import Xsp2_Co1



def mm_conjugate_2B(g, check=True, ntrials=20, verbose = 0):
    """Conjugte a 2B invlution in the monster into the center of G_x0

    Given a 2B involution :math:`g` in the monster group, the function 
    tries to find an element :math:`h` with :math:`h^{-1} g h = z`,
    where :math:`z` is te central involution in the subgroup
    :math:`G_{x0}` of the monster.
    """
    if not isinstance(g, MMGroupWord):
        err = "Object must be element of the monster of type MMGroupWord"
        raise TypeError(err)
    m = g.group
    m_1, z = m(), m(('x', 0x1000))
    in_G_x = g.in_G_x0()
    if (verbose or check) and g*g != m_1:
        err = "Element is not an involution in the monster group"
        raise ValueError(err)
    if in_G_x:
        return Xsp2_Co1(g).conjugate_2B_involution(m)
    if verbose:
        print("Function mm_conjugate_2B: conjugate g to the centre")
        print("g=", g)
    for i in range(ntrials):
        if verbose:
             print("\nmm_conjugate_2B trial", i)
        s = m.rand_mm(1 + max(3, (i >> 2))) if i else m_1
        x = g**s
        o, y = (x * z).half_order(60)
        if verbose:
           print("x = g**s; s =", s)
           print("x * z has order o =", (x * z).order())
        if o == 0 or o & 1 or y is None:
            continue
        # Here o is even and y = (x * z)**(o/2) is an involution,
        # and y commutes with x and with z. Thus y is in G_x0.
        # Next we represent y as an element of G_x0.
        y.in_G_x0()
        if verbose:
            assert (x*z)**(o>>1) == y
            assert y**2 == m_1
            print("y = (x*z)**(o/2) = ", y)
            print("y has characters",  y.chi_G_x0())
        # Try to conjugate y to the central element z.
        # Continue the loop if this fails.
        try:
            h1 = Xsp2_Co1(y).conjugate_2B_involution(m)
        except ValueError:
            continue
        # Now y**h1 = z.
        if verbose:
            print("h1=", h1)
            y_conj = y**h1
            y_conj.in_G_x0()
            assert y_conj == z, y_conj
        x1 = x**h1
        # x1 commutes with y**h1 = z; so x1 is in G_x0.
        # Next we represent x1 as an element of G_x0
        x1.in_G_x0()
        if verbose:
            print("x1=", x1)
            print("x1 has characters",  x1.chi_G_x0())
        # Try to conjugate x1 to the central element z.
        # Raise ValueError the loop if this fails.
        h2 = Xsp2_Co1(x1).conjugate_2B_involution(m)
        # Now x1**h2 = x**(h1*h2) = g**(s*h1*h2) = z.
        # So we may return s*h1*h2
        t = (s*h1*h2).reduce()
        t.in_G_x0()
        t.reduce()
        if verbose:
            print("Result found:")
            print(t)
            print("Function mm_conjugate_2B terminated successfully.\n")
        return t
    err = "Conjugation of element to central involution failed"
    raise ValueError(err)
