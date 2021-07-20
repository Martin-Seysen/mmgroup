import numpy as np

from mmgroup import MMGroup, MMGroupWord
from mmgroup.mm_order import check_mm_order
from mmgroup.structures.xsp2_co1 import Xsp2_Co1



def mm_conjugate_2B(g, check=True, ntrials=3):
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
    m_1, m_z = m(), m(('x', 0x1000))
    in_G_x = g.in_G_x0()
    if check and g*g != m_1:
        err = "Element is not an involution in the monster group"
        raise ValueError(err)
    if in_G_x:
        return Xsp2_Co1(g).conjugate_2B_involution(m)
    for i in range(ntrials):
        s = m.rand_G_x0() if i else m_1
        x = g**s
        o, y = (x * z).half_order(60)
        if o == 0 or o & 1 or y is None:
            continue
        # Here o is even and y = (x * z)**(o/2) is an involution,
        # and y commutes with x and with z. Thus y is in G_x0.
        # Next we represent y as an element of G_x0.
        y.inG_x0()
        # Try to conjugate y to the central element z.
        # Continue the loop if this fails.
        try:
            h1 = Xsp2_Co1(y).conjugate_2B_involution(m)
        except ValueError:
            continue
        # Now y**h1 = z.
        x1 = x**h1
        # x1 commutes with y**h1 = z; so x1 is in G_x0.
        # Next we represent x1 as an element of G_x0
        x1.in_G_x0()
        # Try to conjugate x1 to the central element z.
        # Raise ValueError the loop if this fails.
        h2 = Xsp2_Co1(x1).conjugate_2B_involution(m)
        # Now x1**h2 = x**(h1*h2) = g**(s*h1*h2) = z.
        # So we may return s*h1*h2
        t = (s*h1*h2).reduce()
        t.in_G_x0()
        return t
    err = "Conjugation of element to central involution failed"
    raise ValueError(err)
