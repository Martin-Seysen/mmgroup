import numpy as np
import datetime
import time

from mmgroup import MMGroup, MM0, MM
from mmgroup.structures.mm_order import check_mm_order
from mmgroup.structures.xsp2_co1 import Xsp2_Co1



def mm_conjugate_involution(g, check=True, ntrials=20, verbose = 0):
    """Conjugte an involution in the monster to a standard element

    Given an involution :math:`g` in the monster group, the function 
    tries to find an element :math:`h` with :math:`h^{-1} g h = z`,
    where :math:`z` a standard representative of the involution class
    as in method ``conjugate_involution`` of class ``Xsp2_Co1``.
    """
    if not isinstance(g, (MM,MM0)):
        err = "Object must be element of the monster of type MM"
        raise TypeError(err)
    m = g.group
    m_1, z = m(), m('x', 0x1000)
    in_G_x = g.in_G_x0()
    if (verbose or check) and g*g != m_1:
        err = "Element is not an involution in the monster group"
        raise ValueError(err)
    if in_G_x:
        return Xsp2_Co1(g).conjugate_involution(m)
    if verbose:
        print("Function mm_conjugate_involution: conjugate g to the centre")
        print("g=", g)
    for i in range(ntrials):
        if verbose:
             print("\nmm_conjugate_involution trial", i)
        s = m('r', 1 + max(3, (i >> 2))) if i else m_1
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
        assert y.in_G_x0()
        if verbose:
            assert (x*z)**(o>>1) == y
            assert y**2 == m_1
            print("y = (x*z)**(o/2) = ", y)
            print("y has characters",  y.chi_G_x0())
        # Try to conjugate y to the central element z.
        # Continue the loop if this fails.
        try:
            itype = 3
            itype, h1 = Xsp2_Co1(y).conjugate_involution(m)
            assert itype == 2
        except (ValueError, AssertionError):
            if verbose:
                print("itype should be 2 but found %d" % itype) 
            continue
        # Now y**h1 = z.
        if verbose:
            print("itype=", itype)
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
        itype, h2 = Xsp2_Co1(x1).conjugate_involution(m)
        # Now x1**h2 = x**(h1*h2) = g**(s*h1*h2) = z.
        # So we may return s*h1*h2
        t = (s*h1*h2).reduce()
        t.in_G_x0()
        t.reduce()
        if verbose:
            print("Result found (itype = %s):" % itype)
            print(t)
            print("Function mm_conjugate_involution terminated successfully.\n")
        return itype, t
    err = "Conjugation of element to central involution failed"
    raise ValueError(err)


def reduce_via_power(g, ntrials=20, verbose = 0):
    if verbose:
        print("\nTrying to shorten the element g of the monster:")
        print("g =", g)
    for i in range(ntrials):
        x0 = g.group("r", max(1, (i >> 2)))
        #x0 = g.group("r", "G_x0")
        g1 = g * x0
        o, g2 = g1.half_order()
        if o & 1 or g2 is None:
            continue
        try:
            it, h = mm_conjugate_involution(g2, check=False, ntrials=1, 
                      verbose = 0)
            assert it == 2
        except (ValueError, AssertionError):
            if verbose:
                err = "Conjugation of element to central involution failed"
                print(err + "\n")
            continue
        # Now g2**h = z, so g1**h  is in G_x0
        g1h = g1**h
        g1h.in_G_x0()
        g1_new = g1h**(h**-1)
        if verbose:
            assert g1_new == g1
        g_new = g1_new * x0**-1
        g_new.reduce()
        if verbose:
            print("Shortened value of element g:")
            print(g_new)
            assert g_new == g
            s = "Element of the monster has been shortened successfully"
            print(s + "\n")
        return g_new
    err = "Shortening of an element of the monster failed"
    raise ValueError(err)

        

