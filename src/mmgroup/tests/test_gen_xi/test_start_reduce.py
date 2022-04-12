"""Test some C functions in file ``gen_leech.c``

In this script we test the C functions ``gen_leech2_start_type4``
and ``gen_leech2_start_type4`` in file  ``gen_leech.c``.

We use the terminology defined in
the document *The C interface of the mmgroup project*, 
section *Description of the mmgroup.generators extension*.

The tests in ths script are valid under the assumption that all
tests in the scripts ``test_leech2.py``, ``test_gen_type2``, and
``test_gen_type4`` marked with ``@pytest.mark.gen_xi`` have been 
passed. 
"""


from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, choices #, shuffle, sample
from numbers import Integral
from collections import defaultdict

import numpy as np
import pytest

from mmgroup import XLeech2, Xsp2_Co1, PLoop, GCode, AutPL, Cocode, GcVector

from mmgroup import mat24
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_subtype
from mmgroup.generators import gen_leech2_op_atom
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_start_type4
from mmgroup.generators import gen_leech2_start_type24

# Standard vector in the Leech lattice mod 2 in Leech lattice encoding
# The standard fram \Omega
OMEGA = 0x800000
# The standard type-2 vector \beta
BETA = 0x200 

assert Cocode(BETA).syndrome() == GcVector(0xC)

#####################################################################
# Test function gen_leech2_start_type4
#####################################################################


#####################################################################
# Reference implementation of gen_leech2_start_type4()

def ref_leech2_start_type4(v):
    """Reference implementation for function gen_leech2_start_type4

    This is equivalent to function ``leech2_start_type4```.
    """
    v &= 0xffffff
    t = gen_leech2_subtype(v)
    if (t & 0xf0) != 0x40:
        return min(-1, -(t >> 4))
    if v == OMEGA:
        return 0
    t2 = gen_leech2_subtype(v ^ BETA)
    return t2 if (t2 & 0xf0) == 0x20 else t

#####################################################################
# Implementation of gen_leech2_start_type4() similar to C function


def suboctad_type(octad, w, coc):
    """Return suboctad type.

    Let ``octad`` be an octad, i.e. a bit vector of length 8. Let
    w = 1 (mod 2) if ``octad`` denotes an octad and ``w = 0`` 
    if ``octad`` denotes a complemented octad. Let ``coc`` be
    an even cocode vector in cocode representation.

    The function returns 0 in bit 1 of the return value if the 
    cocode word ``coc`` can be written as a subset of the octad,
    and 1 in bit 1 otherwise.

    the function returns ``1 + w + bit_weight(coc)/2`` in bit 0
    of the return value.

    Then Leech lattice vector ``x_octad + x_coc`` is of subtype 
    0x22 if the return value is zero. Otherwise it is of
    subtype 0x44 (or 0x46) if bit 1 of the return value is 
    (or 1).   
    """
    lsb = mat24.lsbit24(octad)            # ls bit of octad
    syn = mat24.cocode_syndrome(coc, lsb) # cocode syndrome
    wsub = octad & syn == syn             # wsub = 1 if coc is suboctad
    cw = mat24.cocode_weight(coc) >> 1    # cw = cocode_weight(v) / 2
    return 2 * (1 - wsub) + ((w ^ cw ^ 1) & 1)



def leech2_start_type4(v):
    """Return subtype of a Leech lattice frame ``v`` used for reduction

    The function returns the subtype of a vector ``v`` of type 4 in
    the Leech lattice modulo 2. Parameter ``v2`` must be in Leech 
    lattice encoding. The function returns the subtype of ``v`` that
    will be used for reduction in function ``gen_leech2_reduce_type4``.

    In that function we take care of the special case that ``v + v0``
    is of type 2 for a specific short vector ``v0``. 
    
    A simpler (but slower) implementation of thhis function is:

    If ``v ^ v0`` is of type 2 the return the subtype of ``v ^ v0``.
    Otherwise return the subtype of ``v``.

    The function returns 0 if ``v`` is equal to  ``Omega`` and
    a negative value if ``v`` has not type 4.

    This is a refernece implementation for function 
    ``gen_leech2_start_type4()`` in file ``gen_leech.c``.
    """
    if v & 0x7ff800 == 0:
        # Then v or v + Omega is an even cocode element.
        # Return 0 if v == Omega and -1 if v ==  0.
        if v & 0x7fffff == 0:
            return 0 if v & 0x800000 else -1 
        # Let w be the cocode weight. Return -2 if w == 2.
        if mat24.cocode_weight(v) != 4:
            return -2
        # Here v has type 4. Let v23 be the standard type-2 vector.
        # Return 0x20 if v ^ v23 has type 2 and 0x40 otherwise.
        return 0x20 if mat24.cocode_weight(v ^ 0x200) == 2 else 0x40  
    if mat24.scalar_prod(v >> 12, v):
        # Then v has type 3 and we return -3
        return -3
    if v & 0x800:
        # Then the cocode word 'coc' of v is odd.
        coc = (v ^ mat24.ploop_theta(v >> 12)) & 0xfff
        syn = mat24.cocode_syndrome(coc)
        # If 'coc' has weight 1 then v is of type 2 and we return -2.
        if (syn & (syn - 1)) == 0:
            return -2
        # Here v is of type 4. 
        # Return 0x21 if v ^ v23 is of type 2 and 0x43 otherwise.
        if (syn & 0xc) == 0xc and (v & 0x200000) == 0:
            return 0x21
        return 0x43
    # Let w be the weight of Golay code part divided by 4
    w = mat24.gcode_weight(v >> 12) 
    if w == 3:
        # Then the Golay code part of v is a docecad and we return 0x46.
        return 0x46
    # Here the Golay code part of v is a (possibly complemented) octad.
    # Add Omega to v if Golay part is a complemented octad.
    v ^= (w & 4) << 21
    # Put w = 1 if that part is an octad and w = 0 otherwise.   
    w = (w >> 1) & 1

    # Let 'octad' be the octad in the Golay code word in vector rep.
    octad = mat24.gcode_to_vect(v >> 12)
    coc = v ^ mat24.ploop_theta(v >> 12)  # cocode element of v
    # Return -2 if v is of type 2.
    sub = suboctad_type(octad, w, coc)
    if sub == 0:
        return -2
    # Return 0x22 if v ^ v23 is shsort
    if suboctad_type(octad, w, coc ^ 0x200) == 0:
        return 0x22
    # Otherwise return the subtype of v
    return 0x44 if sub & 2 else 0x42


#####################################################################
# Test data for function gen_leech2_start_type4()


def type4_testdata(ntests = 1000):
    """The function yields type-4 vectors in the Leech lattice mod 2

    These type-4 vectors are given in **Leech lattice encoding**.
    """
    testdata = [
        0, 
        0x800,
        0xC,
    ]
    for t in testdata:
         yield t
    for type_ in range(2,4):
        for i in range(5):
            yield XLeech2('r', type_).ord 
    for i in range(24):
         yield (1 << i) + ((randint(1,10) << 24))
         yield (1 << i) ^ 0xffffff
    for i in range(ntests):
        yield XLeech2('r', 4).ord
    for i in range(100):
        v =  XLeech2('r', 2) * XLeech2(BETA)
        if v.type == 4:
            return v.ord

#####################################################################
# Test function gen_leech2_start_type4()


@pytest.mark.gen_xi
def test_start_type4(ntests = 1000, verbose = False):
    """Test function gen_leech2_start_type4()

    The function obtains test vectors from generator function
    type4_testdata(); the number of test vectors given by ``ntests``.

    Each test vector v is a type-4 vector in the Leech lattice mod 2.
    The test function checks that the C function
    ``gen_leech2_start_type4(v)`` returns the same result as both
    functions, ``leech2_start_type4(v)`` and
    ``ref_leech2_start_type4(v)``.          
    """
    for n, v in enumerate(type4_testdata(ntests)):
        if verbose:
            print("Test %d, v = %s, subtype = %s, subtpe(v2) = %s" % (
                n, hex(v), hex(gen_leech2_subtype(v)),
                hex(gen_leech2_subtype(v ^ BETA)) 
        ))
        t_py = leech2_start_type4(v)
        t_ref = ref_leech2_start_type4(v) 
        assert t_py == t_ref, (hex(t_py), hex(t_ref))
        t_c = gen_leech2_start_type4(v)
        assert t_c == t_ref, (hex(t_c), hex(t_ref))
    


#####################################################################
# Test function gen_leech2_start_type24
#####################################################################


def ref_leech2_start_type24(v):
    r"""Reference implementation for function gen_leech2_start_type24

    The function returns the subtype of a vector v of type 2 in
    the Leech lattice modulo, provided that  v + \beta  is of
    type 4. It returns 0 in the special case  v = \beta + \Omega,
    and a negative value if v is not of type 2 or v + \beta
    is not of type 4.
    """
    v &= 0xffffff
    t = gen_leech2_subtype(v)
    if (t & 0xf0) != 0x20:
        return -1
    if (v & 0x7fffff) == 0x200:
        return 0 if v & 0x800000 else -1
    t2 = gen_leech2_subtype(v ^ BETA)
    return t if (t2 & 0xf0) == 0x40 else -1


def type24_testdata(ntests):
    r"""The function yields type-2 vectors in the Leech lattice mod 2

    These type-2 vectors are orthogonal to the standard type-2 vector
    \beta in the real Leech lattice and given in **Leech lattice
    encoding**.

    The function also yields some vectors in the Leech lattice mod 2
    that are not of this shape, for testing illegal input for function
    ``gen_leech2_start_type24``.
    """
    XLeech2_BETA = XLeech2(BETA)
    testdata = [
        0, 
        0x800,
        0xC,
        0x200,
        0x800200,
    ]
    for t in testdata:
         yield t
    for type_ in range(2,5):
        for i in range(5):
            yield XLeech2('r', type_).ord 
    for i in range(24):
         yield (1 << i) + ((randint(1,10) << 24))
         yield (1 << i) ^ 0xffffff
    for i in range(ntests):
        yield XLeech2('r', 2).ord
    for i in range(200):
        v =  XLeech2('r', 2)
        if (v * XLeech2_BETA).type == 4:
            return v.ord



@pytest.mark.gen_xi
def test_start_type24(ntests = 1000, verbose = 0):          
    """Test function gen_leech2_start_type24()

    The function obtains test vectors from generator function
    type24_testdata(); the number of test vectors given by ``ntests``.

    Each test vector v is a type-2 vector in the Leech lattice mod 2.
    The test function checks that the C function 
    ``gen_leech2_start_type24(v)``  returns the same value as
    function ``ref_leech2_start_type4(v)``.          
    """
    for n, v in enumerate(type24_testdata(ntests)):
        if verbose:
            print("Test %d, v = %s, subtype = %s, subtpe(v2) = %s" % (
                n, hex(v & 0xffffff), hex(gen_leech2_subtype(v)),
                hex(gen_leech2_subtype(v ^ BETA)) 
        ))
        #t_py = leech2_start_type24(v)
        t_ref = ref_leech2_start_type24(v) 
        #assert t_py == t_ref, (hex(t_py), hex(t_ref))
        t_c = gen_leech2_start_type24(v)
        assert t_c == t_ref, (hex(t_c), hex(t_ref))
    

