"""Test C function dealing with Leech lattice vectors mod 3 of type 2

In this script we test functions dealing with vectors of type 2 in the 
Leech lattice modulo 3. These functions are implemented in file 
gen_leech3.c and available in the extension mmgroup.generators.

We use the terminology defined in
the document *The C interface of the mmgroup project*, 
section *Description of the mmgroup.generators extension*.
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint, choices #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral
from collections import defaultdict
from multiprocessing import Pool, TimeoutError

import numpy as np
import pytest

from mmgroup import MM0
from mmgroup.mat24 import MAT24_ORDER, ploop_theta
from mmgroup.mat24 import bw24 as mat24_bw24
from mmgroup.generators import gen_leech3to2
from mmgroup.generators import gen_leech3to2_short
from mmgroup.generators import gen_leech2to3_short
from mmgroup.generators import gen_leech3_neg
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_subtype
from mmgroup.generators import gen_leech2_type2

from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
 

from mmgroup.tests.test_gen_xi.test_gen_type4 import str_v3


#####################################################################
# Creating test vectors
#####################################################################




def type2_testvectors():
    """Yield short test vectors in the Leech lattice mod 2

    The function yields pairs ``(v2, vtype)``, where ``v2 ``is a short
    vector (of type 2) in the Leech lattice mod 2 in **Leech lattice
    encoding**, and ``vtype`` is the subtype of ``v``.

    This function uses function from file mm_aux.c for converting
    random indices of the rep 196884x to short vectors.
    """
    TEST_PARS = [
        (300, 2 * 276, 50, 0x20),
        (852, 759 * 64, 500, 0x22),
        (49428, 2048 * 24, 500, 0x21),
    ] 
    for (start, nvectors, ntests, vtype) in TEST_PARS:
         for i in range(ntests):
             v_extern = randint(start, start + nvectors - 1)
             v_sparse = mm_aux_index_extern_to_sparse(v_extern)
             v2 = mm_aux_index_sparse_to_leech2(v_sparse)
             assert v2 != 0, (v_extern, hex(v_sparse), hex(v2))
             #print(v_extern, hex(v_sparse), hex(v2), hex(vtype))
             yield v2, vtype


#####################################################################
# Test mapping of type 2 Leech vectors modulo 3 to vectors modulo 2
#####################################################################


def v3_to_v2(v3):
    """Map type-2 vector from Leech lattice mod 3 to Leech lattice mod 2

    Here parameter ``v3`` is a type-2 vector in the Leech lattice mod 2 
    in Leech lattice mod 3 encoding. 

    The function returns a type-2 vector in the Leech lattice mod 2
    corresponding to ``v3`` in Leech lattice encoding.

    The result is unique. The function returns 0 if ``v3`` is not of
    type 2 in the Leech lattice mod 3
    
    Tis function is a wrapper for the C function ``gen_leech3to2_short`` 
    in file ``gen_leech3.c``. 
    """
    result = gen_leech3to2_short(v3)
    assert result != 0, (str_v3(v3), weight_v3(v3), hex(result))
    result1 = gen_leech3to2(v3)
    assert result1 >> 24 == 2
    assert result1 & 0xffffff == result
    return result

def v2_to_v3(v2):
    """Map type-2 vector from Leech lattice mod 2 to Leech lattice mod 3

    This function reverses the peration of function ``v3_to_v2``.
    The result is unique up to sign only.
    
    Tis function is a wrapper for the C function ``gen_leech2to3_short`` 
    in file ``gen_leech3.c``. 
    """
    result = gen_leech2to3_short(v2)
    assert result != 0, (hex(v2), hex(result))
    return result


@pytest.mark.gen_xi
def test_type3(verbose = 0):
    r"""Test conversion of type-2 vectors 

    The function converts random sshort vectors in the Leech
    lattice mod2 2 to vectors in the Leech lattice mod 3, and
    vice versa, and check that these conversions are consistent. 

    It also checks functions in module ``checks  
    """
    weights = defaultdict(int)
    for ntest, (v2, vtype) in enumerate(type2_testvectors()):
        if verbose:
            print("\nTEST %s" % (ntest+1))
            print("v2 = " , hex(v2))
        assert gen_leech2_subtype(v2) == vtype
        assert gen_leech2_type2(v2) == vtype

        v3 = v2_to_v3(v2) 
        assert v3_to_v2(v3) == v2
        v3n = gen_leech3_neg(v3)
        assert v3_to_v2(v3n) == v2
    
    
    
