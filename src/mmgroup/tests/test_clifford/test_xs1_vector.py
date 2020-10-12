from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix
from mmgroup.structures.xs1_co1 import Xs12_Co1, str_leech3
from mmgroup.structures.xs1_co1 import get_error_pool
from mmgroup.tests.spaces.clifford_space import Space_ZY
from mmgroup.clifford12 import xp2co1_chain_short_3, xp2co1_elem_to_qs
from mmgroup.clifford12 import xp2co1_short_2to3, xp2co1_short_3to2
from mmgroup.clifford12 import xp2co1_find_chain_short_3

MMSpace3 = Space_ZY.mmspace
MMGroup3 = MMSpace3.group



STD_V3  = 0x8000004





#####################################################################
# Convert matrix of type QStateMatrix to a complex matrix
#####################################################################

def create_test_vectors():
   vector_data = [ 
      (0x001, ("B", 2, 3)),
      (0x001, ("C", 15, 19)),
      (0x001, ("T", 643, 51)),
      (0x001, ("T", 199, 7)),
      (0x001, ("X", 0, 23)),
      (0x001, ("X", 137, 23)),
      (0x001, ("X", 1897, 1)),
   ]
   group_data = [
      [('d', 0x124)],
   ]
   for x4096, x24 in vector_data:
       for g in group_data:
           yield x4096, x24, g


def tuple_to_leech3(tag, i0, i1):
    sign = 1
    if tag[0] == '-':
        sign = -1
        tag = tag[1:]
    x2 = MMSpace3.index_to_short_mod2(tag, i0, i1)
    print("x2=", hex(x2))
    x3 = xp2co1_short_2to3(x2)
    print(list(map(hex, get_error_pool(15))))
    print("x3=", hex(x3))
    if sign == -1:
        x3 ^= 0xffffffffffff
    print("x3=", hex(x3))
    return x3
    


def short3_abs_equal(x1, x2):
    x1, x2 = int(x1), int(x2)
    return not bool((x1 ^ (x1 >> 24)  ^ x2 ^ (x2 >> 24)) & 0xffffff)

def map_v3(v, g, expected = None, verbose = 0):
    src = np.zeros(3, dtype = np.uint64)
    dest = np.zeros(3, dtype = np.uint64)
    src[0] = 0x8000004
    dest[0] = g.short3
    src[2] = v if isinstance(v, Integral) else tuple_to_leech3(*v)
    src[1] = xp2co1_find_chain_short_3(src[0], src[2])
    print(str_leech3(dest[0]))
    xp2co1_chain_short_3(xp2co1_elem_to_qs(g._data), src, dest)
    ok = dest[-1] != 0
    if not expected is None:
        ok = ok and short3_abs_equal(dest[-1], expected)
    if verbose or not ok:
        print("Map a short Leech vector v by an element g of Xs2Co1:")
        if not isinstance(v, Integral):
            print("v = %s = %s" % (v, str_leech3(src[2])))
        else:
            print("v = %s" % str_leech3(src[2]))
        print("g =", str(g))
        for i in range(3):
            print(" ",  str_leech3(src[i]), "->", str_leech3(dest[i]))
        for i in range(3):
            v2 =  xp2co1_short_3to2(src[i])       
            v2c = g.qs.pauli_conjugate(v2);
            v3c = xp2co1_short_2to3(v2c)
            print("   -> 0x%06x -> 0x%06x -> %s" %
                (v2, v2c, str_leech3(v3c)))
        print("")
        if not expected == None:
            print("Expected result: %s" % str_leech3(expected))
        print("Obtained result: %s\n" % str_leech3(dest[-1]))
        if not ok:
            ERR = "Mapping of Leech vector failed"
            raise ValueError(ERR)
    return dest[-1]


@pytest.mark.qstate1
def test_vector(verbose = 1):
    FORMAT_REDUCED = False
    for ntest, (x4096, x24, g) in enumerate(create_test_vectors()):
        if verbose:
            print("\nTEST %s" % (ntest+1))
            print("vector =", x4096, x24)
        vm = Space_ZY.unit(x4096, x24)
        if verbose:
            print(vm.as_tuples())
            print(vm)
            vm.dump()
        v3 = vm.as_mmspace_vector() 
        if verbose:
            print(g)
        g3 = MMGroup3(*g)
        if verbose:
            print(g)
        gm = Xs12_Co1(*g)
        if verbose:
            print("g3 = ", gm)
        try:
            wm = vm * gm
            map_v3(x24, gm, wm.short3, verbose = 0)
        except ValueError:
            map_v3(x24, gm, wm.short3, verbose = 1)
            raise
        w3_op =  wm.as_mmspace_vector()
        w3_mult =  v3 * g3
        if verbose:
            print(w3_op)
            print(w3_mult)
        #assert w3_op == w3_mult
            
            
            
            
