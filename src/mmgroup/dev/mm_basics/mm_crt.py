

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
import warnings
from numbers import Integral


import numpy as np

from mmgroup.generate_c import UserDirective, UserFormat

from mmgroup.dev.mm_basics.mm_basics import INT_BITS
from mmgroup.dev.mm_basics.mm_basics import MM_Const
from mmgroup.structures.abstract_rep_space import mod_inv


MODULI = [7, 31, 127, 255]
N = np.prod(MODULI)
NMAX = N//2
INV = dict(zip(MODULI, [N//m * mod_inv(N//m, m) for m in MODULI]))
# So INV[p] = 1 (mod p),  INV[p] = 0 (mod N/p), p = 7, 31, 127, 255

class MM_CrtCombine:
    TAB_7_31 = np.zeros(256, dtype = np.uint32)
    MASK = 0xffffffff
    for i7 in range(8):
        for i31 in range(32):
            a = (i7 * INV[7] + i31 * INV[31] + NMAX) % N
            TAB_7_31[32 * i7 + i31] = a
    TAB_127 = np.zeros(128, dtype = np.uint32)
    for i127 in range(128):
        a = (i127 * INV[127]) % N - N
        TAB_127[i127] = a & 0xffffffff
    TAB_255 = np.zeros(256, dtype = np.uint32)
    for i255 in range(256):
        a = (i255 * INV[255]) % N - N
        TAB_255[i255] = a & 0xffffffff

    # Chinese remaindering of a7, a31, a127, and a255 works as follows:
    # a = TAB_7_31[32 * a7 + a31] + TAB_127[a127] + TAB_255[a255];
    # a += -((a & 0x80000000) >> 31) & N;
    # a += -((a & 0x80000000) >> 31) & N;
    # a -= NMAX;
    # Now a = a7 (mod7), a = a31 (mod 31), a = a127(mod 127),
    # a = a255(mod 255), -N/2 < a < N/2.
        
    tables = {
        "CRT_TAB_7_31" : TAB_7_31,
        "CRT_TAB_127" : TAB_127,
        "CRT_TAB_255" : TAB_255,
        "CRT_NPROD" : N,
    }

    directives = {}



class Tables(MM_CrtCombine):
    def __init__(self, *args, **kwds):
        pass

