from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__
from numbers import Integral
import time
from math import floor
from random import randint, shuffle, sample
from collections import defaultdict
from itertools import chain
from collections.abc import Iterable

import numpy as np
import pytest


if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import MM0, AutPL, PLoop, Cocode, Xsp2_Co1, MMV, XLeech2
from mmgroup.generators import gen_leech2_type
from mmgroup.tests.test_involutions.make_involution_samples import invariant_count_type2
from mmgroup.clifford12 import xsp2co1_elem_find_type4
from mmgroup.clifford12 import xsp2co1_involution_find_type4
from mmgroup.clifford12 import xsp2co1_elem_conj_G_x0_to_Q_x0
from mmgroup.clifford12 import chk_qstate12
from mmgroup.clifford12 import xsp2co1_elem_involution_class
from mmgroup.clifford12 import xsp2co1_involution_invariants
from mmgroup.clifford12 import xsp2co1_leech2_count_type2
from mmgroup.tests.test_involutions.test_involution_invariants import INVOLUTION_SAMPLES



# Character chi_196883 of the monster group taken from the ATLAS 
# for some classes.
# E.g MONSTER_CHARACTERS[4][0] is the character of class 4A,
# MONSTER_CHARACTERS[4][1] is the character of class 4B, etc.
MONSTER_CHARACTERS = {
   1: (196883,),
   2: (4371, 275),
   4: (275, 51, 19, -13),
   8: (35, 11, -1, -5, 3, -1),
}


def get_monster_character(cl):
    order = (cl >> 4) & 0xf
    no = (cl & 0xf) - 1
    return MONSTER_CHARACTERS[order][no]


def Co1_class(g):
    invar = np.zeros(12, dtype = np.uint64)
    n = xsp2co1_involution_invariants(g._data, invar)
    if n in [0,1]:
        return 0
    if n in [8,9]:
        return 1
    if  n == 12:
        invar[0] = 0
        n2 = xsp2co1_leech2_count_type2(invar, 12)
        #print(g, n2)
        return 3 if n2 else 2
    raise ValueError("Cannot compute Co_1 class")



INVOLUTION_MAP = {}


@pytest.mark.involution
def test_xsp2co1_elem_involution_class(verbose = 0):
    global INVOLUTION_MAP
    minus = Xsp2_Co1('x', 0x1000)
    length, classes = 0, set()
    for data, s in INVOLUTION_SAMPLES:
        g = Xsp2_Co1(s)
        cl = xsp2co1_elem_involution_class(g._data)
        if verbose:
             print("%04x" % cl, data)
        classes.add(cl)
        length += 1
        for j in range(100):
            h = g**Xsp2_Co1('r', 'G_x0')
            #print("   ", h)
            assert cl == xsp2co1_elem_involution_class(h._data)
        chi = get_monster_character(cl)
        assert chi == data[1][0], (chi, data[1][0])
        Co1_cls = Co1_class(g)
        assert (cl >> 8) & 0xf == Co1_cls, ( (cl >> 8) & 0xf, Co1_cls)
        cl_neg = xsp2co1_elem_involution_class((g * minus)._data)
        fused = cl == cl_neg
        assert (cl >> 12) & 1 == (not fused), (hex(cl), hex(cl_neg))
        o = g.order() 
        assert o == (cl >> 4) & 0xf
        sqrt_minus1 = bool(cl >> 13)
        if sqrt_minus1:
             assert o & 1 == 0, o
             assert g ** (o//2) == minus 
        elif  o & 1 == 0:
             assert g ** (o//2) != minus    
        INVOLUTION_MAP[cl] = data[1]     
    
    assert len(classes) == length



def int_type(data_list, signed = False):
    assert len(data_list) > 0
    if isinstance(data_list[0], Iterable):
        data_list = list(chain(*(data_list[:])))
    if signed:
        lmax = max([abs(x) for x in data_list])
        prefix = "int"
    else:
        assert min(data_list) >= 0
        lmax = max(data_list)
        prefix = "uint"
    assert lmax < 1 << 32
    t = 32
    for i in (16, 8):
        if lmax < 1 << i:
            t = i
    return prefix + str(t) + "_t"
    


def make_C_table(name, d):
    print("#define LEN_%s_TABLE %d" % (name, len(d)))
    pre_im = sorted(d.keys())
    im = [d[x] for x in pre_im]
    min_pre_im = min(pre_im)
    max_pre_im = max(pre_im)
    if isinstance(im[0], Iterable):
        s = []
        for x in im :
            s.append( "{" + ", ".join(map(str,x)) + "}")
        str_im =  ", ".join(s) 
        #print(str_im)
        dim = "[%d][%s]" % (len(d), len(im[0]))
    else:
        max_im , min_im = max(im), min(im)
        str_im = ", ".join([str(x) for x in im]) 
        dim = "[%d]" % (len(d))
    str_pre_im = ", ".join(map(hex, pre_im)) 
    i_type = int_type(pre_im)
    print("static %s KEYS_%s_TABLE[%d] = {\n %s\n};" 
         % (i_type, name, len(d), str_pre_im))     
    i_type = int_type(im, signed = True)
    print("static %s DATA_%s_TABLE%s = {\n %s\n};" 
         % (i_type, name, dim, str_im))     
        


def display_involution_characters():
    s = """// The following tables have been created automatically
// by executing file ``test_display_characters.py`` 
// in module ``mmgroup.tests.test_involutions``.
"""
    print("\n" + s, end = "")
    I_SMALL = {} 
    I_MM = {}
    I_4096 = {}
    I_98280 = {}
    for inv, chi in INVOLUTION_MAP.items():
         i_small = (inv >> 8) & 0xf
         assert 0 <= i_small < 4
         img_i_small = chi[1:3]
         i_mm = inv & 0xff
         img_i_mm = chi[0]
         i_4096 = inv
         img_i_4096 = chi[3]
         i_98280 = inv
         img_i_98280 = chi[0] - chi[1] - chi[2] * chi[3]

         dicts = [I_98280, I_SMALL, I_MM, I_4096]
         keys = [i_98280, i_small, i_mm, i_4096]
         images = [img_i_98280, img_i_small, img_i_mm, img_i_4096]

         for d, k, im in zip(dicts, keys, images):
             if k in d:
                 assert d[k] == im
             else:
                 d[k] = im
    comments = [
        "Table for character chi_98280",
        "Table for characters chi_24 and chi_299",
        "Table for character chi_196883",
        "Table for character chi_4096",
    ]
    tables = dicts
    names = ["I_CHI_98280", "I_CHI_SMALL", "I_CHI_MM", "I_CHI_4096"]
    for comment, name, table in zip(comments, names, tables):
        print("\n// " + comment)
        make_C_table(name, table)


if __name__ == "__main__":
    test_xsp2co1_elem_involution_class(verbose = 1)
    #print(INVOLUTION_MAP)
    display_involution_characters()


