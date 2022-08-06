
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import sample
import datetime
import time


import pytest

import numpy as np
from multiprocessing import Pool, TimeoutError, cpu_count

from mmgroup.mm_reduce import mm_reduce_M
from mmgroup import MM0, MM
from mmgroup.mm_reduce import GtWord
from mmgroup.generators import gen_leech2_op_word_leech2
from mmgroup.generators import mm_group_invert_word

#####################################################################################
# Test data for wortd shortening
#####################################################################################



def word_shorten_testdata(ntests = 3):
    testdata = [
       "", "pdx", "xt", "xtpx", "ltl",
    ] 
    for s in testdata:
        yield MM0([(tag, 'n') for tag in s])
    for r in range(0, 5):
        for k in range(ntests):
            yield MM0('r', r)



#####################################################################################
# Test word shortening class
#####################################################################################


def display_subwords(text, w):
    print(text)
    for w in w.subwords():
        print("", MM0('a', w.data,), ('t', w.t_exp), ", img_O:",
            hex(w.img_Omega)
        )

OMEGA = 0x800000

def check_subwords(gtw, g_ref=None, is_reduced=False, verbose=0, text=""):
    messages = {
       1: "Inconsistent value of GtWord object!",
       2: "Value of GtWord object differs from expected value!",
       4: "Wrong image of OMEGA in GtWord object!",
       8: "GtWord object is not reduced!",
      16: "Bad atom in GtWord object!",
      32: "Inconsistent length indicator in GtWord object!",
    }
    g1 =  gtw.mmdata(MM0)
    w = gtw.subwords()
    g2 = MM0()
    err = 0
    for x in w:
        g2 *= MM0('a', x.data,) * MM0('t', x.t_exp) 
    err |= bool(g1 != g2) << 0
    ref_error = False
    if g_ref is not None:
        g_ref = MM0(g_ref)
        ref_error = g1 != g_ref
        err |= bool(ref_error) << 1
    for x in w:
         img_o = gen_leech2_op_word_leech2(OMEGA, x.data, len(x.data), 0)
         err |= bool(img_o != x.img_Omega) << 2
         if is_reduced:
            err |= bool(not reduced) << 3
         for atom in x.data:
            if not atom >> 28 in [1,2,3,4,6]:
                err |= 1 << 4
         err |= bool(x.length != len(x.data)) << 5
    if verbose  or err:
        print(text if text else "Subwords of GtWord object")
        for x in w:
            print("", MM0('a', x.data,), ('t', x.t_exp), ", img_O:",
                hex(x.img_Omega)
            )
        print("Value:", g1)
        if ref_error:
            print("Expected:", g_ref)
        if err:
            for key, text in messages.items():
                 if err & key:
                     print(text)
            raise ValueError("Error in GtWord object")
   


    if g_ref is not None:
        g_ref = MM0(g_ref)

@pytest.mark.mmgroup 
def test_shorten_pyx(ntests = 5, verbose = 0):
    print("")
    for i, g in enumerate(word_shorten_testdata(ntests)):
        if verbose:
            print("Test", i+1)
            print("g =", g)
        g_length = GtWord.n_subwords(g) 
        gtw = GtWord(g_length)
        gtw.append(g)
        check_subwords(gtw, g, verbose = verbose, is_reduced = False)
        gi = g.mmdata
        mm_group_invert_word(gi, len(gi))
        #print([hex(x) for x in gi])
        gtwi = GtWord(gi)  
        check_subwords(gtwi, g**-1, verbose = verbose, is_reduced = False)

