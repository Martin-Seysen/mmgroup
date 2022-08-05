
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

#####################################################################################
# Test data for wortd shortening
#####################################################################################



def word_shorten_testdata():
    testdata = [
       "", "pdx", "xt", "xtpx", "ltl",
    ] 
    for s in testdata:
        yield MM0([(tag, 'n') for tag in s])
    for r in range(0, 5):
        for k in range(3):
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


@pytest.mark.mmgroup 
def test_shorten_pyx(verbose = 0):
    print("")
    for i, g in enumerate(word_shorten_testdata()):
        if verbose:
            print("Test", i+1)
            print("g =", g)
        g_length = GtWord.n_subwords(g) 
        gtw = GtWord(g_length)
        gtw.append(g)
        if verbose:
            display_subwords("Subwords of g", gtw)
        g1 = gtw.mmdata(g.__class__)
        ok = g1 == g
        if not ok:
            err = "Error in shortening word g to g1"
            print(err)
            print("g =", g)
            display_subwords("Subwords of g", gtw)
            print("g =", g)
            raise ValueError(err)
        elif verbose:
            print("Shorter:", g1)
            


