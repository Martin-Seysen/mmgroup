
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
    SAMPLE_TAGS = "dxyppllt" * 12
    test_elements = [
       [('l',1), ('l',2), ('t', 1), ('p',2)],
       [('l',1), ('l',2)],
       "M0<l_1*l_2*l_1*x_0cc1h*d_0e1eh*p_210218281*t_1>",
       "M0<x_1383h*p_184462007*l_2*l_1*y_709h*x_198dh*t_2>",
       "M0<y_0b17h*l_1*y_34h*x_1ba5h*d_84ah*p_87380788*l_2*l_2*l_1*y_1c5ch*d_691h*l_2*l_1*x_1029h*d_4b5h*l_1*t_2>",

    ]
    for g in test_elements:
        yield MM0(g)

    testdata = [
       "", "pdx", "xt", "xtpx", "ltl",
    ] 
    for s in testdata:
        yield MM0([(tag, 'n') for tag in s])
    for r in range(0, 5):
        for k in range(ntests):
            yield MM0('r', r)
    for i in range(5*ntests):
        yield MM0([(t, 'n') for t in sample(SAMPLE_TAGS, 30)])



#####################################################################################
# Test word shortening class
#####################################################################################


def display_subwords(text, w):
    print(text)
    fpos, words = w.subwords()
    for i, w in enumerate(words):
        s = "->" if i == fpos else "  "
        print(s, MM0('a', w.data,), ('t', w.t_exp), ", img_O:",
            hex(w.img_Omega)
        )

OMEGA = 0x800000

def check_subwords(gtw, g_ref=None, is_reduced=False, verbose=0, text=""):
    #print("CHECK", verbose, text)
    messages = {
       1: "Inconsistent value of GtWord object!",
       2: "Value of GtWord object differs from expected value!",
       4: "Wrong image of OMEGA in GtWord object!",
       8: "GtWord object is not reduced!",
      16: "Bad atom in GtWord object!",
      32: "Inconsistent length indicator in GtWord object!",
    }
    g1 =  gtw.mmdata(MM0)
    _, w = gtw.subwords()
    g2 = MM0()
    err = 0
    for x in w:
        g2 *= MM0('a', x.data) * MM0('t', x.t_exp) 
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
        text = text if text else "Subwords of GtWord object"
        display_subwords(text, gtw)
        print("Value:", g1)
        if ref_error:
            print("g in:", g_ref)
        if err:
            for key, text in messages.items():
                 if err & key:
                     print(text)
            raise ValueError("Error in GtWord object")
   


def py_load_word(g, check = True, verbose = True):
    g_length = GtWord.n_subwords(g) 
    gtw = GtWord(g_length)
    gtw.append(g)
    vb = verbose > 1
    gtw.seek(1,1)
    if check: 
       check_subwords(gtw, g, verbose = vb, text="start word")
    while not gtw.eof():
       success = gtw.rule_join()
       if not success:
          gtw.seek(1,0)
       if check: 
          check_subwords(gtw, g, verbose = vb, text="intermediate word")
    return gtw


#@pytest.mark.mmm 
@pytest.mark.mmgroup 
def test_shorten_pyx(ntests = 2, verbose = 0):
    for i, g in enumerate(word_shorten_testdata(ntests)):
        if verbose:
            print("\nTest", i+1)
            print("g =", g)
        gtw = py_load_word(g, check = True, verbose = verbose)

        check_subwords(gtw, g, verbose = verbose, text = "result word")
        if verbose:
            print("\nTest inverse of g")
        a = g.mmdata
        mm_group_invert_word(a, len(a))
        gi = MM0('a', a)
        if verbose > 1:
            print("mmdata:", [hex(x) for x in gi.mmdata])
        gtwi = py_load_word(gi, check = True, verbose = verbose)
        check_subwords(gtwi, g**-1, verbose = verbose, 
            text = "result word")


#####################################################################################
# Benchmark shortening in monster group
#####################################################################################


N_MUL_SAMPLES = 16

def mark_MM_element(g):
    def mark_atom(v):
        tag = (v >> 28) & 7
        return "TxE"[tag - 5] if tag >= 5 else "" 
    return "".join([mark_atom(v) for v in g.mmdata])


def make_mul_samples(n = N_MUL_SAMPLES, verbose = 0):
    indices = list(range(n))
    glist = [ MM('r', 'M').reduce() for i in range(n) ]
    if verbose:
        for g in glist:
            print(mark_MM_element(g))
    return indices, glist


def benchmark_shorten(ncases = 32, verbose = 0):
    indices, glist = make_mul_samples(verbose = verbose)
    t_start = time.process_time()
    for i in range(ncases):
        gtw = GtWord(glist[i & 15])  
    t = time.process_time() - t_start
    return t / ncases


#@pytest.mark.mmm 
@pytest.mark.bench 
@pytest.mark.mmgroup 
def test_benchmark_mul(ntests = 10000, verbose = 0):
    print("")
    for i in range(1):
        t = benchmark_shorten(ntests, verbose=verbose) 
        s = "Runtime of word shortening in class  MM, %d tests: %.3f us" 
        print(s % (ntests, 1000000*t))




