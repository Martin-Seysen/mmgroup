
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import sample
import datetime
import time
from random import randint


import pytest

import numpy as np
from multiprocessing import Pool, TimeoutError, cpu_count

from mmgroup import MM0, MM, MM_from_int
from mmgroup.generators import gen_leech2_op_word_leech2
from mmgroup.generators import mm_group_invert_word
from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup import mat24
from mmgroup.mat24 import MAT24_ORDER
from mmgroup import mm_op
from mmgroup.bitfunctions import bitparity


def import_all():
    global mm_reduce_M, gt_word_shorten
    global GtWord
    from mmgroup.mm_reduce import mm_reduce_M, gt_word_shorten
    from mmgroup.mm_reduce import GtWord



#####################################################################################
# Decoding an integer to an element of the Monster
#####################################################################################


COCODE_0   =  0x800  # internal rep of cocode word [0]
COCODE_01  =  0x600  # internal rep of cocode word [0,1]
COCODE_STD =  0x200  # internal rep of cocode word [2,3]
OMEGA = 0x800000

ERR_DECODE_TYPE4 = "Could not decode %s to Leech lattice type-4 vector"
ERR_DECODE_TYPE2 = "Could not decode %s to Leech lattice type-2 vector"
ERR_DECODE_RANGE = "Integer to be decoded out or range"



def expand_23bit_type4(i):
    # Insert a zero bit into i at position k,
    # where k is the position of the lowest bit of i >> 11.
    # Shift the bits in i at positions > k left by one.
    i &= 0x7ffffff
    if (i & 0x7ff800) == 0:
        raise ValueError(ERR_DECODE_TYPE4 % hex(i))
    b = i >> 11 
    b = b & (0 - b)   # Bit k is set in b, other bits cleared
    i = (i & (b - 1)) | ((i & (0 - b)) << 1) 
    # flip bit at position k if parity of i & (i >> 12) is odd
    if bitparity(i & (i >> 12) & 0xfff):
        i ^= b
    # exchange bits 11 and 23
    j = (i ^ (i >> 12)) & 0x800
    i ^= j ^ (j << 12)

    
    type_ = gen_leech2_type(i)
    # Correct a vector of type 2 to a vector of subtype 40 or 48     
    if (type_ == 2):
        j =  (i >> 12) & 0x7ff;
        coc = (mat24.ploop_theta(j) ^ COCODE_0 ^ i) & 0xfff
        if coc == 0:
            assert j != 0, ERR_DECODE_TYPE4 % hex(i)
            return j
        if coc == COCODE_01:
            return j ^ OMEGA
        raise ValueError(ERR_DECODE_TYPE4 % hex(i))
    if (type_ == 4):
        return i
    raise ValueError(ERR_DECODE_TYPE4 % hex(i))
   
        
def expand_17bit_type2(i):
    i &= 0x1ffff
    i = mm_op.mm_aux_index_extern_to_sparse(i);
    assert i > 0, ERR_DECODE_TYPE2 % hex(i)
    i = mm_op.mm_aux_index_sparse_to_leech2(i);
    assert i > 0, ERR_DECODE_TYPE2 % hex(i)
    return i;



def extract_int(n, start, length, tag = 0, verbose = 0):
    value = ((n >> start) & ((1 << length) - 1)) | (tag << 28)
    if verbose:
         print("extract", length, hex(value))
    return value, start + length

def decompress_to_array(n, verbose = 0):
    LEN_A = 80
    a = np.zeros(LEN_A, dtype = np.uint32)
    if n <= 0 or n >= (1 << 255) or (n & ((1 << 64) - 1) == 0):
        raise ValueError(ERR_DECODE_RANGE)
    p, pos_n = extract_int(n, 0, 28, 0, verbose)
    with_t = 0
    if p < MAT24_ORDER:
        a[0], pos_n = extract_int(n, pos_n, 11, 4, verbose)
        a[1], pos_n = extract_int(n, pos_n, 13, 3, verbose)
        a[2], pos_n = extract_int(n, pos_n, 12, 1, verbose)
        a[3] = 0x20000000 + p
        pos_a = 4
    elif p == MAT24_ORDER:
        pos_a = 0
    elif p == MAT24_ORDER + 1:        
        i, pos_n = extract_int(n, pos_n, 1, 5, verbose)
        a[0] = i + 1
        pos_a = 1;
    elif p == MAT24_ORDER + 2:        
        i, pos_n = extract_int(n, pos_n, 17, 0, verbose)
        c = expand_17bit_type2(i)
        if verbose:
            print("extracted type 2:", hex(c))
        s = gen_leech2_reduce_type2(c, a)         
        assert 0 <= s <= 6, ERR_DECODE_TYPE % hex(i)
        mm_group_invert_word(a, s)
        pos_a += s
        with_t = 1
    else:
        raise ValueError(ERR_DECODE_RANGE)

    while True:
        i, pos_n = extract_int(n, pos_n, 23 + with_t, 0, verbose)
        if with_t and i >= 2:
            a[pos_a] = 0x50000001 + (i & 1)
            pos_a += 1
        i >>= with_t
        with_t = 1
        if i < 2:
            return a[:pos_a]        
        c = expand_23bit_type4(i)
        if verbose:
            print("extracted type 4:", hex(c))
        s = gen_leech2_reduce_type4(c, a[pos_a:])
        assert 0 <= s <= 6, ERR_DECODE_TYPE % hex(i)
        mm_group_invert_word(a[pos_a:], s)
        pos_a += s
        assert pos_a + 7 <= LEN_A, ERR_DECODE_RANGE


            
def decompress(n, verbose = 0):
     return MM('a', decompress_to_array(n, verbose))


#####################################################################################
# Test compressing in monster group
#####################################################################################


def compress_testcases(ntests = 10):
    cases = [
        'M<y_0c00h>',
        'M<y_0c8dh*d_394h*p_240807357*t_1>',
    ]
    for s in cases:
        yield MM(s), False
    yield MM(), True
    tags = "yxdptl"
    maxvals = [0x7ff, 0x1fff, 0xffff,  MAT24_ORDER - 1, 2, 2] 
    for tag, maxv in zip(tags,  maxvals):
        yield MM(tag, 1), True
        yield MM(tag, maxv), True
        if maxv > 2:
            yield MM(tag, randint(2, maxv-1)), True
    for i in range(ntests):
        for subgroup in ['N_0', 'G_x0']*3 + ['M']:
            yield MM('r', subgroup), False
        



@pytest.mark.mmgroup 
def test_mm_decompress(ntests = 12, verbose = 0):
    print("Test conversion of Monster elements to integers")
    #expand_23bit_type4(0x45678f)
    #expand_17bit_type2(0x673f)
    for nt, (g, small) in enumerate(compress_testcases(ntests)):
        import_all()
        if verbose:
            print("test %d" % (nt+1))
            print("g=", g)
            gg = g.copy().reduce()
            tmp = GtWord(gg.mmdata).as_int_debug_compress()
            print("debug", [hex(int(x)) for x in tmp])
        n = g.as_int()
        if verbose:
            print("n=", hex(n))
        if verbose or nt < 30:    
            g1 = decompress(n, verbose)
            if verbose:
                print("g1=", g1)
            assert g.reduce() == g1.reduce()
        if small:
            assert n < 1 << 64
        g2 =  MM_from_int(n)     
        assert g == g2

