from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import numpy as np
from random import randint

import pytest

from mmgroup.mm_op import mm_sub_test_prep_pi_64
from mmgroup import mat24





def get_prep_pi64(delta, pi):
    """Obtain component tbl_perm64 from function mm_sub_prep_pi()

    More specifically, function mm_sub_prep_pi() in module 
    mm_tables.c takes a pair (delta pi) indicating a automorphism
    of the Parker loop. That function computes a structure
    containing a component 'tbl_perm64'. This function performs
    the corresponding computation for the input (delta, pi) and
    return component 'tbl_perm64' of the computed structure
    as a numpy array of length 2*759 with dtype = uint32.

    For the interpretation of the returned data, see documentation 
    of function mm_tables.c
    """
    a = np.zeros(759*7, dtype = np.uint32)
    mm_sub_test_prep_pi_64(delta, pi, a)
    return a.reshape(759,7)

def parse_a(a):
    """Interpret data returned by function get_prep_pi64()

    These data describe an automorpism of the Parler loop.

    The function returns a list of tuples (source, sign, fields),
    one for each octde 0,...,759-1. 'source' is the number of
    the octad which is the preimage of octad i. 'sign' is the
    sign of the image of octad 'source'. 'fields' is a 6-tuple
    of 6-bit integers describing the mapping of the suboctads,
    see documentaiton in file mm_tables.c.
    """
    a1 = [None] * 759
    for i, row in enumerate(a):
        idx = row[0] & 0x3ff
        sgn = (row[0] >> 12) & 1
        a1[i] = (idx, sgn, [int(i) for i in row[1:]])
    return a1


def dump_pi64(delta, pi, a, format = 0):
    print("raw data for delta = 0x%x, pi = %d:" % (delta, pi))
    for i in range(759):
         print("%3d: %s" % (i, a[i]))



def map_octad(oct, delta, pi):
    """Map octad through a Parker loop automorphsim to octad.

    Here 'oct' is an octad given in ocatad representation,
    interpreted as a (positive) element of the Parker loop.
    The pair (delta, pi), with delte a cocode element and pi
    the number of a permutation in Mat24, is an automorphism
    of the Parker loop. 

    The function returns a pair (sign, img_octad) representing
    an element of the Parker loop with a given sign. img_octad 
    is the number of the correspondig octad in octad 
    representation.

    The return value is the image of the octad under the
    Parker loop automorphism (delta, pi).
    """
    gcode = mat24.octad_to_gcode(oct)
    perm = mat24.m24num_to_perm(pi)
    m = mat24.perm_to_autpl(delta, perm)
    image = mat24.op_ploop_autpl(gcode, m)
    image_oct = mat24.gcode_to_octad(image)
    sign = (image >> 12) & 1
    return sign, image_oct
 
def suboctad_to_cocode(suboctad, oct):
    """Compute code vector from octad and suboctad

    The octad 'oct' must be given in 'octad' representation and
    the suboctad must be given as in function 
    mat24_suboctad_to_cocode in module mat24_functions.c.
    The function returns the suboctad as a cocode element in 
    'cocode' representation.    
    """
    return  mat24.suboctad_to_cocode(suboctad, oct)
   

def map_cocode(c, pi):
    """Map cocode element with permutation pi

    'c' is a cocode element in 'cocode' representation.
    'pi' is the number of a permutation in the Mathieu group
    Mat24. The function returns the image of 'c' under the 
    permutation 'pi' in 'cocode' representation.
    """
    perm = mat24.m24num_to_perm(pi)
    return mat24.op_cocode_perm(c, perm)


def one_test_prep_pi(delta, pi, verbose = 0):
    a = get_prep_pi64(delta, pi)
    a1 = parse_a(a)
    if verbose:
        print("testing delta = %x, pi = %d" % (delta, pi))
        if verbose > 2:
            dump_pi64(delta, pi, a)
        if verbose > 1:
            for i, l in enumerate(a1): print(i, l)
    for i in range(100):
        oct = randint(0,758)
        preimage = a1[oct]
        pre_oct, pre_sign, pre_suboctads = preimage
        sign, image_oct =  map_octad(pre_oct, delta, pi)
        assert image_oct == oct, (image_oct, oct)
        for j, pre_suboctad in enumerate(pre_suboctads):
            suboctad = 2*2**j - 1
            pre_cocode = suboctad_to_cocode(pre_suboctad, pre_oct)
            cocode = suboctad_to_cocode(suboctad, oct)
            img_cocode = map_cocode(pre_cocode, pi)
            #print(j, pre_suboctad, suboctad, pre_cocode, cocode, pre_oct, oct)
            assert cocode == img_cocode, (i, j, cocode, img_cocode)
        assert sign == pre_sign, (sign, pre_sign)

def pi_testcases():
    testdata = [
        (0, 0),
    ] 
    for x in testdata: yield x  
    for i in range(80):
        yield randint(0, 0xfff), randint(0, mat24.MAT24_ORDER - 1)


def even_testcases():
    testdata = [
        (0, 0),
    ] 
    for x in testdata: yield x  
    for i in range(80):
        yield randint(0, 0xfff) & 0x7ff, randint(0, mat24.MAT24_ORDER - 1)


@pytest.mark.mm_op
@pytest.mark.parametrize("testcases", [pi_testcases])
def test_prep_pi(testcases, verbose = 0):
    e = "\n" if verbose else ""
    print("Test function mm_sub_prep_pi... ", end = e)
    for delta, pi in testcases():
        one_test_prep_pi(delta, pi, verbose = verbose) 
    print("passed")



    