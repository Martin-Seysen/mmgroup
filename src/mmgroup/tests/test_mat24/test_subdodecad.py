from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import pytest

import types
import sys
import re
import os
from operator import __or__, __xor__, __and__
from functools import reduce
import random


 
from mmgroup.bitfunctions import v2, lmap, lrange

from mmgroup.generate_c  import prepend_blanks
from mmgroup.dev.mat24.mat24aux import generate_golay24_decode


def lsbit24(x):
    return v2(x | 0x1000000)


from mmgroup import mat24
from mmgroup.dev.mat24.mat24_ref import Mat24


#########################################################################
# Test completion of heptads and octads
#########################################################################


def random_dodecad(gc):
    """Return random dodecad in gcode rpresentation

    """
    for i in range(1000):
        v = random.randint(0,0xfff)
        if gc.gcode_weight(v) == 3:
            return v
    raise ValueError("Could not find a docecad")


def random_subdodecad(gc, odd=False):
    """Return test sample for function  mat24_cocode_as_subdodecad

    The function returns a test sample ``(d, c, i)`` for function
    ``mat24_cocode_as_subdodecad``. Here ``d`` is a dodecad (in gcode
    representation, ``c`` is a cocode element (in cocode 
    representation), and ``i`` is 24  (if ``odd`` is False) of
    an integer ``0 <= i < 24`` representing a bit poition ``otherwise``.

    Then calling ``mat24_cocode_as_subdodecad(d, c, i)`` should succeed.
    """
    d = random_dodecad(gc)
    if odd:
        v = gc.gcode_to_vect(d)
        for j in range(1000):
             i = random.randint(0,23)
             if v & (1 << i) == 0:
                  return d, random.randint(0, 0xfff), i
    else:
        for j in range(1000):
             c = random.randint(0,0xfff)
             if gc.scalar_prod(d ^ 0x800, c) == 0:
                 return d, c, 24
    raise ValueError("Could not find a docecad test sample")


   


def get_testdata(gc):
    for i in range(400):
        for odd in (0, 1):
            yield random_subdodecad(gc, odd)

              
        
  



@pytest.mark.mat24
def test_subdodecad(gc = mat24, verbose = 0):
    for n, (d, c, i)  in enumerate(get_testdata(gc)):
        ok = True
        try:
            res = gc.cocode_as_subdodecad(c, d, i)
        except ValueError:
            err = "Function cocode_as_subdodecad has rturned an error"
            ok = False
        if ok:
            c_res = gc.vect_to_cocode(res) 
            if c_res != c:
               err = ("Cocode Error: expected: 0x%03x, obtained:  0x%03x" %
                   (c, c_res))
               ok = False
        if ok:
            cc = res if i >= 24 else res & ~(1 << i)
            d_bits =  gc.gcode_to_vect(d)
            if cc & d_bits != cc:
               err = "Error: Result is not a subset of the dodecad" 
               ok = False
        if verbose or not ok:
            print("Test", n)
            d_bits =  gc.gcode_to_vect(d)
            print("Dodecad = 0x%03x, bits = 0x%06x" % (d, d_bits)) 
            c_bits = gc.cocode_syndrome(c, gc.lsbit24(d_bits)) 
            print("Cocode  = 0x%03x, bits = 0x%06x, weight = %d" 
                % (c, c_bits, gc.bw24(c_bits))) 
            print("Bit pos. =", i)
            if ok:
                print("Result = 0x%06x" % res) 
            if not ok:
                raise ValueError(err)




             
