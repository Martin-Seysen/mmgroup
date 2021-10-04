"""Modelling vectors of the representation of the monster

We deal with the rational representation [Seys19] of the monster 
group MM, which is  based on the Conway's construction [Conw85] of 
the monster, modulo various small integers p = 2**n-1, n <= 8.
Here the integer p is called the modulus. See [Conw85], [Seys19], 
for details. We also refert to the documentation of the C modules 
mm_aux.c and mat24_functions.c in this package.
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import sys
import os
from numbers import Integral
from random import randint, sample
import numpy as np
import pytest

from mmgroup.tests.spaces.sparse_mm_space import SparseMmSpace
from mmgroup.tests.spaces.sparse_mm_space import SparseMmV
from mmgroup.tests.groups.mgroup_n import MGroupNWord 

V = SparseMmV
G = MGroupNWord

#################################################################################
# Test functions
#################################################################################


def direct_testdata():
    """Yield test triples (v, g, v*g)

    with v a vector, g a group element and v*g the expected product of v and g. 
    """
    space = V(7)
    g = G("d", 1) * G("y", 1)
    v = -space("Y", 0, 0)
    yield v*g, g**(-1), v
    g = G("t", 2) 
    v = 2*space('C', 0, 1)
    yield v*g, g**2, v

@pytest.mark.space
def test_sparse_direct():
    for v, g, expected in direct_testdata():
        if v * g != expected:
            print("v =", v)
            print("g =", g)
            print("1/g =", g**(-1))
            #assert g * g**(-1) == g**(-1) * g == g**0
            print("v * g expected =", expected)
            print("v * g obtained =", v * g)
            w = v
            txt = "v"
            for g0 in g:
                w *= g0
                print("%s * %s = %s" % (txt, g0, v))
                txt = " "
            raise ValueError("Bad vector multiplication result")
      



def mmv_data(n_tests, p, v_tags, g_tags, n_comp = 1):
    Vp = V(p)
    def rand_group_word(g_tags):
        tags = sample(g_tags, randint(3,5)) 
        return G([('x','r') for x in tags]) 
    def rand_vector(v_tags, n_comp):
        tags = [sample(v_tags, 1)[0] for i in range(n_comp)]
        return Vp([(x,'r') for x in tags])       
    
    for i in range(n_tests):
        g1 = rand_group_word(g_tags)
        g2 = rand_group_word(g_tags)
        yield rand_vector(v_tags, n_comp), g1, g2

def mmv_testdata(p):
    Vp = V(p)
    def mu(tag, *data):
        return Vp(tag, *data)
    def na(tag, *data):
        return G(tag, *data) 
    #n_group =  mmv_space.group
    testdata = [
        
        (mu('T',0,1) , na('t', 1) , na('t', 1) ),
        (mu('B',0,1) , na('t', 1) , na('d', 1177)),
        (mu('B',0,1) , na('t', 1) , na('p', 2221)),
        (3*mu('B',0,1) + mu('A',1,1), na('t', 1)*na('d', 1177) , na('y', 17)),
        (3*mu('B',0,1) + mu('A',1,1), na('t', 1)*na('p', 2221) , na('y', 17)),
        (3*mu('A', 1, 1), na('t', 1), na('t',1)),
        (2*mu('A', 0, 1), na('t', 1), na('t',1)),
        (mu('Y', 0, 0), na("d", 1), na("y", 1)),
        (mu('Y', 0, 0), na("y", 1), na("d", 1)),
        (mu('X', 0x6d3, 23), na("d", 1), na("d",0x331)),
        (mu('X', 0x6d3, 23), na("p", 7851), na("d",0x331)),
        (mu('X', 0x6d3, 23), na("d", 1), na("p", 197652)),
        (mu('X', 0x6d3, 23), na("p", 7851), na("p", 197652)),
    ]
    for d in testdata:
         yield d
    v_tags, g_tags = "ABCTXYZ", "dpxytl"
    for d in mmv_data(100, p,  v_tags, g_tags, 3):
         yield d
    print(".")
    v_tags, g_tags = "ABCTXYZ", "dpxytl"
    for d in mmv_data(3, p,  v_tags, g_tags, 1000):
         yield d
    print(".")



def short_testdata(mmv_space):
    #v_tags, g_tags = "ABCTXYZ", "dpxytl"
    v_tags, g_tags = "ABCTXYZ", "l"    
    v_tags, g_tags = "YZ", "l"    
    yield from mmv_data(5, mmv_space, v_tags, g_tags, n_comp = 1)

def l_testdata(mmv_space):
    v = mmv_space(("Y",0x400,0))
    group =  mmv_space.group
    yield v, group(("l",1)), group()

TEST_PAR = [
   # (3, l_testdata),
   # (7, short_testdata),
    (7, mmv_testdata)
]

@pytest.mark.space
@pytest.mark.parametrize("p, mmv_testdata", TEST_PAR )
def test_sparse_space(p, mmv_testdata):
    mmv = V(p)

    v = mmv('X', 0x6d3, 23)
    w =  mmv('Y', 0x126, 20)
    v += mmv('X', 0x126, 20) / 4
    maxlen = 0
    print( -17 * v / 8 + (w >> 16))
    for i, (v, g1, g2) in enumerate(mmv_testdata(p)):
        al = (v * g1) * g2
        ar = v * (g1 * g2)
        maxlen = max(maxlen, len(al))
        if al != ar:
            print("\nTest computing v * g1 * g2")
            print("v = ", v)
            print("g1 = ", g1)
            print("g2 = ", g2)
            print("v*g1 = ", v*g1)
            print("g1*g2 = ", g1*g2)
            print("(v*g1)*g2 = ", al)
            print("v*(g1*g2) = ", ar)
            print("Diff = ", al - ar)
            raise ValueError("Error in computing v * g1 * g2")
        g3 = (g1 * g2)**(-1)
        v_back = al * g3
        if v != v_back:
            print("\nTest computing v * g1 * g2 * g3 with g1 * g2 * g3 = 1")
            print("v = ", v)
            print("g1 = ", g1)
            print("g2 = ", g2)
            print("g3 = ", g3)
            print("g1*g2*g3 = ", g1*g2*g3)
            print("v*g1*g2 = ", al)
            print("v*g1*g2*g3 = ", v_back)
            print("Diff = ", v_back - v)
            raise ValueError("Error in computing v * g1 * g2 * g3")
            
        if i < 3:
            print("v = ", v)
            print("g1 =", g1)
            print("g2 =", g2)
            print("v*g1*g2 = ", al)
            print("")
            
    print("Max vector length is", maxlen)
    
      