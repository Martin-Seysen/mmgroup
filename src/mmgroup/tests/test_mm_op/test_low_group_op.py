from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint
from collections import defaultdict

import pytest


from mmgroup.tests.test_mm_op.test_group_op import test_op, test_rand_op
from mmgroup.tests.test_mm_op.test_group_op import one_test_rand_op
from mmgroup.mm_space import characteristics
from mmgroup.tests.spaces.spaces import MMTestSpace
from mmgroup.structures.mm0_group import MM0


PRIMES = characteristics()

########################################################################
# Low level multiplication
########################################################################


""" Low level tests

The following function mul_group() is used as a low-level group
operation function in the test functions test_op() and test_rand_op().
This allows testing the main operations of the monster group on a 
lower level than the test in module test_group_op.
"""



def mul_atom(v1, tag, i, v2):
    data1 = v1.data
    data2 = v2.data
    mm = v1.ops
    if tag == 'd':
        mm.op_delta(data1, i, data2) 
    elif tag == 'p':
        mm.op_pi(data1, 0, i, data2)
    elif tag == 't':
        mm.op_t(data1,  i, data2)
    elif tag == 'l':
        mm.op_xi(data1,  i, data2)
    elif tag == 'x':
        mm.op_xy(data1, 0, i, 0, data2)
    elif tag == 'y':
        mm.op_xy(data1, i, 0, 0, data2)
    else:
        raise TypeError("Bad tag %s in monster operation" % str(t))
 

def mul_group(v, g):
    v = v.copy()
    v1, v2 = v, v * 0
    for tag, i in g.as_tuples():
        mul_atom(v1, tag, i, v2)
        v1, v2 = v2, v1
    if v != v1:
       v.data = v1.data[:]
    v.check()
    return v
   

@pytest.mark.mm_op
def test_low_group_op():
    print("")
    test_op(f_mul = mul_group, verbose = 0)
    test_rand_op(f_mul = mul_group, verbose = 0)



########################################################################
# Low level test for operation delta * pi
########################################################################

"Here we test the function op_pi() implementig a group operation."

def f_mul_delta_pi(v, g):
    t = g.as_tuples()
    d = defaultdict(int)
    d.update(t)
    v1  = v.copy()  # empty vector of same shape as v
    v.ops.op_pi(v.data, d['d'], d['p'], v1.data)
    return v1





@pytest.mark.mm_op
def test_op_delta_pi(verbose = 0):
    print("Testing group operation delta * pi")
    for p in PRIMES:
        space = MMTestSpace(p)
        group = space.group
        for i in range(50):
            g = MM0( [("d",'r'), ("p",'r')] )
            basis = "D" * 3 + "ABC" * 10 + "TXYZ" * 20 
            v = space('R')
            one_test_rand_op(v, g, basis, f_mul_delta_pi, verbose)
    print("passed")
            

########################################################################
# Low level test for operation y * x * delta
########################################################################

"Here we test the function op_xy() implementig a group operation."

def f_mul_yx(v, g):
    t = g.as_tuples()
    #print("ttt", dict(t))
    d = defaultdict(int)
    d.update(t)
    v1  = v.copy()  # empty vector of same shape as v
    v.ops.op_xy(v.data, d['y'], d['x'], d['d'], v1.data)
    return v1


@pytest.mark.mm_op
def test_op_yx(verbose = 0):
    print("Testing group operation yx")
    for p in PRIMES:
        space = MMTestSpace(p)
        group = space.group
        for i in range(50):
            g = group( [("y","r"), ("x","r"), ("d","r")] )
            basis = "D" * 3 + "ABC" * 10 + "TXYZ" * 20 
            v = space('R')
            one_test_rand_op(v, g, basis, f_mul_yx, verbose)
    print("passed")



def f_mul_omega(v, g):
    t = g.as_tuples()
    if len(t):
        assert t[0][0] == "x" and t[0][1] & ~0x1800 == 0
        x = t[0][1]
    else:
        x = 0
    v1 = v.copy()
    v.ops.op_omega(v1.data, x)
    return v1

@pytest.mark.mm_op
def test_op_omega(verbose = 0):
    print("Testing group operation x_Omega")
    for p in PRIMES:
        space = MMTestSpace(p)
        group = space.group
        for i in range(5):
            for d in (0, 0x800, 0x1000, 0x1800):
                g = group([("x", d)])
                basis = "D" * 3 + "ABC" * 10 + "TXYZ" * 20 
                v = space('R')
                one_test_rand_op(v, g, basis, f_mul_omega, verbose)
    print("passed")
