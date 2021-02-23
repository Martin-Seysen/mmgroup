from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint

import pytest


from mmgroup.tests.test_mm.test_group_op import test_op, test_rand_op
from mmgroup.tests.test_mm.test_group_op import one_test_rand_op
from mmgroup.mm_space import characteristics
from mmgroup.tests.spaces.spaces import MMTestSpace


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
    space = v1.space
    if tag == 'd':
        space.mm.op_delta(data1, i, data2) 
    elif tag == 'p':
        space.mm.op_pi(data1, 0, i, data2)
    elif tag == 't':
        space.mm.op_t(data1,  i, data2)
    elif tag == 'l':
        space.mm.op_xi(data1,  i, data2)
    elif tag == 'x':
        space.mm.op_xy(data1, 0, i, 0, data2)
    elif tag == 'y':
        space.mm.op_xy(data1, i, 0, 0, data2)
    else:
        raise TypeError("Bad tag %s in monster operation" % str(t))
 

def mul_group(v, g):
    v = v.copy()
    v1, v2 = v, v.space.zero()
    for tag, i in g.as_tuples():
        mul_atom(v1, tag, i, v2)
        v1, v2 = v2, v1
    if v != v1:
       v.data = v1.data[:]
    v.check()
    return v
   

@pytest.mark.mm
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
    assert t[0][0] == "d" and t[1][0] == "p" 
    v1 = v.space.zero()
    v.space.mm.op_pi(v.data, t[0][1], t[1][1], v1.data)
    return v1





@pytest.mark.mm
def test_op_delta_pi(verbose = 0):
    print("Testing group operation delta * pi")
    for p in PRIMES:
        space = MMTestSpace(p)
        group = space.group
        for i in range(50):
            g = group.word( ("d",), ("p",) )
            basis = "D" * 3 + "ABC" * 10 + "TXYZ" * 20 
            one_test_rand_op(space, g, basis, f_mul_delta_pi, verbose)
    print("passed")
            

########################################################################
# Low level test for operation y * x * delta
########################################################################

"Here we test the function op_xy() implementig a group operation."

def f_mul_yx(v, g):
    t = g.as_tuples()
    assert t[0][0] == "y" and t[1][0] == "x" and t[2][0] == "d" 
    y, x, d = t[0][1], t[1][1], t[2][1]  
    v1 = v.space.zero()
    v.space.mm.op_xy(v.data, y, x, d, v1.data)
    return v1


@pytest.mark.mm
def test_op_yx(verbose = 0):
    print("Testing group operation yx")
    for p in PRIMES:
        space = MMTestSpace(p)
        group = space.group
        for i in range(200):
            g = group.word( ("y",), ("x",), ("d",) )
            basis = "D" * 3 + "ABC" * 10 + "TXYZ" * 20 
            one_test_rand_op(space, g, basis, f_mul_yx, verbose)
    print("passed")



def f_mul_omega(v, g):
    t = g.as_tuples()
    assert t[0][0] == "x" and t[0][1] & ~0x1800 == 0
    v1 = v.copy()
    v.space.mm.op_omega(v1.data, t[0][1])
    return v1

@pytest.mark.mm
def test_op_omega(verbose = 0):
    print("Testing group operation x_Omega")
    for p in PRIMES:
        space = MMTestSpace(p)
        group = space.group
        for i in range(5):
            for d in (0, 0x800, 0x1000, 0x1800):
                g = group.word( ("x", d) )
                basis = "D" * 3 + "ABC" * 10 + "TXYZ" * 20 
                one_test_rand_op(space, g, basis, f_mul_omega, verbose)
    print("passed")
