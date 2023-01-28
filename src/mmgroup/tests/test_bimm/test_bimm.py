r"""This module implements tests for the Bimonster.

Class ``BiMM`` in this module implements an element of the Bimonster
:math:`\mathbb{M} \wr 2` as described in the documentation of this
application. Let :math:`\mbox{IncP3}` be the Coxeter group as in
that documentation. Function ``P3_BiMM`` maps a word of generators
of :math:`\mbox{IncP3}` into the Bimonster. The generators of
:math:`\mbox{IncP3}` correspond to the points and lines  of the
projective plane :math:`\mbox{P3}` over the field :math:`\mathbb{F}_3`.
A point or a line  of :math:`\mbox{P3}` is implemented as an
instance of class ``P3_node`` in module ``inc_p3``.

There is also a natural mapping from the automorphism group of
:math:`\mbox{P3}` into the Bimonster compatible with the mapping
from the Coxeter group into the Bimonster. Function
``AutP3_BiMM`` computes that mapping. An automorphism  of
:math:`\mbox{P3}` is implemented as an instance of class ``AutP3``
in module ``inc_p3``. For background see  :cite:`Nor02`.

"""
import os
import sys
from numbers import Integral
from random import randint, sample, choice

import numpy as np
import pytest


from mmgroup import  MM


from mmgroup.bimm import P3_node, P3_is_collinear
from mmgroup.bimm import P3_incidence as inc
from mmgroup.bimm import  AutP3, AutP3_MM, BiMM, P3_BiMM, AutP3_BiMM
from mmgroup.bimm.p3_to_mm import PointP3



#####################################################################
# Tests
#####################################################################



def P3_BiMM_slow(pl):
    result = BiMM()
    for x in pl:
        result *= P3_BiMM(x)
    return result


def do_test_P3_BIMM_cases(max_length):
    LL = list(range(26)) * 5
    for i in range(max_length):
        pl = sample(LL, i)
        assert P3_BiMM(pl) == P3_BiMM_slow(pl)

def do_test_P3_BiMM(multiprocessing = False):
    print("Testing function P3_BiMM")
    NTESTS = 5
    MAX_LENGTH = 10
    if multiprocessing:
        from multiprocessing import Pool
        with Pool() as pool:
            pool.map(do_test_P3_BIMM_cases, [MAX_LENGTH] * NTESTS)
        pool.join()
    else:
        for x in range(NTESTS):
            do_test_P3_BIMM_cases(MAX_LENGTH)
    print("passed")
  

def BiMM_Coxeter_Exp(x1, x2):
    mi, ma = min(x1,x2), max(x1, x2)
    assert 0 <= mi <= ma < 26
    if mi < 13 and ma >= 13 and (mi + ma) % 13 in (0,1,3,9):
         return 3
    return 2 if mi != ma else 1
          
    
def do_test_Coxeter_orders_for_one_node(x):
    for y in range(26):
        o = (P3_BiMM(x) * P3_BiMM(y)).order()
        e = BiMM_Coxeter_Exp(x, y)         
        assert e == o, (x, y, e, o)    

def do_test_Coxeter_orders(multiprocessing = False):
    print("Testing  orders of Coxeter generators of BiMM")
    if multiprocessing:
        from multiprocessing import Pool
        with Pool() as pool:
            pool.map(do_test_Coxeter_orders_for_one_node, range(26))
        pool.join()
    else:
        for x in range(26):
            do_test_Coxeter_orders_for_one_node(x)
    print("passed")

     
def do_test_AutP3_BiMM(ntests=10, verbose = 0):
    print("Test embedding of AutP3 into the Bimonster")
    for i in range(ntests):
        g_autp3 = AutP3('r')
        g = AutP3_BiMM(g_autp3)
        for j in range(4):
            node_p3 =  P3_node(randint(0,25))
            node = P3_BiMM(node_p3)
            img_p3 = node_p3 * g_autp3
            img = P3_BiMM(img_p3)
            assert img == node ** g
    print("passed")



def random_hexagon():
    while 1:
        points = sample(range(13), 3)
        if not P3_is_collinear(points):
            p1, p2, p3 = points
            return p1, inc(p1, p2), p2, inc(p2, p3), p3, inc(p3, p1)


def check_hexagon_relation(u, v, w, x, y, z):
    r"""Check relations in a hexagon in the Bimonster.

    Here ``(u, v, w, x, y, z)`` must be a hexagon in the 
    incidence graph of P3. We take the relation for that
    hexagon from [Iva99], Theorem 8.2.2.
    """
    a = P3_BiMM([u, x, v, y, w, z])
    assert a.order() == 4

def do_test_hexagon_relations(ntests = 10):
    print("Testing some (random) hexagon relations in P3")
    for i in range(ntests):
        check_hexagon_relation(*random_hexagon())
    print("passed")


def do_test_spider_relation():
    print("Testing the spider relation in Y_555")
    spider = P3_BiMM('a,b1,c1,a,b2,c2,a,b3,c3')
    assert spider.order() == 10, spider.order()
    other_spider = ['a', 'b1', 'c1', 'a', 'b2', 'c2', 'a', 'b3', 'c3']
    assert P3_BiMM(other_spider * 10) == BiMM(1)
    print("passed")
    


@pytest.mark.bimm
def test_all():
    do_test_P3_BiMM(multiprocessing = True)
    do_test_AutP3_BiMM()
    do_test_hexagon_relations() 
    do_test_spider_relation()


@pytest.mark.bimm
@pytest.mark.slow
def test_slow():
    do_test_Coxeter_orders(multiprocessing = True)


if __name__ == "__main__":
    test_all()


  