
import os
import sys
from numbers import Integral
from random import randint, sample, choice

import numpy as np
import pytest


from mmgroup import XLeech2, Cocode, PLoop, MM
from mmgroup.clifford12 import xsp2co1_elem_from_mapping

from mmgroup.bimm import inc_p3
from mmgroup.bimm import P3_incidences
from mmgroup.bimm.inc_p3 import P3_point_set_type
from mmgroup.bimm import AutP3, AutP3_MM
from mmgroup.bimm import Norton_generators

from mmgroup.bimm.p3_to_mm import PointP3, StarP3, compute_StarP3
from mmgroup.bimm.p3_to_mm import Precomputed_AutP3


#####################################################################
# Testing the 'Points' and 'Stars' 
#####################################################################
   
    

def do_test_P_Pstar():
    from mmgroup.bitfunctions import bitparity
    def not_commuting(a, b):
       return  (a**(-1))**b * a  != XLeech2() 
       
    assert StarP3(range(13)) == XLeech2()
    P0_list = [PointP3([0,i]) for i in range(13)]
    Pstar_list = [StarP3(i) for i in range(13)]
    for i, P in enumerate(P0_list):
        if i: assert XLeech2(P).type == 2
        i_bitlist = 1 ^ (1 << i)
        for istar, Pst in enumerate(Pstar_list):
            istar_bitlist = 1 << istar
            ncomm = bitparity(i_bitlist & istar_bitlist)
            assert ncomm == not_commuting(P, Pst), (i, istar, ncomm, P, Pst)

        for P2 in P0_list:
            assert not_commuting(P, P2) == 0

    for i in range(13):
        compute_StarP3(i, check = True)

    for i, Pst in enumerate(Pstar_list):
        if i: assert XLeech2(Pst).type == 4, i
        for Pst2 in Pstar_list:
            assert not_commuting(Pst, Pst2) == 0





def do_test_P_sets(N = 1000, verbose = False):
    types = {}
    star_types = {}
    s ="""Test that the types of the products of the points and of the stars
(as elements of 2**{1+24}) are consistent with the geometry of P3."""
    print(s)
    for i in range(2,13,2):
        for n in range(N):
            s = sample(range(13), i)
            set_type = P3_point_set_type(s)
            type_ = PointP3(s).type
            if set_type in types:
                 assert types[set_type] == type_
            else:
                 if (verbose): print(set_type, type_)
                 types[set_type] = type_ 

    for i in range(1, 13):
        for n in range(N):
            s = sample(range(13), i)
            set_type = P3_point_set_type(s)
            star_type = XLeech2(StarP3(s)).type
            if set_type in star_types:
                 assert star_types[set_type] == star_type
            else:
                 if verbose: print("*", set_type, star_type)
                 star_types[set_type] = star_type

    print("Test passed")
   

#####################################################################
# Test stuff
#####################################################################


def do_test_construction_P3(ntests = 500, verbose = 0):
    print("Test construction of AutP3")
    for i in range(10):
        #print("Test", i+1)
        p = AutP3('r', zip([0,1,3],[0,1,9]))
        #print(p)

    p = AutP3(zip(range(13,25), range(14,26)))
    #print(p)
    #print(p.order())


    orders = np.zeros(14, dtype = np.uint32)
    AutP3_ONE = AutP3()
    for i in range(ntests):
        if verbose:
            print("Test", i+1)
        p = AutP3('r')
        p1 = AutP3(zip(p.point_map(), range(13)))
        assert p*p1 == AutP3_ONE
        orders[p.order()] += 1
        m = AutP3_MM(p)
    print('  Orders of elements of AutP3 [1,...,13], %d samples:' % ntests)
    print('    %s' % orders[1:])
    print('  Good orders:', Precomputed_AutP3.good_orders)
    print('  Difficult orders:', list(Precomputed_AutP3.bad_orders.keys()))
    print('  %d elements of AutP3 have been split with %d trials.' %
           (Precomputed_AutP3.n_splits, Precomputed_AutP3.n_split_trials))
    print('  Number of stored elements of AutP3:', 
            Precomputed_AutP3.num_MM - 1)
    print("Test construction of AutP3 passed")
   



def do_test_random_relations():
    N = 200
    print("Testing %d random relations in the rep of AutP3 in G_x0" % N)
    for i in range(200):
        g1 = AutP3('r')
        g2 = AutP3('r')
        g3 = g1 * g2
        mg1 = AutP3_MM(g1)
        mg2 = AutP3_MM(g2)
        mg3 = AutP3_MM(g3)
        assert mg3 == mg1 * mg2
    print("Test passed")





#####################################################################
# Test of the module
#####################################################################


@pytest.mark.bimm
def test_all():
    do_test_P_Pstar()
    do_test_P_sets()
    do_test_construction_P3(500)
    do_test_random_relations()
    generators = Norton_generators(check = True)
    NAMES = "stuvx"
    print("Generators in Norton's construction of the Monster:")
    for name, gen in zip(NAMES, generators):
        print("%s = %s" % (name, gen))




