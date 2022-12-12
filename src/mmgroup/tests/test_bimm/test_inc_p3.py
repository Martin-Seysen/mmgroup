from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from math import log
from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.bimm import P3_node, P3_incidences, P3_incidence
from mmgroup.bimm import P3_remaining_nodes, AutP3
from mmgroup.bimm.inc_p3 import p3_list


#####################################################################
# Tests
#####################################################################


def show_Y555():
    def name(v):
        return P3_node(v).name()
    print("Nodes in the Y_555 graph")
    print("a:", name("a"))
    A = lambda i,j : "bcdef"[i] + str(j)
    for j in range(1,4):
        d = [A(i, j) + ": " + name(A(i, j)) for i in range(5)]
        print(", ".join(d))
    print("More nodes in the plane with names as in the ATLAS")
    print("f:", name("f"))
    F = lambda i,j : "agz"[i] + str(j)
    for j in range(1,4):
        d = [F(i, j) + ": " + name(F(i, j)) for i in range(3)]
        print(", ".join(d))


def generate_inc_P3_test_cases_with_Y555_names(verbose = 0):
    """Yield expected incidence relations from a table in the ATLAS

    The function returns pairs ``(p, L)``, where ``p`` is the number
    of a node in ``P3``, and ``L`` is the list of the numbers of
    the nodes incident with node ``p``.

    The ATLAS :cite:`Atlas` defines names for all 26 nodes of the
    projective plane ``P3``; and it states the incidence relations
    between these nodes in a table in the description of the Monster.
    The list ``REL`` define inside this function is almost a
    a direct copy of that table in ATLAS.

    This generator function yields pairs ``(p,L)`` of incidence
    relations as described above. Therefore it uses the table ``REL``
    and the mapping from the ``Y_555`` notation of nodes (used in
    the ATLAS) to the numbers of the nodes used in this class.
    This mapping is described in the documentation of this
    application and and implemented in class ``P3_node``
    """
    REL = [ # This is a essentially copy of a table in the ATLAS
      'a,bi,bj,bk,f',   'zi,ai,cj,ck,ei',
      'ai,zi,bi,fj,fk', 'bi,a,ai,ci,gi',
      'ci,zj,zk,bi,di', 'di,ci,ei,gj,gk',
      'ei,zi,di,fi,f',  'fi,aj,ak,ei,gi',
      'gi,bi,dj,dk,fi', 'f,a,ei,ej,ek',
      # 'ai,bi,bi,fj,fk', # This is an entry that should fail
    ]
    PERMS = [ # list of lists of permutations of n+1 numbers 1,2,...
        [[1]],
        [[1,2],[2,1]],
        [[1,2,3], [1,3,2], [2,1,3], [2,3,1], [3,1,2], [3,2,1]]
    ]
    for r in REL:
        for i, a in enumerate("ijk"):
            r = r.replace(a, '{%d}' % i)
        n = max(int(x) for x in r if x.isdigit())
        for pi in PERMS[i]:
            data = p3_list(r.format(*pi))
            if verbose:
                 print(r.format(*pi), 'maps to', data)
            yield data[0], data[1:]


def do_test_inc_P3_Y555_names(verbose = 0):
    """Check the incidence of ``P3`` stated in the ATLAS.

    Function ``generate_inc_P3_test_cases_with_Y555_names`` generates
    incidence relations in ``P3`` obtained from a table in the
    ATLAS. Function ``test_inc_P3_Y555_names`` simply tests these
    relations.
    """
    print("Test incidence relations in P3 with names taken from ATLAS")
    for p, inc in generate_inc_P3_test_cases_with_Y555_names(verbose):
        inc_p = set([x._ord for x in P3_incidences(p)])
        assert inc_p == set(inc), (p, inc, P3_incidences(p))
    print("passed")


@pytest.mark.bimm
def test_all_inc_p3():
    print("Test classes P3_node and AutP3")
    show_Y555()   
    print("The graph Y_555 has a vertex e.g.", P3_node("f2"))     
    a = AutP3("r", "b1:b2 , d1:d2")
    assert P3_node("c1") * a  ==  P3_node("c2")
    c = AutP3("c1:c2, e1:e2, c2:c3, e2:e3, c3:c1, e3:e1 ")
    assert c.order() == 3, c.order()
    b = AutP3("b1:b2, d1:d2, b2:b3, d2:d3")
    #print(b, P3_node("b3")*b)
    assert b.order() == 3, b.order()
    do_test_inc_P3_Y555_names()
    print("Test passed")






         


