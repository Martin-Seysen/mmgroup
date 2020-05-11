
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import sys
import os
import collections
import re
import warnings
from random import randint
import threading
import datetime
import pytest

from mmgroup.tests.spaces.sparse_mm_space import SparseMmSpace


NJOBS = 5
GROUP_N_ONLY = 0

good_mm_orders = set(range(1,37)) | set(range(38,43)) | set(range(44,49))
good_mm_orders.update([50, 51, 52, 54, 55, 56, 57, 59, 60, 62, 66, 
    68, 69, 70, 71, 78, 84, 87, 88, 92, 93, 94, 95, 104, 105, 110, 119])

max_mmm_order = max(good_mm_orders) 


p = 241
rep_mm = SparseMmSpace(p)
grp = rep_mm.group


def do_test_mm_order(v,  m, verbose = 0):
    v1, n = v.copy(), 0
    while n <= max_mmm_order:
        v1, n = v1 * m, n+1
        if verbose: print("\r", n, end = "")
        if v1 == v:
             return n
    return None


v_tags = "ABCTTXXXYYYZZZ"

def rand_v():
    v = rep_mm.zero()
    for s in v_tags:
        v +=  rep_mm.unit(s)
    return v
    
    
def rand_m_monomial():
    m = grp.neutral()
    for s in "pxy":
        m *=  grp.rand_word(s)
    if not GROUP_N_ONLY :
        m *= grp.rand_word('t', randint(1,2))
    m *= grp.atom('l', randint(1,2))   
    return m

def rand_m():
    m = grp.neutral()
    for i in range(4):
        m *=  rand_m_monomial()   
    return m

def do_test_order(display = True):
    v, m = rand_v(), rand_m()
    if display: 
        print("\rv =", v)
        print("\rm =", m)
    order = do_test_mm_order(v,  m, display)
    ok = order in good_mm_orders
    st = "ok" if ok else "error"
    if display: 
        print("\rorder is", order, ",", st)
    s =  "v = " + str(v) + "\nm = " + str(m) 
    s += "\norder = " + str(order) + ", " + st + "\n"
    return ok, s



@pytest.mark.space
@pytest.mark.slow
@pytest.mark.very_slow
def test_orders(ntests = NJOBS):
    nerrors = 0
    for i in range(ntests):
        print("Test %d," % (i+1), datetime.datetime.now())
        ok, _ =  do_test_order()
        nerrors += not ok
    print("%d tests, %d errors, " % (ntests, nerrors),
            datetime.datetime.now())
    if nerrors:
        raise ValueError("Illegal order of monster group element found")






