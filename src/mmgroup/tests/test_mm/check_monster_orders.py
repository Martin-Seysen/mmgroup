
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import randint
import datetime
import time

from multiprocessing import Pool, TimeoutError
from collections import defaultdict
from scipy.stats import chisquare

from mmgroup.mm_space import MMSpace
from mmgroup.mm import INT_BITS


################################################################
# Check that monster group elements have coorect orders
################################################################

p = 3
space = MMSpace(p)
group = space.group


good_mm_orders = set(range(1,37)) | set(range(38,43)) | set(range(44,49))
good_mm_orders.update([50, 51, 52, 54, 55, 56, 57, 59, 60, 62, 66, 
    68, 69, 70, 71, 78, 84, 87, 88, 92, 93, 94, 95, 104, 105, 110, 119])

max_mmm_order = max(good_mm_orders) 



def one_test_mm_order(v,  m, verbose = 0):
    v1, n = v.copy(), 0
    while n <= max_mmm_order:
        v1, n = v1 * m, n+1
        if v1 == v:
             return n
    return None


def rand_v():
    return  space.rand_vector()

    
def rand_m_monomial():
    m = group.neutral()
    for s in "pxy":
        m *=  group.rand_word(s)
    m *= group.rand_word('t', randint(1,2))
    m *= group.atom('l', randint(1,2))   
    return m

def rand_m(n_entries = 4):
    m = group.neutral()
    for i in range(n_entries):
        m *=  rand_m_monomial()   
    return m

def random_test_order(n_entries = 4, display = True):
    v, m = rand_v(), rand_m()
    order = one_test_mm_order(v,  m, display)
    ok = order in good_mm_orders
    st = "ok" if ok else "error"
    if display: 
        print("\rorder is", order, ",", st)
    s =  "\nm = " + str(m) 
    s += "\norder = " + str(order) + ", " + st + "\n"
    return ok, order, s




def check_mm_orders(ntests, display = True):
    print("\nTesting orders of elements of the monster group")
    nerrors = 0
    order_sum = 0
    start_time = datetime.datetime.now()
    print(start_time)
    t_start = time.process_time()
    for i in range(ntests):
        t = time.process_time()
        if display:
            print("Test %d, CPU time = %.3f s" % (i+1, t) )
        ok, order, _ =  random_test_order(display = display)
        nerrors += not ok
        order_sum += order
    t = time.process_time() - t_start
    print("started: ", start_time)
    print("finished:", datetime.datetime.now())
    print("CPU time = %.3f s, per test: %.3f ms" % (t, 1000*t/ntests))
    print("CPU time per standard operation: %.5f ms" % (1000.0*t/order_sum))
    print("%d tests, %d errors, " % (ntests, nerrors))
    if nerrors:
        raise ValueError("Error in orders of monster group elements")




################################################################
# Chisquare test of orders of monster group elements
################################################################

MM_WORD_SIZE = 20  # No of elementary operations to construct
                   # an element of the monster


MIN_CHISQU = 500   # Min No of cases for chisquare test

class ChisquareOrder:
    self_centralizeing_classes = [ # from the ATLAS
        71, 71, 78, 78, 78, 84, 84, 84, 87, 87, 88, 88, 92, 92, 93, 93,
        94, 94, 95, 95, 104, 104, 105, 110, 119, 119
    ]
    probabilities = defaultdict(float)
    for n in self_centralizeing_classes:
        probabilities[n] += 1.0/n
    probabilities[0] = 1.0 - sum(probabilities.values())
    positions = sorted(probabilities)
    f_exp = [None] * len(positions)
    for i, order in enumerate(positions): 
        f_exp[i] = probabilities[order]
    min_order = min(self_centralizeing_classes)
    chisquare_ = chisquare

    def __init__(self, p = p):
        self.obtained = defaultdict(int)
        self.p = p
        self.total = 0
        self.order_sum = 0
        self.errors = 0 
        self.word_size = MM_WORD_SIZE

    def add(self, ok, order):
        assert order in good_mm_orders
        key = order if order >= self.min_order else 0
        self.obtained[key] += 1
        self.total += 1
        self.order_sum += order
        self.errors += not ok

    def chisquare(self):
        f_obt = [self.obtained[key] for key in self.positions]
        sum_obt = sum(f_obt)
        f_exp = [sum_obt * x for x in self.f_exp]
        chisq, p = chisquare(f_obt, f_exp = f_exp)
        return chisq, p

    def is_ok(self):
        if self.errors:
            return False
        if self.total < MIN_CHISQU:
            return True
        _, prob = self.chisquare()
        return prob > 1.0e-6

    def show_result(self):
        description = (
"""Chisquare test of distribution of orders >= %d in the monster M,
%d degrees of freedom, characteristic p = %d, %d-bit C
random element of MM built from %d factors,
%d tests, %d MM operations, %d errors.
"""     )
        s = description % (
            self.min_order, 
            len(self.f_exp) - 1,
            self.p, INT_BITS,  self.word_size,
            self.total, self.order_sum, self.errors
        )
        if self.errors == 0 and self.total >= MIN_CHISQU:
            st = "\nChisquare test statistics = %.3f, p = %.4f\n"
            chisq, p = self.chisquare()
            s += st % (chisq, p)
        return s
 



def one_test_order(args):
    v, m = args
    order = one_test_mm_order(v,  m)
    ok = order in good_mm_orders
    return ok, order



def get_test_values(ntests):
    for i in range(ntests):
         yield rand_v(), rand_m(MM_WORD_SIZE)



def statistics_chisqu_orders(results, start_time = None):
    if not start_time is None:
        end_time = datetime.datetime.now()
    chisq = ChisquareOrder()  
    for i, (ok, order) in enumerate(results):
        st = "ok" if ok else "error"
        chisq.add(ok, order)

    print("\n" + chisq.show_result())

    if not start_time is None: 
        ntests, order_sum = chisq.total, chisq.order_sum
        diff_time = end_time - start_time
        t = diff_time.total_seconds() 
        print("started: ", start_time)
        print("finished:", end_time)
        print("time = %.3f s, per test: %.3f ms" % (t, 1000*t/ntests))
        print("time per standard operation: %.5f ms" % (1000.0*t/order_sum))

    return chisq.is_ok()
    

def check_chisqu_orders(ntests, nprocesses = 1, verbose = False):
    start_time = datetime.datetime.now()
    header = "\nChisquare test of distribution of orders in the monster M,"
    print(header)
    print("%d tests, %d processes" % (ntests, nprocesses))
    print("started: ", start_time)
    testvalues = get_test_values(ntests)
    if nprocesses > 1:
        with Pool(processes = nprocesses) as pool:
            results = pool.map(one_test_order, testvalues)
        pool.join()
    else:
        results_ = map(one_test_order, testvalues)
        results = []
        for i, x in enumerate(results_):
            ok, order = x
            if verbose:
                print("Test %d, order = %3d, %s" % (i+1, order, ok) )
            else:
                print("\r %d     " % i, end = "")
            results.append(x)               
    return statistics_chisqu_orders(results, start_time) 




