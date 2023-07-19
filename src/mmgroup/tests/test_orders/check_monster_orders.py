
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

import sys
import os
from random import randint
import datetime
import time

from multiprocessing import Pool, TimeoutError
from collections import defaultdict


from mmgroup import MM0, MMV
from mmgroup.mm_space import MMSpace
from mmgroup.mm_op import INT_BITS
from mmgroup.tests.chisquare import chisquare

################################################################
# Class and character for the monster information taken from GAP
################################################################


#The following information has been obtained from the GAP package: 
GAP_INFO = """
gap> t := CharacterTable("M");  #! The character table of the Monster group
CharacterTable( "M" )
gap> ClassNames(t, "ATLAS");  #! Classes of the Monster in ATLAS notatation 
[ "1A", "2A", "2B", "3A", "3B", "3C", "4A", "4B", "4C", "4D", "5A", "5B", "6A", "6B", "6C", "6D", "6E", "6F", "7A",
  "7B", "8A", "8B", "8C", "8D", "8E", "8F", "9A", "9B", "10A", "10B", "10C", "10D", "10E", "11A", "12A", "12B",
  "12C", "12D", "12E", "12F", "12G", "12H", "12I", "12J", "13A", "13B", "14A", "14B", "14C", "15A", "15B", "15C",
  "15D", "16A", "16B", "16C", "17A", "18A", "18B", "18C", "18D", "18E", "19A", "20A", "20B", "20C", "20D", "20E",
  "20F", "21A", "21B", "21C", "21D", "22A", "22B", "23A", "23B", "24A", "24B", "24C", "24D", "24E", "24F", "24G",
  "24H", "24I", "24J", "25A", "26A", "26B", "27A", "27B", "28A", "28B", "28C", "28D", "29A", "30A", "30B", "30C",
  "30D", "30E", "30F", "30G", "31A", "31B", "32A", "32B", "33A", "33B", "34A", "35A", "35B", "36A", "36B", "36C",
  "36D", "38A", "39A", "39B", "39C", "39D", "40A", "40B", "40C", "40D", "41A", "42A", "42B", "42C", "42D", "44A",
  "44B", "45A", "46A", "46B", "46C", "46D", "47A", "47B", "48A", "50A", "51A", "52A", "52B", "54A", "55A", "56A",
  "56B", "56C", "57A", "59A", "59B", "60A", "60B", "60C", "60D", "60E", "60F", "62A", "62B", "66A", "66B", "68A",
  "69A", "69B", "70A", "70B", "71A", "71B", "78A", "78B", "78C", "84A", "84B", "84C", "87A", "87B", "88A", "88B",
  "92A", "92B", "93A", "93B", "94A", "94B", "95A", "95B", "104A", "104B", "105A", "110A", "119A", "119B" ]
gap> Irr(t)[2];  #! Character of degree 196883
Character( CharacterTable( "M" ), [ 196883, 4371, 275, 782, 53, -1, 275, 51, 19, -13, 133, 8, 78, 77, 14, -3, 5, -1,
  50, 1, 35, 11, -1, -5, 3, -1, 26, -1, 21, 5, -4, 20, 0, 16, 14, 5, 6, -1, -2, 5, -3, 13, 1, -1, 11, -2, 10, 2, 9,
  7, -2, 8, -1, 3, -1, 7, 6, -3, 6, 2, -1, 5, 5, 5, 1, 0, -3, 2, 4, 5, -2, -1, 4, 4, 0, 3, 3, 2, 2, -1, -2, -1, -1,
  -1, 1, 3, -1, 3, 3, 2, 2, 2, 2, 2, -2, 1, 2, 2, 3, -1, 2, -1, 2, 0, 2, 2, 1, 1, -2, 1, 2, 0, 1, 2, -1, 0, 1,
  1, 2, -1, 1, 1, -1, 1, 0, 0, 1, 1, 0, -1, 0, 0, 0, 1, -1, -1, 1, 1, 0, 0, 0, 1, 0, -1, 0, 0, 1, 0, -1, -1, -1, 0,
  0, 1, -1, 0, -2, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, -1, -1, -2, -1, -1, -1, 0, 0, -1, -1, -1, -1, 0,
  0, 0, 0, -1, -1, 0, -1, -1, -1 ] )
gap> SizesCentralizers(t);  #! Sizes of the centralizers of the classes
[ 808017424794512875886459904961710757005754368000000000, 8309562962452852382355161088000000,
  139511839126336328171520000, 3765617127571985163878400, 1429615077540249600, 272237831663616000,
  8317584273309696000, 26489012826931200, 48704929136640, 8244323942400, 1365154560000000, 94500000000,
  774741019852800, 2690072985600, 481579499520, 130606940160, 1612431360, 278691840, 28212710400, 84707280,
  792723456, 778567680, 143769600, 23592960, 12582912, 3096576, 56687040, 2834352, 887040000, 18432000, 12000000,
  6048000, 480000, 1045440, 119439360, 22394880, 17418240, 1161216, 884736, 483840, 373248, 276480, 82944, 23040,
  73008, 52728, 1128960, 150528, 35280, 2721600, 145800, 10800, 9000, 12288, 8192, 8192, 2856, 34992, 23328, 15552,
  3888, 3888, 1140, 76800, 28800, 24000, 19200, 1200, 960, 52920, 6174, 3528, 504, 2640, 2112, 552, 552, 6912, 4608,
  3456, 2304, 1152, 864, 864, 576, 384, 288, 250, 624, 312, 486, 243, 4704, 2688, 896, 168, 87, 10800, 7200, 2880,
  1800, 360, 240, 240, 186, 186, 128, 128, 594, 396, 136, 2100, 70, 1296, 648, 216, 72, 76, 702, 117, 78, 78, 400,
  320, 80, 80, 41, 504, 504, 168, 126, 352, 352, 135, 184, 184, 92, 92, 94, 94, 96, 50, 51, 104, 52, 54, 110, 112,
  56, 56, 57, 59, 59, 360, 240, 120, 120, 60, 60, 62, 62, 132, 66, 68, 69, 69, 140, 70, 71, 71, 78, 78, 78, 84, 84,
  84, 87, 87, 88, 88, 92, 92, 93, 93, 94, 94, 95, 95, 104, 104, 105, 110, 119, 119 ]
"""

def find_table(name):
    """Return table in GAP_INFO after the comment starting with 'name'"""
    s = GAP_INFO[GAP_INFO.find("#! " + name):]
    copen, cclose = s.find("["), s.find("]")
    return eval(s[copen:cclose+1])

ClassNames = find_table("Classes")
ClassOrders = [int(s[:-1]) for s in ClassNames]
CharacterValues = find_table("Character")
SizesCentralizers = find_table("Sizes of the centralizers")

assert len(ClassNames) == len(CharacterValues) == len(SizesCentralizers)

################################################################
# Check that monster group elements have coorect orders
################################################################

p = 3
space = MMV(3)
group = MM0



good_mm_orders = set(ClassOrders)

max_mm_order = max(good_mm_orders) 



def one_test_mm_order(v,  m, verbose = 0):
    v = v.copy()
    v1, n = v.copy(), 0
    while n <= max_mm_order:
        v1, n = v1 * m, n+1
        if v1 == v:
             return n
    return None


def rand_v():
    return  space('R')

    

def rand_m(n_entries = 4):
    return group('r', n_entries)

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
        if ok:
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


MIN_CHISQU = 560   # Min No of cases for chisquare test

class ChisquareOrder:
    probabilities = defaultdict(float)
    orders = set(ClassOrders)
    good_orders = set()
    for order, csize in zip(ClassOrders, SizesCentralizers):
        probabilities[order] += 1.0/csize
        if probabilities[order] >= 1.0/111:
            good_orders.add(order)
    max_small = max(orders - good_orders)
    for x in orders:
        if x <= max_small:
             del probabilities[x]
    min_order = min(probabilities)
    probabilities[0] = 1.0 - sum(probabilities.values())
    chisquare_ = chisquare

    def __init__(self, p = p):
        self.obtained = defaultdict(int)
        self.p = p
        self.total = 0
        self.order_sum = 0
        self.errors = 0 
        self.word_size = MM_WORD_SIZE

    def add(self, ok, order):
        ok = ok and order in self.orders
        if ok:
             key = order if order >= self.min_order else 0
             self.obtained[key] += 1
             self.total += 1
             self.order_sum += order
        self.errors += not ok

    def chisquare(self):
        f_obt = [self.obtained[key] for key in self.probabilities]
        sum_obt = sum(f_obt)
        f_exp = [sum_obt * self.probabilities[key] 
            for key in self.probabilities]
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
            len(self.probabilities) - 1,
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
    verbose = 1
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




