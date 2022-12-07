import os
import sys
from numbers import Integral
from random import sample, choice

import numpy as np

if not r"." in sys.path:
    sys.path.append(r".")
import inc_p3

from mmgroup import  MM
from mmgroup.structures.abstract_group import AbstractGroup
from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.abstract_group import singleton

from inc_p3 import incidence, find_collinear_points, p3_list
from p3_to_mm import PointP3, StarP3
from p3_to_mm import AutP3
from p3_to_mm import Norton_generators

#####################################################################
# class BiMM
#####################################################################


def gcd(a, b):
    """Return greatest common divisor of ``a`` and ``b``"""
    assert a > 0 and b > 0
    while b > 0: a, b = b, a % b
    return a



class BiMM(AbstractGroupWord):
    """This class models an element of the BiMonster.

    Yet to be documented.
    """
    __slots__ = "m1", "m2", "alpha"
    group_name = "BiMM"
    group = None       # will be set to StdBiMMGroup later

    def __init__(self, m1 = None, m2 = None, alpha = 0):
        if isinstance(m1, BiMM):
            self.m1 = m1.m1
            self.m2 = m1.m2
            self.alpha = m1.alpha
        else:
            self.m1 = MM(m1)
            self.m2 = MM(m2)
            self.alpha = alpha & 1
            
    def orders(self):
        """Return orders(!) of element of the group BiMM"""
        if self.alpha:
             a, par = self * self, 2
        else:
             a, par = self, 1
        return a.m1.order(), a.m2.order(), par

    def order(self):
        o1, o2, s = self.orders()
        return s * o1 * o2 // gcd(o1, o2)

 
 
@singleton
class BiMMGroup(AbstractGroup):
    word_type = BiMM   
    conversions = {}

    def __init__(self):
        super(BiMMGroup, self).__init__()


    def __call__(*args, **kwds):
        raise TypeError("Class AutP3Group object is not callable")

    def atom(self, tag = None, data = None):
        err = "Class AutPlGroup has no attribute 'atom'"
        raise AttributeError(err)

    @staticmethod
    def _imul(g1, g2):
        if g1.alpha:
            m1, m2 = g1.m1 * g2.m2, g1.m2 * g2.m1
        else:
            m1, m2 = g1.m1 * g2.m1, g1.m2 * g2.m2
        return BiMM(m1, m2, g1.alpha ^ g2.alpha) 

    @staticmethod
    def _invert(g1):
        if g1.alpha:
            m1, m2 = g1.m2, g1.m1
        else:
            m1, m2 = g1.m1, g1.m2
        return BiMM(m1**(-1), m2**(-1), g1.alpha )

    def copy_word(self, g1):
        return BiMM(g1.m1, g1.m2, g1.alpha)

    def _equal_words(self, g1, g2):
        return g1.m1 == g2.m1 and g1.m2 == g2.m2 and g1.alpha == g2.alpha

    def str_word(self, g):
        """Convert group atom g to a string

        """
        s1, s2 = str(g.m1), str(g.m2)
        al = ',*' if g.alpha else ''
        return "BiMM(%s, %s, %s)" % (s1, s2, al)


StdBiMMGroup = BiMMGroup()   # This is the only instance of AutP3Group

BiMM.group = StdBiMMGroup

 

#####################################################################
# Enter points and lines into the BiMonster
#####################################################################

      
def BiMM_Coxeter_Exp(x1, x2):
    mi, ma = min(x1,x2), max(x1, x2)
    assert 0 <= mi <= ma < 26
    if mi < 13 and ma >= 13 and (mi + ma) % 13 in (0,1,3,9):
         return 3
    return 2 if mi != ma else 1


PL_VALUES = None
PL_DATA = None
MAXLEN_PL_DATA = 0
ALPHA = BiMM([],[], 1)

def points_lines_list():
    global PL_VALUES, PL_DATA, MAXLEN_PL_DATA
    PL_VALUES = []
    s, t, u, v, x =  Norton_generators()
    for i in range(13):
        P_i = PointP3(list(range(13)) + [i])
        bm = BiMM(P_i, P_i, 1)
        PL_VALUES.append(bm)
    for i in range(13):
        L_i = (v*x) ** (u**i)
        bm = ALPHA * BiMM(L_i, L_i ** (-1))
        PL_VALUES.append(bm)
    data = [[list(bm.m1.mmdata),list(bm.m2.mmdata)] for bm in PL_VALUES]
    MAXLEN_PL_DATA = max(max(len(m1), len(m2)) for m1, m2 in data)
    #print(maxlen)
    PL_DATA = np.zeros((26, 2, MAXLEN_PL_DATA), dtype = np.uint32)
    for i in range(26):
        for j in range(2):
            d = data[i][j]
            PL_DATA[i, j, :len(d)] = d
points_lines_list()

 


def P3_BiMM(pl):
    a = np.zeros((2, len(pl), MAXLEN_PL_DATA), dtype = np.uint32)
    pl = p3_list(pl)
    len_pl = len(pl)
    for i in range(0, len_pl & -2, 2):
        k = pl[i]
        a[0,i] = PL_DATA[k,0]
        a[1,i] = PL_DATA[k,1]
        k = pl[i+1]
        a[0,i+1] = PL_DATA[k,1]
        a[1,i+1] = PL_DATA[k,0]
    al = len_pl & 1
    if al:
        k = pl[-1]
        a[0,-1] = PL_DATA[k,0]
        a[1,-1] = PL_DATA[k,1]
    return BiMM(MM('a', a[0].ravel()), MM('a', a[1].ravel()), al)
   



#####################################################################
# Tests
#####################################################################



def P3_BiMM_slow(pl):
    result = BiMM()
    for x in pl:
        result *= PL_VALUES[x]
    return result


def test_P3_BiMM():
    print("Testing function P3_BiMM")
    LL = list(range(26)) * 5
    for i in range(20):
        for j in range(3):
            pl = sample(LL, 3)
            assert P3_BiMM(pl) == P3_BiMM_slow(pl)
    print("passed")
  

          
    
def test_Coxeter_orders_for_one_node(x):
    for y in range(26):
        o = (PL_VALUES[x] * PL_VALUES[y]).order()
        e = BiMM_Coxeter_Exp(x, y)         
        assert e == o, (x, y, e, o)    

def test_Coxeter_orders(multiprocessing = False):
    print("Testing  orders of Coxeter generators of BiMM")
    if multiprocessing:
        from multiprocessing import Pool
        with Pool() as pool:
            pool.map(test_Coxeter_orders_for_one_node, range(26))
        pool.join()
    else:
        for x in range(26):
            test_Coxeter_orders_for_one_node(x)
    print("passed")

     




def random_hexagon():
    inc = incidence
    while 1:
        points = sample(range(13), 3)
        if not find_collinear_points(points):
            p1, p2, p3 = points
            return p1, inc(p1, p2), p2, inc(p2, p3), p3, inc(p3, p1)


def check_hexagon_relation(u, v, w, x, y, z):
    """Check relations in a hexagon in the BiMonster.

    Here ``(u, v, w, x, y, z)`` must be a hexagon in the 
    incidence graph of P3. We take the relation for that
    hexagon from [Iva99], Theorem 8.2.2.
    """
    a = P3_BiMM([u, x, v, y, w, z])
    assert a.order() == 4

def test_hexagon_relations(ntests = 10):
    print("Testing some (random) hexagon relations in P3")
    for i in range(ntests):
        check_hexagon_relation(*random_hexagon())
    print("passed")


def test_spider_relation():
    print("Testing the sipder relation in Y_555")
    spider = P3_BiMM('a,b1,c1,a,b2,c2,a,b3,c3')
    assert spider.order() == 10, spider.order()
    print("passed")
    


def test_all():
    test_P3_BiMM()
    test_Coxeter_orders(multiprocessing = True)
    test_hexagon_relations() 
    test_spider_relation()


if __name__ == "__main__":
    test_all()


  