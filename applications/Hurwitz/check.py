import sys
import os
import time
from datetime import datetime
from collections import OrderedDict, defaultdict
from random import sample, randint, shuffle
from multiprocessing import Pool

import numpy as np

sys.path.append(r".")
from find_generators import  MM


########################################################################
# Statistics of orders of element of a subgroup of the Monster
########################################################################


class MonsterOrderStatistics():
    def __init__(self):
        self.n = 0
        self.orders = defaultdict(int)
        self.verifiers = {}
        self.verifier_found = False

    def add_order(self, order, r = None):
        self.n += 1
        self.orders[order] += 1
        if order in [41, 59, 71] and not self.verifier_found:
            self.verifiers[order] = r
            self.verifier_found = True
        if order == 94:
            self.verifiers[order] = r

    def is_Monster(self):
        high_orders = sum(self.orders[x] for x in [41, 59, 71])
        return high_orders > 0 and self.orders[94] > 0

    def is_probably_not_Monster(self):
        n = self.n
        if n < 50: return False
        low_orders = sum(self.orders[x] for x in range(20)) / n 
        high_orders = sum(self.orders[x] for x in range(50,120)) / n
        return low_orders > 0.2 and high_orders < 0.2

    def get_verifiers(self):
        return self.verifiers if len(self.verifiers) >= 2 else {}


########################################################################
# Verifying set of generators of a Hurwitz subgroup of the Monster
########################################################################


def mm_product(mm_list):
    """Return product of Monster group elements of a list

    Here ``mm_list`` must be a list (or a tuple) of instances of
    class ``mmgroup.MM``. The function returns the product of these
    elements as an instance of class ``mmgroup.MM``.
    """
    a = np.concatenate([g.mmdata for g in mm_list])
    return MM('a', a).reduce()

class HurwitzVerifyer:
    z = MM('x', 0x1000)   # central involution of G_x0
    N = 64                # No of factors of random element of Monster
    ORDER_MSG = "order %3d: g_sample is 0x%016x"
    NEGATIVE_MSG = "Could not prove that the Monster is a Hurwitz group"
    POSITIVE_MSG = """Presentation of the Monster as a Hurwitz group. Let
  z  = %s,
  g3 = %s,
  ge = %s,
  g  = ge**(-1) * g3 * g.
Then z has order 2, g has order 3, and z*g has order 7.
z and g generate the Monster.
"""


    
    def __init__(self, g3, ge):
        z = self.z               # central involution z of G_x0 
        self.g3 = g3 = MM(g3)    # input g3, should have order 3 in MM 
        self.ge = ge =  MM(ge)   # input ge
        self.g = g = g3**ge      # element g = g3**ge of order 3 in MM
        self.zg = z * g          # element g_0 = z * g
        self.zgi = z * g**(-1)   # element g_1 = z * g**(-1)
        self.zg_list = [self.zg, self.zgi]  # The list [g_0, g_1]
        self._dict4 = {}         # map (b3,..b0)_2 to g_b0 * ... * g_b3
        self._dict8 = {}         # map (b7,..b0)_2 to g_b0 * ... * g_b7
        self.is_Monster = False  # True if z and g generate the Monster MM
        self.verifiers = {}      # dict {o:x}, with x a binary number
                                 # (b63,...,b0)_2 such that 
                                 # g_b0 * ... * g_b63 has order o

    def generators(self):
        return self.z, self.g
        

    def check_orders(self):
        assert self.z.order() == 2
        assert self.g.order() == 3
        assert self.zg.order() == 7

    def dict4(self, n):
        try:
            return self._dict4[n] 
        except KeyError:
            assert 0 <= n < 16 
            factors = [self.zg_list[(n >> i) & 1] for i in range(4)]
            self._dict4[n] = mm_product(factors)
            return self._dict4[n]
    
    def dict8(self, n):
        try:
            return self._dict8[n] 
        except KeyError:
            assert 0 <= n < 256 
            factors = [self.dict4((n >> i) & 15) for i in (0,4)]
            self._dict8[n] = mm_product(factors)
            return self._dict8[n]
 
    def product(self, n):
        assert 0 < n < 1 << self.N
        indices = [(n >> i) & 255 for i in range(0, self.N, 8)]
        factors =  [self.dict8(x) for x in indices]
        return mm_product(factors)

    def order_rand(self):
        r = randint(0, (1 << self.N) - 1)
        h = self.product(r)
        return self.N, r, h.order()


    def find_order_rand(self, n = 1000, verbose = 0):
        statistics = MonsterOrderStatistics()
        mask = 0
        self.is_Monster = False
        for i in range(n):
            _, r, order = self.order_rand()
            statistics.add_order(order, r)
            if verbose:
                print(self.ORDER_MSG % (order, r))
            if statistics.is_Monster():
                self.verifiers = statistics.get_verifiers()
                self.is_Monster = True
                break
            if statistics.is_probably_not_Monster():
                break
        self.verify_trials = statistics.n
        return self.is_Monster

    def verify(self, n_trials = 1000, verbose = 0):
        self.check_orders()
        return self.find_order_rand(n_trials, verbose)


    def display_verification(self):
        if not self.is_Monster:
            print(self.NEGATIVE_MSG)
            return
        print(self.POSITIVE_MSG  % (self.z, self.g3, self.ge))
        ORDER_VERIFY_MSG = "Order verfications (%d orders tested):"
        print (ORDER_VERIFY_MSG % self.verify_trials)
        for order, r in self.verifiers.items():
            print(self.ORDER_MSG % (order, r))
            
        
########################################################################
# Some generating a subgroup of the Monster
########################################################################


g3 = 'M0<y_74eh*x_1289h*d_29ah*p_12584189*l_2*p_2391360>'
ge = 'M0<y_51bh*x_9c9h*d_0f93h*p_233471242*t_1*p_169474141*l_2*t_2*p_83125779*l_2*t_1*p_20529902*l_2*t_2*p_201168724*l_1*t_1*p_64766107*l_2*t_2*p_13175458*l_1*t_1*p_77830155*l_1>'



g3 = 'M0<y_3h*x_0a94h*d_552h*p_33396464*l_2*p_3281280>'
ge = 'M0<y_491h*x_1c11h*d_0edeh*p_198726816*t_2*p_129314519*l_2*t_1*p_167991521*l_2*t_1*p_196142570*l_1*t_2*p_121538749*l_1*t_2*p_201007290*l_2*t_2*p_79665729*l_1*t_1*p_42715090*l_1>'

g3 = 'M0<y_784h*x_1416h*d_0bach*p_32531246*l_1*p_6547200>'
ge = 'M0<y_37ah*x_0a36h*d_0cdeh*p_187160018*t_2*p_127556155*l_2*t_2*p_167461608*l_1*t_1*p_217416341*l_1*t_1*p_100596922*l_1*t_2*p_242019761*l_2*t_1*p_6660783*l_2*t_1*p_198436069*l_2>'

g3 = 'M0<y_738h*x_16f9h*d_0fb8h*p_33444865*l_1*p_3302400>'
ge = 'M0<y_37h*x_1fceh*d_0dbfh*p_24749007*t_2*p_199649670*l_1*t_1*p_3497192*l_1*t_1*p_183457350*l_2*t_2*p_211869746*l_1*t_1*p_231080318*l_1*t_2*p_83239144*l_1*t_1*p_24788002*l_1>'

g3 = 'M0<y_2d5h*x_8f6h*d_79dh*p_31127249*l_2*p_27840>'
ge = 'M0<y_245h*x_6f5h*d_919h*p_190113051*t_1*p_138040333*l_1*t_2*p_12036469*l_2*t_1*p_20261480*l_2*t_2*p_109245691*l_1*t_1*p_86954362*l_2*t_2*p_110098494*l_1*t_1*p_123914620*l_1>'

g3 = 'M0<y_0bah*x_1896h*d_0bf0h*p_2340037*l_2*p_1168320>'
ge = 'M0<y_60fh*x_19c5h*d_0b78h*p_211232352*t_2*p_77937845*l_2*t_1*p_180496119*l_1*t_2*p_183535467*l_1*t_1*p_70180627*l_1*t_2*p_180641970*l_2*t_1*p_49210458*l_1*t_1*p_173613137*l_2>'


g3 = 'M0<x_25ah*d_48bh*p_932725*l_2*p_299520>'
ge = 'M0<y_515h*x_880h*d_0c11h*p_178102840*t_2*p_224426777*l_1*t_1*p_234052114*l_2*t_2*p_80553737*l_1*t_1*p_207786986*l_2*t_1*p_170832908*l_1*t_2*p_242391017*l_1*t_2*p_111166208*l_1>'


########################################################################
# Some samples generating the Monster
########################################################################


g3 = 'M0<y_0dah*x_0c8bh*d_468h*p_45301176*l_2*p_79432320>'
ge = 'M0<y_751h*x_1b4bh*d_0e29h*p_4777506*t_1*p_208261539*l_1*t_2*p_50724252*l_2*t_2*p_92719177*l_1*t_2*p_125129694*l_1*t_2*p_69104434*l_2*t_2*p_218887867*l_1*t_2*p_211978497*l_1>'

########################################################################
# Main program
########################################################################






if __name__ == "__main__":
    h = HurwitzVerifyer(g3, ge)
    h.check_orders()
    ok = h.verify(verbose = 1)
    h.display_verification()
    z, g = h.generators()
    if ok:
        print("Order of commutator [z, g] is", commutator(z,g).order())





