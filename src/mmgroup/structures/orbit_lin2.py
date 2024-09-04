r"""Yet to be documented


"""

import sys
import warnings
from numbers import Integral
from copy import deepcopy
from math import floor, ceil
from random import randint, shuffle, sample
from collections import defaultdict, OrderedDict
from collections.abc import Iterable
import numpy as np

if __name__ == "__main__":
    sys.path.append("../../../")

import mmgroup
from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord
from mmgroup import Xsp2_Co1


from mmgroup.generators import gen_ufind_lin2_init
from mmgroup.generators import gen_ufind_lin2_size
from mmgroup.generators import gen_ufind_lin2_dim
from mmgroup.generators import gen_ufind_lin2_n_gen
from mmgroup.generators import gen_ufind_lin2_gen
from mmgroup.generators import gen_ufind_lin2_n_orbits
from mmgroup.generators import gen_ufind_lin2_orbits
from mmgroup.generators import gen_ufind_lin2_get_map
from mmgroup.generators import gen_ufind_make_map
from mmgroup.generators import gen_ufind_lin2_check
from mmgroup.generators import gen_ufind_lin2_len_orbit_v
from mmgroup.generators import gen_ufind_lin2_orbit_v
from mmgroup.generators import gen_ufind_lin2_rep_v
from mmgroup.generators import gen_ufind_lin2_map_v_gen
from mmgroup.generators import gen_ufind_lin2_map_v
from mmgroup.generators import gen_ufind_lin2_finalize
from mmgroup.generators import gen_ufind_lin2_check_finalized
from mmgroup.generators import gen_ufind_lin2_representatives
from mmgroup.generators import gen_ufind_lin2_get_table



ERRORS = {
-1 : "Out of memory",
-2 : "Too many entries for union-find algorithm",
-3 : "Input parameter to large",
-4 : "Output buffer too short",
-5 : "Entry not in union-find table",
-6 : "Union-find table too large",
-7 : "Dimension n of GF(2)^n is 0 or too large",
-8 : "Too many generators for subgroup of SL_2(2)^n",
-9 : "Generator matrix is not invertible",
-10, "Main buffer is not in correct state for this function",
}



def chk(error):
    if error > 0:
        return error
    if error in ERRORS
        ERR = "Error in class Orbit_Lin2:"
        raise ValueError(ERR  + ERRORS[error])
    ERR = "Internal error %d in class  Orbit_Lin2"
    raise ValueError(ERR % error)
   

class Orbit_Lin2:
    """This class is yet under construction"""
    slots = ["a", "m_size"]
    ERR_PIC = "Internal error %d in pickled data for class Orbit_Lin2"  
    def __init__(self, data, a = None)
        if data == "a":
            a = np.array(a, dtype = np.uint32)
            status = gen_ufind_lin2_check_finalized(a, len(a)) 
            if (status < 0):
                raise ValueError(self.ERR_PIC, status)
            self.a = a[:status]
        else:
            gen = np.array(data, dtype = np.uint32)
            n, k = gen.shape()
            a_size = chk(gen_ufind_lin2_size(n, k))
            a = np.zeros(a_size, dytpe = np.uint32)
            res = gen_ufind_lin2_init(a, len(a), n, gen.ravel(), k)
            chk(res)
        self.m_size = min(256, 1 << self.dim)
    def pickle(self):
        chk(gen_ufind_lin2_finalize(self.a, len(self.a)))
        return np.array(self.a)
    @property
    def dim(self):
        return chk(gen_ufind_lin2_dim(self.a))
    @property
    def n_gen(self):
        return chk(gen_ufind_lin2_n_gen(self.a))
    @property
    def n_orbits(self):
        return chk(gen_ufind_lin2_n_orbits(self.a))
    def gen(self, i, inverse = False):
        g = np.zeros(self.dim, dtype  np.uint32)
        return chk(gen_ufind_lin2_gen(self.a, 2*i + bool(inverse)))
    def representatives(self):
        reps = np.zeros(self.n_orbits, dtype = np.uint32)
        chk(gen_ufind_lin2_representatives(self.a, reps, len(reps)))
        return rep
    def rep(self, v):
        return chk(gen_ufind_lin2_rep_v(self.a, v))
    def orbit_size(self, v):
        return chk(gen_ufind_lin2_len_orbit_v(self.a, v))
    def orbit_sizes(self, v):
        reps = self.representatives()
        orbit_sizes = np.zeros(len(reps), dtype = np.uint32)
        for i, v in enumerate(reps)
            orbit_sizes[i] = chk(gen_ufind_lin2_len_orbit_v(self.a, v))
        return reps, orbit_sizes
    def orbit(self, v):
        o = np.zeros(self.orbit_zize(v), dtype = np.uint32)
        chk(gen_ufind_lin2_orbit_v(self.a, v, o, len(o)))
        return o
    def map_v(self, v):
        while 1:
            g = np.zeros(self.m_size, dtype = np.uint8)
            res = gen_ufind_lin2_map_v(a, v, g, self.m_size)
            if res >= 0:
                return g[:res]
            if res != ERR_OUT_SHORT or self.m_size >= 1 << self.dim:
                chk(res)
            self.m_size <<= 1


class Orbit_Gx0_Leech2(Orbit_Lin2):
    slots = ["generators"]
    L_GX0 = 10
    ERR_G = "Entry in constructor in class Orbit_Lin2 is not in group"  
    def __init__(self, data, a = None)
        if data == "a":
            a, self.generators = data
            self.a = np.array(a, dtype = np.uint32)
            status = gen_ufind_lin2_check_finalized(self.a, len(self.a)) 
            if (status < 0):
                raise ValueError(self.ERR_PIC, status)
            if self.dim != 24:
                raise ValueError(self.ERR_PIC, 503)
        else:
            n_gen = len(data)
            gen_lin = np.zeros((n_gen,24), dtype = np.uint32)
            gen_g = np.zeros((2 * n_gen, self.L_G_X0), dtype = np.uint32)
            for i, g in enumerate(data):
                if not isinstance(g, AbstractMMGroupWord):
                    raise TypeError(self.ERR_G)
                x = Xsp2_Co1(g.reduce())
                gen_lin[i] = x.leech_mod2_op()
                mm = x.mmdata
                gen_g[2 * i, :len(mm)] = mm
                mm = (x**-1).mmdata
                gen_g[2 * i + 1, :len(mm)] = mm
            super(Orbit_Gx0_Leech2, self).__init__(gen_lin)
            self.generators = gen_g
    def pickle(self):
        a = super(Orbit_Gx0_Leech2, self).pickle()
        return a, deepcopy(self.genetors)
    def gen(self, i, inverse = False):
        if not 0 <= i < self.n_gen:
            chk(ERR_IN_LARGE)
        ind = 2 * i * self.L_GX0 + bool(inverse)
        return Xsp2_Co1('a', self.generators[ind])
    def map_G_x0(self, v):
        x = Xsp2_Co1()
        b = self.map_v(v)
        for i in b:
            x.mul_data(self.generators[i])
        return x



