r"""Yet to be documented


"""

import sys
import warnings
from numbers import Integral
from typing import Callable
from copy import deepcopy
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
from mmgroup.generators import gen_ufind_lin2_add
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
from mmgroup.generators import gen_ufind_lin2_orbit_lengths


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
-10: "Object is not in correct state for this method",
}


ERR_GEN_UFIND_LIN2_GEN   = -8


def chk(error):
    if error >= 0:
        return error
    if error in ERRORS:
        ERR = "Error in class Orbit_Lin2:"
        raise ValueError(ERR  + ERRORS[error])
    ERR = "Internal error %d in class  Orbit_Lin2"
    raise ValueError(ERR % error)
   

class Orbit_Lin2:
    r"""Model orbits of a group acting on the vector space GF(2)^n

    Let :math:`G` be a any group and let :math:`H` be a  group of
    :math:`n \times n` bit matrices acting on the vector space
    :math:`V = \mbox{GF}_2^n` by right multiplication. Furthermore,
    let :math:`\rho` be a homomporphism from :math:`G` to :math:`H`.
    An instance of this class models the action of the group
    :math:`G` as a permutation group on the vector space :math:`V`
    via the homomorphism :math:`\rho`.

    In the constructor of this class we just have to enter a callable
    function ``map`` that maps an element of  :math:`G`  to the
    :math:`n \times n` bit matrix :math:`\rho(G)`. Here a bit vector
    in :math:`\V` is implemented as an integer :math:`v`, where the bit
    of valence :math:`2^i` in the binary representation of :math:`v`
    is bit :math:`i` of the vector. An :math:`n \times n`bit matrix
    is implemented as an array of :math:`n` integers, with each entry
    corresponding to a row vector of the matrix. The only requirement
    for elements of :math:`G` is that they can be multiplied and
    raised to the power of -1 (for inversion), 0, and 1.

    A generator of the group can be added with method
    ``add_generator``. Then the dimension :math:`n` of :math:`V`
    is obtained automatically by applying the function ``map`` given
    in the constructor to the first group element added by method
    ``add_generator``. Here  ``1 <= n <= 24`` must hold. So we may e.g.
    compute orbits in the Leech lattice mod 2 under the action of the
    Conway group  :math:`\mbox{Co}_1`.

    There are methods for obtaining the orbits of :math:`V` under the
    action of :math:`G`, and for finding an element of :math:`G` that
    maps an element of :math:`V` to a given element in its orbit.

    After adding all generators, and before retriving information
    about orbits, a set of Schreier vectors is computed as described
    in :math:`HE05`, Section 4.1. This can be done manually via
    method ``finalize``. It is done automatically by calling any
    methods obtaining infomration about orbits. After computing
    the Schreier vectors, no more generator can be added.

    Todo: document pickle!!!
    """
    ERR_PIC = "Internal error %d in pickled data for class Orbit_Lin2"
    MAX_N_GEN = 127
    def __init__(self, map = None):
        if map is None:
            map = lambda x : x
        if isinstance(map, Callable):
            self.a = None
            self.gen = []
            self.map_gen = None
            self.map = map
            self.pickle_gen = lambda x : x
            self.unpickle_gen = lambda x : x
        else:
            a, gen, map_gen, map_, pickle, unpickle = map
            l_a = chk(gen_ufind_lin2_check_finalized(a, len(a)))
            self.a = a[:l_a]
            self.gen = unpickle(gen)
            self.map_gen = map_gen
            self.map = map_
            self.pickle_gen = pickle
            self.unpickle_gen = unpickle
    def finalize(self):
        if self.a is None:
            ERR = "No generators present in Orbit_Lin2 object"
            raise ValueError(ERR)
        chk(gen_ufind_lin2_finalize(self.a))
        try:
            self.a.flags.writeable = False
            self.map_gen.flags.writeable = False
        except:
            pass
        self.gen_inverse = [x ** -1 for x in self.gen]
        self.gen_neutral = self.gen[0] ** 0
    def add_generator(self, g, use = True):
        rep = np.array(self.map(g), dtype = np.uint32)
        if self.a is None:
            dim = len(rep)
            max_n_gen = 32 if dim <= 16 else self.MAX_N_GEN
            a_size = chk(gen_ufind_lin2_size(dim, max_n_gen))
            self.a = np.zeros(a_size, dtype = np.uint32)
            chk(gen_ufind_lin2_init(self.a, len(self.a), dim, max_n_gen))
            self.map_gen = np.zeros(max_n_gen, dtype = np.uint32)
            self.m_size = min(256, 1 << self.dim)
        num_gen = len(self.gen)
        num_gen_a = gen_ufind_lin2_n_gen(self.a)
        status = gen_ufind_lin2_add(self.a, rep, len(rep), use)
        dim = chk(gen_ufind_lin2_dim(self.a))
        if len(rep) != dim:
            ERR = "Representations of generators are incompatible"
            raise valueError(ERR)
        if status == ERR_GEN_UFIND_LIN2_GEN:
            newsize = chk(gen_ufind_lin2_size(dim, self.MAX_N_GEN))
            np.pad(self.a, (0, newsize - len(self.a)))
            chk(gen_ufind_lin2_pad(self.a, newsize, self.MAX_N_GEN))
            np.pad(self.map_gen, (0, self.MAX_N_GEN - len(self.map_gen)))
            status = gen_ufind_lin2_add(self.a, rep, len(rep), use)
        chk(status)
        self.gen.append(g)
        if status > 0:
            self.map_gen[num_gen_a] = num_gen
    def pickle(self):
        r"""Save the information stored in the object

        The function returns an (opaque) array of 32-bit integers
        containing the complete information stored in this object.
        This may may e.g. be written to a file.

        if ``data`` is an array returned by this method then 
        calling the constructor of this class with
        ``Orbit_Lin2('a', data)`` creates a clone of the current
        object.
        """
        self.finalize()
        gen = self.pickle_gen(self.gen)
        return (self.a, gen, self.map_gen, self.map,
            self.pickle_gen, self.unpickle_gen)
    @property
    def dim(self):
        r"""Return dimension ``n`` of the vector space over ``GF(2)``"""
        return chk(gen_ufind_lin2_dim(self.a))
    @property
    def n_gen(self):
        r"""Return number of generators of the group"""
        return chk(gen_ufind_lin2_n_gen(self.a))
    def n_orbits(self):
        r"""Return number of orbits on the vector space

        The function returns the number of orbits on the vector
        space under the action of the group.
        """
        return chk(gen_ufind_lin2_n_orbits(self.a))
    def gen(self, i, inverse = False):
        r"""Return the i-th generator of the group

        The generator is returned as an ``n`` times ``n`` bit matrix
        in an array of ``n`` 32-bit integers a desribed in the
        constructor of this class.

        If parameter ``inverse`` is ``True`` then the inverse of
        the generator is returned.
        """
        g = np.zeros(self.dim, dtype = np.uint32)
        return chk(gen_ufind_lin2_gen(self.a, 2*i + bool(inverse)))
    def representatives(self):
        r"""Return representatives of orbits on the vector space

        The method returns a pair ``(reps, sizes)`` of arrays of 32-bit
        integers of equal length. Array ``reps`` stores an ordered list
        of representatives of the orbits on the vector space under the
        action of the group. Thes represetatives are vectors encoded as
        integers as decribed in the constructor of of this class. Array
        ``sizes`` stores the sizes of these orbits. Both arrays have
        length ``self.n_orbits()``.
        """
        reps = np.zeros(self.n_orbits(), dtype = np.uint32)
        chk(gen_ufind_lin2_representatives(self.a, reps, len(reps)))
        sizes = np.zeros(self.n_orbits(), dtype = np.uint32)
        chk(gen_ufind_lin2_orbit_lengths(self.a, sizes, len(sizes))) 
        return reps, sizes
    def orbit_rep(self, v):
        r"""Return a representative of an orbit on the vector space

        Given a vector ``v`` in the vector space, the function returns
        the representative of the orbit of ``v`` under the the action
        of the group. Vector ``v`` and the retrurned vactor are
        encoded as integers as encoded as decribed in the constructor
        of this class.
        """
        return chk(gen_ufind_lin2_rep_v(self.a, v))
    def orbit_size(self, v):
        return chk(gen_ufind_lin2_len_orbit_v(self.a, v))
    def orbit(self, v):
        o = np.zeros(self.orbit_size(v), dtype = np.uint32)
        chk(gen_ufind_lin2_orbit_v(self.a, v, o, len(o)))
        return o
    def _map_v_word_a(self, v, image = None):
        self.finalize()
        while 1:
            g = np.zeros(self.m_size, dtype = np.uint8)
            res = gen_ufind_lin2_map_v(a, v, g, self.m_size)
            if res >= 0:
                g =  g[:res]
                break
            if res != ERR_OUT_SHORT or self.m_size >= 1 << self.dim:
                chk(res)
            self.m_size <<= 1
        if image is None:
            return g
        if self.rep(v) != self.rep(image):
            ERR = "Preimage and image vector are not in the same orbit"
            raise ValueError(ERR)
        g1 = self.map_v(image)
        return np.concatenate((g, np.flip(g1) ^ 1))

    def map_v_word_G(self, v, image = None):
        w0 = self._map_v_word_a(v, image = None)
        return [(self.map_gen(x) >> 1, -1 ** x) for x in w0]
    def map_v_G(self, v, image = None):
        w0 = self._map_v_word_a(v, image = None)
        g = self.gen_neutral
        gg = [self.gen, gen_inverse]
        for w in w0:
            g *= gg[w & 1, w >> 1]
        return g

