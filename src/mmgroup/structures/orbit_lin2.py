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
from mmgroup.generators import gen_ufind_lin2_n_max_gen
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
from mmgroup.clifford12 import bitmatrix64_vmul

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

    A list of generators of the group can be passed with parameter
    ``generators`` of the constructor. A generator may also be added
    with method ``add_generator``. The dimension :math:`n` of :math:`V`
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
    def __init__(self, map = None, generators = []):
        if map is None:
            map = lambda x : x
        if isinstance(map, Callable):
            self.a = None
            self.gen = []
            self.map_gen = None
            self.map = map
            for g in generators:
                self.add_generator(g)
        else:
            # Abuse arguments of constructor for passing a pickled object
            a, gen = map   # pickled data
            self.a = a
            self.map, unpickle = generators  # pickled functions
            if self.a is None:
                self.gen = []
                self.map_gen = None
            else:
                len_a = chk(gen_ufind_lin2_check_finalized(a, len(a)))
                self.map_gen = a[len_a:]
                self.gen = gen if unpickle is None else unpickle(gen)
        self.m_size = 32
    def finalize(self):
        if self.a is None:
            ERR = "No generators present in Orbit_Lin2 object"
            raise ValueError(ERR)
        chk(gen_ufind_lin2_finalize(self.a))
        try:
            #self.a.flags.writeable = False
            #self.map_gen.flags.writeable = False
            pass
        except:
            pass
        self.gen_inverse = [x ** -1 for x in self.gen]
        self.gen_neutral = self.gen[0] ** 0
        return self
    def _size_a_buf(self, n_gen, dim):
        """Size of main buffer for given No of generators an dimension"""
        return chk(gen_ufind_lin2_size(dim, n_gen)) + n_gen
    def _init_a_buf(self, dim = None):
        """Initialize main buffer for a given dimension"""
        assert self.a is None
        max_n_gen = min(4 * dim, self.MAX_N_GEN)
        a_size = self._size_a_buf(max_n_gen, dim)
        self.a = np.zeros(a_size, dtype = np.uint32)
        chk(gen_ufind_lin2_init(self.a, len(self.a), dim, max_n_gen))
        self.map_gen = np.zeros(max_n_gen, dtype = np.uint32)
    def _extend_a_buf(self, n_more = 1):
        """Extend main buffers for more generators

        The function tries to extend the buffers so that they can store
        up to ``n_more`` more generators. It raises ``ValueError``
        if this is not possible. The function returns the current
        number of generators.
        """
        assert self.a is not None
        n_gen = chk(gen_ufind_lin2_n_gen(self.a))
        n_max_gen = chk(gen_ufind_lin2_n_max_gen(self.a))
        new_n_gen = n_gen + n_more
        if new_n_gen > n_max_gen:
            new_n_max_gen = max(2 * n_max_gen, new_n_gen)
            new_n_max_gen = min(new_n_max_gen, self.MAX_N_GEN)
            if new_n_max_gen < new_n_gen:
                chk(ERR_GEN_UFIND_LIN2_GEN)
            new_a_size = self._size_a_buf(new_n_max_gen, self.dim)
            self.a = np.pad(self.a, (0, new_a_size - len(self.a)))
            chk(gen_ufind_lin2_pad(self.a, new_a_size, new_n_max_gen))
            pad_map_gen = (0, new_n_max_gen - n_max_gen)
            self.map_gen = np.pad(self.map_gen, pad_map_gen)
        return n_gen
    def add_generator(self, g, use = True):
        """Add generator ``g`` to the object

        Parameter ``use`` should usually be ``True``

        Todo: document parameter ``use``
        """
        rep = np.array(self.map(g), dtype = np.uint32)
        if self.a is None:
            self._init_a_buf(len(rep))
        n_gen = self._extend_a_buf(1)
        if len(rep) != chk(gen_ufind_lin2_dim(self.a)):
            ERR = "Representations of generators are incompatible"
            raise valueError(ERR)
        status = chk(gen_ufind_lin2_add(self.a, rep, len(rep), use))
        if status > 0:
            self.map_gen[n_gen] = len(self.gen)
        self.gen.append(g)
    def pickle(self, f_pickle_gen = None, f_unpickle_gen = None):
        r"""Save the information stored in the object

        The function returns an pair ``(data, functions)``. Here
        ``data`` is a tuple of data that may be written to a file
        or a shelve. ``functions`` is a tuple of functions.

        Caling the constructor ``Orbit_Lin2(data, functions)`` creates
        a copy of the current object. Before copying the object it is
        finalized if any generators have been added.

        Argument ``f_pickle_gen`` if an optional argument for pickling
        the list of generators of the group. If this is set then
        argument ``f_unpickle_gen`` should be a function that reverses
        the effect of functzion ``f_pickle_gen``.
        """
        if self.a is None:
            return (None, None), (self.map, f_unpickle_gen)
        self.finalize()
        len_a = chk(gen_ufind_lin2_check_finalized(self.a, len(self.a)))
        self.a[len_a : len_a + len(self.map_gen)] = self.map_gen
        len_a += len(self.map_gen)
        g = self.gen if f_pickle_gen is None else f_pickle_gen(self.gen)
        return (self.a[:len_a], g), (self.map, f_unpickle_gen)
    @property
    def dim(self):
        r"""Return dimension ``n`` of the vector space over ``GF(2)``"""
        return chk(gen_ufind_lin2_dim(self.a))
    @property
    def n_gen(self):
        r"""Return number of generators of the group

        Deprecated!!!!!
        """
        return chk(gen_ufind_lin2_n_gen(self.a))
    def n_orbits(self):
        r"""Return number of orbits on the vector space

        The function returns the number of orbits on the vector
        space under the action of the group.
        """
        return chk(gen_ufind_lin2_n_orbits(self.a))
    def gen(self, i, inverse = False):
        r"""Return the i-th generator of the group

        Deprecated!!!!

        The generator is returned as an ``n`` times ``n`` bit matrix
        in an array of ``n`` 32-bit integers a desribed in the
        constructor of this class.

        If parameter ``inverse`` is ``True`` then the inverse of
        the generator is returned.
        """
        g = np.zeros(self.dim, dtype = np.uint32)
        return chk(gen_ufind_lin2_gen(self.a, 2*i + bool(inverse)))
    def generators(self):
        """Return list of generators of the group"""
        return self.gen
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
        """Return the size of th orbit of the vector ``v``"""
        return chk(gen_ufind_lin2_len_orbit_v(self.a, v))
    def orbit(self, v):
        o = np.zeros(self.orbit_size(v), dtype = np.uint32)
        chk(gen_ufind_lin2_orbit_v(self.a, v, o, len(o)))
        return o
    def _map_v_word_a(self, v, image = None):
        self.finalize()
        while 1:
            g = np.zeros(self.m_size, dtype = np.uint8)
            res = gen_ufind_lin2_map_v(self.a, v, g, self.m_size)
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
        gg = [self.gen, self.gen_inverse]
        for w in w0:
            g *= gg[w & 1][w >> 1]
        return g

    def mul_v_g(self, v, g):
        rep = np.array(self.map(g), dtype = np.uint64)
        return bitmatrix64_vmul(v, rep, self.dim)

    def rand_stabilizer(self, v, size, r, n, history = False):
        if history:
             raise NotImplementedError("History isnot implemented")
        rg = Random_Subgroup(self.gen, r, n, history)
        gen = []
        for i in range(size):
            g = rg()
            img_v = self.mul_v_g(v, g)
            g1 = self.map_v_word_G(img_v, v)
            gen.append(g * g1)

class Random_Subgroup:
    r"""Generate random elements in a group given by generators

    We follow :math:`HE05`, Section 3.2. The constructor of this
    class corresponds to Algorithm ``PrInitialize`` in that section;
    and the main method ``rand`` of this class corresponds to
    Algorithm ``PrRandom`` in that section. Method ``rand``
    returns a ramdom element of the group.

    The constructor takes a set ``generators`` of generators of the
    group. Parameters ``r`` and ``n`` specify the number of group
    elements kept internally and the number of initialization
    steps, repectively, as in Algorithm  ``PrInitialize``.

    According math:`HE05`, ``r`` should be greater than 10, and
    ``n`` should be much greater than 50.

    If parameter ``history`` is ``True`` then we record the history
    of the random generator for future use. In principle this allows
    to represent the generated random elements as words in the
    original generators.
    """
    def __init__(self, generators, r, n, history = False):
        generators = list(generators)
        self.n_gen = len(generators)
        q, r = divmod(max(r, self.n_gen), self.n_gen)
        neutral = generators[0] ** 0
        self.x = [neutral] + generators * q  + generators[:r]
        self.history = None
        if history:
            assert len(self.x) < 0x8000
            self.nsteps = 0
            self.history = np.zeros((512, 2), dtype = np.uint16)
        for i in range(n):
            self.step()

    def rand(self):
        """Return a random element of the group"""
        rmax = 2 * len(self.x) - 1
        s = randint(2, rmax)
        t = randint(2, rmax - 2)
        if (t >= s):
            t += 2
        if self.history is not None:
            if self.nsteps >= len(self.history):
                self.history = np.pad(self.history,
                    ((0, self.nsteps), (0,0)))
            self.history[nsteps] = s, t
            nsteps += 1
        e = -1 ** t
        s1, t1 = s >> 1, t >> 1
        if s & 1:
            self.x[s1] = self.x[s1] * self.x[t1] ** e
            self.x[0] = self.x[0] * self.x[s1]
        else:
            self.x[s1] = self.x[t1] ** e * self.x[s1]
            self.x[0] = self.x[s1] * self.x[0]
        return self.x[0]
