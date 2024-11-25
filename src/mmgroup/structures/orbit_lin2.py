r"""We support abtract groups acting on certain sets.

.. warning::

  The functions in this module are yet experimental and subject to
  change in the future. We recommend to contact the author  before
  using them in research projects.

Let :math:`G` be a finite group and :math:`V` be a set. An action
of :math:`G` on :math:`V` is a mapping :math:`\rho` from :math:`G`
into the symmetric permutation group of the set :math:`V`. The
current version support actions on sets :math:`V` with a certain
structure only. Here :math:`V` may be

- A vector space :math:`V = \mbox{GF}_2^n, n \leq 24`, with  :math:`G`
  acting as a group of linear transformations on :math:`V`.
  Such an action is modelled by class ``Orbit_Lin2``.

- An elementary Abelian 2 group :math:`V` of structure
  :math:`2^n, n \leq 32`, with :math:`G` acting on :math:`V` via a
  homomorphism from :math:`G` to :math:`V`.
  Such an action is modelled by class ``Orbit_Elem2``.

In both cases a group is given by a set of generators
:math:`(g_{0} \ldots g_{k-1})` and a mapping :math:`\rho` that
maps each of these generators to an element of :math:`V`. The only
requirements for the elements :math:`g_{i}` are that they support
multiplication (for the group operation) and exponentiation with
an integer exponent (e.g. for computing inverses an the neutral
element). The mapping :math:`\rho` must be implemented as a
function taking a single group element :math:`g_{i}` and returning
:math:`\rho(g_{i})` as an element of :math:`V`. The encoding of
elemnts of :math:`V` is is described in the doccumentation of the
corresponding class.

The methods of these classes allow the computation of orbits of
:math:`V` under the action of :math:`G`. We also may compute
Schreier vectors on an orbits of :math:`V` (under the action of
:math:`G`), as described in :cite:`HE05`, Section 4.1. If :math:`G`
is acting on an elementary Abelian 2 group :math:`V`, we may
compute the kernel of that action instead.

The two types of actions mentioned above are desigend to
compute in the managable 2-local subgroups of the Monster, especially
in the groups :math:`G_{x0}` of structure :math:`2^{1+24}.\mbox{Co}_1`,
and in the group :math:`N_{0}` of structure
:math:`2^{2+11+22}.(\mbox{M}_{24} \times \mbox{S}_{3})`.

In principle the functionality in this module can be used for
computing the order of a subgroup of such a managable subgroup,
or to perform a constructive membership test for such a subgroup.
The current implementation should be considered as
*under construction* and may be extended in future versions.

Class ``Random_Subgroup`` in this module implements a simple
random generator for generating pseudorandom elements of a
group given by generators.
"""

import sys
import warnings
from numbers import Integral
from typing import Callable
from copy import deepcopy
from random import randint, shuffle, sample
from collections import defaultdict, OrderedDict
from collections.abc import Iterable
from operator import __mul__
from functools import reduce

import numpy as np

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
from mmgroup.generators import gen_ufind_lin2_compressed_size
from mmgroup.generators import gen_ufind_lin2_compress
from mmgroup.clifford12 import bitmatrix64_vmul
from mmgroup.clifford12 import bitmatrix64_echelon_h
from mmgroup.clifford12 import bitmatrix64_t
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.clifford12 import leech2matrix_subspace_eqn

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
-11: "Duplicate entry in orbit list",
}



ERR_GEN_UFIND_LIN2_DIM   = -7
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
    r"""Model orbits of a group acting on the vector space :math:`GF(2)^n`

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
    in :math:`V` is implemented as an integer :math:`v`, where the bit
    of valence :math:`2^i` in the binary representation of :math:`v`
    is bit :math:`i` of the vector. An :math:`n \times n` bit matrix
    is implemented as an array of :math:`n` integers, with each entry
    corresponding to a row vector of the matrix. The only requirement
    for elements of :math:`G` is that they can be multiplied and
    raised to the power of -1 (for inversion), 0, and 1.

    A list of generators of the group :math:`G` can be passed with
    parameter ``generators`` of the constructor. A generator may also
    be added with method ``add_generator``. The dimension :math:`n`
    of :math:`V` is obtained automatically by applying the function
    ``map`` given in the constructor to the first generator of the
    group. Here  ``1 <= n <= 24`` must hold. So we may e.g. compute
    orbits in the Leech lattice mod 2 under the action of the
    Conway group  :math:`\mbox{Co}_1`.

    There are methods for obtaining the orbits of :math:`V` under the
    action of :math:`G`, and for finding an element of :math:`G` that
    maps an element of :math:`V` to a given element in its orbit.

    After adding all generators, and before retriving any information
    about orbits, a set of Schreier vectors is computed (and stored
    inside an instance of this class) as described in :cite:`HE05`,
    Section 4.1. This can be done manually via method ``finalize``. It
    is done automatically by calling any methods obtaining infomration
    about orbits. After computing the Schreier vectors, no more
    generators can be added.
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
    @property
    def dim(self):
        r"""Return dimension ``n`` of the vector space over ``GF(2)``"""
        return chk(gen_ufind_lin2_dim(self.a))
    def finalize(self):
        if self.a is None:
            ERR = "No generators present in Orbit_Lin2 object"
            raise ValueError(ERR)
        chk(gen_ufind_lin2_finalize(self.a))
        self.gen_inverse = [x ** -1 for x in self.gen]
        self.gen_neutral = self.gen[0] ** 0
        return self
    def add_generator(self, g, use = True):
        r"""Add generator ``g`` to the object

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
        r"""Save the information stored in the object for pickling

        The function returns an pair ``(data, functions)``. Here
        ``data`` is a tuple of data that may be written to a file or
        a shelve; and ``functions`` is a tuple of functions used for
        reconstructing the object from ``data``. This process is called
        *pickling*. Note that this method finalizes an object before
        pickling it.

        Calling the constructor ``Orbit_Lin2(data, functions)`` creates
        a copy of the pickled object.

        Argument ``f_pickle_gen`` if an optional function used for
        pickling a list of generators of the group. If this is set
        then argument ``f_unpickle_gen`` should be a function that
        reverses the effect of function ``f_pickle_gen``.

        Note that pickling a function (e.g. for storing it in a shelve)
        in python is sometimes difficult. If you create and pickle an
        instance of this class without any generators then that
        instance is not finalized, and a pair ``(data, functions)``
        is returned as usual. In this case component ``functions`` of
        that pair is not data dependent and may be used for unpickling
        similar pickled objects.
        """
        if self.a is None:
            return (None, None), (self.map, f_unpickle_gen)
        self.finalize()
        len_a = chk(gen_ufind_lin2_check_finalized(self.a, len(self.a)))
        self.a[len_a : len_a + len(self.map_gen)] = self.map_gen
        len_a += len(self.map_gen)
        g = self.gen if f_pickle_gen is None else f_pickle_gen(self.gen)
        return (self.a[:len_a], g), (self.map, f_unpickle_gen)
    def n_orbits(self):
        r"""Return number of orbits on the vector space

        The function returns the number of orbits on the vector
        space under the action of the group.
        """
        return chk(gen_ufind_lin2_n_orbits(self.a))
    def generators(self):
        r"""Return the list of the generators of the group"""
        return self.gen
    def representatives(self):
        r"""Return representatives of orbits on the vector space

        The method returns a pair ``(reps, sizes)`` of arrays of 32-bit
        integers of equal length. Array ``reps`` stores an ordered list
        of representatives of the orbits on the vector space under the
        action of the group. These represetatives are vectors encoded as
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
        of the group. Vector ``v`` and the returned representative are
        encoded as integers as encoded as decribed in the constructor
        of this class.
        """
        return chk(gen_ufind_lin2_rep_v(self.a, v))
    def orbit_size(self, v):
        r"""Return the size of the orbit of the vector ``v``"""
        return chk(gen_ufind_lin2_len_orbit_v(self.a, v))
    def orbit(self, v):
        r"""Return orbit of the vector ``v`` as an array of integers"""
        o = np.zeros(self.orbit_size(v), dtype = np.uint32)
        chk(gen_ufind_lin2_orbit_v(self.a, v, o, len(o)))
        return o
    def map_v_word_G(self, v, img = None):
        r"""Yet to be documented!"""
        w0 = self._map_v_word_a(v, img)
        return [(self.map_gen[x >> 1], (-1) ** int(x)) for x in w0]
    def map_v_G(self, v, img = None):
        r"""Find group element transforming a vector inside an orbit

        let ``v`` and ``img``  be vectors of the representation in the
        same orbit under the action of the group, encoded as integers.
        The function returns a group element mapping ``v`` to ``img``.
        Parameter ``img`` defaults to the standard representative of
        the orbit.

        Internally, a Schreier vector is used for computing the
        group element.
        """
        w0 = self._map_v_word_a(v, img)
        g = self.gen_neutral
        gg = [self.gen, self.gen_inverse]
        for w in w0:
            g *= gg[w & 1][w >> 1]
        return g
    def mul_v_g(self, v, g):
        r"""Return product of vector ``v`` with group element ``g``"""
        rep = np.array(self.map(g), dtype = np.uint64)
        return bitmatrix64_vmul(v, rep, self.dim)

    def rand_stabilizer(self, v, size, r, n, history = False):
        r"""Return generators for stabilizer of vector ``v``

        Details are yet to be documented!
        """
        if history:
             raise NotImplementedError("History isnot implemented")
        rg = Random_Subgroup(self.gen, r, n, history)
        gen = []
        for i in range(size):
            g = rg.rand()
            img_v = self.mul_v_g(v, g)
            g1 = self.map_v_G(img_v, v)
            gen.append(g * g1)
        return gen
    def _map_v_word_a(self, v, img = None):
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
        if img is None:
            return g
        if self.orbit_rep(v) != self.orbit_rep(img):
            ERR = "Preimage and image vector are not in the same orbit"
            raise ValueError(ERR)
        g1 = self._map_v_word_a(img)
        return np.concatenate((g, np.flip(g1) ^ 1))
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
    def compress(self, orbits):
        r"""Yet to be documented!!!

        """
        n_gen = chk(gen_ufind_lin2_n_gen(self.a))
        o = np.array(orbits, dtype = np.uint32)
        l_c = chk(gen_ufind_lin2_compressed_size(self.a, o, len(o)))
        c = np.zeros(l_c + n_gen, dtype = np.uint32)
        l_c = chk(gen_ufind_lin2_compress(self.a, o, len(o), c, l_c))
        c[l_c : l_c + n_gen] = self.map_gen[:n_gen]
        pickled = (c[:l_c + n_gen], self.gen), (self.map, None)
        return Orbit_Lin2(*pickled)

class Orbit_Elem2:
    r"""Model orbits of a group acting on an elementary Abelian 2 group

    Let :math:`G` be a any group given by generators and let :math:`V`
    be an elementary Abelian 2 group. Furthermore, let :math:`\rho` be
    a homomporphism from :math:`G` to :math:`V`. An instance of this
    class models the mapping from :math:`G` to :math:`V` via the
    homomorphism :math:`\rho`. A function ``map`` mapping a generator
    :math:`g` of the group :math:`G` to the element :math:`\rho(g)` of
    :math:`V` must be given in the constructor for  this class.
    Elements of :math:`V` are implemented as unsigend 32-bit integers.
    Here an integer :math:`v` stands for a bit vector, where the bit
    of valence :math:`2^i` in the binary representation of :math:`v`
    is bit :math:`i` of the bit vector.

    A list of generators of the group :math:`G` can be passed with
    parameter ``generators`` of the constructor. The dimension
    :math:`n` of :math:`V` is obtained automatically by applying
    the function ``map`` given in the constructor to the all
    generators of the group.

    More details are yet to be documented!
    """
    NROWS = 32
    NCOLS = 32
    def __init__(self, map = None, generators = []):
        if map is None:
            map = lambda x : x
        if isinstance(map, Callable):
            self.gen = []
            self.m = np.zeros(self.NROWS, dtype = np.uint64)
            self.map_gen = np.zeros(self.NROWS + 1, dtype = np.uint8)
            self.map = map
            self._exp = 0
            self.solver = None
            for g in generators:
                self.add_generator(g)
        else:
            ERR = "Bad mapping function for class Orbit_Elem2"
            raise TypeError(ERR)
    @property
    def exp(self):
        r"""Return exponent of elementary Abelian 2 group"""
        return self._exp
    def add_generator(self, g):
        v = int(self.map(g))
        if (v & -0x100000000):
            chk(ERR_GEN_UFIND_LIN2_DIM)
        status = leech2matrix_add_eqn(self.m, self.exp, self.NCOLS, v)
        if status < 0 or status > 2:
            chk(-1001)
        self.map_gen[self._exp] = len(self.gen)
        self._exp += status
        self.gen.append(g)
    def generators(self):
        r"""Return the list of the generators of the group"""
        return self.gen
    def _map_v_word_a(self, v, img = 0):
        v ^= img
        w = leech2matrix_subspace_eqn(self.m, self.exp, self.NCOLS, v)
        if w < 0:
            ERR = "Preimage and image vector are not in the same orbit"
            raise ValueError(ERR)
        return w
    def map_v_word_G(self, v, img = 0):
        r"""Yet to be documented!"""
        w = self._map_v_word_a(v, img)
        return [(self.gen[k], 1) for i, k in enumerate(self.map_gen)
            if (w >> i) & 1]
    def map_v_G(self, v, img = 0):
        r"""Find group element transforming a vector inside an orbit

        let ``v`` and ``img``  be vectors of the representation in the
        same orbit under the action of the group, encoded as integers.
        The function returns a group element mapping ``v`` to ``img``.
        Parameter ``img`` defaults to the standard representative of
        the orbit.
        """
        w = self._map_v_word_a(v, img)
        if w == 0:
            return self.generators()[0] ** 0
        data = [self.gen[k] for i, k in enumerate(self.map_gen)
             if (w >> i) & 1]
        return reduce(__mul__, data)
    def rand_kernel(self, size, r, n, history = False):
        r"""Return generators for the kernel of the homomorphism :math:`\rho`

        Details are yet to be documented!
        """
        if history:
             raise NotImplementedError("History isnot implemented")
        rg = Random_Subgroup(self.gen, r, n, history)
        gen = []
        for i in range(size):
            g = rg.rand()
            v = self.map(g)
            g1 = self.map_v_G(v)
            gen.append(g * g1)
        return gen


class Random_Subgroup:
    r"""Generate random elements in a group given by generators

    We follow :cite:`HE05`, Section 3.2. The constructor of this
    class corresponds to Algorithm ``PrInitialize`` in that section;
    and the main method ``rand`` of this class corresponds to
    Algorithm ``PrRandom`` in that section. Method ``rand``
    returns a random element of the group.

    The constructor takes a set ``generators`` of generators of the
    group. Parameters ``r`` and ``n`` specify the number of group
    elements kept internally and the number of initialization
    steps, repectively, as in Algorithm  ``PrInitialize``.

    According :cite:`HE05`, ``r`` should be greater than 10, and
    ``n`` should be much greater than 50.

    If parameter ``history`` is ``True`` then we record the history
    of the random generator for future use. In principle this allows
    to represent the generated random elements as words in the
    original generators; but at present this feature is not supported.
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
            self.rand()

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
