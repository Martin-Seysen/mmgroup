"""Auxilary class an functions for testing the union-find algorithm.

"""

from collections import defaultdict
from numbers import Integral
from copy import deepcopy
from random import randint
import numpy as np

from mmgroup.generators import gen_ufind_lin2_dim
from mmgroup.generators import gen_ufind_lin2_gen
from mmgroup.generators import gen_ufind_lin2_n_orbits
from mmgroup.generators import gen_ufind_lin2_orbits
from mmgroup.generators import gen_ufind_lin2_n_gen
from mmgroup.generators import gen_ufind_lin2_representatives
from mmgroup.generators import gen_ufind_lin2_map_v_gen


from mmgroup.clifford12 import bitmatrix64_vmul
from mmgroup.clifford12 import bitmatrix64_mul
from mmgroup.clifford12 import bitmatrix64_inv



from mmgroup.general import Orbit_Lin2
a = None

#####################################################################
# Check result of a C function
#####################################################################


def chk(res, expected = None):
    """Check result 'res' of a C function

    If ``expected`` is None, the result ``res`` must be a nonnegative
    number. Otherwise ``res`` must be equal to ``expected``.
    """
    ok = (expected is None and res >= 0) or  res == expected
    if ok:
        return res
    if not ok:
        if isinstance(a, np.ndarray):
            a0 = int(a[0])
            a0 = a0 if a0 < 0x80000000 else a0 - 0x100000000
            print("Dump of array a", a[:6], ", a[0] =",  a0)
        if expected is None:
            err = "Result of C function is %d" % res
        else:
            err = "Result obtained: %d, expected: %d" % (res, expected)
        raise ValueError(err)



#####################################################################
# Bit matrix operations
#####################################################################

def vmatmul(v, m, n):
    """Multiply bit vector v with n times n bit matrix m"""
    w = 0
    for i in range(n):
        if (v >> i) & 1:
            w ^= m[i]
    return w & ((1 << n) - 1)

def is_inv(m, mi, n):
    """Return True if mi is the inverse of the n times n bit matrix m"""
    acc = 0
    for i in range(n):
        acc |= vmatmul(m[i], mi, n) ^ (1 << i)
    return acc == 0


def rand_bitmatrix(n):
    """Return an invertible n-dimensional random bit matrix"""
    assert 1 <= n <= 24
    mask = (1 << n) - 1
    while(1):
        a = np.array([randint(0, mask) for i in range(n)],
            dtype = np.uint32)
        a1 = np.array(a, dtype = np.uint64)
        try:
            bitmatrix64_inv(a1)
            return a
        except:
            continue



#####################################################################
# BitMatrix
#####################################################################



class BitMatrix:
    r"""Model an affine operation on a vector space.

    The constructor takes two parameters ``a``, ``b``, with ``a``
    :math:`n \times n` bit matrix given as an array of ``n`` bit
    vectors, and ``b`` a bit vector. Bit vectors are encoded as
    integers as usual. A pair ``(a, b)`` encodes the mapping
    :math:`v \mapsto v \cdot a + b`. ``b`` defaults to 0.

    Instances of this class may be multiplied or raised to the
    power of an integer. They may be left multipled with an
    integer denoting a bit vector.

    The length of an instance is the dimension of matrix ``a``.
    """ 
    def __init__(self, a, b = 0):
        self.dim = len(a)
        assert 1 <= self.dim <= 24
        self.m = np.array(a, dtype = np.uint32)
        self.m = np.append(self.m, np.uint32(b & 0xffffffff))
        self.m &= self.mask
        self.m[self.dim] |= 1 << self.dim
    @property
    def matrix(self):
        """Return linear part as numpy array"""
        return self.m[:self.dim]
    @property
    def const(self):
        """Return constant part"""
        return int(self.m[self.dim] & self.mask)
    @property 
    def mask(self):
        return (1 << self.dim) - 1
    @property 
    def _m64(self):
       """Return self.m as numpy array of dtype np.uint64"""
       return np.array(self.m, dtype = np.uint64)
    def __mul__(self, other):
        if isinstance(other, BitMatrix):
            a = bitmatrix64_mul(self._m64, other._m64)
            b = int(a[self.dim]) & self.mask
            return BitMatrix(a[:self.dim], b)
        return NotImplemented
    def __pow__(self, e):
        if e > 1:
            sqr = self ** (e // 2)
            return sqr * sqr * self if e & 1 else sqr * sqr
        if e == 1:
            return self
        if e == 0:
            return BitMatrix([1 << i for i in range(self.dim)])
        if e == -1:
            a = bitmatrix64_inv(self._m64)
            b = int(a[self.dim]) & self.mask
            return BitMatrix(a[:self.dim], b)
        if e < -1:
            return (self ** -1) ** e
        return NotImplemented
    def __rmul__(self, other):
        if isinstance(other, Integral):
            b = other | (1 << self.dim)
            res = bitmatrix64_vmul(b, self._m64, self.dim + 1)
            return res & self.mask
        return NotImplemented
    def __len__(self):
        return self.dim
    def __eq__(self, other):
        assert isinstance(other, BitMatrix)
        n = self.dim
        if other.dim != n:
            return False
        return (self.m == other.m).all()
    def map_lin2(self):
        """Map operation to pair suitable for class Orbit_Lin2"""
        return self.matrix, self.const
    def __str__(self):
        try:
            s = "<BitMatrix " + str(self.matrix)
            if self.const:
               s += " + " + hex(self.const)
            return s + ">"
        except:
            return str(self)
    __repr__ = __str__


#####################################################################
# Function constructing instances of class BitMatrix
#####################################################################


def PermBitMatrix(perm, b = 0):
    return BitMatrix([1 << i for i in perm], b)


def RandBitmatrix(dim, affine):
    a = rand_bitmatrix(dim)
    b = randint(0, (1 << dim) - 1) if affine else 0
    return BitMatrix(a, b)

def BitMatrix_from_Orbit_Lin2(a, i):
    """Return i-th generator from orbit array as Bitmatrix

    This is a low-level function
    """
    if isinstance(a, Orbit_Lin2):
        a = a.a
    n = chk(gen_ufind_lin2_dim(a))
    m = np.zeros(n + 1, dtype = np.uint32)
    chk(gen_ufind_lin2_gen(a, i, m, n + 1))
    return BitMatrix(m[:n], m[n])


#####################################################################
# Functions dealing with orbits
#####################################################################


def v_mul_g(a, v, g):
    """Multiply vector v with group word g stored in orbit array a

    This is a low-level function
    """
    if isinstance(a, Orbit_Lin2):
        a = a.a
    n = chk(gen_ufind_lin2_dim(a))
    m = np.zeros(n + 1, dtype = np.uint32)
    for i in g:
        chk(gen_ufind_lin2_gen(a, i, m, n))
        v = vmatmul(v, m, n) ^ m[n]
    return v


def rand_pair_orbit(a):
    """Return pair of random vectors in the same orbit

    Here ``a`` must be an instance of class Orbit_Lin2.
    This will fail for a compressed orbit array ``a``,
    and for an array where all orbits have length 1.
    """
    assert isinstance(a, Orbit_Lin2)
    mask = (1 << a.dim) - 1
    d = {}
    for _ in range(max(100, mask)):
        v = randint(0, mask)
        rep = a.orbit_rep(v)
        if rep in d and d[rep] != v:
            return v, d[rep]
        d[rep] = v
    raise ValueError("No pair in same orbit found")
        

#####################################################################
# Test the Schreier vector in an orbit array
#####################################################################


def do_test_scheierv(orbit, verbose = 1):
    """Test Schreier vector in instance ``orbit`` of ``Orbit_Lin2``"""
    ERR_NREP = "Schreier vectors ends at a non-representative "
    ERR_ILLEG = "Illegal Schreier vector entry at "
    ERR_CYCLE = "Schreier vector has a cycle containing "
    def test_v(v):
        """Return None if test done or a new vector to be tested""" 
        if v in tested:
             return None
        g = chk(gen_ufind_lin2_map_v_gen(a, v))
        if g >= 2 * n_gen:
            if g == 0xfe:
                if v in reps:
                    tested.add(v)
                    return None
                raise ValueError(ERR_NREP + hex(v))  
            else:
                raise ValueError(ERR_ILLEG + hex(v))
        tested.add(v) 
        return v * gen[g]

    assert isinstance(orbit, Orbit_Lin2)
    a = orbit.a
    if a[0] != 2:
        raise ValueError("Orbit has bad status", a[0])
    n = orbit.dim
    n_gen = len(orbit.generators())
    if verbose:
        print("\nCkecking Schreier vector, dim = %d, %d generators" %
             (n, n_gen))
    gen = [BitMatrix_from_Orbit_Lin2(a, i) for i in range(2 * n_gen)]
    reps = orbit.representatives()[0]
    tested = set()
    for v in range(1 << n):
        v1 = v
        for i in range(1 << n):
            v1 = test_v(v1)
            if v1 is None:
                break
        if i >= (1 << n) - 1:
            raise ValueError(ERR_CYCLE + hex(v))
    if verbose:
        print("passed")  

#####################################################################
# Test the orbit decomposition in an orbit array
#####################################################################



def lin2_orbits_raw(a, check_orbits = True, verbose = 0):
    r"""Store orbits of a given orbit array ``a``

    Here ``a`` should be be an instance of class ``Orbit_Lin2``.
    The function returns a (sorted) list of lists that defines the
    partition of the vectors into orbits.

    If ``check_orbits`` is True the correctness of the returned result
    is checked. This may take a long time.
    """
    if isinstance(a, Orbit_Lin2):
        a = a.a
    l_data = 1 << chk(gen_ufind_lin2_dim(a))
    n_sets = gen_ufind_lin2_n_orbits(a)
    data = np.zeros(l_data, dtype = np.uint32)
    indices = np.zeros(n_sets + 1, dtype = np.uint32)
    status = gen_ufind_lin2_orbits(
        a, data, l_data, indices, n_sets + 1)
    assert status >= 0, (2, status, l_data, n_sets)
    assert indices[n_sets] == l_data, "Bad orbit partition"
    assert set(data) == set(range(l_data)), "Bad orbit partition"
    result = []
    for i in range(n_sets):
        result.append(data[indices[i] : indices[i+1]])
    if verbose:
        print("Orbit array contains %d orbits:" % len(result))
        for data in result:
            print(data)
    if check_orbits:
        gen = [BitMatrix_from_Orbit_Lin2(a, i) for i in range(
            0, 2*gen_ufind_lin2_n_gen(a), 2)]
        for orbit in result:
            s, s_new = set(), set([orbit[0]])
            while len(s_new):
                v = s_new.pop()
                s.add(v)
                for g in gen:
                    v_new = v * g
                    if not (v_new in (s | s_new)):
                        s_new.add(v_new)
            assert s == set(orbit), "Orbits are wrong"

    reps = [o[0] for o in result]
    a_rep = np.zeros(l_data, dtype = np.uint32)
    l = gen_ufind_lin2_representatives(a, a_rep, len(a_rep))
    assert (reps == a_rep[:l]).all()
  
    return result




#####################################################################
# Test the orbit decomposition in an orbit array
#####################################################################


def do_test_orbit_array(g, g_list = None, chk = 0, verbose = 0):
    r"""Test orbits and Schreier vector for orbit array ``g``

    ``g`` must be an instance of class ``Orbit_Lin2``.

    if ``g_list`` is a list of instances of class ``BitMatrix``
    the we check if this is equal to the list of generators
    stored in ``g``. Otherwise ``g_list`` must be None.

    If ``check_orbits`` is True the the correctness of the orbit
    decomposition is checked. This may take a long time.
    """
    # General test stuff
    assert isinstance(g, Orbit_Lin2) 
    dim = g.dim
    if g_list:
        generators = g.generators()
        assert len(g_list) == len(generators)
        for i, gen in enumerate(generators):
            gen_i = gen**(-1)
            assert gen == g_list[i]
            assert gen * gen_i == gen**0
            assert BitMatrix_from_Orbit_Lin2(g, 2*i) == gen
            assert BitMatrix_from_Orbit_Lin2(g, 2*i+1) == gen_i
    g1 = deepcopy(g)
    if chk:
        lin2_orbits_raw(g1, True, verbose)
    g1._finalize()
    do_test_scheierv(g1, verbose = verbose)

    mask = (1 << dim) - 1
    for i in range(3):
        v1, v2 = rand_pair_orbit(g1)
        g_t = g1.map_v_G(v1, v2)
        assert v1 * g_t == v2 ,(i, v1, v2)




