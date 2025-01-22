"""Auxilary class an functions for testing the union-find algorithm.

"""

from collections import defaultdict
from numbers import Integral
from random import randint
import numpy as np

from mmgroup.generators import gen_ufind_lin2_dim
from mmgroup.generators import gen_ufind_lin2_gen

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


