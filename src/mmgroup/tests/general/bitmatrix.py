"""Test C functions for the union-find algorithm.

In this module we test the C functions in module ``gen_ufind_lin2.c``
using their python wrappers in the ``mmgroup.generators`` extension.
In that file we define an *orbit array* containing orbits of the
action of a group on a linear space GF(2)^n. We also test the python
class ``Orbit_Lin2`` for wrapping such an orbit array. On the way we
also test the C functions in module ``gen_ufind_lin2.c`` implementing
a union-find algorithm used by the functions in ``gen_ufind_lin2.c``.

For testing the C functions we deal with the affine group H of
structure 2^8.L_3(2) acting naturally as a permutation group on the
affine space A8 corresponding to the linear space GF(2)^3. A8 has
size 8. Let P := 2**A8 be the power set of A8. P has a natural
structure an a vector space GF(2)^8. For our tests we compute the
orbits on the space P under the action of H using the functions in
module ``gen_ufind_lin2.c`` and also the python class ``Orbit_Lin2``.
These orbits are well known, so that that we can check the
correctness of the functions to be tested. We generate H as a
group of permutation matrices acting the vector space P.
"""

from collections import defaultdict
from numbers import Integral
import numpy as np
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
    """Multiply vector v with group word g stored in orbit array a"""
    n = chk(gen_ufind_lin2_dim(a))
    m = np.zeros(n, dtype = np.uint32)
    for i in g:
        chk(gen_ufind_lin2_gen(a, i, m, n))
        v = vmatmul(v, m, n)
    return v

class BitMatrix:
    """Model an affine operation on a vector space.

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

