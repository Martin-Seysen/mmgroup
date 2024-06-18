"""Implements fast arithmetic of bit vectors

A bit vector b = (b[0],...,b[n-1]) , b[i] in (0,1) is coded as the
integer sum( [ b[i] * 2**i  for i in range(n) ] ). The is a bit 
inconvenient to  read for a human, but with this encoding functions 
operating on bit vectors can be defined and implemented in a much more 
consistent way. Function bin() also allows to print the bits of an 
integer in reverse order.  The first index of a vector or matrix is 
always zero, as usual in Python or C. 

Generally speaking, bit vectors are mainly used as modules with a
group operating on them. Elements of such groups may be matrices
or permutations on bit vectors. Here bit vectors are considered
as row vectors, and bit matrices or permutation always operate
on the right. Thus  x * g1 * g2, where x is a bit vector and g1,
g2 are elements of a group, means that first g1 is applied to v
and then g2 is applied to the bit vector v * g1.

Bit matrices are represented as lists of integers, where the i-th 
element of the list encodes the i-th row vector of the matrix. 
Function bit_mat_mul() multiplies an arbitrary number of
matrices. The first factor may also be a bitvector. In the same 
way as an integer in binary representation can be preceeded by
an arbitrary number of of zero bits (without changing its value),
we assume that a list representing a matrix can be followed by
an arbitrary number of zero entries so that multiplcation of
bit matrices is always possible.

A permutation p is also represented as a list of integers.
Here p[i] = k means that the bit at position i is mapped to
position k. More strictly, we will interpret a permutation p
represented by a list [ p[0], p[1],...,p][n-1] ] with
p[i] = k as a bit matrix b with b[i,j]  = 1 if j==k and
b[i,j] = 0 if j != k. Function  bit_perm_mul() multiplies a
bit vector or a permutation with an arbitrary number of
permutations. Function  bit_perm_2_mat() transforms a
permutation into an equivalent bit matrix.

""" 
from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint
import types
from operator import __or__ , __xor__
from functools import reduce
from collections.abc import Iterable
from numbers import Integral

import numpy
from numpy import array, zeros, uint8, uint16, uint32, uint64, int64


def isinteger(x):
    """Deprecated, use numbers.Integral() instead!"""
    return  isinstance(x, (Integral, numpy.integer))

def lmap(*args):
    """For compatibility with map() in python 2"""
    return list(map(*args))


def lrange(*args):
    """For compatibility with range() in python 2"""
    return list(range(*args))


def bitlen(x):
    """returns the length of the binary representation of an integer x >= 0"""
    assert x >= 0, "bitlen(x) requires nonegative integer x"
    try:
        return (int(x)).bit_length()
    except:  # for python versions < 2.7 or 3.1:
        l, m, sh = 0, 0x10000, 16
        while x >= m:   
            #x, l = x>>sh, l+sh
            m, sh = m << sh, sh << 1
        #assert m == 1 << sh
        while sh:
            if x >= m: x, l = x>>sh, l+sh
            sh >>= 1
            m = m >> sh  
        while x: x, l = x>>1, l+1
        return l



def hibit(x):
   """return the position of the highest bit of an integer"""
   if x==0: raise ZeroDivisionError("Zero has no highest bit")
   return bitlen(x) - 1





def v2(*args):
   """Return the 2-adic value of an integer or of a set of integers.

   For an integer ``i = o * 2**k``, ``o`` odd, we have ``v2(i) = k``.

   If ``x`` is iterable or a numpy array then ``v2(x)`` is the 
   minimum of the values ``v2(i)`` for all entries of ``i`` of ``x``.

   If several arguments ``x_1,...,x_n`` are given then 
   ``v2(x_1,...,x_n)`` is the minimum of all values ``v2(x_i)``.

   The function raises ZeroDivisionError if all integers occuring as 
   arguments (or contained in an argument) are zero.
   """
   x = 0
   for a in args:
       if isinstance(a, Integral):
           x |= int(a) 
       elif isinstance(a, numpy.ndarray):
           try:
               x |= int(reduce(__or__, a.reshape(-1), 0))
           except TypeError:
               if issubclass(a.dtype.type, np.integer):
                    # workaround for the case a.dtype == numpy.uint64
                    x |= reduce(__or__, map(int,a.reshape(-1)), 0)
               else:
                   raise
       elif isinstance(a, Iterable):
           x |= int(reduce(__or__, a, 0))
       else:
          err = "Function v2() undefined for a %s object"
          raise TypeError(err % type(a))
   l = (x & -x).bit_length() - 1
   if l >= 0:
       return l
   raise ZeroDivisionError("Function v2(x) is undefined for x = 0")



def bits2list(x):
   """returns the list of positions of one-bits of integer x

   List of bit positions is sorted in ascending order
   """
   l = bitlen(x)
   return [ i for i in range(l) if (x >> i) & 1 ]


def list2bits(lst):
   """return number x with bit i of x  set for all i in the list 'lst' """
   return reduce(__or__, [ 1<<i for i in lst ], 0) 





def list_xor(ListList):
   """Given a list of lists, the list of  XOR sums of the cartesian product are returned"""
   if len(ListList) <= 1:
       return ListList[0] if len(ListList) else []
   d = len(ListList) >>1
   l1, l2 = ListXor(ListList[:d]), ListXor(ListList[d:])
   return [int(x) ^ int(y) for x in l1 for y in l2]





def iter_bitweight(k,n):
    """Iterate through all integers 0 <= x < n with bit weight k"""
    if k > 0:
        x = (1 << k) - 1  # smallest integer of bit weight k
        while x < n:
            yield x
            # use Gosper's hack to get next integer x of bit weight k
            c = x & -x
            r = x + c
            x = (((x ^ r) >> 2) // c) | r
    elif n > 0:
        yield 0

def iter_ascending_bitweight(n):
    """Iterate through all integers 0 <= x < n by ascending bit weight"""
    if n > 0:
        yield 0
        k = 1
        while 1 << k <= n:
            for y in iter_bitweight(k,n):
                yield y  
            k = k + 1




def bw24(x):
    """Return the bit weight of the lowest 24 bits of x"""
    x = (x & 0x555555) + ((x & 0xaaaaaa) >> 1)
    x = (x & 0x333333) + ((x & 0xcccccc) >> 2)
    x = (x + (x >> 4)) & 0xf0f0f
    return (x + (x >> 8) + (x >> 16)) & 0x1f 


def bitweight(x):
   """return the bit weight of the  word 'x'. """
   assert x >= 0, "bitweight(x) requires integer x >= 0"
   w = 0
   while x:
        w += bw24(x)
        x >>= 24
   return w 

def bitparity(x):
   """return the bit parity of the  word 'x'. """
   assert x >= 0, "bitparity(x) requires integer x >= 0"
   while x > 0xffffffff:
       x = (x >> 32) ^ (x & 0xffffffff)
   x ^= x >> 16
   x ^= x >> 8
   x ^= x >> 4
   return (0x6996 >> (x & 15)) & 1




def reverse24(x, k=24):
    """reverses order of bits 0,...,k-1 of an integer. Default is k=24.

    If any bits apart form bit 0,...,k-1 are set, an execption is raised.
    """
    if x & -(1 << k):
        raise ValueError("Too high bits are set for bit reversal")
    return sum( ((x >> (k-i-1)) & 1) << i  for  i in range(k) )



def bin(x, l=24, d=8, reverse=False):
   """returns the binary representation of an 'l'-bit integer 'x' as a string

   Leading zeros are printed, blanks are inserted between 'd'-bit blocks.
   'd' may also be a list of bit positions where to insert blanks. If
   'reverse' is True, the string starts with bit 0 (which is fine for bit
   vectors), otherwise it starts with bit l-1 (as usual for integers). 
   """
   s = ""
   if not isinstance(d, list):  d=range(d,l,d)
   if reverse:
      for i in range(l):
         if x & (1<<i): s += '1'
         else: s += '0'
         if i+1 in d: s += ' '
   else:
      for i in range(l-1,-1,-1):
         if x & (1<<i): s += '1'
         else: s += '0'
         if i in d: s += ' '
   return s



def binbv(x, l=24, d=8):
   """shorthand for bin(x, l, d, reverse=True)"""
   return bin(x, l, d, reverse=True)

      


def pivot_binary_high(basis):
   """Pivot over a basis of integers considered as bit vectors

   Here 'basis' list a list of integers corresponding to bit vectors,
   and the elements of the list 'basis' are considered as the rows of a 
   bit matrix defining a basis of a linear space in a vector space 
   over {0,1}.

   The function performs Gaussian elimination to change the basis to 
   a standard form. Here pivoting over bit columns is done, so that
   at the end each pivoted bit column contains exactly one nonzero bit.
   Redundant rows of the original basis are removed.

   Bit columns are not exchanged. The function returns a pair
   (basis, columns), where 'basis' is the modified basis and 'columns'
   is a list of integers containing the index of the pivoted column
   used in the corresponding row of the basis. Here index i corresponds
   to the bit column with valence i.    

   The highest possible bit column is always taken for pivoting.  
   """
   basis, columns = list(basis)[:], []
   for i in range(len(basis)):
       m = max(basis[i:])
       j = basis[i:].index(m) + i
       try:
           mi = hibit(m) 
       except:
           del basis[i:]
           break
       else:
           mask = 1 << mi  
           columns.append(mi)
           basis[i], basis[j] = m, basis[i]
           for j in range(0, len(basis)):
               if j != i and int(basis[j]) & mask:
                    basis[j] ^= m
   return basis, columns   







def pivot_binary_low(basis):
   """Pivot over a basis of integers considered as bit vectors

   Here 'basis' list a list of integers corresponding to bit vectors,
   and the elements of the list 'basis' are considered as the rows of a 
   bit matrix defining a basis of a linear space in a vector space 
   over {0,1}.

   The function performs Gaussian elimination to change the basis to 
   a standard form. Here pivoting over bit columns is done, so that
   at the end each pivoted bit column contains exactly one nonzero bit.
   Redundant rows of the original basis are removed.

   Bit columns are not exchanged. The function returns a pair
   (basis, columns), where 'basis' is the modified basis and 'columns'
   is a list of integers containing the index of the pivoted column
   used in the corresponding row of the basis. Here index i corresponds
   to the bit column with valence i.    

   The lowest possible bit column is always taken for pivoting.  
   """
   basis, columns = basis[:], []
   for i in range(len(basis)):
       m = reduce(__or__, basis[i:], 0)
       try:
           mi = v2(m) 
       except:
           del basis[i:]
           break
       else:
           mask = 1 << mi  
           j = [x & mask for x in basis[i:]].index(mask) + i
           columns.append(mi)
           m = basis[j]
           basis[i], basis[j] = m, basis[i]
           for j in range(0, len(basis)):
               if j != i and basis[j] & mask: 
                   basis[j] ^= m
   return basis, columns   









_maxbits = {uint8:8, uint16:16, uint32:32}

def lin_table(lst, dtype = uint32, t0 = 0):
    """return a linear table t with t[0]=t0, t[2**i]=lst(i).

    For all entries i,j of t we have t[i^j] ^ t[0] = t[i] ^ t[i].
    t[0] is set to the value t0 (default is 0).
    
    The table ist returned as a numpy array of type 'dtype',
    with dtype = uint8, unit16 or uint 32.
    """
    if not 0 <= (reduce(__or__,lst,0) | t0) <= 1 << _maxbits[dtype]:
        print( hex(reduce(__or__,lst,0) | t0) )
        raise ValueError("entries too large for table")
    a = numpy.zeros(1 << len(lst), dtype = dtype)
    a[0] = t0
    for i, x in enumerate(lst):
        a[1<<i:2<<i].fill(x)
        a[1<<i:2<<i] ^= a[0:1<<i]
    return a


def bit_mat_iszero(a):
    """return 0 if a is a zero matrix"""
    if isinteger(a):
         return a == 0
    return sum(lmap(abs,a)) == 0

def bit_mat_mul(a, *args):
    """ Multiplies 'a' with an arbitrary number of bit matrices.

    Parameter 'a' may be a bit vector or a list representing a bit matrix.
    Representation of bit vectors and matrices see comment on module.
    The result returned has the same shape as parameter 'a'. 
    """
    if len(args) == 1:
        b = args[0]
        if isinteger(a):
            return reduce(
                   __xor__, [x for i,x in enumerate(b) if a & (1<<i)], 0 )
        else:
            return [ reduce(
                  __xor__, [x for i,x in enumerate(b) if y & (1<<i)], 0 )
               for y in a ]
    if len(args) == 0:
        return a
    return  bit_mat_mul(bit_mat_mul(a, args[0]), *args[1:])



def bit_mat_transpose(a, ncolumns=None):
    """return the  transpose of bit Matrix a.
 
    Optionally, the number of columns of a may be given.

    """
    ncols = bitlen(reduce(__or__, a, 0))
    if ncolumns is None: 
        ncolumns = ncols 
    elif ncols > ncolumns:
         raise ValueError( "Bad entry in bit matrix" )
    t = [0] * ncolumns
    for i, l in enumerate(a):
        for j in range(ncolumns):
            if (l >> j) & 1:
                 t[j] ^= 1 << i
    return t





def bit_mat_inverse(a):
    """return the inverse of bit Matrix a."""
    msg = "Bit matrix is not invertible"
    ncols = bitlen(reduce(__or__, a, 0))
    if ncols != len(a):
         raise ZeroDivisionError( msg )
    perm = [None] *  ncols
    hicol = 1 << ncols
    ah = [ int(x) + (hicol << i) for i,x in enumerate(a) ]
    for i in range(ncols):
        piv = v2(ah[i])
        if piv >= ncols:
            raise ZeroDivisionError( msg )
        perm[piv] = i
        msk = 1 << piv
        for j in range(ncols):
            if ah[j] & msk and i != j:
                ah[j] ^= ah[i]
    return [ ah[perm[i]] >> ncols  for i in range(ncols) ]


 

def bit_mat_det(a):
    """return the determinant of bit Matrix a."""
    try:
        bit_mat_inverse(a)
        return 1
    except ZeroDivisionError:
        return 0


def bit_mat_orthogonal_complement(a, ncolumns=None):
    """return an orthognal complement of bit Matrix a.

    Optionally, the number of columns of a may be given.
    Here matrix 'a' is considered a matrix of row vectors generating a
    linear subspace W of V = GF(2)**ncolumns. The function returns a
    matrix of row vetors generating the orthognal complement of W with
    respect to the standard Euclidean inner product in V.      
    """
    # This can be improved using the ideas in bit_mat_inverse !!
    if ncolumns == None:
         ncolumns = bitlen(reduce(__or__, a, 0))
    piv_a, piv_cols = pivot_binary_low(a)
    for i in range(ncolumns):
        if not i in piv_cols:
            piv_a.append(1 << i)
    ao = bit_mat_transpose(bit_mat_inverse(piv_a), ncolumns)
    return ao[len(piv_cols):]


def bit_mat_rank(a):
    """return the rank of a bit marix a"""
    return len(pivot_binary_high(a)[0])


def bit_mat_random(nrows, ncolumns):
    """return a random bit matrix with given number of rows and columns"""
    a, rmax = [],  ( 1 << ncolumns) - 1
    for i in range(nrows):
         a.append(randint(0, rmax))
    return a


def bit_mat_basis_spanned(basis, dtype = None):
    """return an array of the bitvectors spanned by a basis

    The function returns an array A of length 2**len(basis) containing
    the bit vectors spanned by the basis.

    Entry A[sum(b_i * 2*i)] is sum(b_i * basis[i]), b[i] = 0,1.

    If dtype is not None, a numpy array of the corresponding dtype
    is returned.
    """
    sp = [0]
    for x in basis:
        sp += [ y ^ x for y in sp ] 
    if not dtype is None:
        sp = numpy.array(sp, dtype = dtype)
    return sp



      

                        

def bit_perm_mul(a, *args):
    """ Multiplies 'a' with an arbitrary number of bit permutations.

    Parameter 'a'  may be a bit vector or a list representing a permutation.
    Representation of bit vectors and permutations see comment on module.
    The result returned has the same shape as parameter 'a'. 
    """
    if len(args) == 1:
        b = args[0]
        if isinteger(a):
            return reduce(
                   __xor__, [1<<x for i,x in enumerate(b) if a & (1<<i)], 0 )
        else:
            return [ b[x]  for  x in a ]
    if len(args) == 0:
        return a
    return  bit_perm_mul(a, *(args[:-2] + ([args[-1][x] for x in args[-2] ],)))  
 
   
    raise NotImplementedError("")


def bit_perm_2_mat(a):
    """transforms a permutation into a bit matrix.

    Representation of permutations and bit matrices see comment 
    on module.  
    """
    if isinteger(a): return a
    return [1 << x for x in a]



def rand_perm(n):
    """returns a list representing a random permutation of n elements"""
    if n == 1: return [0]
    else:
       a = randint(0, n-1)
       p1 = [x for x in range(a,n)+range(a)]
       p2 = [n-1] + rand_perm(n-1)
       return bit_perm_mul(p1,p2)
        


def unnumpy(obj):
    """Change the stange numpy scalars in an object to python integers

    This change is done recursively in tuples or lists.
    """
    if isinstance(obj, list):
         return [unnumpy(x) for x in obj]
    if isinstance(obj, tuple):
         return tuple([unnumpy(x) for x in obj])
    if isinstance(obj, Integral):
         return int(obj)
    return obj



try:
    from math import comb as binomial
except:
    def binomial(n, k):
        """Return binomial coefficient ``n`` over ``k``"""
        assert isinstance(k, Integral) 
        assert isinstance(n, Integral) and n >= 0
        if k + k > n:
            k = n - k
        if k > 0:
            x = 1
            for i in range(1, k+1):
                x = (x * (n + 1 - i)) // i
            return x
        return 1 if k == 0 else 0


