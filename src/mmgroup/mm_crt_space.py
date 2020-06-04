r"""We deal with the 196884-dimensional representation of the monster.

Here class ``MMSpaceCRT`` is an analogue of class ``MMSpace`` that
allows a limited set of operations on the real 
``196884``-dimensional representation of the monster with vectors of
a small norm. Typical examples of such vectors are unit vectors, or
so called *axes*, as described in :cite:`Con85`.

"""
# References in the __docstr__ see file docs/source/references.bib


from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from numbers import Integral
import warnings
from collections import defaultdict
import math 




from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
from mmgroup.structures.mm_space_indices import tuple_to_sparse
from mmgroup.structures.mm_space_indices import numeric_index_to_sparse
from mmgroup.mm_group  import MMGroup, MMGroupWord
from mmgroup.mm_space  import MMSpace, MMSpaceVector, mm_wrapper
from mmgroup.mm_space  import standard_mm_group

from mmgroup.mm import mm_vector, mm_aux_random_mmv
from mmgroup.mm import mm_aux_zero_mmv, mm_aux_reduce_mmv
from mmgroup.mm import mm_aux_mmv_to_sparse
from mmgroup.mm import mm_aux_mmv_set_sparse
from mmgroup.mm import mm_crt_combine, mm_crt_check_v2
from mmgroup.mm import mm_crt_check_g
from mmgroup.mm import mm_crt_norm_int32


_mm_op, _mm_compare = {},  {}
for p in (7, 31, 127, 255):
    _mm_op[p] = mm_wrapper(p).op_word
    _mm_compare[p] = mm_wrapper(p).op_compare

PRECISION = math.log(7 * 31 * 127 * 255) / math.log(2.0) - 4

_index_dict = {
    "A": (      0,   24,  32,  24 ),
    "B": (    768,   24,  32,  24 ),
    "C": (   1536,   24,  32,  24 ),
    "T": (   2304,  759,  64,  64 ),
    "X": (  50880, 2048,  32,  24 ),
    "Z": ( 116416, 2048,  32,  24 ),
    "Y": ( 181952, 2048,  32,  24 ),
}

######################################################################
# Auxiliary functions
######################################################################

def _tuples_to_sparse_dict(*tuples):
    d = defaultdict(int)
    for data in tuples:
        scalar = 1
        if isinstance(data[0], Integral):
            scalar, data = data[0],  data[1:]
        if not data[0] in "ABCTXZYIDE":
            err = "Illegal tag %s in tuple for class MMSpaceCRT"
            raise err % tag            
        for t in tuple_to_sparse(255, *data):
            scalar1, t =  t & 0xff, t & 0xffffff00
            scalar1 = scalar1 if scalar1 < 128 else scalar1 - 255
            d[t] += scalar * scalar1
    return d


def _norm_sparse_dict(d):
    norm = 0
    for i, scalar in d.items():
        factor = 1
        if i & 0xE000000 == 0x2000000:
            i0, i1 = (i >> 14) & 0x7ff, (i >> 8) & 0x3f
            factor += i0 != i1
        norm += factor * scalar * scalar
    return norm
    

def _err_tag(*args, **kwds):
    err = "Bad entry in monster group element"
    raise ValueError(err)



def _iter_group(g):
    start, data, len_ = 0, g.data, len(g.data)
    while start < len_:
        tag = data[start] & 0x70000000
        if  tag >= 0x50000000:
            yield  data[start : start+1], True
            start += 1
        else:
            end_ = start
            while end_ < len_ and data[end_] & 0x70000000 < 0x50000000:
               end_ += 1
            yield data[start : end_], False
            start = end_ 



######################################################################
# Modelling a vector of the 196884-dimensional rep of the monster
######################################################################

class MMSpaceVectorCRT(AbstractMmRepVector):
    """Models a vector in a space of type ``MMSpaceCRT``.

    Such a vector should be constructed by calling an instance ``V``
    of class ``MMSpaceCRT`` which models a real representation of
    the monster group. Calculations in this space are exact in
    fixed-point arithmetic with a precsion of about %.2f bits.    

    ValueError is raised in case of overflow or underflow.
    
    The functionality of this class is a subset of the functionality 
    of class ``MMSpaceVector``. See class ``MMSpaceCRT`` for details.

    A vector may also be reduced modulo ``p = 7, 31, 127, 255`` with
    the modulo operator ``%%``. Then the result is a vector in the 
    vector space ``MMSpace(p)``.
     
    :var space:
        This attribute contains the space to which the vector
        belongs. That space is an instance of class |MMSpaceCRT|.

    .. warning::
       The constructor of this class is not for public use! You
       may call an instance ``V`` of class  |MMSpaceCRT| for
       constructing vectors in the real representation space
       ``V`` of the monster group.
    """ % PRECISION

    def __init__(self, space):
        self.space = space
        self.data = {}
        for p in (7, 31, 127, 255):
            self.data[p] = mm_vector(p)
        self.data_int = np.zeros(247488, dtype = np.int32)
        self.shift = space.shift
        self.expanded = False

    def check(self):
        """Check if the vector is correct

        Raise ValueError if the vector is errorneous.
        """
        return True


    def expand(self):
        """Expand data of vector ``v`` with CRT

        In the array ``v.data_int`` the ``i``-th entry is computed 
        form the ``i``-th entries of the vector modulo 7, 31, 127,
        and 255. 
        """
        if not self.expanded:
            mm_crt_combine(self.data[7], self.data[31], 
                self.data[127], self.data[255], self.data_int)           
            self.expanded = True

    def __mod__(self, p):
        """Return the vector modulo ``p``. 

        ``p``  must be in (7, 31, 127, 255). Actually, we divide 
        the vector by the *scaling factor* befor reducing it 
        modulo ``p``.

        This methos is mainly for tesing.
        """ 
        assert p in (7, 31, 127, 255)
        v = MMSpace(p)()
        np.copyto(v.data, self.data[p]) 
        return v

    def v2(self): 
        """Return the ``2``-adic value of the vector.
 
        If "a, b" are odd integers and ``k`` is an integer then
        the  ``2``-adic value  ``a * 2**k / b`` is ``k``.

        The  ``2``-adic value of a vector is the minimum of the
         ``2``-adic values of its entries, ignoring zero entries.

        The function raises ZeroDivisionError if ``v ==  0``
        """
        v2 = mm_crt_check_v2(self.data[7], self.data[31], 
                self.data[127], self.data[255])
        if v2 >= 24:
            err = "The zero vector has infinite 2-adic value"
            raise ZeroDivisionError(err)
        return v2 - self.space.shift
 
    def inorm(self):
        """Return a scaled norm of the vector as an integer

        For a vector ``v`` we have 

            ``v.fnorm() = v.inorm() * v.factor**2``.

        Where ``v.fnorm()`` is the real norm of ``v``.
        """
        self.expand()
        return mm_crt_norm_int32(self.data_int)

    def fnorm(self):
        """Return norm of vector as a floating point number.

        The norm of a vector in the representation of the monster
        is the squared sum of is entries. Here the squares of all
        entries with index ``("A", i0, i1)`` must be doubled in
        case ``i0 != i1``.
 
        The returned norm is exact.
        """
        return self.norm() * self.factor**2

    @property
    def factor(self):
        """This is the *scaling factor* of the vector

        Since floating-point arithmetic is imprecise, entries of the 
        vector are returned as (signed 32-bit) integers. Each integer
        entry of the vector must be multiplied with the 
        *scaling factor*, which is always a negative power of two.
        """
        return 1.0 / (1 << self.shift)


######################################################################
# class MMSpace
######################################################################



class MMSpaceCRT(MMSpace):
    """Models a ``196884``-dimensional representation of the monster group 

    This class models a real representation of the monster group 
    with fixed-point arithmetic. Calcuations are done by combining
    vectors modulo ``p = 7, 31, 127, 255`` with Chinese remaindering.
    This way we achive about %.2f bit precision.

    The construction of a vector in this space and the computation
    with such vectors works in the same way as in class |MMSpace|.
    But there are som limitaions:

      * Vectors may be constructed as in class |MMSpace|, but
        the arguments of the constructor may be tuples only. 
        Here randomized scalars are illegal.

      * The only operations allowed for vectors are copying, 
        multiplication with a group element, and testing for 
        equality. Vector addition and scalar multiplication
        are illegal.

      * A vector may be reduced modulo one of the primes
        ``p = 7, 31, 127, 255`` using the modulo operatior ``%%``.
        The result is a vector in the space ``MMSpace(p)``.

      * Changing entries of a vector via item assignment is
        illegal. 

      * Entries and subarrays of a vector may be obtained as
        in class |MMSpace|. Here subarrays are returned as
        integers or numpy  arrays with ``dtype = np.int32``.
        On output, each number is multiplied with ``2**k``
        in order to obtain an integer. Here ``k`` is the
        argument given in the cnstructor by parameter ``shift``.
        
        You may multiply such an integer value with the attribute
        ``factor`` of the space for obtaining an exact floating
        point value. But floating point data may lead to imprecise
        results when processing the with a computer algebra system.

      
    :var shift:
        If shift is ``k`` then a vector ``v`` can be represented if 
        all entries of ``v`` are integer multiples of ``2**(-k)``.  
        The absolute value of an entry may be at most 
        ``3513772 / 2**k``. In case ``shift == 20`` (default)
        this is about  ``3.35``.
 
    :var group:
        Contains the group operating on the vector
        space by right multiplication. 
        That group must be an instance of class |MMGroup|.

    Caution!
    
        The norm of an input vector must not exceed 
        ``1.234 * 10**13 * 4.0 ** (-k)``. Otherwise *overflow*
        occurs. In case ``k = 20`` (default) this quantity is about 
        ``11.229``. Here ``k`` is given by parameter ``shift``. 

        The norm of a vector is the sum of the squares of its 
        entries, where for all entries with index ``("A", i0, i1)`` ,
        ``i0 != i1``, the doubled square of the entry must be 
        taken instead.

    Caution!

        When multiplying a vector ``v`` with a generator of the group 
        with tag ``t`` or ``l``, and ``v`` has an entry that is not 
        an integer multiple of ``8 * 2**(-k)``, then a 
        *precision error* may occur. We raise ``ValueError`` in 
        case of a precision error.
      

    """ % PRECISION
    vector_type = MMSpaceVectorCRT
    _max_norm = (7*31*125*255)**2 // 4

    def __init__(self, shift = 20, group = None):
        """Create a 196884-dimensional representation of the monster

        All calculations are done modulo the odd number p
        """
        if not 3 <= shift <= 21:
            raise ValueError("Bad shift factor for class MMSpaceCRT") 
        self.shift = shift
        self.factor_ = 1.0 / (1 << shift)
        if group is None:
            group = standard_mm_group 
        assert isinstance(group, MMGroup) 
        super(MMSpaceCRT, self).__init__(255, group)

    #######################################################################
    # Conversion to a list of tuples 
    #######################################################################

    def _not_supported(self, *args, **kwds):
        err = "Method not supported in class MMSpaceCRT"
        raise NotImplementedError(err)

    as_tuples = _not_supported



    #######################################################################
    # Creating vectors 
    #######################################################################

    def zero(self):
        """Return the zero vector"""
        err = "Cannot create zero vector in space of class MMSpaceCRT"
        return MMSpaceVectorCRT(self)

    def copy_vector(self, v1):
        assert v1.space == self
        v = MMSpaceVectorCRT(self)
        for p in (7, 31, 127, 255):
            np.copyto(v.data[p], v1.data[p])
        if v1.expanded:
            np.copyto(v.data_int, v1.data_int)
        v.expanded = v1.expanded
        return v

    def from_tuples(self, *tuples):
        d = _tuples_to_sparse_dict(*tuples)
        norm = (4 << self.shift) * _norm_sparse_dict(d)
        if norm > self._max_norm:
            err = "Overflow in vector in class MMSpaceCRT"
        v =  MMSpaceVectorCRT(self)
        sh = self.shift
        for p in (7, 31, 127, 255):
            ind = [index + (val << sh) % p for index, val in d.items()]
            ind = np.array(ind, dtype = np.uint32)
            mm_aux_mmv_set_sparse(p, v.data[p], ind, len(ind)) 
        v.expanded = False
        v.shift = self.shift 
        return v        
      
    def __call__(self, *tuples):
       return self.from_tuples(*tuples)      

    def unit(self,  *args):
        """Return a unit vector of the vector space

        Constructing a unit vector without any arguments should
        return the zero vector.
        """
        return self.from_tuples(args)      


    def rand_uniform(self, *args, **kwds):
        err = "Cannot create random vector in space of class MMSpaceCRT"
        raise NotImplementedError(err)

    rand_vector = rand_uniform

    parse = _not_supported

    #######################################################################
    # Obtaining and setting components via sparse vectors
    #######################################################################

    def getitems_sparse(self, *args, **kwds):
        err = "Sparse representation not supported in class MMSpaceCRT"
        raise NotImplementedError(err)

    additems_sparse = getitems_sparse

    setitems_sparse = getitems_sparse


    #######################################################################
    # Conversion from and to to sparse representation 
    #######################################################################

    as_sparse = getitems_sparse


    #######################################################################
    # Vector operations 
    #######################################################################


    def iadd(self, v1, v2):
        err = "Vector addition not supported in class MMSpaceCRT"
        raise NotImplementedError(err)
 
    def imul_scalar(self, v1, a):
        err = "Scalar multiplication not supported in class MMSpaceCRT"
        raise NotImplementedError(err)
           
    #######################################################################
    # Group operation 
    #######################################################################

    def imul_group_word(self, v1, g):
        """Return product v1 * g of vector v1 and group word g.

        v1 may be destroyed.

        This method is called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.
        """
        work =  mm_vector(255)
        assert v1.space == self
        g = self.group(g)
        assert isinstance(g, MMGroupWord)
        g.reduce() 
        for g_data, check in _iter_group(g):
            if (check and mm_crt_check_g(g_data, v1.data[7], 
                    v1.data[31], v1.data[127], v1.data[255])):
                err = "Underflow at group operation in class MMSpaceCRT"
                raise ValueError(err)
            for p in (7, 31, 127, 255):
                _mm_op[p](v1.data[p], g_data, len(g_data), 1, work)
        v1.expanded = False
        return v1  



    vector_mul_exp = _not_supported

    #######################################################################
    # Checking equality
    #######################################################################

    def equal_vectors(self, v1, v2):
        """Return True iff vectors v1 and v2 are equal 

        This method is called for elements v1 and v2 of the space
        'self' only.
        """
        diff = 0
        for p in (7, 31, 127, 255):
            diff |= _mm_compare[p](v1.data[p], v2.data[p])
        return not diff

    #######################################################################
    # Conversion from and to byte format
    #######################################################################

    as_bytes =  _not_supported

    from_bytes = _not_supported
        
    #######################################################################
    #  Checking and reducing a vector
    #######################################################################

    def check(self, v1):
        return True
 
    def reduce(self, v1):
        return v1

    #######################################################################
    # getitem and setitem
    #######################################################################

    
    def vector_get_item(self, v, index):
        v.expand()
        if not isinstance(index, tuple):
            index = (index,)
        try:
            start, max0, stride0, max1 = _index_dict[index[0]]
            a = v.data_int[start : start + max0 * stride0]
            a = np.reshape(a, (max0, stride0), order='C')[:,:max1]
            return a.__getitem__(index[1:])
        except KeyError:
            if index[0] == "D":
                a =  v.data_int[0:300:25]
                return a.__getitem__(index[1:])
            if index[0] == "E":
                shape, indices = numeric_index_to_sparse(
                    0, "E", *index[1:])
                a = a.__getitem__(indices)
                return a if len(shape) else int(a[0])
            err = "Bad index for vector in space of type MMSpacCRT"
            raise TypeError(err)


    def vector_set_item(*args, **kwd):
        err = "Item assigment not supported in space of type MMSpaceCRT"
        raise NotImplementedError(err) 
 

    #######################################################################
    # Formatting a vector 
    #######################################################################


    def str_vector(self, v1):
        return "<vector in space of type MMSpaceCRT>" 

 
    #######################################################################
    # Properties
    #######################################################################

    @property
    def factor(self):
        """Factor of type ``float`` for vector entries.

        When reading entries from a vector then integers or arrays
        of integers are returned. Here each entry must be multiplied
        with this ``factor`` to obtain an exact floating point value
        of that entry.
        """
        return self._factor       



    @property
    def max_norm(self):
        """Maximum feasible norm of an input vector.

        The norm of a vector is the sum of the squares of its 
        entries, where for all entries with index ``("A", i0, i1)`` ,
        ``i0 != i1``, the doubled square of the entry must be 
        taken instead.
        """
        return (self._max_norm * self._factor)**2        


