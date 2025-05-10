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
import re
import numpy as np
from numbers import Integral
import warnings
from collections import defaultdict, OrderedDict
import math 
from functools import partial




from mmgroup.structures.abstract_group import singleton
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
from mmgroup.structures.mm_space_indices import tuple_to_sparse
from mmgroup.structures.mm_space_indices import numeric_index_to_sparse
from mmgroup.structures.mm_space_indices import sparse_from_indices
from mmgroup.structures.abstract_group import AbstractGroupWord
from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord
from mmgroup.mm_space  import MMSpace, MMVector
from mmgroup.mm_space  import standard_mm_group
from mmgroup.structures.parse_atoms import AtomDict 
from mmgroup.structures.parse_atoms import eval_atom_expression, ihex 

from mmgroup.mm_op import mm_aux_mmv_size
from mmgroup.mm_op import mm_vector, mm_aux_random_mmv
from mmgroup.mm_op import mm_aux_zero_mmv, mm_aux_reduce_mmv
from mmgroup.mm_op import mm_aux_mmv_to_sparse
from mmgroup.mm_op import mm_aux_mmv_set_sparse
from mmgroup.mm_op import mm_crt_combine, mm_crt_check_v2
from mmgroup.mm_op import mm_crt_combine_bytes
from mmgroup.mm_op import mm_crt_check_g
from mmgroup.mm_op import mm_crt_norm_int32
from mmgroup.generators import mm_group_n_clear
from mmgroup.generators import mm_group_n_mul_word_scan
from mmgroup.generators import mm_group_n_to_word
from mmgroup.mm_op import mm_op_word
from mmgroup.mm_op import mm_op_compare



PRECISION = math.log(7 * 31 * 127 * 255) / math.log(2.0) - 4



######################################################################
# Auxiliary class vsparse representing a sparse vector
######################################################################


ERR_CRT_TYPE = "Connot construct MMVectorCRT object from type '%s'"


class vsparse:
    def __init__(self, *data):
        self.d = defaultdict(int)
        if len(data) == 0 or not data[0]:
            return
        scalar = 1
        if isinstance(data[0], vsparse):
            self.d.update(data[0].d)
            return
        if isinstance(data[0], Integral):
            scalar, data = int(data[0]),  data[1:]
        if scalar == 0 or len(data) == 0:
            return
        if isinstance(data[0], str):
            if len(data[0]) != 1 or not data[0] in "ABCTXZYIDE0":
                err = "Illegal tag '%s' in tuple for class MMSpaceCRT"
                raise ValueError(err %  data[0]) 
        else:
            raise TypeError(ERR_CRT_TYPE % type(data[0]))          
        for t in tuple_to_sparse(255, *data):
            t = int(t)
            scalar1, t =  t & 0xff, t & 0xffffff00
            scalar1 = scalar1 if scalar1 < 128 else scalar1 - 255
            self.d[t] += scalar * scalar1

    def __imul__(self, other):
        assert isinstance(other, Integral)
        if other == 0:
            self.d.clear()
        else:
            for tag in self.d.keys():
                 self.d[tag] *= other
        return self

    def __mul__(self, other):
        return vsparse(self).__imul__(other)

    __rmul__ = __mul__

    def __neg__(self):
        return self.__mul__(-1) 

    def __pos__(self):
        return self  
 
    def __iadd__(self, other):
        if other == 0:
            return
        assert isinstance(other, vsparse)
        for tag, value in other.d.items():
            self.d[tag] += value
        return self

    def __add__(self, other):
        return vsparse(self).__iadd__(other)

    def __isub__(self, other):
        if other == 0:
            return
        assert isinstance(other, vsparse)
        for tag, value in other.d.items():
            self.d[tag] -= value
        return self

    def __sub__(self, other):
        return vsparse(self).__isub__(other)

    def reduce(self):
        for tag, value in self.d.items():
            if value == 0:
                del self.d[tag]
        return self
        
    def norm(self):
        norm = 0
        self.reduce()
        for tag, value in self.d.items():
            factor = 1
            if tag & 0xE000000 == 0x2000000:
                i0, i1 = (tag >> 14) & 0x7ff, (tag >> 8) & 0x3f
                factor += i0 != i1
            norm += factor * value * value
        return norm    

    def sparse_array(self, p, shift = 0):
        assert p & 1 and p < 256
        self.reduce()
        a = np.zeros(len(self.d), dtype = np.uint32)
        for i, (tag, value) in enumerate(self.d.items()):
            a[i] = tag + (value << shift) % p 
        return a    
  

######################################################################
# Convert string to instance of class vsparse
######################################################################
 
FRAME = re.compile(r"^([A-Za-z_])+\<([0-9]+;)?(.+)\>$") 


def vsparse_from_str(s):
    string = s
    m = FRAME.match(s)
    if m:
        _, p_str, string = m[1], m[2], m[3]
        if p_str:
            err = "Vector is defined module an integer only"
            raise ValueError(err)
    f = AtomDict(vsparse)
    return eval_atom_expression(string, f)
    
 

######################################################################
# Auxiliary functions
######################################################################



def _err_tag(*args, **kwds):
    err = "Bad entry in monster group element"
    raise ValueError(err)



def _iter_group(g):
    start, data, len_ = 0, g.mmdata, len(g.mmdata)
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


def compress_data(data):
    out = np.zeros(196884, dtype = np.int32)
    out[0:24] = data[0:768:33]
    k = 24
    for i in range(1,24):
        out[k:k+i] = data[32*i:32*i+i]
        out[276+k:276+k+i] = data[768+32*i:768+32*i+i]
        out[552+k:552+k+i] = data[1536+32*i:1536+32*i+i]
        k += i
    assert k == 300
    out[852:49428] = data[2304:50880]
    t = data[50880:247488].reshape((3*2048,32))[:,:24]
    out[49428:196884] = t.ravel()
    assert len(out) == 196884
    return out
    

######################################################################
# If check_MMVectorCRT is True then we check a vector of type
# MMVectorCRT for overflow or underflow after each operation
######################################################################

check_MMVectorCRT = True 

######################################################################
# Modelling a vector of the 196884-dimensional rep of the monster
######################################################################

MODULUS = 7*31*125*255
MAX_CRT_NORM = MODULUS**2 // 4
ERR_OVERFLOW = "Overflow in class MMVectorCRT"
ERR_UNDERFLOW = "Underflow in class MMVectorCRT"
ERR_UNDERFLOW_G = "Underflow at group operation in class MMVectorCRT"



class MMVectorCRT(AbstractMmRepVector):
    """Models a vector in a space of type ``MMSpaceCRT``.

    Such a vector should be constructed by calling an instance ``V``
    of class ``MMSpaceCRT`` which models a real representation of
    the monster group. Calculations in this space are exact in
    fixed-point arithmetic with a precision of about %.2f bits.    

    ValueError is raised in case of overflow or underflow.
    
    The functionality of this class is a subset of the functionality 
    of class ``MMVector``. See class ``MMSpaceCRT`` for details.

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
    p = MODULUS
    MAX_CRT_NORM = MAX_CRT_NORM

    def __init__(self, shift, tag = 0, i0 = None, i1 = None):
        self.shift = shift
        self.factor = 1.0 / (1 << shift)
        if not 3 <= self.shift <= 21:
            raise ValueError("Bad shift factor for class MMVectorCRT") 
        self.data = OrderedDict()
        self.data_int = np.zeros(247488, dtype = np.int32)
        self.expanded = False
        d = vsparse()
        if isinstance(tag, MMVectorCRT): 
            for p in (7, 31, 127, 255):
                self.data[p] = tag.data[p].copy()
            self.shl(self.shift - tag.shift)
            return
        elif isinstance(tag, Integral) and not tag:
            for p in (7, 31, 127, 255):
                self.data[p] = MMVector(p)
            return
        elif isinstance(tag, str) and len(tag) == 1:
            d += vsparse(tag, i0, i1)
        elif isinstance(tag, str) and tag == "Axis":
            g = None
            if isinstance(i0, AbstractMMGroupWord):
               i0, g = MMSpace._mm_element_to_axis(i0)
               i1 = None
            for p in (7, 31, 127, 255):
                self.data[p]  = MMVector(p, tag, i0, i1) << self.shift
            if g is not None:
                self *= g
            return
        elif isinstance(tag, str):
            d += vsparse_from_str(tag)
        elif isinstance(tag, list):
            for x in tag:
                if isinstance(x,tuple):
                    d += vsparse(*x)
                elif isinstance(x, str):
                    d += vsparse_from_str(x)
        else:
            err = "Connot construct MMVectorCRT object from type '%s'"
            raise ValueError(err % type(tag))
        for p in (7, 31, 127, 255):
            self.data[p] = v = MMVector(p)
            ind = d.sparse_array(p, shift)
            mm_aux_mmv_set_sparse(p, v.data, ind, len(ind)) 
        if check_MMVectorCRT:
            self.expand()
            if self._inorm > MAX_CRT_NORM:
                raise ValueError(ERR_OVERFLOW)

    def expand(self):
        """Expand data of vector ``v`` with CRT

        In the array ``v.data_int`` the ``i``-th entry is computed 
        form the ``i``-th entries of the vector modulo 7, 31, 127,
        and 255. 
        """
        if not self.expanded:
            v2 = mm_crt_combine(self.data[7].data, 
                self.data[31].data, self.data[127].data, 
                self.data[255].data, self.data_int)  
            self._v2 = v2 - self.shift if v2 < 24 else 24        
            self._inorm = mm_crt_norm_int32(self.data_int)
            self.expanded = True


    def shl(self, sh):
        assert isinstance(sh, Integral)
        if sh == 0:
            return self
        if sh > 0:
            if check_MMVectorCRT:
                self.expand()
                if sh > 24 or self._inorm * 4**sh > MAX_CRT_NORM:
                    raise ValueError(ERR_OVERFLOW)
                self.data_int <<= sh  
                self._inorm *= 4**sh            
            for d in self.data.values():
                d <<= sh
            self._v2 += sh            
        elif sh < 0:
            nsh = -sh
            if check_MMVectorCRT:
                self.expand()
                if nsh > self.v2:
                    raise ValueError(ERR_UNDERFLOW)
                self.data_int >>= nsh           
                self._inorm *= 4**sh            
            for d in self.data.values:
                d >>= nsh
            self._v2 += sh            
        return self     

    def check(self):
        """Check if the vector is correct

        Raise ValueError if the vector is erroneous.
        """
        return True


    def __ilshift__(self, other):
        return self.shl(other)

    def __lshift__(self, other):
        return self.copy().shl(other)

    def __irshift__(self, other):
        return self.shl(-other)

    def __lrhift__(self, other):
        return self.copy().shil(-other)


    def __mod__(self, p):
        """Return the vector modulo ``p``. 

        ``p``  must be in (7, 31, 127, 255). Actually, we divide 
        the vector by the *scaling factor* before reducing it 
        modulo ``p``.

        This method is mainly for testing.
        """ 
        if  p in (7, 31, 127, 255):
            return MMVector(p, self.data[p]) >> self.shift
        elif p in (3, 15):
            v0 = self % 255
            return MMVector(p, v0)
        elif isinstance(p, Integral):
            err = "Cannot reduce MMVectorCRT object modulo %d"
            raise ValueError(err % p)
        else:
            err = "Modulus for reducing MMVectorCRT object must be int"
            raise TypeError(err)

    @property
    def v2(self): 
        """Return the ``2``-adic value of the vector.
 
        If "a, b" are odd integers and ``k`` is an integer then
        the  ``2``-adic value  ``a * 2**k / b`` is ``k``.

        The  ``2``-adic value of a vector is the minimum of the
         ``2``-adic values of its entries, ignoring zero entries.

        The function raises ZeroDivisionError if ``v ==  0``
        """
        v2 = self._v2
        if self._v2 >= 24:
            err = "The zero vector has infinite 2-adic value"
            raise ZeroDivisionError(err)
        return self._v2
 
    @property
    def inorm(self):
        """Return a scaled norm of the vector as an integer

        For a vector ``v`` we have 

            ``v.fnorm() = v.inorm() * v.factor**2``.

        Where ``v.fnorm()`` is the real norm of ``v``.
        """
        self.expand()
        return self._inorm



    def norm(self):
        """Return norm of vector as a floating point number.

        The norm of a vector in the representation of the monster
        is the squared sum of is entries. Here the squares of all
        entries with index ``("A", i0, i1)`` must be doubled in
        case ``i0 != i1``.
 
        The returned norm is exact.
        """
        return self.inorm * self.factor**2

######################################################################
# class MMSpace
######################################################################


@singleton
class MMSpaceCRT(AbstractMmRepSpace):
    """Models a ``196884``-dimensional representation of the monster group 

    This class models a real representation of the monster group 
    with fixed-point arithmetic. Calculations are done by combining
    vectors modulo ``p = 7, 31, 127, 255`` with Chinese remaindering.
    This way we achieve about %.2f bit precision.

    The construction of a vector in this space and the computation
    with such vectors works in the same way as in class |MMSpace|.
    But there are som limitations:

      * Vectors may be constructed as in class |MMSpace|, but
        the arguments of the constructor may be tuples only. 
        Here randomized scalars are illegal.

      * The only operations allowed for vectors are copying, 
        multiplication with a group element, and testing for 
        equality. Vector addition and scalar multiplication
        are illegal.

      * A vector may be reduced modulo one of the primes
        ``p = 7, 31, 127, 255`` using the modulo operator ``%%``.
        The result is a vector in the space ``MMSpace(p)``.

      * Changing entries of a vector via item assignment is
        illegal. 

      * Entries and subarrays of a vector may be obtained as
        in class |MMSpace|. Here subarrays are returned as
        integers or numpy  arrays with ``dtype = np.int32``.
        On output, each number is multiplied with ``2**k``
        in order to obtain an integer. Here ``k`` is the
        argument given in the constructor by parameter ``shift``.
        
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
    vector_type = MMVectorCRT
    _max_norm = (7*31*125*255)**2 // 4

    def __init__(self):
        """Create a 196884-dimensional representation of the monster

        All calculations are done modulo the odd number p
        """
        super(MMSpaceCRT, self).__init__()

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

    def zero(self, shift):
        """Return the zero vector"""
        err = "Cannot create zero vector in space of class MMSpaceCRT"
        return MMVectorCRT(shift)

    def copy_vector(self, v1):
        assert v1.space == self
        v = MMVectorCRT(v1.shift, 0)
        for p in (7, 31, 127, 255):
            np.copyto(v.data[p].data, v1.data[p].data)
        if v1.expanded:
            np.copyto(v.data_int, v1.data_int)
            v._v2 = v1._v2        
            v._inorm = v1._inorm
        v.expanded = v1.expanded
        return v

      
    def __call__(self, *args):
       return self.from_tuples(*args)      


    def set_rand_uniform(self, *args, **kwds):
        err = "Cannot create random vector in space of class MMSpaceCRT"
        raise NotImplementedError(err)


    parse = _not_supported

    #######################################################################
    # Obtaining and setting components via sparse vectors
    #######################################################################


    def getitems_sparse(self, v, a_sparse):
        raise NotImplementedError(err)

    def additems_sparse(self, *args, **kwds):
        err = "Sparse representation not supported in class MMSpaceCRT"
        raise NotImplementedError(err)

    setitems_sparse = additems_sparse


    #######################################################################
    # Conversion from and to to sparse representation 
    #######################################################################

    as_sparse = additems_sparse


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

    def _imul_word(self, v1, g_word, buf):
        if len(g_word):
            vectors = list(v1.data.values())
            if check_MMVectorCRT:
                if mm_crt_check_g(g_word[0], *[v.data for v in vectors]):
                    print("MMVectorCRT: shift = %d, v2 = %d, tag = %s" %
                        (v1.shift, v1.v2, ((g_word[0] >> 28) & 7)))
                    raise ValueError(ERR_UNDERFLOW_G)
            for v in vectors:
                mm_op_word(v.p, v.data, g_word, len(g_word), 1, buf)
        

    def imul_group_word(self, v1, g):
        """Return product v1 * g of vector v1 and group word g.

        v1 is replaced by v1 * g.

        This method is called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.
        """
        assert isinstance(g, AbstractGroupWord) and g.group.is_mmgroup 
        a = g.mmdata
        nn = np.zeros(5, dtype = np.uint32)
        nnw = np.zeros(5, dtype = np.uint32)
        buf = np.zeros(mm_aux_mmv_size(255), dtype = np.uint64)
        while len(a):
            mm_group_n_clear(nn)
            i = mm_group_n_mul_word_scan(nn, a, len(a))
            length =  mm_group_n_to_word(nn, nnw)
            self._imul_word(v1, nnw[:length], buf)
            a = a[i:]
            if len(a):
                self._imul_word(v1, a[:1], buf)    
                a = a[1:]
        del buf 
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
        v1.expand()
        v2.expand()
        v1_2, v2_2 = v1._v2, v2._v2
        if v1_2 == v2_2 == 24:
             return True   # then v1 = v1 = zero
        if v1_2 != v2_2:
             return False
        data1, data2 = v1.data_int, v2.data_int
        if v1.shift > v2.shift:
             data1 = data1 >> (v1.shift - v2.shift)
        if v2.shift > v1.shift:
             data2 = data2 >> (v2.shift - v1.shift)
        return (data1 == data2).all()

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

    
    def vector_get_item(self, v, item):
        assert v.space == self
        if not isinstance(item, tuple):
            item = (item,) 
        shape, a_sparse = sparse_from_indices(255, *item)
        d = {}
        for p in (7, 31, 127, 255):
            vdata = v.data[p]
            asp = a_sparse.copy()
            vdata.space.getitems_sparse(vdata, asp)        
            d[p] = (asp & p).astype(np.uint8)
        l = len(d[7])
        a = np.zeros(l, dtype = np.int32)
        mm_crt_combine_bytes(d[7], d[31], d[127], d[255], l, a)
        a = np.array(a, dtype = float)
        af = np.array(a * v.factor,  dtype = float)
        return af.reshape(shape) if len(shape) else float(af[0])
        


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




StdMMSpaceCRT = MMSpaceCRT()
MMVectorCRT.space = StdMMSpaceCRT


def MMV_CRT(shift):
    return partial(MMVectorCRT, shift)
