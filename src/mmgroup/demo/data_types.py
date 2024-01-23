"""Data structures for the reduction algorithm for the Monster group

The classes defined in this module are for demonstration only.

Do not impart any of these classes in a life application!
"""

from numbers import Integral
import numpy as np

from mmgroup import mat24, XLeech2
from mmgroup.structures.parse_atoms import ihex
from mmgroup.structures.construct_mm import iter_mm       
from mmgroup.structures.construct_mm import load_group_name     
from mmgroup.structures.abstract_group import singleton
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord
from mmgroup.structures.abstract_mm_group import AbstractMMGroup

from mmgroup.generators import gen_leech2_type
from mmgroup.mm_op import mm_vector, mm_op_word
from mmgroup.mm_op import mm_op_compare
from mmgroup.mm_op import mm_op_store_axis
from mmgroup.mm_op import mm_op_vector_add, mm_op_scalar_mul
from mmgroup.mm_reduce import mm_order_load_vector


######################################################################
# Vectors in Leech lattice mod 2
######################################################################

class Leech2:
    r"""Models an vector in the Leech lattice mod 2

    Such vectors may be added. A vector may be multiplied with an
    element of the subgroup :math:`G_{x0}` of the monster.

    Calling the constructor with the string 'Omega' or 'beta' returns
    the vector :math:`\lambda_\Omega` or :math:`\lambda_{\beta}`
    defined in :cite:`Seysen22`, Section 4.1 or 7.1, respectively.

    Calling the constructor with an integer argument returns the 
    element of the Leech lattice mod 2 with that number. Here the
    numbering is as in class :py:class:`~mmgroup.XLeech2`, setting
    the sign bit to zero.
    """
    STD_VECTORS = {'beta' : 0x200, 'Omega' : 0x800000}
    def __init__(self, value):
        if isinstance(value, Integral):
            self.value = int(value) & 0xffffff
        elif value in self.STD_VECTORS:
            self.value = self.STD_VECTORS[value]
        elif isinstance(value, (Leech2, XLeech2)):
            self.value = value.value & 0xffffff
        else:
            ERR = "Cannot construct class Leech2 object from input"
            raise ValueError(ERR)
    
    def __add__(self, other):
        assert isinstance(other, Leech2)
        return Leech2(self.value ^ other.value)

    def __mul__(self, other):
        assert isinstance(other, AbstractMMGroupWord)
        data = other.mmdata
        return Leech2(gen_leech2_op_word(self.value, data, len(data)))

    def __eq__(self, other):
        assert isinstance(other, Leech2)
        return (self.value ^ other.value) & 0xffffff == 0

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def type(self):
        r"""Return type of vector in Leech lattice mod 2.

        The type of a vector is the halved norm of a shortest
        preimage of the vector in the Leech lattice.
        That type is equal to 0, 2, 3, or 4.
        """
        return gen_leech2_type(self.value)

    def str(self):
        v = self.value
        x = (v >> 12) & 0xfff
        d = (mat24.ploop_theta(v >> 12) ^ v) & 0xfff
        return "Leech2<x_%s*d_%s>" % (ihex(x, 3), ihex(d, 3))
    __repr__  = str
    
  

######################################################################
# Elements of the Monster group
######################################################################



class Mm(AbstractMMGroupWord):
    r"""Models an element of the Monster group

    Here it would be tempting to use the class :py:class:`~mmgroup.MM`
    defined in the *API reference* of the **mmgroup** package for 
    computing in the Monster instead. Note that methods of class
    :py:class:`~mmgroup.MM` may apply the reduction algorithm whenever
    appropriate. Hence it would be cheating to demonstrate the reduction
    algorithm with instances of that class.  
    
    So we use this class **Mm** for demonstrating the reduction algorithm
    for the Monster.

    The constructor of the class **Mm** performs the same
    action as the constructor of class :py:class:`~mmgroup.MM`.
    The generators of the Monster used
    here are discussed in :cite:`Seysen22`, Section 3 - 5.

    But  multiplication of instances of class **Mm** is simply a
    concatenation of words, without any attempt to reduce the result.

    We use the following special cases of calling the constructor:

    Calling **Mm('r', k)** constructs a random word in the generators
    of the Monster containing **k** triality elements
    :math:`\tau^{\pm 1}`, for an integer **k**.

    Calling the constructor with the string **'negate_beta'** constructs
    an element of the normal subgroup :math:`Q_{x0}` of the subgroup
    :math:`G_{x0}` of the Monster that exchanges the axes :math:`v^+`
    and :math:`v^-` defined in :cite:`Seysen22`, Section 7.2.
    """
    __slots__ =  "length", "_data"
    ERR_ITER = "A monster group element g is not iterable. Use g.mmdata instead"
    MIN_LEN = 16
    def __init__(self,  tag = None, atom = None, *args, **kwds):
        self._data = np.zeros(self.MIN_LEN, dtype = np.uint32)
        if tag is None:
            self.length = 0
        elif tag == 'negate_beta':
            self._setdata(Mm('x', 0x200).mmdata) 
        else:
            atoms = iter_mm(self.group, tag, atom)
            self._setdata(np.fromiter(atoms, dtype = np.uint32)) 
                  
    def _extend(self, length):
        if length > len(self._data):
            length = max(length, 3*len(self._data) >> 1)
            self._data = np.resize(self._data, length)
              
    @property
    def mmdata(self):
        """Return the internal representation of the group element"""
        return np.copy(self._data[:self.length])

    @property
    def count_triality_elements(self):
        r"""Length of a word representing an element of the Monster

        The function returns the number of the triality elements
        :math:`\tau^{\pm 1}` contained in the word representing the
        group element. This will be used as an indication for the
        length of the word. This is similar to (although not equal to)
        the definition of the length of a word in :cite:`Wilson13`.
        """
        return sum([1 for x in self.mmdata if (x >> 28) & 7 == 5])

    def __bool__(self):
        return True
        
    def __len__(self):
        raise TypeError(self.ERR_ITER)

    def _setdata(self, data):
        self.length = len_ = len(data)
        self._extend(len_)
        self._data[:len_] = data
       
    def __getitem__(self,i):
        raise TypeError(self.ERR_ITER)
 
    def is_reduced(self):
        return True




@singleton
class MmGroup(AbstractMMGroup):
    r"""An instance of this class cintains the opration for class Mm
    """
    word_type = Mm
    group_name = "Mm"

    def __init__(self):
        super(MmGroup, self).__init__()

    def reduce(self, g1):
        return g1          
        
    def _imul(self, g1, g2):
        l1, l2 = g1.length, g2.length
        g1._extend(l1 + l2)
        g1._data[l1 : l1 + l2] = g2._data[:l2]
        g1.length = l1 + l2
        return g1

    def _invert(self, g1):
        w = self.word_type()
        w._setdata(np.flip(g1.mmdata) ^ 0x80000000)
        return w

    def copy_word(self, g1):
        result = self.word_type()
        result._setdata(g1.mmdata)
        return result

    def _equal_words(self, g1, g2):
        raise ValueError("Don't know if monster group elements are equal")


StdMmGroup = MmGroup()
Mm.group = StdMmGroup
load_group_name(StdMmGroup, "M0")

######################################################################
# Vectors in the representation of the Monster mod 15
######################################################################

  

class MmV15(AbstractMmRepVector):
    """Models a vector in the representation of the Monster mod 15

    The constructor constructs a vector in the representation of
    the Monster mod 15 depending on the parameter **vector_name**
    passed as an argument.
    If **vector_name** is the string 'v+' or 'v-' then we construct
    the vector :math:`v^+` or :math:`v^-` defined in :cite:`Seysen22`,
    Section 7.2, respectively.

    If **vector_name** is the string 'v1'  then the we construct a
    precomputed vector :math:`v_1` satisfying the requirements
    in :cite:`Seysen22`, Section 6. The computation of such a
    vector :math:`v_1` is discussed in :cite:`Seysen22`, Appendix B.

    Two such vectors in this class may be checked for equality. A
    vector may be multiplied with an element of the Monster group.
    Other operations are not supported.
    """
    p = 15
    __slots__ =  "data"
    def __init__(self, vector_name = 0):
        self.data = mm_vector(15)
        if not vector_name:
            return
        if vector_name == 'v+':
            # Vector v+ in Leech lattice encoding
            V_PLUS = 0x200
            mm_op_store_axis(15, V_PLUS, self.data)
        elif vector_name == 'v-':
            # Vector v+ in Leech lattice encoding
            V_MINUS = 0x1000200
            mm_op_store_axis(15, V_MINUS, self.data)
        elif vector_name == 'v1':
            mm_order_load_vector(self.data)
        else:
            ERR = "Don't know a reduction vector with name %s"
            raise ValueError(ERR % vector_name)


@singleton
class MmV15Space(AbstractMmRepSpace):
    r"""Models a ``196884``-dimensional rep of the monster group 

    This class contains a collection of functions for manipulating
    vectors in the representation :math:`\rho_{15}` of the Monster.
    Such vectors are instances of class MmV15.
    Most of these function are used implicitly in the operators
    applied to these vectors.
    """
    vector_type = MmV15
    space_name = "MmV15"

    def __init__(self):
        """Create 196884-dimensional representation of Monster mod 15
        """
        pass

    def zero(self, p=15):
        """Return the zero vector"""
        return MmV15(0)

    def copy_vector(self, v1):
        assert v1.space == self
        v = MmV15(0)
        np.copyto(v.data, v1.data)
        return v

    def iadd(self, v1, v2):
        mm_op_vector_add(15, v1.data, v2.data)
        return v1

    def imul_scalar(self, v1, a):
        mm_op_scalar_mul(15, a % 15, v1.data)
        return v1
 
    def imul_group_word(self, v1, g):
        """Return product v1 * g of vector v1 and group word g."""
        work =  mm_vector(15)
        if isinstance(g, AbstractMMGroupWord):
            a = g.mmdata
            mm_op_word(15, v1.data, a, len(a), 1, work)
            return v1
        err = "Multiplicator for vector must be in Monster group"   
        raise TypeError(err) 
 
    def equal_vectors(self, v1, v2):
        """Return True iff vectors v1 and v2 are equal"""
        return not mm_op_compare(15, v1.data, v2.data) 

    def str_vector(self, v1):
        return("<Vector in class MmV15>")


StdMmV15Space = MmV15Space()
MmV15.space = StdMmV15Space







    
