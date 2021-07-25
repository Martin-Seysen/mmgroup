r"""We deal with the 196884-dimensional representation of the monster.

The purpose of the ``mmgroup`` package is to give the user the 
capability to compute in the ``196884``-dimensional representation 
:math:`\rho_p` of the monster group :math:`\mathbb{M}`. Here
:math:`\rho` is the rational representation  ``196884_x`` of the 
monster constructed in in :cite:`Con85`, and :math:`\rho_p`
is obtained from  :math:`\rho` by taking all matrix coefficients
modulo a small odd number ``p``.  This package supports the 
cases ``p = 3, 7, 15, 31, 127, 255``.

We are a bit sloppy, calling :math:`\rho_p` a vector space also
in cases ``p = 15, 255``. The user may select these values for
obtaining representations modulo ``5`` and ``17``, respectively.

For calculating in :math:`\rho_p` the user must first create and 
instance of class |MMGroup| that models the group  :math:`\mathbb{M}`
as described in section :ref:`mmgroup-label`. Then the user must 
create an instance of class |MMSpace| that models a representation 
space for that group modulo ``p``, with  ``p`` as above. By 
construction of :math:`\rho_p` is is natural to label a basis vector 
by a tuple ``(tag, i0, i1)`` with integers ``i0, i1``, as described 
in the table below.

If ``V`` is an instance of class  |MMSpace| then the user may call
``V`` with a variable number of arguments to create a vector in
``V``. Here each argument evaluates to a vector in ``V`` and
the sum of these vectors is returned. Such a vector is an instance
of class |MMSpaceVector|. An argument passed to ``V`` may be a 
tuple ``(factor, tag, i0, i1)``. Here ``tag`` is a single 
capital letter and indices ``i0``, ``i1`` are unsigned integers; 
these three quantities describe a basis vector of the space ``V``.
``factor`` is an optional integer (default is ``1``) describing
a factor to be multiplied by the basis vector. Valid tags and
indices for a basis vector are described in the table below.

Vector addition, subtraction and scalar multiplication can be done with
operators ``+``, ``-`` and ``*`` as usual. Elements of the monster group
operate on vectors by right multiplication. Only vectors in the same 
space can be added. An instance ``MM`` of the monster group (of class
|MMGroup|) must be given in the constructor of a vector space ``V``. 
Then only elements of that instance ``MM`` of the monster group can 
act on vectors in ``V`` by right multiplication.

The following code example generates an instance ``MM`` of the
monster group operating on a vector space  ``V`` of integers
modulo ``3``. It creates an element ``g`` of ``MM`` and a vector
``v`` in ``V`` and displays the result ``g * v`` of the 
operation of ``g`` on ``v``. Note that a basis vector
``(tag, i0, i1)`` is displayed in the form ``tag_i0_i1``.
For a displayed index ``i0`` or ``i1`` the suffix ``h`` means 
hexadecimal notation.


.. code-block:: python

  >>> from mmgroup import MMGroup
  >>> from mmgroup import MMSpace
  >>> # Create an instance MM of the monster group
  >>> MM = MMGroup()
  >>> # Create representation space V for MM (modulo 3)
  >>> V = MMSpace(3, MM)
  >>> # Create an element g of MM 
  >>> g = MM(('d', 0x123), ('p', 217821225))
  >>> # Create a vector v in V 
  >>> v = V(('T', 712, 13), (2, 'X', 0x345, 13))
  >>> # Let g operate on v by right multiplication
  >>> print(v * g)
  v<-T_444_0fh+X_244h_11>
  

The following table shows the tags available for describing basis
vectors. The column ``Vector`` in that table shows the basis vector  
in the notation of :cite:`Seysen20`. The column ``Tuple`` shows the 
basis vector as a tuple ``(tag, i0, i1)``.

.. _table-vector-tags:

.. table:: Basis vectors of the representation of the monster
  :widths: 10 15 75

  ============================ ============= =================================
  Vector                       Tuple         Remarks
  ============================ ============= =================================
  :math:`(ij)_1`               ``('A',i,j)`` ``0 <= i,j < 24``; we have
                                             :math:`(ij)_1 = (ji)_1`.

  :math:`X_{ij}`               ``('B',i,j)`` ``0 <= i,j < 24, i != j``; we 
                                             have :math:`X_{ij} = X_{ji}`.  

  :math:`X^+_{ij}`             ``('C',i,j)`` ``0 <= i,j < 24, i != j``; we 
                                             have :math:`X^+_{ij} = X^+_{ji}`.

  :math:`X^+_{d \cdot \delta}` ``('T',o,s)`` ``o = x.octad``, 
                                             ``s = x.suboctad``,
                                             ``x = SubOctad(d, delta)``,

                                             :math:`d \in \mathcal{P}`,
                                             :math:`\delta \in \mathcal{C^*}`,
                                             :math:`d` an octad, 
                                             :math:`\delta \subset d`,
                                             :math:`\delta` even. See class
                                             |SubOctad|.

                                             We have ``0 <= o < 759``,
                                             ``0 <= s < 64``. 

  :math:`X^+_{d \cdot i}`      ``('X',d,i)`` ``0 <= d < 2048``, 
                                             ``0 <= i < 24``,

                                             ``d`` represents the element
                                             ``PLoop(d)`` of the Parker
                                             loop :math:`\mathcal{P}`

  :math:`d^+ \otimes_1 i`      ``('Z',d,i)`` ``d`` and ``i``  as in case              
                                             ``('X',d,i)``

  :math:`d^- \otimes_1 i`      ``('Y',d,i)`` ``d`` and ``i``  as in case              
                                             ``('X',d,i)``
  ============================ ============= =================================

**Remarks**

  * The space :math:`\rho_p` is an orthogonal space. 
    Two basis vectors are orthogonal except when equal or opposite.
    The basis vectors :math:`(ij)_1` have norm ``2`` in case
    :math:`i \neq j`. All other basis vectors have norm ``1``. 

  * In the tuple ``('T',o,s)`` the integers ``o`` and ``s`` describe
    the instance ``SubOctad(o, s)`` of class |SubOctad|. Here 
    ``o`` and ``s`` may be anything that is accepted as input for
    ``SubOctad(o, s)``

  * In the tuples ``('X',d,i)``, ``('Y',d,i)``, ``('Z',d,i)``,
    the input ``d`` may be an instance of class |PLoop| or 
    anything such that ``PLoop(d)`` is an instance of class |PLoop|.
    For an instance ``d`` of class |PLoop| we have:

      * ``('X',d,i) == -('X',-d,i) ==  ('X',~d,i)`` ,
      * ``('Y',d,i) == -('Y',-d,i) == -('Y',~d,i)`` ,
      * ``('Z',d,i) == -('Z',-d,i) ==  ('Z',~d,i)`` .

    The value of ``('X',~d,i)`` has been set by convention;
    the other equations are forced.

  * The basis vector :math:`d^+ \otimes_1 i` in :cite:`Seysen20` is 
    equal to the basis vector 
    :math:`(-1)^{|d/4|} \cdot d^+ \otimes_1 i` in :cite:`Con85`. This 
    modification simplifies the formulas for the operation of the 
    element :math:`\xi` of  :math:`\mathbb{M}` on :math:`\rho_p`.

  * A first or second index following a tag may be ``'r'``, indicating 
    a random integral index. Any omitted index is interpreted as
    ``'r'``.  For tags ``'A', 'B', 'C', 'T'`` the second index should 
    be ``'r'`` or omitted if the first index is ``'r'``, since in these
    cases the two indices are not sufficiently independent.

  * In order to generate a random multiple of a basis vector, the
    component ``factor`` of a tuple can be set to one of the 
    following values:

      * ``'u'``: equivalent to ``1`` 
      * ``'s'``: a factor ``1`` or ``-1`` selected at random
      * ``'n'``: a random nonzerodivisor modulo ``p`` 
      * ``'r'``: a random integer modulo ``p`` 

"""
# References in the __docstr__ see file docs/source/references.bib


from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import sys
import os
import numpy as np
from numbers import Integral
from random import randint
import warnings
import time
from timeit import default_timer
from importlib import import_module






from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
from mmgroup.mm_group  import MMGroup, MMGroupWord, MM

from mmgroup.generators import rand_get_seed
from mmgroup.mm import mm_vector, mm_aux_random_mmv
from mmgroup.mm import mm_aux_zero_mmv, mm_aux_reduce_mmv
from mmgroup.mm import mm_aux_mmv_to_sparse, mm_aux_mmv_add_sparse
from mmgroup.mm import mm_aux_mmv_extract_sparse
from mmgroup.mm import mm_aux_mmv_set_sparse, mm_aux_mmv_add_sparse
from mmgroup.mm import mm_aux_bytes_to_mmv, mm_aux_mmv_to_bytes
from mmgroup.mm import mm_aux_get_mmv, mm_aux_put_mmv
#from mmgroup.mm import mm_rng_make_seed
from mmgroup.mm import INT_BITS, PROTECT_OVERFLOW
from mmgroup.mm import mm_aux_check_mmv
from mmgroup.mm import mm_aux_get_mmv1
from mmgroup.mm import mm_aux_index_extern_to_sparse
from mmgroup.mm import mm_aux_index_sparse_to_extern
from mmgroup.mm import mm_aux_index_sparse_to_leech
from mmgroup.mm import mm_aux_index_sparse_to_leech2

uint_mmv = np.uint32 if INT_BITS == 32 else np.uint64
#standard_seed = mm_rng_make_seed()
standard_mm_group = MM


TAGS = " ABCTXZY"


######################################################################
# Importing a C wrapper for a specific characteristic 'p'
######################################################################

mm_op = {}
all_characteristics_found = None

def mm_wrapper(p):
    """Return the module dealing with 'characteristic' p"""
    try:
        return mm_op[p]
    except KeyError:
        mm_op[p] = import_module('mmgroup.mm%d' % p)
        return mm_op[p] 

def characteristics():
    """Return list of all *characteristics* ``p`` supported"""
    global all_characteristics_found
    if all_characteristics_found is None:
        for k in range(2,9):
            p = (1 << k) - 1
            if not p in mm_op:
                try:
                    mm_op[p] = import_module('mmgroup.mm%d' % p)
                except ModuleNotFoundError:
                    pass
        all_characteristics_found = sorted(mm_op.keys())
    return all_characteristics_found

######################################################################
# Modelling a vector of the 196884-dimensional rep of the monster
######################################################################

class MMSpaceVector(AbstractMmRepVector):
    """Models a vector in a space of type |MMSpace|.

    Such a vector should be constructed by calling an instance ``V``
    of class |MMSpace| which models a specific representation of
    the monster group. See class |MMSpace| for details.

    Addition and subtraction of vectors in the same space work as 
    usual. This is also the case for scalar multiplication with 
    integers. Elements of the moster group operating on ``V`` act 
    on a vector in ``V`` by right multiplication.

    An entry of a vector ``v`` in the space ``V`` may be adressed
    as ``v[tag, i0, i1]``, where ``tag`` is a single letter in
    the string ``ABCTXYZD`` and ``i0`` and ``i1`` are integers,
    see table :ref:`table-vector-tags` and class |MMSpace| for
    details. Getting and setting entries of a vector is as in
    the ``numpy`` package. Here ``i0`` and ``i1`` may also be 
    slices of integers in the same way as in ``numpy`` arrays. 
    Depending on the ``tag``, the value ``i0`` or ``i1`` may also 
    be a single instance of the appropriate class, as described 
    in the remarks after table :ref:`table-vector-tags` .

    The entries of vector ``v``  also have a linear order as 
    described in method ``tuple_to_index`` of class |MMSpace|.
    Here ``v[i]`` is the ``i``-th entry in that order.  Here
    ``i`` may also be a slice of integers in the same way as
    in a one-dimensional ``numpy`` array.

    The internal representation of a vector ``v`` in this class
    is not part of the public interface. Use ``v.as_bytes()`` to 
    convert ``v`` to an one-dimensional array of ``8``-bit 
    integers.
    Method ``from_bytes`` of class |MMSpace| constructs a
    vector ``v`` as an instance of class |MMSpaceVector| from
    a one-dimensional array of integers of length ``196884``.

    :var space:
        This attribute contains the space to which the vector
        belongs. That space is an instance of class |MMSpace|.

    .. warning::
       The constructor of this class is not for public use! You
       may call an instance ``V`` of class  |MMSpace| for
       constructing vectors in the modular representation space
       ``V`` of the monster group.
    """
    __slots__ = "space", "data"
    def __init__(self, space):
        self.space = space
        self.data = mm_vector(self.space.p)

    def check(self):
        """Check if the vector is correct

        Raise ValueError if the vector is errorneous.
        """
        self.space.check(self)

       
    def mul_exp(self, g, e, break_g = False):
        """Multiply the vector with ``g ** e`` inplace

        The vector is changed and the changed vector is returned.
        
        Afterwards, the vector contains an attribute ``last_timing``
        containing the run time of this operation in seconds.
        This is useful for benchmarking.

        By default, we try to simplify the expression ``g ** e`` 
        before multiplying it with the vector. If ``break_g`` is 
        set,  we always do ``abs(e)`` muliplications with ``g`` or
        with its inverse.
        """
        return self.space.vector_mul_exp(self, g, e, break_g)
 
        

######################################################################
# class MMSpace
######################################################################



class MMSpace(AbstractMmRepSpace):
    r"""Models a ``196884``-dimensional representation of the monster group 

    Any such representation is a vector space which has a small odd
    characteristic ``p``. This means that all coordinates of a vector 
    in this  space are taken modulo ``p``.  We are a bit sloppy here, 
    allowing also some  composite  odd integers ``p``. 

    The package contains highly optimized individual programs for
    each characteristic ``p``. Function ``characteristics()`` in
    this package returns the list of legal values for ``p``. 

    :var p:
        Contains the characteristic of the vector 
        space, which is a small odd integer. Function 
        ``mmgroup.characteristics()`` returns the list
        of legal values ``p``.

    :var group:
        Contains the group operating on the vector
        space by right multiplication. 
        That group is an instance of class |MMGroup|. This 
        defaults to the predefined instance ``mmgroup.MM``
        of class |MMGroup|. 


    Thus the instruction ``V = MMSpace(3, MM)``, with ``MM`` an instance of
    class |MMGroup| creates a representation space ``V`` of characteristic
    ``3`` on which the monster group ``MM`` operates.

    See function ``MMS`` for creating a standard representation space.

    The preferred way to construct a vector in the vector space ``V``
    is to call ``V`` with a  variable number of arguments. The each 
    argument is evaluated to a vector in ``V`` according to the 
    following table, and the sum of all these vectors is returned.

    
    .. table:: Legal types for constructing a vector
      :widths: 25 75

      ====================== ============================================= 
      type                   Evaluates to                               
      ====================== ============================================= 
      ``(f,tag, i0, i1)``    ``tag`` must be one of the letters           
                             ``ABCTXYZ``; ``i0, i1`` must be integers. 

                             The tuple ``(tag, i0, i1)`` denotes a basis 
                             vector as in the table above. 

                             ``f`` is an optional integer factor 
                             for the basis vector, default is ``1``.                       
      ---------------------- --------------------------------------------- 
      ``(f,'D', i0)``        Shorthand for ``(f,'A', i0, i0)``   
      ---------------------- --------------------------------------------- 
      ``(f,'I', i0, i1)``    Shorthand for ``f`` times the sum of the
                             basis vectors 

                             ``('A', i0, i0)`` + ``('A', i1, i1)`` - 
                             ``('A', i0, i1)`` - ``2 ('B', i0, i1)``,

                             see remark below. 
      ---------------------- --------------------------------------------- 
      ``(f,'E', i)``         This is ``f`` times the basis vector with
                             *linear* index ``i``. 

                             See method
                             ``tuple_to_index`` for converting a tuple to         
                             a linear index.              
      ---------------------- --------------------------------------------- 
       class |MMSpaceVector| A deep copy of the given vector is returned.            
      ---------------------- --------------------------------------------- 
       ``str``               For an vector ``v`` in ``V`` we have      
                             ``V(str(v)) == v``. 

                             This is helpful for rereading printed 
                             vectors.       
      ====================== ============================================= 

    Remarks

    The meaning of tuples with tags ``ABCTXYZ`` is explained in
    table :ref:`table-vector-tags`. That table also describes
    the legal values for the indices ``i0``, ``i1``.

    The centralizer of a vector labelled by ``('I', i0, i1)`` in
    the monster has structure :math:`2 \cdot B`, where :math:`B` is
    the Baby Monster, see :cite:`Asc86`, :cite:`Con85`, :cite:`Iva09`
    for background. According to our selected basis and cocycle of
    the Golay Code, it is probably easiest to compute in the 
    centralizer of ``('I', 3, 2)``.
    """
    vector_type = MMSpaceVector

    check_errors = {
        -1: "Bad input value p",
        -2: "A one bit outside a field has been found",
        -3: "A subfield has an illegal nonzero entry at index >= 24",
        -4: "Illegal nonzero diagonal entry", 
        -5: "Symmetric part of vector is not symmetric",
    }

    extend_modulus = {
        3:3, 7:7, 15:15, 31:31, 63:63, 127:127, 255:255,
        5:15, 9:63, 17:255,
    }

    tag_indices = { # tag: (start, nrows, ncolumns, stride)
       "A": (     0,   24, 24, 32),
       "B": (   768,   24, 24, 32),
       "C": (  1536,   24, 24, 32),
       "T": (  2304,  759, 64, 32),
       "X": ( 50880, 2048, 24, 32),
       "Z": (116416, 2048, 24, 32),
       "Y": (181952, 2048, 24, 32),
    }

    def __init__(self, p, group = None):
        """Create a 196884-dimensional representation of the monster

        All calculations are done modulo the odd number p
        """
        try: 
            self.p = self.extend_modulus[p]
            self.p_original = p
        except KeyError:
            raise KeyError("Modulus %s not supported for monster group" % p)
        if group is None:
            group = standard_mm_group 
        assert isinstance(group, MMGroup) 
        super(MMSpace, self).__init__(self.p, group)
        mm = mm_wrapper(p)  # fetch mm_op[p]
        #self.mm = mm_wrapper(self.p)
                            # if we do this, we cannot pickle vectors
        self.MMV_INTS = mm_op[self.p].MMV_INTS 
        self.op_vector_add = mm_op[self.p].op_vector_add
        self.op_scalar_mul = mm_op[self.p].op_scalar_mul
        self.op_word = mm_op[self.p].op_word
        self.op_compare = mm_op[self.p].op_compare
        del mm

    @property
    def mm(self):
        """Return module object mmgroup.mm<p> for characteristic p"""
        return mm_op[self.p]


    #######################################################################
    # Creating vectors 
    #######################################################################

    def zero(self):
        """Return the zero vector"""
        return MMSpaceVector(self)

    def copy_vector(self, v1):
        assert v1.space == self
        v = MMSpaceVector(self)
        np.copyto(v.data, v1.data)
        return v

    def rand_uniform(self, seed = None):
        """Return a uniform distributed random vector.

        ``seed`` is a seed for the random generator. The current version 
        supporst the default seed only. Here some random data taken from 
        the operating system and from the clock are entered into the seed.
        """
        v = MMSpaceVector(self)
        seed = rand_get_seed(seed)
        mm_aux_random_mmv(self.p, v.data, seed) 
        return v

    #######################################################################
    # Obtaining and setting components via sparse vectors
    #######################################################################

    def getitems_sparse(self, v1, a_indices):
        """Get items from vector v1

        Here we assert that v1 is a vector of this vector space and
        that 'a_indices' is a one-dimensional numpy array of type
        np.uint32, containing the coordinates to be read from v1.
 
        The function must add the corresponding coordinate to each
        entry of the array 'sparse_items'. All coordinates must be
        nonnegative and < 256. 

        A zero entry in the array 'a_indices' is ignored.
        """
        assert v1.space == self
        if len(a_indices):
            mm_aux_mmv_extract_sparse(self.p, v1.data, a_indices,
                len(a_indices))
        return a_indices 

    def additems_sparse(self, v, a_indices):
        """Add a vector in sparse representation to vector v.

        This method takes a numpy array 'a_indices' of integers of dtype 
        numpy.uint32 containing the description of a vector v2 in sparse 
        representation. It computes 

             v  =  v + v2 .

        Here vector v is a standard vector in this space.
        """
        if len(a_indices):
            mm_aux_mmv_add_sparse(self.p, a_indices, len(a_indices),
                v.data)
        return v

    def setitems_sparse(self, v, a_indices):
        """Setro selected components of a vector 

        Arguments 'v' and 'a_indices' are as in method getitems_sparse().
        Here the coordinates of vector 'v' described by 'a_indices' are 
        set to the values given in 'a_indices'. 
        The array 'a_indices' is not changed.
        """
        assert v.space == self
        if len(a_indices):
            mm_aux_mmv_set_sparse(self.p, v.data, a_indices, 
                len(a_indices))
        return v


    #######################################################################
    # Conversion from and to to sparse representation 
    #######################################################################

    def as_sparse(self, v1):
        """Yet to be documented!!

        """
        sp = np.zeros(196884, dtype = np.uint32)
        length = mm_aux_mmv_to_sparse(self.p, v1.data, sp)
        return sp[:length]


    #######################################################################
    # Vector operations 
    #######################################################################


    def iadd(self, v1, v2):
        self.op_vector_add(v1.data, v2.data)
        return v1
 
    def imul_scalar(self, v1, a):
        self.op_scalar_mul(a % self.p, v1.data)
        return v1
           
    #######################################################################
    # Group operation 
    #######################################################################

    def imul_group_word(self, v1, g):
        """Return product v1 * g of vector v1 and group word g.

        v1 may be destroyed.

        This method is called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.
        """
        work =  mm_vector(self.p)
        assert v1.space == self
        g = self.group(g)
        assert isinstance(g, MMGroupWord) 
        self.op_word(v1.data, g._data, g.length, 1, work)
        return v1  


    def vector_mul_exp(self, v1, g, e, break_g = False):
        """Compute product v1 * g**e of vector v1 and group word g.

        Here v1 is a vector in this space, e is an integer, g is a 
        group element, and  v1 is replaced by v1 * g**e.

        This method should be  called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.

        If break_g is True, each factor g is multiplied with v1
        separately. Otherwise, the expression  g**e  may be
        optimized. This option is mainly for benchmarking.

        After applying this function to vecter v1, the vector
        v1 has an attribute v1.last_timing containing the run
        time of the C part of this operation in seconds.
        """
        work = mm_vector(self.p)
        assert v1.space == self
        assert -1 << 31 < e < 1 << 31
        g = self.group(g)
        assert isinstance(g, MMGroupWord) 
        length = g.length
        if break_g:
            g._extend(length + 1)
            g._data[length] = 0x70000000
            length += 1
        t_start = time.process_time()
        #t_start = default_timer()
        self.op_word(v1.data, g._data, length, e, work)
        v1.last_timing = time.process_time() - t_start
        #v1.last_timing = default_timer() - t_start
        return v1  


    #######################################################################
    # Checking equality
    #######################################################################

    def equal_vectors(self, v1, v2):
        """Return True iff vectors v1 and v2 are equal 

        This method is called for elements v1 and v2 of the space
        'self' only.
        """
        if self.p == self.p_original:
            return not self.op_compare(v1.data, v2.data) 
        return self.as_bytes(v1) == self.as_bytes(v2)

    #######################################################################
    # Conversion from and to byte format
    #######################################################################

    def as_bytes(self, v1):
        """Return vector 'self' as a byte array

        The result is a numpy array with dtype = uint8 and
        shape = (196884,).
        """
        b = np.zeros(196884, dtype = np.uint8)
        mm_aux_mmv_to_bytes(self.p, v1.data, b)
        if self.p != self.p_original:
            b %= self.p_original
        return b

    def from_bytes(self, b):
        """Construct a vector from a byte array

        Here ``b`` is an array-like object representing a 
        one-dimensional array of ``196884`` integers as in the
        ``numpy`` package. These integers are taken modulo the
        characteristic ``p`` of the space. 

        Array ``b`` represents a vector in this space in linear order
        as described method ``tuple_to_index``. 

        The function returns the vector given by the array ``b`` in
        this vector space as an instance of class |MMSpaceVector|.
        """
        b = np.array(b, dtype = np.int32)
        if len(b.shape) != 1:
            raise TypeError("Bad shape of byte data vector")
        if len(b) != 196884:
            raise TypeError("Bad length of byte data vector")
        b = np.array(b % self.p, dtype = np.uint8)
        v =self.zero()
        mm_aux_bytes_to_mmv(self.p, b, v.data)
        return v

 
        
    #######################################################################
    #  Checking and reducing a vector
    #######################################################################

    def check(self, v1):
        """Check the vector 'self'.

        Raise ValueError if an error is found in vector 'self'.
        """
        if len(v1.data) != self.MMV_INTS + 1:
            err = "MM vector has wrong length"
            raise MemoryError(err)   
        if v1.data[-1] != PROTECT_OVERFLOW:
            err = "Buffer overflow in MM vector detected"
            raise MemoryError(err)
        result = mm_aux_check_mmv(self.p, v1.data)
        #print("check result is", result)
        if not result:
            return
        try:
            err = self.check_errors[result]
        except:
            err = "Unknown error %d in MM vector" % result
        print("\n%s!\n" % err)
        raise ValueError("Error in MM vector")
 
    def reduce(self, v1):
        """Convert vector v1 to a unique reduced form"""
        if self.p == self.p_original:
            mm_aux_reduce_mmv(self.p, v1.data)
        else:
            b = np.zeros(196884, dtype = np.uint8)
            mm_aux_mmv_to_bytes(self.p, v1.data, b)
            b %= self.p_original
            mm_aux_bytes_to_mmv(self.p, b, v1.data)
        return v1

    #######################################################################
    # Conversion between tags and indices
    #######################################################################

    @classmethod
    def tuple_to_index(cls, tag, i0 = -1, i1 = -1):
        r"""Convert tuple ``(tag, i0, i1)`` to a linear index

        By construction, it is most natural to use a tuple 
        ``(tag, i0, i1)`` for indexing an entry of a vector in
        :math:`\rho_p`. There is also a linear order of these
        entries. For an index given by a tuple ``(tag, i0, i1)``
        this method returns the linear index ``0 <= i < 196884``
        corresponding to this tuple. Legal tuples are listed
        in the following table:

        .. table:: Conditions for indices ``i0``, ``i1``
          :widths: 25 35 45 

          ================= ==================== ==================== 
          Tag               Condition for ``i0`` Condition for ``i1``
          ================= ==================== ====================
          ``'A'``           ``0 <= i0 < 24``     ``0 <= i1 < 24``
          ``'B', 'C'``      ``0 <= i0 < 24``     ``0 <= i1 < 24``, 
                                                 ``i1 != i0``
          ``'T'``           ``0 <= i0 < 759``    ``0 <= i1 < 64``
          ``'X', 'Y', 'Z'`` ``0 <= i0 < 2048``   ``0 <= i1 < 24``
          ================= ==================== ==================== 

        For ``tag = 'A', 'B', 'C'``, we also have 
        ``(tag, i0, i1) == (tag, i1, i0)``.

        The linear order of the tuples ``(tag, i0, i1)`` is as follows:

        offset ``0``:

          | ``(A,0,0), (A,1,1),  ..., (A,23,23)``

        offset ``24``:
 
          | ``(A,1,0),`` 
          | ``(A,2,0),  (A,2,1)``
          | . . .
          | ``(A,23,0),  (A,23,1), ...,  (A,23,22)``

        offset ``300``:
 
          Tuples ``(B,i0,i1)``, same order as 
          tuples ``(A,i0,i1)`` for ``i0 > i1``

        offset ``576``:

          Tuples ``(C,i0,i1)``, same order as  
          tuples ``(A,i0,i1)`` for ``i0 > i1``
           
        offset ``852``:
 
          | ``(T,  0,0),  (T,  0,1), ...,  (T,  0,63)`` 
          | ``(T,  1,0),  (T,  1,1), ...,  (T,  1,63)`` 
          | . . .
          | ``(T,758,0),  (T,758,1), ...,  (T,758,63)`` 

        offset ``49428``:

          | ``(X,   0,0),  (X,   0,1), ...,  (X,   0,23)`` 
          | ``(X,   1,0),  (X,   1,1), ...,  (X,   1,23)`` 
          | . . .
          | ``(X,2047,0),  (X,2047,1), ...,  (X,2047,23)`` 

        offset ``98580``:

          Tuples ``(Z,i0,i1)``, same order as corresponding  
          tuples ``(X,i0,i1)``
       
        offset ``147732``:

          Tuples ``(Y,i0,i1)``, same order as corresponding  
          tuples ``(X,i0,i1)``


        Remarks:

        The tuple ``('D', i0)`` is accepted as a shorthand for 
        ``('A', i0,i0)``.

        A tuple ``('E', i0)`` means a linear index ``i0``,  i.e.
        this method returns ``i0`` on input ``('E', i0)``, for
        ``0 <= i0 < 196884``.

        If ``tag`` is an instance of class |MMSpaceVector|, which is
        a nonzero multiple of a basis vector, then the linear index 
        corresponding to that basis vector is returned.        
        """
        i = 0
        if isinstance(tag, str) and len(tag) == 1:
            t = TAGS.find(tag)
            if  t >= 1 and 0 <= i0 < 2048 and 0 <= i1 < 64:
                i = (t << 25) + (i0 << 14) + (i1 << 8) 
            elif  tag == "E" and 0 <= i0 < 196884:
                return i0
            elif tag == "D" and 0 <= i0 < 24:
                i = (1 << 25) + (i0 << 14) + (i0 << 8) 
        elif isinstance(tag, MMSpaceVector):
            sp = tag.as_sparse()
            if len(sp) == 1:
                i = sp[0] 
            else:
                err = "MM vector is not multiple of basis vector"
                raise ValueError(err)
        else:
            raise TypeError("Cannot convert object to MM vector index")
        i_ext = mm_aux_index_sparse_to_extern(i)
        if 0 <= i_ext < 196884:
            return i_ext
        err = "Could not convert tuple with tag %s to MM vector index"
        raise ValueError(err % tag)

 
    @classmethod
    def index_to_tuple(cls, index):
        """Convert linear index to tuple ``(tag, i0, i1)``

        This method reverses the effect of method ``tuple_to_index``.
        Given a linear index ``0 <= i < 196884`` for a basis
        vector, the function returns that index as a tuple
        ``(tag, i0, i1)`` with ``tags`` one letter of the
        string ``"ABCTXYZ"`` and integers ``i0``, ``i1``.

        See method ``tuple_to_index`` for details.

        If ``index`` is an instance of class |MMSpaceVector|, which is
        a nonzero multiple of a basis vector, then the index of that 
        basis vector is taken.        
        """
        if isinstance(index, Integral): 
            if  0 <= i < 196884:
                i = mm_aux_index_extern_to_sparse(index)
                return TAGS[i >> 25], (i >> 14) & 0x7ff, (i >> 8) & 0x3f
            else:
                raise ValueError("MM vector index out of range")
        elif isinstance(index, MMSpaceVector):
            sp = index.as_sparse()
            if len(sp) == 1:
                i = sp[0]
                return TAGS[i >> 25], (i >> 14) & 0x7ff, (i >> 8) & 0x3f
            else:
                raise ValueError("MM vector is not multiple of basis vector")
        else:    
            raise TypeError("Cannot convert object to MM index tuple")


    #######################################################################
    # Conversion to short Leech lattice vector
    #######################################################################


    @staticmethod
    def index_to_sparse(tag, i0 = -1, i1 = -1):
        r"""Auxiliary method for index_to_short

        Convert a tagged tuple ``(tag, i0, i1)`` to a sparse 
        index. That tuple must refer to an index describing a
        short Leech lattice vector.
        See method ``index_to_short`` for details. 
        """
        i = 0x0
        if isinstance(tag, Integral) and 300 <= tag < 98580:
            i = mm_aux_index_extern_to_sparse(tag)
        elif isinstance(tag, str) and len(tag) == 1:
            t = TAGS.find(tag)
            if  t >= 1 and 0 <= i0 < 2048 and 0 <= i1 < 64:
                i = (t << 25) + (i0 << 14) + (i1 << 8) 
            elif  tag == "E" and 300 <= i0 < 98580:
                i = mm_aux_index_extern_to_sparse(i0)
            elif tag == "D" and 0 <= i0 < 24:
                i = (1 << 25) + (i0 << 14) + (i0 << 8) 
        elif isinstance(tag, MMSpaceVector):
            sp = tag.as_sparse()
            if len(sp) == 1:
                i = sp[0]
        else:    
            err = "Cannot convert object to short Leech lattice vector"
            raise TypeError(err)
        return i

    @staticmethod
    def index_to_short(tag, i0 = -1, i1 = -1):
        r"""Convert index to a short Leech lattice vector

        If ``tag`` is an integer, this is interpreted a linear index 
        for a basis vector in the representation :math:`\rho_p` as in 
        method ``tuple_to_index``. Otherwise the tuple ``(tag, i0, i1)`` 
        is interpreted a standard index for such a basis vector.
        If ``tag`` is an instance of class |MMSpaceVector|, which is
        a nonzero multiple of a basis vector, then that basis vector 
        is taken.
 
        Some but not all of these basis vectors correspond to short
        vector in the Leech lattice up to sign. If this is the
        case then the function returns a short vector of the Leech
        lattice as a numpy array of signed ``32``-bit integers of 
        length ``24``. The norm of the returned vector is ``32``.

        The function raises ValueError if the basis vector in 
        :math:`\rho_p` does not correspond to a short vector.

        The sign of the short vector is undefined, but two calls
        referring to the same basis vector return the same short
        vector.   
        """
        i = MMSpace.index_to_sparse(tag, i0, i1)
        v = np.zeros(24, dtype = np.int32)
        if mm_aux_index_sparse_to_leech(i, v) == 0:
            return v          
        err = "Vector does not map to short Leech lattice vector"
        raise ValueError(err)

    @staticmethod
    def index_to_short_mod2(tag, i0 = -1, i1 = -1):
        r"""Convert index to a short Leech lattice vector modulo 2

        The tuple ``(tag, i0, i1)`` is interpreted as a basis vector
        as in method ``index_to_short``.
 
        Some but not all of these basis vectors correspond to short
        vector in the Leech lattice up to sign. If this is the
        case then the function returns a short vector of the Leech
        lattice modulo 2 as an integer. 

        The function raises ValueError if the basis vector in 
        :math:`\rho_p` does not correspond to a short vector.
        """
        i = MMSpace.index_to_sparse(tag, i0, i1)
        i2 = mm_aux_index_sparse_to_leech2(i)
        if i2 != 0:
            return i2          
        err = "Vector does not map to short Leech lattice vector"
        raise ValueError(err)


_dict_MMS = {}


def MMS(p):
   """Return a standard space ``MMSpace(p)`` 

   The function returns an instance of class |MMSpace| of
   characteristic ``p``. The group operating on this space is
   the group ``mmgroup.MM``, which is an instance of class
   |MMGroup|.

   Different calls to this function with the same parameter
   ``p`` return exactly the same object.
   """
   global _dict_MMS
   try:
       return _dict_MMS[p]
   except KeyError:
       _dict_MMS[p] = MMSpace(p)
       return _dict_MMS[p]



