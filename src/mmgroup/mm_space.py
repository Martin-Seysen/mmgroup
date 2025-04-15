r"""We deal with the 196884-dimensional representation of the monster.

The monster :math:`\mathbb{M}` has a 196884-dimensional 
representation :math:`\rho` with coefficients in 
:math:`\mathbb{Z}[\frac{1}{2}]`, see :cite:`Con85`. 
Representation  :math:`\rho` is called 
:math:`196884_x` in  :cite:`Con85`. We obtain a modular 
representation :math:`\rho_p` of :math:`\mathbb{M}` by 
reducing the matrix coefficients of the representation 
:math:`\rho` modulo an odd number :math:`p`. 
The representation  :math:`\rho_p` can be implemented very
efficiently if :math:`p + 1` is a small power of two.
The current version of the *mmgroup* package supports the 
representations :math:`\rho_p` of :math:`\mathbb{M}` for 
``p = 3, 7, 15, 31, 127, 255``. In the sequel 
we are a bit sloppy, calling :math:`\rho_p` a vector space also
in cases ``p = 15, 255``. The user may select these values for
obtaining representations modulo ``5`` and ``17``, respectively.


One purpose of the ``mmgroup`` package is to give the user the 
capability to compute in the ``196884``-dimensional representation 
:math:`\rho_p` of the monster group :math:`\mathbb{M}`.

For calculating in :math:`\rho_p` the user may create an instance 
of class |MM| that models an element ``g`` of the monster group  
:math:`\mathbb{M}` as described in section :ref:`mmgroup-label`. 
Then the user may create an instance ``v`` of class |MMVector| 
that models a vector in :math:`\rho_p`,  with ``p`` as above. 
The element ``g`` of :math:`\mathbb{M}` acts on the 
vector ``v`` by right multiplication.




Creating basis vectors of the representation space :math:`\rho_p`
.................................................................

By construction of :math:`\rho_p` it is natural to label a basis 
vector of :math:`\rho_p` by a tuple ``(tag, i0, i1)`` with integers 
``i0, i1``, as  described in the next section.  Here ``tag`` is a 
single capital letter and indices ``i0``, ``i1`` are unsigned integers; 
these three quantities describe a basis vector of the space ``V``.

To create a basis vector of  :math:`\rho_p` corresponding to the
tuple ``(tag, i0, i1)`` the user may call the constructor
``MMVector(p, tag, i0, i1)``.

For a fixed characteristic ``p`` the function ``MMV(p)`` returns 
an object corresponding to the vector space  :math:`\rho_p`.
Thus the sequence

 .. code-block:: python

  >>> from mmgroup import MMV  
  >>> V = MMV(p)  
  >>> v = V(tag, i0, i1)

creates the same basis vector ``v`` as the 
statement ``v = MMVector(p, tag, i0, i1)``. 




Description of the basis vectors of the representation space
.............................................................


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

  :math:`X^+_{d \cdot \delta}` ``('T',o,s)`` ``o`` and ``s`` can be obtained 
                                             as follows. Let,  
                                             ``x = SubOctad(d, delta)``. Then
                                             ``x.vector_tuples()`` returns
                                             the tuple ``(1, 'T', o, s)``.

                                             :math:`d \in \mathcal{P}`,
                                             :math:`\delta \in \mathcal{C^*}`,
                                             :math:`d` an octad, 
                                             :math:`\delta \subset d`,
                                             :math:`\delta` even. 

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
    the instance ``SubOctad(o, s)`` of class |XLeech2|. Here 
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
    a random integral index. Any omitted index after a given index 
    ``'r'`` is interpreted as ``'r'``.  For tags ``'A', 'B', 'C', 'T'`` 
    the second index should  be ``'r'`` or omitted if the first index 
    is ``'r'``, since in these cases the two indices are not 
    sufficiently independent.


Representing a vector as a list of tuples
.........................................

For creating an arbitrary 
(but yet sparse) vector in :math:`\rho_p` the user may call
``MMVector(p, data_list)``, where ``data_list`` is a list of
tuples ``(factor, tag, i0, i1)``. Here ``(tag, i0, i1)``
describes a unit vector and the (optional) integer ``factor`` is 
a scalar factor for that unit vector. Then the constructor returns 
the linear combination of the unit vectors given in the list.

In order to generate a random multiple of a basis vector, inside a 
list of tuples, the (optional) component ``factor`` of a tuple can 
be set to one of the following values:

      * ``'u'``: equivalent to ``1`` 
      * ``'s'``: a factor ``1`` or ``-1`` selected at random
      * ``'n'``: a random nonzerodivisor modulo ``p`` 
      * ``'r'``: a random integer modulo ``p`` 


Operations on vector in the representation :math:`\rho_p` 
.........................................................

Vector addition, subtraction and scalar multiplication can be done with
operators ``+``, ``-`` and ``*`` as usual. Elements of the monster group
operate on vectors by right multiplication. Only vectors modulo
the same characteristic ``p``  can be added. 

The following code example generates an element ``g`` of the
monster group and a vector ``v`` in the vector space  :math:`\rho_3` 
and displays the result ``g * v`` of the operation of ``g`` on ``v``. 
Note that a basis vector ``(tag, i0, i1)`` is displayed in the 
form ``tag_i0_i1``. For a displayed index ``i0`` or ``i1`` the 
suffix ``h`` means hexadecimal notation.


.. code-block:: python

  >>> # An instance of class MM represents an element of the monster
  >>> from mmgroup import MM
  >>> # Function MMV is used to create a representation space
  >>> from mmgroup import MMV
  >>> # Create representation space V for the monster (modulo 3)
  >>> V = MMV(3)
  >>> # Create an element g of the monster group
  >>> g = MM([('d', 0x123), ('p', 217821225)])
  >>> # Create a vector v in the representation space V 
  >>> v = V([('T', 712, 13), (2, 'X', 0x345, 13)])
  >>> # Let g operate on v by right multiplication
  >>> print(v * g)
  MV<3;-T_444_0fh+X_244h_11>
  

Special tags for creating vectors in the representation :math:`\rho_p` 
.......................................................................

Apart from the standard tags ``A``, ``B``, ``C``, ``T``, ``X``, ``Y``,
and ``Z``, the constructor of class |MMVector| accepts a variety of 
special tags. Details are given in the following list:

.. table:: Special tags for constructing a vector
   :widths: 25 75

   ========================= ============================================== 
   tag, i0, i1               Evaluates to                               
   ========================= ============================================== 
   ``('D', i0)``             Shorthand for ``('A', i0, i0)``   
   ------------------------- ---------------------------------------------- 
   ``('I', i0, i1)``         Shorthand for the sum of the basis vectors 

                             ``('A', i0, i0)`` + ``('A', i1, i1)`` - 
                             ``('A', i0, i1)`` - ``2 ('B', i0, i1)``,

                             see remark below. 
   ------------------------- ---------------------------------------------- 
   ``('J', i0, i1)``         Shorthand for the sum of the basis vectors 

                             ``('A', i0, i0)`` + ``('A', i1, i1)`` - 
                             ``('A', i0, i1)`` + ``2 ('B', i0, i1)``,

                             see remark below. 
   ------------------------- ---------------------------------------------- 
   ``('U')``                 Tag ``'U'`` (without any further parameters)
                             stands for the sum 
                             
                             ``('A', 0, 0)`` + ``('A', 1, 1)`` + ...
                             +  ``('A', 23, 23)``.
   ------------------------- ---------------------------------------------- 
   ``('E', i)``              This is the basis vector with
                             *linear* index ``i``. 

                             See the following subsection for details.              
   ------------------------- ---------------------------------------------- 
   ``('S', data)``           Here ``data`` is an array-like object (as
                             defined in the numpy package) that encodes
                             a vector in sparse representation as an 
                             array of unsigned 32-bit integers. 

                             This is for internal use. For details see
                             the following subsection: 

                             Sparse representation of vectors 
                             in :math:`\rho_p`.            
   ------------------------- ---------------------------------------------- 
   ``('V', data)``           Here ``data`` is an array like object (as
                             defined in the numpy package) that encodes  
                             a one-dimensional array of integers of
                             length 196884. The order of the entries
                             is as in the next section.              
   ------------------------- ---------------------------------------------- 
   ``('R')``                 Tag ``'R'`` (without any further parameters)
                             stands for the generation of a uniform
                             distributed random vector.
   ------------------------- ---------------------------------------------- 
   ``(i, v)``, ``i`` integer Here ``i`` is an integer and ``v`` is a
                             vector, i.e.an instance of class
                             |MMVector|. Then the vector ``i * v``
                             is generated. This is useful for  
                             extending the modulus ``p`` of a 
                             vector, see  remark below.
   ========================= ============================================== 

Remarks

The vectors labelled by ``('I', i0, i1)`` and ``('J', i0, i1)`` are 
axes of the elements :math:`x_\delta` and :math:`x_{-1} x_\delta` of 
the monster, respectively, see  :cite:`Con85`. The centralizers of 
these elements and also of their axes in the monster have structure 
:math:`2 \cdot B`, where :math:`B` is the Baby Monster, see 
:cite:`Asc86`, :cite:`Con85`, :cite:`Iva09` for background. E.g. 
the axes  :math:`v^+` and :math:`v^-` in :cite:`Seysen22` may be 
obtained as ``MMVector(15, 'I', 2, 3)`` and
``MMVector(15, 'J', 2, 3)``, respectively.

The modulus of a vector can be extended as follows. Assume that 
``v`` is a vector in :math:`\rho_3`, given as an instance of class
|MMVector|. Then ``w = 5 * v`` is a well-defined vector in 
:math:`\rho_{15}`. Vector ``w`` can be constructed as follows:
``V15 = MMV(15); w = V15(5, v)``. Note that the construction
``w = V15(5*v)`` leads to an error, since ``5*v`` is defined
modulo ``3``, but not modulo ``15``. 


Linear order of the basis vectors
..................................

By construction of the representation :math:`\rho_p`, it is most natural 
to use a tuple ``(tag, i0, i1)`` for indexing a basis vector of 
:math:`\rho_p`. There is also a linear order of these entries. Such a 
linear order is required e.g. for expressing a vector in :math:`\rho_p` 
as an array of integers  modulo :math:`p`. Therefore an index given by a 
tuple ``(tag, i0, i1)`` is mapped to linear index ``0 <= i < 196884``. 
Legal standard tuples are listed in the  following table:

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





Sparse representation of vectors in :math:`\rho_p` 
....................................................

Internally, a vector :math:`\rho_p` is sometimes given in *sparse*
representation. This is useful for sparse vectors containing many
zero entries. In *sparse* representation the vector is given as
an array of unsigned 32-bit integers. Here each entry of the 
array encodes a multiple of the basis vector. Such an entry it 
interpreted as follows:

.. table:: Bit fields in an entry of a sparse representation
          :widths: 28 24 24 24 24 
 
          ========= =========== =========== ========== ========== 
          Component ``tag``     ``i0``      ``i1``     ``factor``
          Bits      ``27 - 25`` ``24 - 14`` ``13 - 8`` ``7 - 0``
          ========= =========== =========== ========== ==========

This corresponds to the tuple ``(factor, tag, i0, i1)``
which denotes a multiple of a basis vector as described in
subsection *Representing a vector as a list of tuples*. 

Tags are mapped to integers as follows:

.. table:: Mapping of integers to tags
          :widths: 14 14 14 14 14 15 15

          ======= ======= ======= ======= ======= ======= ======= 
          ``'A'`` ``'B'`` ``'C'`` ``'T'`` ``'X'`` ``'Z'`` ``'Y'``
           ``1``   ``2``   ``3``   ``4``   ``5``   ``6``   ``7``
          ======= ======= ======= ======= ======= ======= ======= 


Other data types accepted as tags for vectors in :math:`\rho_p` 
................................................................

Some data types are accepted as tags and interpreted as described in the
following table. In this case parameters `i0`, `i1` after a tag must not
be set.
    
    .. table:: Legal types for constructing a vector
      :widths: 20 80

      ====================== ============================================== 
      type                   Evaluates to                               
      ====================== ============================================== 
      List of tuples         This evaluates to a vector as described in
                             subsection 
                             *Representing a vector as a list of tuples*.

                             An entry of such a list may also be an 
                             instance of class |MMVector| or a string. 
      ---------------------- ---------------------------------------------- 
       class |MMVector|      A deep copy of the given vector is returned.            
      ---------------------- ---------------------------------------------- 
       class |XLeech2|       If ``tag`` is of type |XLeech2| then the
                             (possibly negative) basis vector 
                             corresponding to that ``tag`` 
                             in :math:`Q_{x0}` is created.
      ---------------------- ---------------------------------------------- 
       ``str``               For an vector ``v`` in ``V`` we have      
                             ``V(str(v)) == v``. 

                             This is helpful for rereading printed 
                             vectors.       
      ---------------------- ----------------------------------------------
      ``('Axis', g)``        If a pair containing the string ``'Axis'`` and
                             a 2A involution ``g`` in the Monster group is 
                             given then we contruct the axis corresponding
                             to that involution (with norm 8) as
                             described in :cite:`Con85`. Here ``g``
                             should be an instance of class |XLeech2| or
                             |MM|. A integer or string ``g`` describes the
                             element ``MM('q', g)``.
      ====================== ============================================== 



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
from functools import partial





from mmgroup.structures.abstract_group import singleton
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepVector
from mmgroup.structures.abstract_mm_rep_space import AbstractMmRepSpace
from mmgroup.structures.abstract_mm_group import AbstractMMGroupWord
from mmgroup.structures.mm0_group  import MM0Group, MM0
from mmgroup.structures.abstract_mm_group import is_mmgroup_word
from mmgroup.structures.abstract_mm_rep_space import add_vector

from mmgroup.generators import rand_get_seed
from mmgroup.mm_op import mm_vector, mm_aux_random_mmv
from mmgroup.mm_op import mm_aux_zero_mmv, mm_aux_reduce_mmv
from mmgroup.mm_op import mm_aux_mmv_to_sparse, mm_aux_mmv_add_sparse
from mmgroup.mm_op import mm_aux_mmv_extract_sparse
from mmgroup.mm_op import mm_aux_mmv_set_sparse, mm_aux_mmv_add_sparse
from mmgroup.mm_op import mm_aux_bytes_to_mmv, mm_aux_mmv_to_bytes
from mmgroup.mm_op import INT_BITS, PROTECT_OVERFLOW
from mmgroup.mm_op import mm_aux_check_mmv
from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_index_sparse_to_extern
from mmgroup.mm_op import mm_aux_index_sparse_to_leech
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
from mmgroup.mm_op import mm_aux_index_leech2_to_sparse
from mmgroup.mm_op import mm_aux_hash
from mmgroup.mm_op import mm_aux_mmv_size

from mmgroup import mm_op
from mmgroup.mm_op import mm_op_copy, mm_op_compare
from mmgroup.mm_op import mm_op_checkzero, mm_op_vector_add
from mmgroup.mm_op import mm_op_scalar_mul, mm_op_compare_mod_q
from mmgroup.mm_op import mm_op_store_axis
from mmgroup.mm_op import mm_op_pi, mm_op_xy, mm_op_omega
from mmgroup.mm_op import mm_op_t_A, mm_op_word
from mmgroup.mm_op import mm_op_word_tag_A, mm_op_word_ABC
from mmgroup.mm_op import mm_op_eval_A
from mmgroup.mm_op import mm_op_eval_X_count_abs
from mmgroup.mm_op import mm_op_scalprod

uint_mmv = np.uint32 if INT_BITS == 32 else np.uint64
#standard_seed = mm_rng_make_seed()
standard_mm_group = MM0


TAGS = " ABCTXZY"


_characteristics = None

def characteristics():
    """Return list of all *characteristics* ``p`` supported"""
    global _characteristics
    if _characteristics is not None:
        return _characteristics
    _characteristics = mm_op.characteristics()
    return _characteristics

######################################################################
# Imports on demand
######################################################################

import_pending = True

def complete_import():
    global XLeech2 
    global std_q_element
    from mmgroup import XLeech2
    from mmgroup.structures.construct_mm import std_q_element
    import_pending = False    


import_mm_reduce_pending = True

def complete_import_mm_reduce():
    global mm_reduce_2A_axis_type
    from mmgroup.mm_reduce import mm_reduce_2A_axis_type
    import_mm_reduce_pending = False



#####################################################################
# Modelling a vector of the 196884-dimensional rep of the monster
######################################################################


    

class MMVector(AbstractMmRepVector):
    r"""Models a vector in a representation of the monster group.

    The construction of vectors in the representations :math:`\rho_p`
    of the monster group (using e.g. the constructor of this class)
    is documented at the  beginning of section :ref:`mmrep-label`.
    There a basis vector of that representation is described as a
    tuple ``(tag, i0, i1)``, where ``tag`` is a single capital letter
    and ``i0``, ``i1`` are integers. So we may invoke the constructor
    of this class as follows:

    :param p:   Modulus ``p`` of the representation  :math:`\rho_p`. 

    :param tag: Tag of the unit vector to be constructed. In the
                standard case this is a single capital letter
                describing the type of the basis vector.

    :param i0:  First index in the tuple describing the basis vector.
    
    :param i2:  Second index in the tuple describing the basis vector.
    


    Linear operations can be used for for creating linear combination 
    of unit vectors. This is just the simplest case of the  
    construction of a vector in  :math:`\rho_p`. Section
    :ref:`mmrep-label` contains a large number of further 
    possibilities for creating such vectors.

    Addition and subtraction of vectors in the same space work as 
    usual. This is also the case for scalar multiplication with 
    integers. An element ``g`` of the monster group (given as
    an instance of class |MM|)  operates on a vector ``v`` (given
    as an instance of this class) by right multiplication.

    An entry of a vector ``v`` (given as an instance of this class)
    may be addressed as ``v[tag, i0, i1]``, where ``tag`` is a single 
    letter in the string ``ABCTXYZD`` and ``i0`` and ``i1`` are 
    integers, such that the tuple ``(tag, i0, i1)`` describes a 
    basis vector. Getting and setting entries of a vector is as in
    the ``numpy`` package. Here ``i0`` and ``i1`` may also be 
    slices of integers in the same way as in ``numpy`` arrays. 
    Depending on the ``tag``, the value ``i0`` or ``i1`` may also 
    be a single instance of the appropriate class, as described 
    in the remarks after table :ref:`table-vector-tags` .

    The entries of vector ``v``  also have a linear order as 
    described at the beginning of section  :ref:`mmrep-label`.
    Here ``v['E', i]`` is the ``i``-th entry in that order.  Index
    ``i`` may also be a slice of integers in the same way as
    in a one-dimensional ``numpy`` array.

    An instance ``x`` of class ``XLeech2`` correponding to a short
    vector in the Leech lattice is mapped to a (possibly negated)
    basis vector of :math:`\rho_p`. In this case ``v[x]`` is
    the co-ordinate of ``v`` with respect to that basis vector.

    The internal representation of a vector ``v`` in this class
    is not part of the public interface. Use ``v['E']`` to 
    convert ``v`` to an one-dimensional array of ``8``-bit 
    integers of length 196884.
    """
    __slots__ = "p", "data"
    #group = MM
    def __init__(self, p, tag = 0, i0 = None, i1 = None):
        if p not in characteristics():
             if isinstance(p, Integral):
                 err = "Bad characteristic p = %d for class MMVector"
                 raise ValueError(err % p)
             else:
                 err = "Characteristic p for class MMVector must be integer"
                 raise TypeError(err)
        self.p = p
        self.data = mm_vector(p)
        add_vector(self, tag, i0, i1)

    def check(self):
        """Check if the vector is correct

        Raise ValueError if the vector is erroneous.
        """
        self.space.check(self)

       
    def mul_exp(self, g, e = 1, break_g = False):
        r"""Multiply the vector with ``g ** e`` inplace

        Here ``g`` is an element of the monster group represented
        as an instance of class |MM| and ``e`` is an integer.
        The vector is updated and the updated vector is returned.
        
        Afterwards, the vector contains an attribute ``last_timing``
        containing the run time of this operation in seconds.
        This is useful for benchmarking.

        By default, we try to simplify the expression ``g ** e`` 
        before multiplying it with the vector. If ``break_g`` is 
        set,  we always do ``abs(e)`` multiplications with ``g`` or
        with its inverse.
        """
        return self.space.vector_mul_exp(self, g, e, break_g)
 
    def eval_A(self, v2, e = 0):
        r"""Internal method, not for public use

        The part of this vector with tag 'A' corresponds to  a
        symmetric 24 times 24 matrix :math:`A`. 

        Let :math:`v_2` be a short vector in the Leech lattice 
        modulo 2 given by  parameter ``v2``, where ``v2`` is an 
        instance of class |XLeech2|. If ``v2`` is an integer then 
        this is converted ``XLeech2(v2)``.

        Then in the Leech lattice shortest preimage :math:`v` of
        :math:`v_2` is determined up to sign and :math:`v A v^\top`
        is determined uniquely.

        In case ``e = 0`` (default) the function returns 
        :math:`v_2 A v_2^\top`. Otherwise the function returns
        :math:`v_2 (A \tau^e) v_2^\top`, where :math:`\tau` is
        the triality element in the monster group.
 
        The current version supports vectors modulo ``p = 15`` only.

        The short Leech lattice vector :math:`v_2` (of norm 4) is 
        scaled to norm 32 as usual, when :math:`v_2` is given in 
        integer coordinates.

        The function raises ValueError if :math:`v_2`  is not
        a short Leech lattice vector or if ``p != 15``.
        """
        if import_pending:
            complete_import()
        if isinstance(v2, XLeech2):
            v2 = v2.value
        v1 = np.zeros(24*4, dtype = np.uint64)
        mm_op_t_A(self.p, self.data, e % 3, v1)
        if self.p != 15:
            err = "Method eval_A is implemented for vectors mod 15 only"
            raise ValueError(err)
        res = mm_op_eval_A(self.p, v1, v2)
        if res < 0:
            err = "Method eval_A failed on vector in rep of monster"
            raise ValueError(err)
        return res
        

    def count_short(self):
        r"""Count certain entries of the vector

        The function counts the absolute values of all entries 
        of the monomial part of the vector.
   
        Here the monomial part consists of the entries with
        tags 'B', 'C', 'T', 'X'. These entries correspond to 
        the short vectors of the Leech lattice.

        The function returns a tuple of length ``(p + 1)/2``, where
        ``p`` is the characteristic of the vector. Entry ``i`` of
        that tuple contains the number of entries of the monomial 
        part of the vector with absolute value ``i``. 

        The returned tuple is invariant under operation of the 
        subgroup :math:`G_{x0}` of the monster group.

        Caution: The function is implemented for characteristic
        ``p = 15`` only!
        """
        #if import_pending:
        #    complete_import()
        if self.p != 15:
            err = "Method supported for characteristic p = 15 only"
            raise ValueError(err)
        a = np.zeros(8, dtype = np.uint32)
        mm_op_eval_X_count_abs(self.p, self.data, a)
        return tuple(a)

    def axis_type(self, e = 0):
        r"""Return axis type if this vector is a 2A axis 

        If this vector is a 2A axis then the function computes the 
        type of the 2A axis. Each 2A axis corresponds uniquely to 
        a 2A involution :math:`x` in the monster. Let :math:`z` be 
        the central involution in the subgroup :math:`G_{x0}` of 
        the monster. Then the type of the 2A axis ``v`` is the class 
        of the product :math:`xz` in the monster. Such a class must 
        be 2A, 2B, 4A, 4B, 4C, 6A, 6C, 6F, 8B, 10A, 10B, or 12C.
   
        In case ``e = 0`` (default) the function returns the type 
        of the axis as a string if this vector is a 2A axis. It 
        may return ``None`` if this vector is not a 2A  axis.

        In case ``e != 0`` the function replaces this vector
        :math:`v` by :math:`v \cdot \tau^e`, where :math:`\tau`
        is the triality element in the monster group.

        The function works for a vector space modulo ``p = 15`` 
        only. It raises ValueError in case ``p != 15`` 

        Caution:

        This is a quick disambiguation of the type of a 2A axis. The 
        function may return any axis type if this is not  a 2A axis.
        """
        if self.p != 15:
            err = "Method supported for characteristic p = 15 only"
            raise ValueError(err)
        e, v = e % 3, self.data
        if e:
            v = np.zeros(24*4, dtype = np.uint64)
            mm_op_t_A(self.p, self.data, e, v)     
        if import_mm_reduce_pending:
            complete_import_mm_reduce()
        t = mm_reduce_2A_axis_type(v)
        if t == 0:
            return None
        s = str((t >> 28) & 15) + "?ABCDEFGHIJKLMNO"[(t >> 24) & 15]
        return s

    def hash(self, tags = None):
        r"""Return a hash value of the vector

        If parameter ``tag`` is set then it must be a string containing
        some of the letters ``ABCTXZY``. Then the hash value is computed
        over entries with a tag contained in that string only.
        """
        skip = 0
        if tags is not None:
            for i, t in emumerate("ABCTXZY"):
                if t not in tags:
                     skip |= 1 << i
        return int(mm_aux_hash(self.p, self.data, skip))



    def __mod__(self, p):
        """Return the vector modulo ``p``. 

        The function returns a vector object of class ``MMVector``
        if reduction modulo ``p`` is possible.
        """ 
        if  p == self.p:
            v = MMVector(self.p, 0)
            np.copyto(v.data, self.data)
            return v
        elif isinstance(p, Integral) and p > 0 and self.p % p == 0:
            return MMVector(p, self)
        elif isinstance(p, Integral):
            err = "Cannot reduce MMVector object modulo %d"
            raise ValueError(err % p)
        else:
            err = "Modulus for reducing MMVector object must be int"
            raise TypeError(err)

            

######################################################################
# class MMSpace
######################################################################


@singleton
class MMSpace(AbstractMmRepSpace):
    r"""Models a ``196884``-dimensional representation of the monster group 

    This class contains a collection of functions for manipulating
    vectors in the representation :math:`\rho_p` of the monster group.
    Such vectors are instances of class |MMVector|.
    Most of these function are used implicitly in the operators
    applied to these vectors.

    Some of these function may be helpful for the user; so we document 
    them in the *mmgroup API reference*.    
    """
    vector_type = MMVector
    space_name = "MV"

    check_errors = {
        -1: "Bad input value p",
        -2: "A one bit outside a field has been found",
        -3: "A subfield has an illegal nonzero entry at index >= 24",
        -4: "Illegal nonzero diagonal entry", 
        -5: "Symmetric part of vector is not symmetric",
    }



    def __init__(self):
        """Create a 196884-dimensional representation of the monster

        All calculations are done modulo the odd number p
        """
        pass

    @property
    def mm(self, p):
        """Return module object mmgroup.mm<p> for characteristic p"""
        characteristics()
        return mm_op[p]


    #######################################################################
    # Creating vectors 
    #######################################################################

    def zero(self, p):
        """Return the zero vector"""
        return MMVector(p, 0)

    def copy_vector(self, v1):
        assert v1.space == self
        v = MMVector(v1.p, 0)
        np.copyto(v.data, v1.data)
        return v

    def set_rand_uniform(self, v1, seed = None):
        """Return a uniform distributed random vector.

        ``seed`` is a seed for the random generator. The current version 
        supports the default seed only. Here some random data taken from 
        the operating system and from the clock are entered into the seed.
        """
        seed = rand_get_seed() if seed is None else seed
        mm_aux_random_mmv(v1.p, v1.data, seed) 
        return v1

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
        if len(a_indices):
            mm_aux_mmv_extract_sparse(v1.p, v1.data, a_indices,
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
            mm_aux_mmv_add_sparse(v.p, a_indices, len(a_indices),
                v.data)
        return v

    def setitems_sparse(self, v, a_indices):
        """Set selected components of a vector 

        Arguments 'v' and 'a_indices' are as in method getitems_sparse().
        Here the coordinates of vector 'v' described by 'a_indices' are 
        set to the values given in 'a_indices'. 
        The array 'a_indices' is not changed.
        """
        if len(a_indices):
            mm_aux_mmv_set_sparse(v.p, v.data, a_indices, 
                len(a_indices))
        return v




    #######################################################################
    # Conversion from and to to sparse representation 
    #######################################################################

    def as_sparse(self, v1):
        """Yet to be documented!!

        """
        sp = np.zeros(196884, dtype = np.uint32)
        length = mm_aux_mmv_to_sparse(v1.p, v1.data, sp)
        return sp[:length]


    #######################################################################
    # Vector operations 
    #######################################################################


    def iadd(self, v1, v2):
        if v1.p == v2.p:
            mm_op_vector_add(v1.p, v1.data, v2.data)
            return v1
        else:
            err = "Cannot add vectors modulo differnt numbers"
            raise ValueError(err)
 
    def imul_scalar(self, v1, a):
        mm_op_scalar_mul(v1.p, a % v1.p, v1.data)
        return v1
 

    #######################################################################
    # Support for axes 
    #######################################################################



          
    @classmethod
    def _mm_element_to_axis(cls, mm):
        if import_pending:
            complete_import()
        try:
            inv_type, h = MM0(mm).conjugate_involution()
            assert inv_type == 1
        except:
            err = "Group element after tag 'Axis' must be a 2A involution"
            raise ValueError(err)
        return "v+", h**-1



    @classmethod
    def _make_axis(cls, p, i0, *_):
        if import_pending:
            complete_import()
        if isinstance(i0, str):
            try:
                i0 = std_q_element('q', i0)
            except:
                err = "Unknown string after tag 'Axis'"
                raise ValueError(err)
        elif isinstance(i0, XLeech2):
            i0 = i0.ord
        if isinstance(i0, Integral):
            v = cls.vector_type(p, 0)
            mm_op_store_axis(p, i0, v.data)
            return v
        elif isinstance(i0, AbstractMMGroupWord):
            std_axis, g = cls._mm_element_to_axis(i0)
            return cls._make_axis(p, std_axis) * g
        else:
            err = "Illegal type %s after tag 'Axis'"
            raise VaueError(err % type(i0))



    #######################################################################
    # Group operation 
    #######################################################################

    def imul_group_word(self, v1, g):
        """Return product v1 * g of vector v1 and group word g.

        v1 may be destroyed.

        This method is called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.
        """
        work =  mm_vector(v1.p)
        if isinstance(g, AbstractMMGroupWord):
            a = g.mmdata
            mm_op_word(v1.p, v1.data, a, len(a), 1, work)
            return v1
        err = "Multiplicator for MM vector must be int or in MM group"   
        raise TypeError(err) 
 

    def vector_mul_exp(self, v1, g, e, break_g = False):
        """Compute product v1 * g**e of vector v1 and group word g.

        Here v1 is a vector in this space, e is an integer, g is a 
        group element, and  v1 is replaced by v1 * g**e.

        This method should be  called for elements v1 of the space
        'self' and for elements g of the group 'self.group' only.

        If break_g is True, each factor g is multiplied with v1
        separately. Otherwise, the expression  g**e  may be
        optimized. This option is mainly for benchmarking.

        After applying this function to vector v1, the vector
        v1 has an attribute v1.last_timing containing the run
        time of the C part of this operation in seconds.
        """
        work = mm_vector(v1.p)
        assert v1.space == self
        assert -1 << 31 < e < 1 << 31
        assert isinstance(g, MM0) 
        length = g.length
        if break_g:
            g._extend(length + 1)
            g._data[length] = 0x70000000
            length += 1
        t_start = time.perf_counter()
        #t_start = default_timer()
        mm_op_word(v1.p, v1.data, g._data, length, e, work)
        v1.last_timing = time.perf_counter() - t_start
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
        if v1.p == v2.p:
            return not mm_op_compare(v1.p, v1.data, v2.data) 
        return False

    #######################################################################
    # Conversion from and to byte format
    #######################################################################

    def as_bytes(self, v1):
        """Return vector 'self' as a byte array

        The result is a numpy array with dtype = uint8 and
        shape = (196884,).
        """
        b = np.zeros(196884, dtype = np.uint8)
        mm_aux_mmv_to_bytes(v1.p, v1.data, b)
        return b

    def from_bytes(self, p, b):
        """Construct a vector from a byte array

        Here ``b`` is an array-like object representing a 
        one-dimensional array of ``196884`` integers as in the
        ``numpy`` package. These integers are taken modulo the
        characteristic ``p`` of the space. 

        Array ``b`` represents a vector in this space in linear order
        as described method ``tuple_to_index``. 

        The function returns the vector given by the array ``b`` in
        this vector space as an instance of class |MMVector|.
        """
        b = np.array(b, dtype = np.int32)
        if len(b.shape) != 1:
            raise TypeError("Bad shape of byte data vector")
        if len(b) != 196884:
            raise TypeError("Bad length of byte data vector")
        b = np.array(b % p, dtype = np.uint8)
        v =self.zero(p)
        mm_aux_bytes_to_mmv(p, b, v.data)
        return v

 
        
    #######################################################################
    #  Checking and reducing a vector
    #######################################################################

    def check(self, v1):
        """Check the vector 'self'.

        Raise ValueError if an error is found in vector 'self'.
        """
        if len(v1.data) != mm_aux_mmv_size(v1.p) + 1:
            err = "MM vector has wrong length"
            raise MemoryError(err)   
        if v1.data[-1] != PROTECT_OVERFLOW:
            err = "Buffer overflow in MM vector detected"
            raise MemoryError(err)
        result = mm_aux_check_mmv(v1.p, v1.data)
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
        mm_aux_reduce_mmv(v1.p, v1.data)
        return v1

    #######################################################################
    # Conversion between tags and indices
    #######################################################################

    @staticmethod
    def index_to_sparse(tag, i0 = -1, i1 = -1):
        r"""Convert tuple ``(tag, i0, i1)`` to a sparse index

        The method converts an index referring to a basis
        vector in the representation of the Monster to a linear 
        index.

        Input parameters are as in method ``index_to_linear``, but
        the function returns the sparse index corresponding to the
        input parameters instead of the linear index.
        """
        if isinstance(tag, str) and len(tag) == 1:
            t = TAGS.find(tag)
            if  t >= 1 and 0 <= i0 < 2048 and 0 <= i1 < 64:
                return  (t << 25) + (i0 << 14) + (i1 << 8)
            elif  tag == "E" and 0 <= i0 < 196884:
                return mm_aux_index_extern_to_sparse(i0)
            elif tag == "D" and 0 <= i0 < 24:
                return  (1 << 25) + (i0 << 14) + (i0 << 8)
            elif tag == "S":
                if i0 >= 0:
                    return  i0 & 0x0fffff00
            else:
                raise ValueError("Cannot convert tuple to MM vector index")
        elif isinstance(tag, Integral):
            if  0 <= tag < 196884:
                return mm_aux_index_extern_to_sparse(tag)
            else:
                raise ValueError("MM vector index out of range")
        elif isinstance(tag, MMVector):
            sp = tag.as_sparse()
            if len(sp) == 1:
                return sp[0] & 0xfffff00
            else:
                err = "MM vector is not multiple of basis vector"
                raise ValueError(err)
        elif import_pending:
            try:
                complete_import()  
                if isinstance(tag, XLeech2):
                    i = mm_aux_index_leech2_to_sparse(tag.ord & 0xffffff)
                    if i > 0:
                        return i                   
                    err = "Vector in Leech lattice mod 2 is not short"
                    raise ValueError(err)
            except:
                pass
        ERR = "Cannot convert %s object to MM vector index"
        raise TypeError(ERR % type(tag))



    @classmethod
    def index_to_linear(cls, tag, i0 = -1, i1 = -1):
        r"""Convert an index to a linear index

        The method converts an index referring to a basis
        vector in the representation of the Monster to a linear 
        index. Starndard tuples as in the constuctor of a vector
        in that representation are accepted. Furthermore, the
        following tags or tuples are accepted as input:

        The tuple ``('D', i0)`` is accepted as a shorthand for 
        ``('A', i0,i0)``.

        A tuple ``('S', i0)`` means an index ``i0`` in *sparse*
        representation.

        A tuple ``('E', i0)`` means a linear index ``i0``, for
        ``0 <= i0 < 196884``.

        If ``tag`` is an integer ``i0`` then this is equivalent
        to an input  ``('E', i0)``.

        If ``tag`` is an instance of class |MMVector|, which is
        a nonzero multiple of a basis vector, then the index 
        corresponding to that basis vector is taken as input.

        If ``tag`` is an instance of class |XLeech2| encoding a
        short vector in the Leech lattice mod 2 then the index
        corresponding the that short vector is taken as input.

        A sign (possibly encoded in the input data) is ignored.         
        """
        i = cls.index_to_sparse(tag, i0 , i1)
        i_ext = mm_aux_index_sparse_to_extern(i)
        if 0 <= i_ext < 196884:
            return i_ext
        err = "Could not convert tuple with tag %s to MM vector index"
        raise ValueError(err % tag)


    @classmethod
    def tuple_to_index(cls, tag, i0 = -1, i1 = -1):
        r"""Deprecated, equivalent to method ``index_to_linear``"""
        W = "Method tuple_to_index() of class MMSpace is deprecated. "
        "Use method index_to_linear() instead!" 
        warnings.warn(W, UserWarning)
        return cls.index_to_linear(tag, i0 = -1, i1 = -1)
 
    @classmethod
    def index_to_tuple(cls, tag, i0 = -1, i1 = -1):
        """Convert an index to tuple ``(tag, i0, i1)``

        The function converts an index referring to an entry of
        a vector in the representation of the Monster to tuple of
        shape  ``(tag, i0, i1)`` with ``tags`` one letter of the
        string ``"ABCTXYZ"`` and integers ``i0``, ``i1``.

        Input parameters are as in method ``index_to_linear``.       
        """
        i = cls.index_to_sparse(tag, i0 , i1)
        return TAGS[i >> 25], (i >> 14) & 0x7ff, (i >> 8) & 0x3f


    #######################################################################
    # Conversion to short Leech lattice vector
    #######################################################################



    @staticmethod
    def index_to_short(tag, i0 = -1, i1 = -1):
        r"""Convert index to a short Leech lattice vector
 
        The function converts an index referring to an entry of
        a vector in the representation of the Monster to a
        short Leech lattice vector. 

        Input parameters are as in method ``index_to_linear``.       

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

        The function converts an index referring to an entry of
        a vector in the representation of the Monster to a short
        vector in Leech lattice mod 2. The result is returned as
        an 24-bit integer as described in the documentation of
        class |XLeech2|.

        Input parameters are as in method ``index_to_linear``.       
 
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


StdMMSpace = MMSpace()
MMVector.space = StdMMSpace



def MMV(p):
    r"""Return an object corresponding to the space :math:`\rho_p`

    Here characteristic ``p`` is as in the constructor of
    class |MM|. Thus the sequence

    .. code-block:: python

       >>> V = MMV(p)  
       >>> v = V(tag, i0, i1)

    returns the same vector in :math:`\rho_p` as the 
    statement ``MMVector(p, tag, i0, i1)``. 

    Function ``characteristics`` in module ``mmgroup`` returns
    the list of legal values for the characteristic ``p``.

    Technically, ``MMV(p)`` is the partial application of
    the constructor of class |MMVector| to the number ``p``. 
    """
    return partial(MMVector, p)


def mmv_scalprod(v1, v2):
    r"""Return scalar product of two vectors

    The function returns the scalar product of the vectors :math:`v_1`
    and :math:`v_2`. These two vectors must be instances of
    class |MMVector|; and they must be vectors in the same
    representation :math:`\rho_p`.
    """
    if not (isinstance(v1, MMVector) and isinstance(v2, MMVector)):
        err = "Arguments of mmv_scalprod() must be of type MMVector"
        raise TypeError(err)
    p = v1.p
    if (p != v2.p):
        err = "Arguments of mmv_scalprod() must have same modulus p"
        raise ValueError(err)
    res = mm_op_scalprod(p, v1.data, v2.data)
    assert res >= 0, "Error in function mm_op_scalprod()"
    return res
    

def order_vector():
    r"""Return the precomputed order vector

    The function returns the precomputed order vector :math:`v_1`
    (which is an element of the representation  :math:`\rho_{15}`
    of the monster) as an instance of class |MMVector|.

    The properties required for such a vector :math:`v_1`
    are defined in :cite:`Seysen22`.
    """
    v = MMVector(15)
    try:
        from mmgroup.mm_reduce import mm_order_load_vector
        mm_order_load_vector(v.data)
    except (ImportError, ModuleNotFoundError):     
        from mmgroup.dev.mm_reduce.py_mm_order import ov
        mm_op_copy(15, ov.order_vector.data, v.data)
    return v


