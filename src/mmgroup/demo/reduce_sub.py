# Subfunctions for the reduction algorithm for the Monster group
"""
The Python module **mmgroup.demo.reduce_sub** contains auxiliary
functions required for the reduction on an element of the Monster
group.

Most functions in this module correspond to functions defined in
:cite:`Seysen22`. Cross references are given in the documentation
of a function.

Implementations of the functions in this module are usually just Python
wrappers for the corresponding C functions in the **mmgroup** package.
"""

from numbers import Integral
import numpy as np

from mmgroup.demo import Leech2, Mm,  MmV15

from mmgroup.generators import gen_leech3to2
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.generators import gen_leech2_reduce_type2
from mmgroup.generators import gen_leech2_reduce_type2_ortho
from mmgroup.generators import mm_group_n_scan_word
from mmgroup.clifford12 import leech2_matrix_basis
from mmgroup.clifford12 import leech2_matrix_expand
from mmgroup.clifford12 import leech2_matrix_radical
from mmgroup.mm_op import mm_op_norm_A, mm_op_eval_A_rank_mod3
from mmgroup.mm_op import mm_op_eval_A, mm_op_eval_X_find_abs
from mmgroup.mm_op import mm_op_t_A
from mmgroup.mm_reduce import mm_reduce_2A_axis_type
from mmgroup.mm_reduce import mm_order_check_in_Gx0



def mat15_norm(v):
    r"""Compute norm of a matrix related to a vector in :math:`\rho_{15}`

    :param v: vector in representation :math:`\rho_{15}` of the Monster
    :type v: class MmV15

    Let :math:`A` be the 24 times 24 matrix corresponding to part
    :math:`300_x` of vector **v**. The function returns the norm
    of the matrix :math:`A` (reduced mod 15).

    :return: Norm of matrix :math:`A` (mod 15)
    :rtype: int 

    This corresponds to computing :math:`||M(v)|| \pmod{15}` in
    :cite:`Seysen22`, Section 8.2.
    """
    # Implementation is technical in order to be fast
    assert isinstance (v, MmV15) and v.p == 15
    return mm_op_norm_A(15, v.data)


def mat15_apply(v, l2):
    r"""Apply matrix related to a vector in :math:`\rho_{15}` to Leech lattice vector

    :param v: vector in representation :math:`\rho_{15}` of the Monster
    :type v: class MmV15
    :param l2: vector in the Leech lattice mod 2 of type 2
    :type l2: class Leech2

    Let :math:`A` be the 24 times 24 matrix corresponding to part
    :math:`300_x` of vector **v**.

    The function returns :math:`l_2 A l_2^\top` (reduced mod 15),
    where :math:`l_2` is any preimage of **l2** in the Leech
    lattice.

    :return:  :math:`l_2 A l_2^\top \pmod{15}` 
    :rtype: int 

    This corresponds to computing :math:`M(v, l_2) \pmod{15}` in
    :cite:`Seysen22`, Section 8.3. 
    """
    # Implementation is technical in order to be fast
    assert isinstance (v, MmV15) and v.p == 15
    assert isinstance(l2, Leech2)
    return mm_op_eval_A(15, v.data, l2.value)


def mat15_rank_3(v, k):
    r"""Compute rank of a matrix related to a vector in :math:`\rho_{15}` 

    :param v: vector in representation :math:`\rho_{15}` of the Monster
    :type v: class MmV15
    :param k: scalar factor for the 24 times 24 unit matrix
    :type k: int

    Let :math:`A` be the 24 times 24 matrix corresponding to part
    :math:`300_x` of vector **v**. Let :math:`A_k` be the matrix
    :math:`A` minus **k**  times the unit matrix, with coefficients
    taken modulo 3. Matrix :math:`A_k` has a natural interpretation
    as a matrix acting on the Leech lattice (mod 3).

    :return:  tuple(**r, l2**) with **r** the rank of the matrix :math:`A_k`
    :rtype: tuple(int, class Leech2)

    If the kernel of :math:`A_k` is
    one-dimensional and its preimage in the Leech lattice contains
    a nonzero vector **l3** of type at most 4 then we let **l2** be 
    the (unique) vector in the Leech lattice mod 2 with
    **l2 = l3** (mod 2). Otherwise we put **l2 = None**.

    This corresponds to computing the pair 
    :math:`(\mbox{rank}_3(v,k), l_2)`, as defined in
    :cite:`Seysen22`, Section 8.3,
    with :math:`l_2 \in \mbox{ker}_3(v,k)`
    if the conditions mentioned above are satisfied.
    """
    # Implementation is technical in order to be fast
    assert isinstance (v, MmV15) and v.p == 15
    r, l3 = divmod(mm_op_eval_A_rank_mod3(15, v.data, k), 1 << 48)
    u = gen_leech3to2(l3)
    if r != 23 or u <= 0:
        return r, None
    return r, Leech2(u & 0xffffff)

 
def vect15_S(v, k):
    r"""Compute Leech lattice vectors related to a vector in :math:`\rho_{15}`

    :param v: vector in representation :math:`\rho_{15}` of the Monster
    :type v: class MmV15
    :param k: a number 0 < **k** < 8
    :type k: int

    The basis vectors of part :math:`98280_x` of **v** are in
    one-to-one correspondence with the shortest nonzero vectors of the
    Leech lattice, up to sign. The function returns the list of short
    Leech  lattice vectors such that the corresponding
    co-ordinate of (part :math:`98280_x`  of) **v**  has absolute
    value **k** (modulo 15). 

    :return: List of short Leech lattice vectors as described above
    :rtype: list[class Leech2]

    This corresponds to computing :math:`S_k(v)` in :cite:`Seysen22`,
    Section 8.3.
    """
    # Implementation is technical in order to be fast
    MAXLEN = 1024
    assert isinstance (v, MmV15) and v.p == 15
    assert 0 < k < 8
    a = np.zeros(MAXLEN, dtype = np.uint32)
    length = mm_op_eval_X_find_abs(15, v.data, a, MAXLEN, k, 0)
    assert 0 <= length <= MAXLEN
    return [Leech2(i) for i in a[:length]] 


def leech2_span(l2_list):
    r"""List vectors in a subspace of the Leech lattice modulo 2

    :param l2_list: List of vectors in the Leech lattice mod 2
    :type l2_list: list[class Leech2]

    The function returns the list of all vectors in the Leech
    lattice mod 2 that are in the subspace spanned by the vectors
    in the list **l2_list**. 

    :return: List of short Leech lattice vectors as described above
    :rtype: list[class Leech2]

    This corresponds to computing the set :math:`\mbox{span}(S)`, 
    as defined in :cite:`Seysen22`, Section 8.3, where :math:`S`
    is the set of vectors given by the list **l2_list**.
    """
    # Implementation is technical in order to be fast
    MAXDIM = 10
    basis = np.zeros(24, dtype = np.uint64)
    expanded = np.zeros(1 << MAXDIM, dtype = np.uint32)
    l2_array = np.array([x.value for x in l2_list], dtype = np.uint32)
    dim = leech2_matrix_basis(l2_array, len(l2_array), basis, 24)
    assert 0 <= dim <= MAXDIM;
    length = leech2_matrix_expand(basis, dim, expanded);
    assert 0 <= length <= 1 << MAXDIM
    result = [x for x in expanded[:length]]
    result.sort()
    return [Leech2(x) for x in result]
   

def leech2_rad(l2_list):
    r"""List vectors in a subspace of the Leech lattice modulo 2

    :param l2_list: List of vectors in the Leech lattice mod 2
    :type l2_list: list[class Leech2]

    The function returns the list of all vectors in the Leech
    lattice mod 2 that are in the radical of the subspace spanned
    by the vectors in the list **l2_list**. The radical of a
    subspace is the intersection of that subspace with its
    orthogonal complement.

    :return: List of short Leech lattice vectors as described above
    :rtype: list[class Leech2]

    This corresponds to computing the set :math:`\mbox{rad}(S)`, 
    as defined in :cite:`Seysen22`, Section 8.3, where :math:`S`
    is the set of vectors given by the list **l2_list**.
    """
    # Implementation is technical in order to be fast
    MAXDIM = 10
    basis = np.zeros(24, dtype = np.uint64)
    expanded = np.zeros(1 << MAXDIM, dtype = np.uint32)
    l2_array = np.array([x.value for x in l2_list], dtype = np.uint32)
    dim =  leech2_matrix_radical(l2_array, len(l2_array), basis, 24)
    assert 0 <= dim <= MAXDIM;
    length = leech2_matrix_expand(basis, dim, expanded);
    assert 0 <= length <= 1 << MAXDIM
    result = [x for x in expanded[:length]]
    result.sort()
    return [Leech2(x) for x in result]



def map_type4_to_Omega(l2):
    r"""Map a type-4 vector in the Leech lattice mod 2 to :math:`\lambda_\Omega`

    :param l2: Vector of type 4 in the Leech lattice mod 2
    :type l2: class Leech2

    The function returns an element :math:`g` of the group
    :math:`G_{x0}` that maps the vector **l2** of type 4 in
    the Leech lattice mod 2 to the standard type-4 vector
    :math:`\lambda_\Omega`.  

    :return: group element  :math:`g` mapping **l2**  to :math:`\lambda_\Omega`
    :rtype: class Mm


    An implementation of this function is
    discussed in :cite:`Seysen22`, Appendix A. 
    """
    # Implementation is technical in order to be fast
    assert isinstance(l2, Leech2)
    r = np.zeros(6, dtype = np.uint32)
    len_r = gen_leech2_reduce_type4(l2.value, r)
    assert 0 <= len_r <= 6
    return Mm('a', r[:len_r])


def map_type2_to_standard(l2):
    r"""Map a type-2 vector in Leech lattice mod 2 to :math:`\lambda_\beta`

    :param l2: Vector of type 2 in the Leech lattice mod 2
    :type l2: class Leech2

    The function returns an element :math:`g` of the group
    :math:`G_{x0}` that maps the vector **l2** of type 2 in
    the Leech lattice mod 2 to the standard type-2 vector
    :math:`\lambda_\beta`.

    :return: group element :math:`g` mapping **l2**  to :math:`\lambda_\beta`
    :rtype: class Mm

    An implementation of this function is
    discussed in :cite:`Seysen22`, Appendix C. 
    """
    # Implementation is technical in order to be fast
    assert isinstance(l2, Leech2)
    r = np.zeros(6, dtype = np.uint32)
    len_r = gen_leech2_reduce_type2(l2.value, r)
    assert 0 <= len_r <= 6
    return Mm('a', r[:len_r])


def map_feasible_type2_to_standard(l2):
    r"""Map a feasible type-2 vector in the Leech lattice mod 2 to :math:`\lambda_\beta+\lambda_\Omega`

    :param l2: Feasible vector of type 2 in the Leech lattice mod 2
    :type l2: class Leech2

    The function returns an element :math:`g` of the group
    :math:`G_{x0}` that maps the *feasible* vector **l2** of
    type 2 in the Leech lattice mod 2 to the standard feasible
    type-2 vector :math:`\lambda_\Omega + \lambda_\beta`. 
    Here the term *feasible* is defined in :cite:`Seysen22`,
    Section 9.1.

    :return: group element :math:`g` mapping **l2**  to :math:`\lambda_\beta+\lambda_\Omega`
    :rtype: class Mm

    An implementation of this function is
    discussed in :cite:`Seysen22`, Appendix D. 
    """
    # Implementation is technical in order to be fast
    assert isinstance(l2, Leech2)
    r = np.zeros(6, dtype = np.uint32)
    len_r = gen_leech2_reduce_type2_ortho(l2.value, r)
    assert 0 <= len_r <= 6, len_r
    return Mm('a', r[:len_r])



def find_triality_element_for_axis(v, axis_orbits):
    r"""Try to transform an axis in :math:`\rho_{15}` into a given axis orbit

    :param v: an axis in the representation :math:`\rho_{15}` of the Monster
    :type v: class MmV15
    :param axis_orbits: List of expected types of the transformed axis 
    :type axis_orbits: list[str]

    Let :math:`v` be the axis  given by **v**. The function computes 
    the axes :math:`v \cdot \tau^e` for :math:`e =\pm 1`, where
    :math:`\tau` is the triality element in the Monster.
    
    If possible, the function returns an element
    :math:`\tau^e, e = \pm 1` such that the type of the axis
    :math:`v \cdot \tau^e` is in the list **axis_orbits**
    of axis types. Names of axis types are as in
    :cite:`Seysen22`, Section 8.2.
    
    The function raises **ValueError** if this is not possible.
    It does not change the axis stored in **v**. 

    :return: An element :math:`\tau^e` of the Monster as described above
    :rtype: class Mm

    Using function **axis_orbit** in module **mmgroup.demo.reduce_axis**
    this function can be implemented as follows: 

    .. code-block:: python

       def find_triality_element_for_axis(v, axis_orbits)
           for e in [1, -1]:
               if axis_orbit(v * Mm('t', e)) in axis_orbits:
                   return  Mm('t', e)
           raise ValueError
    
    The implementation used here is much faster, since it computes
    fewer co-ordinates of the transformed axes. 
    """
    # Implementation is technical in order to be fast
    assert isinstance(v, MmV15)
    img = np.copy(v.data[:576])
    for e in [1,2]:
        mm_op_t_A(15, v.data, e, img);
        ax_t = mm_reduce_2A_axis_type(img) >> 24;
        ax_type = str(ax_t // 16) + chr(ax_t % 16 + ord('A') - 1) 
        #print("e=", e, ax_type, axis_orbits)
        if ax_type in axis_orbits:
            return Mm('t', e)
    raise ValueError("Could not simplify axis")


    

def find_in_Nx0(v):
    r"""Identify an element of the subgroup :math:`N_{x0}` of the Monster

    :param v: vector in representation :math:`\rho_{15}` of the Monster
    :type v: class MmV15

    Let  :math:`v_1` be the precomputed vector in :math:`\rho_{15}`
    defined in the constructor of class
    :py:class:`~mmgroup.demo.MmV15`.
    Vector **v** must be an image of :math:`v_1` under an element
    of the subgroup :math:`N_{x0}` of the Monster. Then the function
    returns an element :math:`g` of :math:`N_{x0}` with
    :math:`{\bf v} g = v_1`, if such an element :math:`g` exists.
    Otherwise it raises **ValueError**.
    Note that at most one element :math:`g` of the Monster may satisfy
    that condition.
    
    :return: An element :math:`g` of the Monster as described above
    :rtype: class Mm

    An implementation of this identification of an element of the
    group :math:`N_{x0}` is described in
    :cite:`Seysen22`, Section 6. The precomputed vector :math:`v_1`
    satisfies the conditions stated in the same section.
    """
    # Here the implementation just gets the job done
    assert isinstance(v, MmV15)
    work = np.zeros(2*len(v.data), dtype = v.data.dtype)
    g = np.zeros(11, dtype = np.uint32)
    # In the following function call the argument mode = 9 means:
    #  - the result g maps v to v_1 (and not vice versa)
    #  - preserve input v; this requires a bigger work buffer 
    res = mm_order_check_in_Gx0(v.data, g, 9, work)
    ok = 0 <= res <= len(g) # check success of last function
    # Next we check if the result g is in N_x0
    ok = ok and mm_group_n_scan_word(g, res) == res
    if not ok:
        print(res)
        ERR = "Could not identify element of N_x0 from its image"
        raise ValueError(ERR)
    return Mm('a', g[:res])
    



    
