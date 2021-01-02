import sys
import os
import warnings

from random import randint
from collections import OrderedDict

import numpy as np

import mmgroup
from mmgroup.mm_group import MMGroup, MMGroupWord
from mmgroup.mm_space import MMSpace
from mmgroup import structures
from mmgroup.generators import mm_group_check_word_n
from mmgroup.generators import mm_group_words_equ
from mmgroup.generators import mm_group_n_mul_element
from mmgroup.generators import mm_group_n_reduce_word 
from mmgroup.mm import mm_vector
from mmgroup.mm15 import op_copy as mm_op15_copy
from mmgroup.mm15 import op_compare as mm_op15_compare
from mmgroup.mm15 import op_word as mm_op15_word


###########################################################################
# Creating and saving an order vector for computing the group order 
###########################################################################

MMV15 = MMSpace(15)
MM = MMV15.group





_DIR = os.path.split(structures.__file__)[0]
FILENAME = os.path.join(_DIR, "mm_order_vector.py")

order_vector = None

def factors(n):
    """Return sorted list of prime factors of n, for n <= 120"""
    assert 1 <= n < 121, "Number too large"
    primes = [2,3,5,7]
    factors = []
    for p in primes:
        q, r = divmod(n, p)
        if r == 0:            
            factors.append(p)
            n = q
            while n % p == 0:
                n = n // p
    if n > 1:
        factors.append(n)
    return factors


def rand_mm_element(size = 1):
    """Return a random element ``g`` of the moster group

    Here ``size`` is a measure for the requested word length. 
    of ``g``. 

    ``g`` is returend as an element of the standard instance of
    class  ``mmroup.MMGroup``.
    """    
    group = MM
    m = group.neutral()
    for i in range(size):
        m *= group.rand_word('t', randint(1,2))
        m *= group.atom('l', randint(1,2))   
        for s in "yxdp":
            m *=  group.rand_word(s)
    return m
    


def find_element_of_order(order, p, v = 0, minsize = 1, verbose = 0):
    """Return an element ``g`` the monster of order ``order``.

    We ``g`` that check acts as an element of given ``order`` on a
    vector ``v`` in characteristic ``p``.  So ``g`` has that order
    with a very high probalibility, if ``v`` is selected at random.
    But, theoretically, the order of ``g`` may be a multiple of 
    ``order``. ``v`` should be specified as a list of tuples or set 
    to set to zero for generating a random vector.

    ``g`` is returned as an element of the standard instance of
    class  ``mmroup.MMGroup``.

    Parameter ``minsize`` specifies a measure for requested word 
    length of ``g``.  If no element ``g`` of the requested size
    is found, we look for larger candidates. There is no need 
    for changing the default value of ``minsize``, except if you
    want to make experiments for accelerating special cases.
    """
    if verbose:
        print("Find element of order %s in monster group" % order)
        if (v and len(v) < 6):
            print("using vector: ", v)
    rfactors = reversed(factors(order))
    cofactors = [0] + [order // x for x in rfactors] + [order]
    maxsize = minsize + 10
    sp = MMSpace(p)
    v = sp(*v) if v else sp.rand_vector()
    while(minsize < maxsize):
        for n in range(1000):
            ok = True
            g = rand_mm_element(minsize)
            w = v.copy()
            for i in range(1, len(cofactors)-1):
                w.mul_exp(g, cofactors[i] - cofactors[i-1])
                if v == w:
                    ok = False
                    break
            if ok:
                w.mul_exp(g,  cofactors[-1] - cofactors[-2])
                if v == w:
                    if verbose:
                        msg = "Element found after %d trials"
                        print(msg % n)
                    return g
        minsize += 1
    err = "No element of order %s found in the monster group"
    raise ValueError(err % order)
            
 

def compute_order_vector(s_v71, s_g71, s_v94, s_g94, verbose = 0):
    """Compute vector for testing the order of an element of the monster.

    The function returns a vector ``w`` in the 194884-dimensional
    representation of the monster such that only the neutral
    element of the monster fixes ``w``. We call such a vector 
    an *order vector*.
    
    The function requires the following input:

       vector ``v_71`` and group element ``g_71`` such that 
       ``v_71 * g_71**k`` = ``v_71`` (mod 5) iff ``k`` is a 
       multiple  of 71.

       vector ``v_94`` and group element ``g_94`` such that 
       ``v_94 * g_94**k`` = ``v_94`` (mod 3) iff ``k`` is a 
       multiple  of 94.

    All these arguments must be passed as strings, in the way as 
    instances of classes ``MMGroupWord`` and ``MMSpaceVector``
    are represented as strings. This is a very compact, but yet
    readable notation for group elements and vectors.  

    The function returns a vector w with

       ``w = sum  v71 * g_71**i  (mod 5),``      with  71 terms,

       ``w = sum  v94 * g_94**(2*i)  (mod 3),``  with  47 terms.

    The conditions for ``v71, g71, v94, g94`` given above are 
    checked. we also check ``w != 0 (mod 3)`` and 
    ``w != 0 (mod 5)``, ignoring all entries of ``w`` with 
    tag ``A``.

    The vector ``w`` is returned as an element of the space given 
    by class ``MMSpace(15)``.

    The function raises an exception if any of the checks given
    above fails.

    According to  :cite:`LPWW98`, section 7, the vector ``w`` is 
    fixed by the  neutral element of the monster group only. Note 
    that vector ``w`` is obtained by Chinese remaindering two 
    vectors in characteristic ``3`` and ``5``, so that the 
    operation of an element ``g`` of the monster on ``w`` actually 
    corresponds to the operation on the two vectors ``w mod 3`` 
    and ``w mod 5``.

    Inputs ``s_v71, s_g71, s_v94, s_g94`` can be found by function
    ``find_order_vector``. Searching for such inputs may take a 
    long time. Function ``write_order_vector`` writes these data 
    to a python file. 
    """
    if verbose:
         print("Computing order vector from input:")
         print("g71 =", s_g71)
         print("g94 =", s_g94)
         print("v71 =", s_v71)
         print("v94 =", s_v94)
    sp = MMSpace(15)
    group = sp.group
    v71 = sp(s_v71)
    g71 = group(s_g71)
    v94 = sp(s_v94)
    g94 = group(s_g94)

    assert 3*v94 == sp.zero()
    g47 = g94 ** 2
    w94 = v94.copy()
    w = v94.copy()
    for i in range(1, 47):
        w94 *= g47
        w += w94
        if i == 1:
            assert v94 != w94
        if i == 23:
            assert v94 != w94  * g94
    assert v94 == w94 * g47
       
    assert 5*v71 == sp.zero()
    w71 = v71.copy()
    w += v71.copy()
    for i in range(1, 71):
        w71 *= g71
        w += w71
        if i == 1:
            assert v71 != w71
    assert v71 == w71 * g71

    w3 = 5 * w
    assert np.count_nonzero(w3["Z",:500]), "g47"
    assert w3 * g47 == w3

    w5 = 3 * w
    assert np.count_nonzero(w5["Z",:500]), "g71"
    assert w5 * g71 == w5

    if verbose:
         print("Order vector computed successfully")
    return w
   

def find_one_order_vector(verbose = 0):
    """Compute suitable inputs for function ``compute_order_vector``

    The function returns a pair ``(d, v)``.

    Here ``v`` is an order vector as defined in function
    ``compute_order_vector``. ``d`` is an ordered dictionary that 
    maps the strings ``"v71", "g71", "v94", "g94"`` (in that order) 
    to suitable input values for function ``compute_order_vector``.

    When calling function ``compute_order_vector()`` with these 
    inputs, we would obtain the the vector ``v``. 
    
    The function may fail with a small probability, so we have
    to apply it repeatedly.
    """
    v71 =  [(3,"Y",0,0)]
    v94 =  [
        ("Y",0,0),
        ("X",randint(0,2047),randint(0,23)),
        ("Z",randint(0,2047),randint(0,23)),
        ("T",randint(0,758),randint(0,63)),
    ]
    g71 = find_element_of_order(71, 15, v71, verbose=verbose) 
    g94 = find_element_of_order(94, 3, v94, verbose=verbose) 
    sp = MMSpace(15)
    v71 = sp(*v71)
    v94 = 5 * sp(*v94)
    names = ["v71", "g71", "v94", "g94"]
    data = [v71, g71, v94, g94]
    d = OrderedDict(zip(names, map(str,data)))
    v = compute_order_vector(*(d[s] for s in names), verbose = verbose)
    return d, v


def find_order_vector(verbose = 0):
    """Repeated application of function find_one_order_vector()

    Function find_one_order_vector() is repeated until it returns
    a correct result.

    Parameters and return value are as in function 
    find_one_order_vector(). 
    """
    for i in range(100):
        try:
            return find_one_order_vector(verbose)
        except:
            pass
    raise
       

def write_order_vector(d, verbose = 0):
    """Write data returned by function ``find_order_vector`` to file

    Function ``find_order_vector`` returns a pair ``(d, v)`` 
    where ``d`` is a dictionary containing data required for a
    fast recomputation of the *order vector* ``v``.

    This function writes the entries ``key : value`` of 
    dictionary ``d`` to a python file in the form::

      key = value

    We may import this python file for obtaining the data 
    required for recomputing ``v``. 
    """
    if verbose:
        print("Writing order vector to file:\n%s" % FILENAME)
    with open(FILENAME, "wt") as f:
        print("""# This is an automatically generated python file.
# It contains data required by module mmgroup.structures.mm_order.
# Recomputing these data may take a long time!

""", file = f)
        for item in d.items():
             print('%s = "%s"' % item, file = f)
    if verbose:
        print("Order vector written to file.")


###########################################################################
# Obtaining the actual order vector
###########################################################################


def get_order_vector(recompute = False, verbose = 0):
    """Return an *order vector*

    The function returns an order vector suitable for computing
    the order of an element of the monster group. See function
    ``compute_order_vector`` for the definition of an order vector. 

    Data for fast computation of an order vector is read from a 
    file if present. Otherwise this data is recalculated and 
    written to that file.

    The order vector is also stored in memory for reuse.

    The order vector is returned as an element of the vector space
    ``MMSpace(15)``.

    If ``recompute`` is True then the data for computing the order
    vector is always recomputed.
    """
    global order_vector
    if recompute:
        order_vector = None
        try:
            os.remove(FILENAME)
        except:
            pass
    if not order_vector is None:
        return order_vector
    try:
        if verbose:
            print("Load order vector data from file")
        from mmgroup.structures.mm_order_vector import g71, v71, g94, v94
        order_vector = compute_order_vector(v71, g71, v94, g94, verbose)
        return order_vector
    except:
        if verbose:
            "Generate data for order vector"
        data, order_vector = find_order_vector(verbose)
        try:
            write_order_vector(data, verbose)
            return order_vector
        except:
            raise
            warn = """Cannot write vector for calculating group order to file
%s.
Recomputing such a vector may take a long time!!
"""
            warnings.warn(warn % FILENAME, UserWarning)
            return order_vector



###########################################################################
# Computing the order of an element of the monster
###########################################################################


def check_mm_order(g, max_order = 119, mode = 0):
    """Return order of monster group element ``g``.

    if ``order(g) < max_order`` return ``order(g)``; else return ``0``.

    If mode is ``0`` (default) we first check if ``g`` is in the 
    subgroup ``N_0 of`` the monster. If this is the case the we check 
    the order of ``g``  by calculating in ``N_0``.

    Othewise we compute the minimum ``i`` such that
    ``v * g**i == v`` for the *order vector* ``v`. 
    """
    assert isinstance(g, MMGroupWord)
    g.reduce()
    if mode == 0:
        n0 = np.zeros(5, dtype = np.uint32)
        status = mm_group_check_word_n(g._data, g.length, n0)
        if status == 0:
            return 1
        if status == 1:
            n1 = np.copy(n0)
            for i in range(2, max_order+1):
                mm_group_n_mul_element(n1, n0)
                if not mm_group_n_reduce_word(n1):
                    return i
            return 0

    v = get_order_vector().data
    w = mm_vector(15)
    work = mm_vector(15)
    mm_op15_copy(v, w)
    for i in range(1, max_order+1):
        mm_op15_word(w, g._data, g.length, 1, work)
        if not mm_op15_compare(v, w):
            return i
    return 0

    

###########################################################################
# Check equality of two elements of the monster
###########################################################################
 

def check_mm_equal(g1, g2, mode = 0):
    """Return ``g1 == g2`` for elements ``g1, g2`` of the monster.

    If ``mode == 0`` (default) we first try to check equality inside 
    in the subgroup ``N_0 of`` the monster, which may be considerbly 
    faster. 

    If ``mode != 0`` or this is not possible we check if 
    ``v * g**i == v`` holds for an *order vector* ``v`. 

    We just check the data in ``g1`` and ``g1``, ingnoring
    ``g1.group`` and ``g2.group``.
    """
    assert isinstance(g1, MMGroupWord)
    assert isinstance(g2, MMGroupWord)
    g3 = np.zeros(2 * (g1.length + g2.length) + 1, dtype = np.uint32)
    status = mm_group_words_equ(g1._data, g1.length,
        g2._data, g2.length, g3)
    if status < 2:
        return not status

    v = get_order_vector().data
    w = mm_vector(15)
    work = mm_vector(15)
    mm_op15_copy(v, w)
    mm_op15_word(w, g3, status - 2, 1, work)
    return not mm_op15_compare(v, w)


###########################################################################
# Test
###########################################################################



if __name__ == "__main__":
    get_order_vector(recompute = True, verbose = 1)
    order_vector = None
    get_order_vector(recompute = False, verbose = 1)

