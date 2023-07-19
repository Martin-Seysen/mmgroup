r"""This module computes an order vector in a representation of the monster.

In :cite:`Seysen22` we define a vector :math:`v_1` in the 
representation  :math:`\rho_{15}` of the monster group 
:math:`\mathbb{M}` that is used for recognizing elements of the 
subgroup :math:`G_{x0}` of :math:`\mathbb{M}`. Vector :math:`v_1`
is also used for computing the order of an element of 
:math:`\mathbb{M}`, so we call  :math:`v_1` an **order vector**
here.

:math:`v_1` is constructed from two vectors :math:`v_{71} \in \rho_3`
and :math:`v_{94} \in \rho_5` via Chinese remaindering. This module
contains functions for computing vectors  :math:`v_{71}` and
:math:`v_{94}` with the properties required in :cite:`Seysen22`.

These computations are the most time-consuming part of the
computation of :math:`v_1`; and we use multiprocessing for 
performing these computations. 

All computations in the monster group are done with instances 
of class ``MM0``, but not ``MM``. The reason for this is that 
class ``MM`` requires the existence of an order vector.
"""

import os
from random import randint
from collections import OrderedDict
from multiprocessing import Pool, TimeoutError, Lock

import numpy as np

import mmgroup
from mmgroup import structures, MM0, MMV
from mmgroup.mm_space import MMSpace, MMVector
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
from mmgroup.mm_op import mm_aux_mmv_extract_sparse
from mmgroup.mm_op import mm_op_eval_A_rank_mod3 



#######################################################################
# Parallel search function
#######################################################################


def _repeat_f(repetitions, f, args):
    """Auxiliary function for function par_search()
    
    This function executes ``f(*args)`` repeatedly (at most
    ``repetitions`` times) until the return value is a result 
    different from None.  It returns the pair ``(result, n)``, 
    where ``result`` is the result of the last execution of 
    ``f``, and ``n`` is the number of executions of `f``.
    
    """
    for i in range(repetitions):
        result = f(*args)
        if result is not None:
            return i+1, result
    return repetitions, None         


def par_search(trials, f, *args, **kwds):
    """Repeat executing ``f(*args)`` until execution is successful 
    
    Function ``par_search`` executes ``f(*args)`` repeatedly until 
    a result different from ``None`` is returned. Then function 
    ``par_search`` returns that result.
    
    Execution of ``f`` is repeated (at least) ``trials`` times.
    If ``f`` returns ``None`` in all cases then ``ValueError`` is
    raised. 
    
    Execution of function ``f`` runs in parallel. The number of
    parallel processes is equal to ``os.cpu_count()``. One can 
    use the keyword argument ``processes`` to adjust the number 
    of processes.
    
    In each round function ``f`` is executed ``chunksize`` times,
    where ``chunksize`` can be given by a keyword argument. By
    default we make an educated guess for a good chunk size.
    If ``chunksize`` is too small then there might be too much 
    overhead for setting up the parallelization. If ``chunksize``
    is too large then one must always wait for at least
    ``chunksize`` executions of ``f``, even if the probability
    of success is much larger than ``1 / chunksize``.    
    
    Function ``f`` usually calls a random generator for obtaining
    test values. That random generator should return different 
    random  data when running in different processes. Here a 
    random generator like ``random.randint()`` is ok.

    Preferably, the arguments in ``args`` should be simple objects
    like numbers or strings or lists of such simple objects.  You
    must be aware of the fact that each process works on its own
    copy of the input data.

    If the keyword argument ``verbose`` is True then the 
    function displays an indication of its activity.    
    """
    import math
    trials = max(1, math.ceil(trials))
    processes = kwds.get("processes", 0)
    if not processes:
        processes = os.cpu_count()
    else:
        processes = max(1, min(processes, os.cpu_count()))
    chunksize = kwds.get("chunksize", 0)
    if not chunksize:
        chunksize = max(1, math.ceil(math.sqrt(trials)))
    chunksize = max(1, min(chunksize, trials))
    repetitions = math.ceil(chunksize / processes)
    chunksize = repetitions * processes
    n_iter = math.ceil(trials / chunksize)
    n = 0
    verbose = bool(kwds.get("verbose", False))
    arg_list = [(repetitions, f, args)] * processes
    for i in range(n_iter):        
        with Pool(processes = processes) as pool:
            results = pool.starmap(_repeat_f, arg_list, chunksize = 1)
        pool.join()
        for i, result in results:
            if result is not None:
                n += processes * i
                s = "\rSuitable candidate found after about %d tests"
                if verbose:
                    print(s % n)
                return result
        n += chunksize
        if verbose:
            print("\r%5d candidates tested   " % n, end = "")
    if verbose:
        print("\n")
    err = "No suitable candidate found"
    raise ValueError(err) 


#######################################################################
# Find an element of the monster of a given order
#######################################################################

MMV3 = MMV(3)
MMV15 = MMV(15)


def get_factors(n):
    """Return sorted list of factors of n, for n <= 120"""
    assert 1 <= n < 121, "Number too large"
    return [i for i in range(1, n+1) if n % i == 0]


def rand_mm_element(size = 1):
    """Return a random element ``g`` of the monster group

    Here ``size`` is a measure for the requested word length. 
    of ``g``. The returned random element is not uniform
    distributed. It is just sufficiently random, so that
    elements of a high order occur sufficiently often

    ``g`` is returned as an instance of class ``MM0``.
    """    
    data = []
    for i in range(size):
        data += [('t','n'), ('l','n')]
        data += [(s,'r') for s in "yxdp"]
    return MM0(data)

    
def rand_elem_of_order(factors, elem_size):   
    r"""Return random element of the monster group of a given order

    The function creates a random element ``g`` of the monster by 
    calling ``g = rand_mm_element(elem_size)``. It returns ``g`` as 
    a string if ``g`` has order ``o``, where ``o`` is  a multiple of 
    ``factors[-1]`` and ``None`` otherwise. Here ``factors`` must be
    the list of all factors of the value ``factors[-1]`` in natural
    order.

    The order is checked by applying a random vector ``v`` in the 
    representation :math:`\rho_3` and by multiplying ``v`` with 
    powers of ``g``. 
    """    
    g = rand_mm_element(elem_size)
    v = MMV3('R')
    w = v.copy()
    for i in range(1, len(factors)-1):
        w.mul_exp(g, factors[i] - factors[i-1])
        if v == w:
            return None
    w.mul_exp(g,  factors[-1] - factors[-2])
    return str(g) if v == w else None
    
    

def find_element_of_order(order, minsize = 1, verbose = 0):
    """Return an element ``g`` the monster of order ``order``.

    We ``g`` that check acts as an element of given ``order`` on a
    vector ``v`` in characteristic ``p``.  So ``g`` has that order
    with a very high probability, if ``v`` is selected at random.
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
        print("Find element g%s of order %s in monster group" 
            % (order, order))
    factors = [0] +  get_factors(order) # [order // x for x in rfactors] + [order]
    if verbose:
        print("Factors =" , factors[1:])
    maxsize = minsize + 10
    for elem_size in range(minsize, maxsize):
        try:
            g = par_search(1000, rand_elem_of_order, factors, 
                   elem_size, chunksize = max(70, order),
                   verbose = verbose
                )
            if verbose:
                 print("Found g%d =" % order, g)
            return g
        except ValueError:
            continue
    err = "No element of order %s found in the monster group"
    raise ValueError(err % order)
 

#######################################################################
# Find a vector stabilizing by an element of order n
#######################################################################

def stabilizer_vector(v, g, n):
    """Compute a vector stabilized by an element of the monster

    Le ``g`` be an element of the monster group of order ``n`` and 
    ``v`` a vector in a representation of the monster. We return the
    vector ``sum(v * g**i for i  in range(n))`` which is stabilized
    by ``g``. We always return ``None`` if that sum is 0 or a 
    multiple of the 1 element in the representation space. The 
    last condition is checked with a fast crude check only.  
    """
    vg = v.copy()
    w = v.copy()
    for i in range(1, n):
        vg *= g 
        w += vg
    assert v == vg * g
    if (w['B'] == 0).all():
        return None
    return w

#######################################################################
# Obtaining a vector mod 3 stable under an element of order 71
#######################################################################




    
def gA_from_type4(v_type4):
    r"""Compute a group element reducing a type-4 Leech lattice vector 
    
    Let ``v_type_4`` be a vector of type 4 in the Leech lattice mod 2
    in **Leech lattice encoding**.  The function returns an element
    ``g`` of the subgroup :math:`G_{x0}` of the monster that maps
    `v_type_4`` to the standard frame :math:`\Omega`. The group
    element ``g`` is returned as a string.   
    """
    g_data = np.zeros(10, dtype = np.uint32)
    len_g = gen_leech2_reduce_type4(v_type4, g_data)
    gA = MM0('a', g_data[:len_g])
    return str(gA)

def make_v71_sample(g71):
    r"""Compute a vector stabilized by a group element of order 71
    
    Let ``g71`` be an element of the monster group of order 71. The
    function generates  a random vector ``v71`` in the representation
    of the monster modulo 3 and computes the vector
    ``w = sum(v71 * g71**i for i in range(71))``. ``w`` is
    stabilized by ``g71``. 
    
    The function returns triple containing ``v71`` if ``w`` is 
    non-trivial and ``w`` satisfies an additional property
    described below. Otherwise the function returns ``None``. The 
    function succeeds with probability about 1/700; so it must be 
    repeated several times.
    
    The part of the vector ``w`` labelled with the tag ``A`` 
    corresponds to a symmetric bilinear form over the Leech lattice
    modulo 3. The function succeeds if that bilinear form corresponding
    to ``w`` has a one-dimensional eigenspace with the eigenvalue 
    ``diag`` containing a type-4 eigenvector in the Leech lattice. 
    
    In case of success we return the triple ``(v71, gA, diag)``. Here 
    ``gA`` is an element  of the subgroup :math:`G_{x0}` of the 
    monster group such that the bilinear form on the Leech lattice 
    corresponding  to `w * gA`` has a type-4 eigenvector in the 
    standard frame math:`\Omega`. 

    The values ``v71`` and ``gA`` are returned as strings. ``diag``
    is returned as an integer.        
    """
    r1 = [("s", x, "r") for x in "XTXTX"]
    g71 = MM0(g71)
    v71 = MMV3(r1)
    w71 = stabilizer_vector(v71, g71, 71)
    if not w71:
        return None
    for diag in range(3):
        v3 = mm_op_eval_A_rank_mod3(3, w71.data, diag) & 0xffffffffffff
        if v3:
            v_type4 = gen_leech3to2_type4(v3)
            if v_type4:
                return str(v71), gA_from_type4(v_type4), diag               
    return None 


def find_vector_71_mod3(verbose = 0):
    r"""Compute a vector stabilized by a group element of order 71

    The function computes an element ``g71`` of order 71 of the
    monster and a vector ``v71``  in the representation of the 
    monster modulo 3, such that the vector
    ``w = sum(v71 * g71**i for i in range(71))`` is not trivial. 
    Then ``w`` is stabilized by ``g71``. More precisely, the 
    projection of ``w`` onto the 196883-dimensional irreducible
    representation of the monster is not trivial.
    
    The vector ``w`` satisfies the following additional property:
    
    The part of a vector ``w`` labelled with the tag ``A`` always
    corresponds to a symmetric bilinear form over the Leech lattice
    modulo 3. The bilinear form corresponding to the computed vector
    ``w`` has a one-dimensional eigenspace with the eigenvalue 
    ``diag`` containing a type-4 eigenvector in the Leech lattice.

    Such a vector ``w`` satisfies the requirements for the
    vector :math:`v_{71} \in \rho_3` in :cite:`Seysen22`. 
    
    The function also computes an element ``gA`` of the subgroup 
    :math:`G_{x0}` of the monster  group such that the bilinear 
    form on the Leech lattice  corresponding  to ``w * gA`` has a 
    type-4 eigenvector in the standard frame math:`\Omega`. 
    
    The function returns the tuple ``(g71, v71, gA, diag)``. 
    Here ``diag`` is returned as an integer; the other components  
    are returned as strings.
    """    
    s_g71 = find_element_of_order(71, verbose = 1)
    assert isinstance(s_g71, str)
    g71 = MM0(s_g71)
    if verbose:
        print("Find vector stabilized by g71") 
    s_v71, s_gA, diag = par_search(2.0e5, make_v71_sample, s_g71, 
        verbose = 1, chunksize = 300)
    if verbose:
        print("g71 =", s_g71)
        print("v71 =", s_v71, "# (mod 3)")
        print("gA =", s_gA)
        print("kernel_diagonal =", diag, "# (mod 3)" )
    return s_g71, s_v71, s_gA, diag


#######################################################################
# Obtaining a vector mod 15 stable under an element of order 94
#######################################################################




def make_v94_sample(s_g94):
    """Auxiliary function for function ``find_vector_94_mod5``
 
    Given ab element ``g94`` of order 94 in the monster, the function
    computes a vector ``v94``  in the representation of the  monster 
    modulo 15, such that the vector
    ``w = 3 * sum((-1)**i * v94 * g94**i for i in range(94))``
    is not trivial. Then `w`` is stabilized by ``g94**2``, but not 
    by ``g94``.
    
    The function returns the vector ``v94`` in case of  success
    and ``None in case of failure``.    
    """
    g94 = MM0(s_g94)
    r1 = [("n", x, "r") for x in "XTYTZ"]
    v94 = 3 *  MMV15(r1)
    w = stabilizer_vector(v94 - v94 * g94, g94**2, 47)
    return None if w is None else str(v94)



def do_find_vector_v94_mod5(s_g94, verbose = 0):
    """Compute a vector stabilized by a group element of order 47
 
    The function computes and element ``g94`` of order 94 in the
    monster and a vector ``v94``  in the representation of the
    monster modulo 15, such that the vector
    ``w = 3 * sum((-1)**i * v94 * g94**i for i in range(94))``
    is not trivial. Then `w`` is stabilized by ``g94**2``, but not 
    by ``g94``.
    
    The function returns the pair ``(g94, v94)`` as a pair of 
    strings.    
    """
    if verbose:
        print("Find vector v94 stabilized by g94**2, but not by g94")
    try:
        s_v94 = par_search(2, make_v94_sample, s_g94,
            chunksize = 1, verbose = verbose)
        if verbose:
            print("v94 =", s_v94, "# (mod 15)")
        return  s_v94
    except ValueError:
        if verbose:
             print("No suitable candidate for v94 found")
        return None




def find_vector_v94_mod5(verbose = 0):
    r"""Compute a vector stabilized by a group element of order 94

    The function computes an element ``g94`` of order 94 of the
    monster,  and a (sparse) vector ``v94``  in the representation 
    of the monster modulo 15, such that the vector
    ``w = 3 * sum((-1)**i * v94 * g94**i for i in range(94))``
    is not trivial. Then `w`` is stabilized by ``g94**2``, but not 
    by ``g94``.
    
    Such a vector ``w`` satisfies the requirements for the
    vector :math:`v_{94} \in \rho_5` in :cite:`Seysen22`. 

    The function returns the the pair ``(g94, v94)`` in case of 
    success and a pair ``(None, None)`` in case of failure.  In 
    of success, both objects, ``g94`` and ``v94``, are returned
    as strings. The function fails with a very small probability.  
    """
    s_g94 = find_element_of_order(94, verbose = verbose)
    if s_g94 is None:
        return None, None
    s_v94 = do_find_vector_v94_mod5(s_g94, verbose = verbose)
    return s_g94, s_v94



















