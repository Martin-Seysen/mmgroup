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
import sys
from random import randint
from collections import OrderedDict, defaultdict
from multiprocessing import Pool, TimeoutError, Lock

import numpy as np


if __name__ == "__main__":
    sys.path.append(os.path.join("..", "..", ".."))

import mmgroup
from mmgroup import structures, MM0, MMV
from mmgroup.mm_space import MMSpace, MMVector
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.mm_op import mm_aux_index_sparse_to_leech2
from mmgroup.mm_op import mm_aux_mmv_extract_sparse
from mmgroup.mm_op import mm_op_eval_A_rank_mod3 


import_pending = True



def complete_import():
    """Internal function of this module

    If you need any of the objects declared above, do the following:

    if import_pending:
        complete_import()
    """
    global import_pending, mat24
    from mmgroup import mat24
    import_pending = False

    



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
# Check if a symmetric 24 times 24 matrix can be hashed nicely
#######################################################################




def hash_mat24(m, verbose = 0):
    """Compute the hash value described in function ``nicely_hashable``

    The function returns the list ``hlist`` of the hash values of the
    24 rows of the 24 times 24 matrix ``a`` whose entries are integers
    mod 3. The hash value of the row is the number of its nonzero
    entries plus 32 times the diagonal entry.
    """
    hlist = []
    if verbose:
        print("Absolute values of matrix A, diagonal entry, count")
    for i in range(24):
        v = m[i] % 3
        diag = v[i]
        s1 = np.sum(v != 0)
        if verbose:
            print([int(x) for x in v], diag, s1)
        s = 32 * diag +  s1
        hlist.append(int(s))
    return hlist

def hash_unique(hlist):
    """Select unique hash values from result of ``hash_mat24()``

    Let ``hlist`` be the list of the 24 hash values computed by
    function ``hash_mat24``. The function returns a dictionary
    mapping these hash values to the corresponding row indices.
    Duplicate entries occuring in that list are dropped.
    """
    num = defaultdict(int)
    index = {}
    for i, h in enumerate(hlist):
        num[h] += 1
        index[h] = i
    hdict = {}
    for h, n in num.items():
        if n == 1:
            hdict[h] = index[h]
    return hdict



def find_umbral_heptad(ilist):
    """Find umbral heptad in a list ``ilist`` of indices

    Given a list ``ilist`` of indices of a Golay code vectors, the
    function tries to find an umbral heptad contained in that list.
    In case of success it computes an element of the Mathieu group
    M_24 that transforms an umbral heptad in that list to the
    standard umbral heptad [0,1,2,3,4,5,8]. Then it returns an
    integer ``x`` so that ``MM0('a', x)`` performs this permutation.
    In case of failure the function returns 0.
    """
    if import_pending:
        complete_import()
    if len(ilist) < 7:
        return 0;
    if len(ilist) != len(set(ilist)):
        raise ValueError("Internal error") 
    ilist = sorted(ilist)[:9]


    v = sum(1 << i for i in ilist)
    bw = mat24.bw24(v)
    while bw > 7:                      # try to delete bits from v
        v_tl = v                       # until v has weight 7
        while v_tl:
            v_hd = v_tl & -v_tl        # get next bit v_hd of v
            v_dec = v ^ v_hd           # v_dec = v \setminus v_hd
            syn = mat24.syndrome(v_dec, 0)
            if syn & v_dec != syn:     # if syndrom of v_dec is not
               v = v_dec               # a subset of v_dec then
               bw -= 1                 # delete bit v_dec from v
               break
            v_tl ^= hd
        if v_tl == 0:                  # abort if no bits v_hd remain
            return 0
    assert mat24.bw24(v) == 7          # v should be an umbral heptad
    special = v & mat24.syndrome(v, 0) # special entry of heptad v
    v &= ~special                      # the other six entries of v
    # Let heptad be the selected heptad as a list of indices
    heptad = mat24.vect_to_list(v, 6) + mat24.vect_to_list(special, 1)
    assert len(heptad) == 7
    pi = mat24.perm_from_heptads(heptad, [0,1,2,3,4,5,8])
    mat24.perm_check(pi)
    a =  mat24.perm_to_m24num(pi)
    return 0x20000000 + a
       

          

def nicely_hashable(a, verbose = 0):
    """Checks if a nice hash value of a matrix can be computed

    Here ``a`` is a 24 times 24 matrix of integers modulo 3 that
    corresponds to part 'A' of of a vector in the Griess alebra.
    For row ``i`` of ``a`` we compute the hash value:

    h[i] = a[i,i] + (number of nonzero entries in row a[i]) .

    So the h[i] are invariant under the subgroup 2^(1+24+11)
    acting on part 'A' of the Griess algebra. They are permuted
    according to the factor group M_24 of the group

    N_x0 = 2^(1+24+11).Mat_24

    acting on 'A'. So the h[i] can be used to detect tha action of
    the factor group Mat_24, privided that the are sufficiently
    disjoint. More specifically, there must be an *umbral heptad*,
    which is a subset of size 7 of the set {0,...,23} of which M_24
    acts, that is not contained in an octad.

    The hash value for any row in that umbral heptad must be unique
    among the hash values of all rows.

    The function returns ``True`` if a suitable heptad can be found.
    """
    #return 1
    hlist = hash_mat24(a)
    hdict = hash_unique(hlist)
    if verbose:
        print("No of good hashes:", len(hdict), ", list of hashes:")
        print(sorted(hlist))
        print("Matrix A:", "kernel diagonal")
        print(a)
    return find_umbral_heptad(hdict.values())



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

    We compute an element ``g`` of the Moster that acts as an element
    of given ``order`` on a random vector ``v`` in characteristic
    ``p``.  So ``g`` has that order with a very high probability.
    But, theoretically, the order of ``g`` may be
    a multiple of  ``order``.

    ``g`` is returned as a string in standard mmgroup notation.

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


####################################################################
## Compute entries of order vector s59 for finding 2-group
####################################################################


def map_y(y_index):
    if import_pending:
        complete_import()
    i, j = (y_index >> 14) & 0x1f, (y_index >> 8) & 0x1f
    vect = (1 << i) + (1 << j)
    gc = mat24.vect_to_cocode(vect)
    assert 0 <= gc < 0x800
    return gc 
    
    
def map_x(x_index):
    v2 = mm_aux_index_sparse_to_leech2(x_index) 
    return ((v2 & 0xfff) << 12) | ((v2 >> 12) & 0xfff)    


def eqn_system(vector, tag_indices, map_f, n, error_msg):
    entries = [vector[index] for index in tag_indices]
    indices = [MMSpace.index_to_sparse(*x) for x in tag_indices]
    matrix= np.zeros(24, dtype = np.uint64)
    rows, cols = 0, n
    out_indices = []
    for (entry, index) in zip(entries, indices):
        if entry == 1:
            eqn = map_f(index)
            new_rows = leech2matrix_add_eqn(matrix, rows, cols, eqn)
            if new_rows:
                out_indices.append(index)
                rows += new_rows
                if rows == n:
                    return out_indices
    raise ValueError("Check %s failed, %d rows found" % (error_msg, rows)) 


def eqn_sign(vector):
    for i in range(100):
        for j in range(24):
            if vector["Z", i, j] == 1:
                return [ MMSpace.index_to_sparse("Z", i, j) ]
    raise ValueError("Check for sign failed") 


def check_v(v):
    r"""Return list of entries uses for obtaining 2-subgroup of G_x0

    Let math:`Q_{x0}` be the subgroup  of structure math:`2^{1+24}`
    of math:`N_{x0}` (of structure math:`2^{1+24+11}.M_{24}`) as
    usual; and let math:`N_2` be the subgroup of math:`N_{x0}` of
    structure math:`2^{1+24+11}`.

    Here input ``v`` must be a vector in the rep of the Monster mod 3.
    The function returns a list ``sp`` of 11+24+1 coordinate positions
    of ``v`` that can be used to determine a ``g`` in math:`N_2` from 
    ``v * g``. The returned list of  coordinate positions is given in
    sparse notation. vector ``v`` has coordinate 1 at all these
    positions.

    The first 11 positions in ``sp`` are at part 'A' of the vector.
    They are independent in the sense that they uniquely determine
    an ``y_d``, with ``d`` in the Parker loop (modulo its centre),
    such that ``y_d`` transforms any ``v * g`` to  ``v * h``,
    for ``g`` in  math:`N_2` and ``h`` math:`Q_{x0}`.

    The next 24 positions in ``sp`` are at parts 'BCTX' of the
    vector. They are independent in the sense that they uniquely
    determine an element ``r`` of the Leech lattice mod 2, such that
    ``x_r`` transforms any ``v * g`` to  ``v * :math:`x_{\pm 1}`,
    for ``g`` in  math:`Q_{x0}`.

    Finally, the last entry of ``sp`` is a coordinate position in
    part 'ZY' that can be used to  obtain the sign
    in  ``v1`` * :math:`x_{\pm 1}`.
    """
    Y_INDICES = [("A",i,j) for i in range(10) for j in range(i+1, 24)]
    X_INDICES = [("B",i,j) for i in range(10) for j in range(i+1, 24)]
    X_INDICES += [("C",0,j)  for j in range(1, 24)]
    BASIS = [0] + [1 << i for i in range(11)]
    X_INDICES += [("X",i,j) for j in range(0, 24)  for i in  BASIS]
    result_y = eqn_system(v, Y_INDICES, map_y, 11, "Y")
    result_x = eqn_system(v, X_INDICES, map_x, 24, "X")
    result_sign = eqn_sign(v)
    result = result_y + result_x + result_sign
    return [int(x) for x in result]
    



#######################################################################
# Obtaining a vector mod 3 stable under an element of order 71
#######################################################################


    
def gA_from_type4(v_type4):
    r"""Compute a group element reducing a type-4 Leech lattice vector 
    
    Let ``v_type_4`` be a vector of type 4 in the Leech lattice mod 2
    in **Leech lattice encoding**.  The function returns an element
    ``g`` of the subgroup :math:`G_{x0}` of the monster that maps
    `v_type_4`` to the standard frame :math:`\Omega`. The group
    element ``g`` is returned as an instance of class ``MM0``.   
    """
    g_data = np.zeros(10, dtype = np.uint32)
    len_g = gen_leech2_reduce_type4(v_type4, g_data)
    return MM0('a', g_data[:len_g])

def make_v_p_sample(p, g_p):
    r"""Compute a vector stabilized by a group element of order p
    
    Let ``g_p`` be an element of the monster group of order ``p``. The
    function generates  a random vector ``v_p`` in the representation
    of the monster modulo 3 and computes the vector
    ``w = sum(v_p * g_p**i for i in range(p))``. ``w`` is
    stabilized by ``g_p``.
    
    The function returns a triple containing ``v_p`` if ``w`` is
    non-trivial and ``w`` satisfies an additional property
    described below. Otherwise the function returns ``None``. The 
    function succeeds with small probability; so it must be 
    repeated several times.
    
    The part of the vector ``w`` labelled with the tag ``A`` 
    corresponds to a symmetric bilinear form over the Leech lattice
    modulo 3. The function succeeds if that bilinear form corresponding
    to ``w`` has a one-dimensional eigenspace with the eigenvalue 0
    containing a type-4 eigenvector in the Leech lattice.
    
    In case of success we return the triple ``(v_p, gA)``. Here
    ``gA`` is an element  of the subgroup :math:`G_{x0}` of the 
    monster group such that the bilinear form on the Leech lattice 
    corresponding  to `w * gA`` has a type-4 eigenvector in the 
    standard frame math:`\Omega`. 

    The values ``v_p`` and ``gA`` are returned as strings.
    """
    r1 = [("s", x, "r") for x in "XTXTX"]
    g_p = MM0(g_p)
    v_p = MMV3(r1)
    w_p = stabilizer_vector(v_p, g_p, p)
    if not w_p:
        return None
    v3 = mm_op_eval_A_rank_mod3(3, w_p.data, 0) & 0xffffffffffff
    if v3:
         v_type4 = gen_leech3to2_type4(v3)
         if v_type4:
            gA = gA_from_type4(v_type4)      
            w_p_all = w_p * gA
            w_p_A = w_p_all['A']
            perm_num = nicely_hashable(w_p_A)            
            if perm_num > 0:
                try:
                    w_p_all = w_p_all * MM0('a', [perm_num])
                    yx = check_v(w_p_all)
                except:
                    #raise
                    return None
                gA1 = gA * MM0('a', [perm_num])
                return str(v_p), str(gA1), yx
    return None 


def find_vector_p_mod3(p, verbose = 0):
    r"""Compute a vector stabilized by a group element of order p

    The function computes an element ``g_p`` of order ``p`` of the
    monster and a vector ``v_p``  in the representation of the
    monster modulo 3, such that the vector
    ``w0 = sum(v_p * g_p**i for i in range(p))`` is not trivial.
    Then ``w0`` is stabilized by ``g_p``. More precisely, the
    projection of ``w0`` onto the 196883-dimensional irreducible
    representation of the Monster is not trivial.

    The function also computes an element ``gA`` of the subgroup 
    :math:`G_{x0}` such that the vector ``w = w0 * gA`` has the
    following additional properties. 
    
    The part ``w['A']`` of a vector ``w`` labelled with tag ``A``
    is always a symmetric bilinear form over the Leech lattice
    modulo 3. The bilinear form ``w['A']`` has a one-dimensional
    kernel spanned by a basis vector of the Leech lattice.

    Furthermore, the hash values for the rows of matrix ``w['A']``
    computed by function ``hash_mat24`` are unique among hash
    values for rows 0,1,2,3,4,5,8.

    Such a vector ``w`` satisfies the requirements for the
    vector :math:`v_{_p} \in \rho_3` in :cite:`Seysen22`.
    
    The function returns the tuple ``(g_p, v_p, gA, yx)``.

    Here ``yx`` is a list of entries of vector ``w`` (in sparse
    notation) that can be used to idetify a ``g`` in math:`G_{x0}`
    from ``v * v``, see function ``check_v``.
    """
    if import_pending:
        complete_import()
    s_g_p = find_element_of_order(p, verbose = 1)
    assert isinstance(s_g_p, str)
    g_p = MM0(s_g_p)
    if verbose:
        print("Find vector stabilized by g_p")
    s_v_p, s_gA, yx = par_search(2.0e5, make_v_p_sample, p, s_g_p,
        verbose = 1, chunksize = 300)
    if verbose:
        print("g_p =", s_g_p)
        print("v_p =", s_v_p, "# (mod 3)")
        print("gA =", s_gA)
    return s_g_p, s_v_p, s_gA, yx



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


def assemble_vector_mod3(p, s_g, s_v, s_gA):
    """Compute order vector from output of ``find_vector_p_mod3``

    The function computes the order vector ``v`` from the input ``p``
    passed to function ``find_vector_p_mod3`` and the results
    ``s_g, s_v, s_gA`` returnd by that function. It returns ``v`` if
    ``v`` has the properties stated in function ``find_vector_p_mod3``
    and fails if this is not the case.
    """
    g = MM0(s_g) 
    v = MMV3(s_v)
    v1 = MMV3(s_v)
    for i in range(1, p):
        v1 *= g
        v += v1
    assert v1 * g == MMV3(s_v)
    v *= MM0(s_gA)
    v3 = mm_op_eval_A_rank_mod3(3, v.data, 0)
    assert v3 >> 48 == 23
    v_type4 = gen_leech3to2_type4(v3)
    assert v_type4 == 0x800000
    hashes = hash_mat24(v["A"])
    num_h = defaultdict(int)
    for h in hashes:
        num_h[h] += 1
    for i in [0,1,2,3,4,5,8]:
        assert num_h[hashes[i]] == 1
    return v

if __name__ == "__main__":
    s_g71, s_v71, s_gA, yx = find_vector_p_mod3(71, verbose = 0)
    print("Group element g of order 71 is:\n" + s_g71)
    print("Seed vector (mod3) is:\n" + s_v71)
    print("Group element for adjusting eigenvector is:\n" + s_gA)
    v = assemble_vector_mod3(71, s_g71, s_v71, s_gA)
    print("Part 'A' of order vector (mod 3)")
    print(v['A']) 
    print("Hash values of rows of part 'A'")
    hashes = hash_mat24(v['A'], verbose = 0)
    for i in range(24):
       print(" %02d" % hashes[i], end = " " * (i==7))
    vhex = np.vectorize(hex)
    print("\nyx  =\n%s" % vhex(np.array(yx)))
    print()














