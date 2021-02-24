import os
from random import randint
from collections import OrderedDict
from multiprocessing import Pool, TimeoutError, Lock

import numpy as np

import mmgroup
from mmgroup import structures
from mmgroup.mat24 import vect_to_cocode
from mmgroup.mm_space import MMSpace
from mmgroup.generators import gen_leech3to2_type4
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.clifford12 import leech3matrix_kernel_vector
from mmgroup.clifford12 import leech3matrix_watermark
from mmgroup.clifford12 import leech2matrix_add_eqn
from mmgroup.mm import mm_aux_index_sparse_to_leech2
from mmgroup.mm import mm_aux_mmv_extract_sparse
from mmgroup.mm_order import stabilizer_vector
from mmgroup.mm_order import make_order_vector
from mmgroup.mm_order import map_y, map_x

_DIR = os.path.split(structures.__file__)[0]
PY_FILENAME = os.path.join(_DIR, "order_vector_data.py")


#######################################################################
# Parallel search function
#######################################################################


def _repeat_f(repetitions, f, args):
    """Auxilieary function for function par_search()
    
    This function executes ``f(*args)`` repeatedly (at most
    ``repetitions`` times) until the return value is a result 
    differnt from None.  It returns the pair ``(result, n)``, 
    where ``result`` is the result of the last execution of 
    ``f``, and ``n`` is the number of ececutions of `f``.
    
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

MMV3 = MMSpace(3)
MMV15 = MMSpace(15)
MM = MMV3.group
assert  MMV15.group == MM

def get_factors(n):
    """Return sorted list of prime factors of n, for n <= 120"""
    assert 1 <= n < 121, "Number too large"
    return [i for i in range(1, n+1) if n % i == 0]


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
    
def rand_elem_of_order(factors, elem_size):   
    g = rand_mm_element(elem_size)
    v = MMV3.rand_vector()
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
# Some precomputed stuff
#######################################################################

USE_PRECOMPUTED = False


if USE_PRECOMPUTED:
    S_g71 = "M<t_1*l_2*y_75h*x_8e1h*d_27eh*p_93090391>"
    S_v71 = "v<-T_695_3bh+T_705_0ch-X_283h_18-X_31ch_6+X_548h_21>" # (mod 3)
    S_gA = "M<(1/p_205578880)*(1/l_1)*(1/p_139991040)*(1/l_2)*(1/p_1392)*(1/l_2)>"
    S_diag = 2



#######################################################################
# Obtaining a vector mod 3 stable under an element of order 71
#######################################################################




    
def gA_from_type4(v_type4):
    """Comupute a group element reducing a type-4 Leech lattice vector 
    
    Let ``v_type_4`` be a vector of type 4 in the Leech lattice mod 2
    in **Leech lattice encoding**.  The function returns an element
    ``g`` of the subgroup :math:`G_{x0}` of the monster that maps
    `v_type_4`` to the standard frame :math:`\Omega`. The group
    element ``g`` is returned as a string.   
    """
    g_data = np.zeros(10, dtype = np.uint32)
    len_g = gen_leech2_reduce_type4(v_type4, g_data)
    gA = MM.from_data(g_data[:len_g])
    return str(gA)

def make_v71_sample(g71):
    """Compute a vector stabilized by a group element of order 71
    
    Let ``g71`` be an element of the monster group of order 71. The
    function generates  a random vector ``v71`` in the representatoon
    of the monster modulo 3 and computes the vector
    ``w = sum(v71 * g71**i for i in range(71))``. ``w`` is
    stabilized by ``g71``. 
    
    The function returns triple containing ``v71`` if ``w`` is 
    non-trivial and ``w`` satisfies an additional property
    described below. Ohterwise the function returns ``None``. The 
    function succeeds with probability about 1/700; so it must be 
    repeated several times.
    
    The part of the vector ``w`` labelled with the tag ``A`` 
    corresponds to a symmetric bilinear form over the Leech lattice
    modulo 3. The function succeeds if that bilinear form correponding
    to ``w`` has a one-dimesional eigenspace with the eigenvalue 
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
    v71 = MMV3(*r1)
    w71 = stabilizer_vector(v71, g71, 71)
    if not w71:
        return None
    for diag in range(3):
        v3 = leech3matrix_kernel_vector(3, w71.data, diag)
        if v3:
            v_type4 = gen_leech3to2_type4(v3)
            if v_type4:
                return str(v71), gA_from_type4(v_type4), diag               
    return None 


def find_vector_71_mod3(verbose = 0):
    """Compute a vector stabilized by a group element of order 71

    The function computes an element ``g71`` of order 71 of the
    monster and a vector ``v71``  in the representatoon of the 
    monster modulo 3, such that the vector
    ``w = sum(v71 * g71**i for i in range(71))`` is not trivial. 
    Then ``w`` is stabilized by ``g71``. 
    
    The vector ``w`` satisfies the following additional property:
    
    The part of a vector ``w`` labelled with the tag ``A`` always
    corresponds to a symmetric bilinear form over the Leech lattice
    modulo 3. The bilinear form correponding to the computed vector
    ``w`` has a one-dimesional eigenspace with the eigenvalue 
    ``diag`` containing a type-4 eigenvector in the Leech lattice. 
    
    The function also computes an element ``gA`` of the subgroup 
    :math:`G_{x0}` of the monster  group such that the bilinear 
    form on the Leech lattice  corresponding  to ``w * gA`` has a 
    type-4 eigenvector in the standard frame math:`\Omega`. 
    
    The function returns the tuple ``(g71, v71, gA, diag)``. 
    Here ``diag`` is returned as an integer; the other components  
    are returned as strings.
    """    
    s_g71 = find_element_of_order(71, verbose = 1)
    g71 = MM(s_g71)
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
    computea a vector ``v94``  in the representatoon of the  monster 
    modulo 15, such that the vector
    ``w = 3 * sum((-1)**i * v94 * g94**i for i in range(94))``
    is not trivial. Then `w`` is stabilized by ``g94**2``, but not 
    by ``g94``.
    
    The function returns the vector ``v94`` in case of  success
    and ``None in case of failure``.    
    """
    g94 = MM(s_g94)
    r1 = [("n", x, "r") for x in "XTYTZ"]
    v94 = 3 *  MMV15(*r1)
    w = stabilizer_vector(v94 - v94 * g94, g94**2, 47)
    return None if w is None else str(v94)



def find_vector_v94_mod5(s_g94, verbose = 0):
    """Compute a vector stabilized by a group element of order 47
 
    The function computes and element ``g94`` of order 94 in the
    monster and a vector ``v94``  in the representatoon of the
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



#######################################################################
# Check that the test vector supports reduction
#######################################################################

Y_INDICES = [("A", i, j) for i in range(2) for j in range(i+1, 24)]
X_INDICES = [("B", i, j) for i in range(2) for j in range(i+1, 24)]
X_INDICES += [("C", 0, j)  for j in range(1, 24)]
BASIS = [0] + [1 << i for i in range(11)]
X_INDICES += [("X", i, j) for j in range(0, 24)  for i in  BASIS]
del BASIS

    
def eqn_system(vector, tag_indices, map_f, n):
    entries = [vector[index] for index in tag_indices]
    indices = [MMV15.index_to_sparse(*x) for x in tag_indices]
    matrix= np.zeros(24, dtype = np.uint64)
    rows, cols = 0, n
    out_indices = []
    for (entry, index) in zip(entries, indices):
        if 1 <= entry < 15:
            eqn = map_f(index)
            new_rows = leech2matrix_add_eqn(matrix, rows, cols, eqn)
            if new_rows:
                out_indices.append(index)
                rows += new_rows
                if rows == n:
                    break
    if rows < n:
        return None
    mask = (1 << n) - 1
    out_matrix = [int(matrix[i]) & mask for i in range(n)] 
    return  out_indices   


def eqn_sign(vector):
    for i in range(2):
        for j in range(24):
            if 1 <= vector["Z", i, j] < 15:
                return [ MMV15.index_to_sparse("Z", i, j) ]
    return None


def augment_v_data(v, data):
    a = np.array(data, dtype = np.uint32)
    mm_aux_mmv_extract_sparse(15, v.data, a, len(a))
    return [x for x in a]


def check_v(v, verbose = 0):
    vB = v['B']
    for p in (3,5):
        if (vB % p == 0).all():
            if verbose:
                print("Vector may be zero (mod %d)" % p)
            return None
    mark = np.zeros(24, dtype = np.uint32)
    if leech3matrix_watermark(15, v.data, mark) < 0:
        if verbose:
            print("Permutation watermarking failed")
        return None
    result_y = eqn_system(v, Y_INDICES, map_y, 11)
    if result_y is None:
        if verbose:
            print("Check Y failed")
        return result
    result_x = eqn_system(v, X_INDICES, map_x, 24)
    if result_x is None:
        if verbose:
            print("Check X failed")
        return result
    result_sign = eqn_sign(v)
    if result_sign is None:
        if verbose:
            print("Check for sign failed")
        return result
    results = [result_y, result_x, result_sign]
    return tuple([augment_v_data(v, data) for data in results])
    



#######################################################################
# Seach for the relevant data
#######################################################################

def str_data(text, data):
    if isinstance(data, list):
        s = "%s = [\n   " % text
        for i, x in enumerate(data):
            s += hex(x) + ","
            s += "\n   " if i % 6 == 5 else " "
        s += "\n]\n"
    elif isinstance(data, str):
        s = '%s = \"%s\"\n' % (text, data)
    elif isinstance(data, int):
        s = '%s = %d\n' % (text, data)        
    else:
        raise TypeError("type " + str(type(data)))
    return s


def find_order_vector(verbose = 1):
    if USE_PRECOMPUTED:
        s_g71, s_v71, s_gA, diag =  S_g71, S_v71, S_gA, S_diag
    else:
        s_g71, s_v71, s_gA, diag = find_vector_71_mod3(verbose)
    for trials in range(200,-1,-1):
        s_g94 = find_element_of_order(94, verbose = verbose)
        s_v94 = find_vector_v94_mod5(s_g94, verbose = verbose)
        if s_v94 is None:
            continue
        v = make_order_vector(s_g71, s_v71, s_gA, diag, s_g94, s_v94)
        tag_data = check_v(v, verbose=verbose)
        if tag_data is not None:
            if verbose:
                print("Tests for v94 passed")
            break
    if not trials:
        err = "No suitable vector in the monster representation found"
        raise ValueError(err) 
    v_data =  s_g71, s_v71, s_gA, diag, s_g94, s_v94
    V_NAMES =  ["S_G71", "S_V71", "S_GA", "DIAG_VA", "S_G94", "S_V94"]
    TAG_NAMES =  ["TAGS_Y", "TAGS_X", "TAG_SIGN"]   
    if verbose:        
        for text, data in zip(TAG_NAMES, tag_data):
            print(str_data(text, data))
    result = OrderedDict(zip(V_NAMES + TAG_NAMES, v_data + tag_data))
    return result


def write_order_vector(result):
    print("Writing file " + PY_FILENAME)
    f = open(PY_FILENAME, "wt")
    for text, data in result.items():
        print(str_data(text, data), file = f)
    f.close()
    
    
if __name__ == "__main__":
    result = find_order_vector(verbose = 1)
    write_order_vector(result)

# AAAAAAAAAA