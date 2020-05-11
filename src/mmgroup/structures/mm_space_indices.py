from numbers import Integral, Number
import warnings
from collections.abc import Iterable
from random import randint
from importlib import import_module

import numpy as np

    
from mmgroup.structures.abstract_rep_space import mod_rand_invertible
from mmgroup.structures.abstract_rep_space import mod_rand_unit

from mmgroup.structures.parse_atoms import ihex 

try:
    from mmgroup.mm import mm_aux_array_extern_to_sparse
except (ImportError, ModuleNotFoundError):
    ERR_MM_LIN = "Linear access to MM vectors not supported"
    warnings.warn(ERR_MM_LIN, UserWarning)   





########################################################################
# Description of tags
########################################################################




# Ofssets for tags in sparse notation
tag_offsets = {
   "A": 0x2000000,
   "B": 0x4000000,
   "C": 0x6000000,
   "T": 0x8000000,
   "X": 0xa000000,
   "Z": 0xc000000,
   "Y": 0xe000000,
}

# Boundaries for the two indices following a standard tag
tag_lengths = {
        "A": (24,24), "B": (24,24), "C": (24,24), "T":(759,64),
        "X": (2048,24), "Y": (2048,24), "Z": (2048,24), "D": (24,24)
}

# Numpy array dtype used for the sparse representation
U32 = np.uint32

# For each tag, the following dictionary contains the letter for 
# the tag and the formatting functions for the following indices. 
str_tag = {
    0: (None, None, None),
    1: ("A", str, str),
    2: ("B", str, str),
    3: ("C", str, str),
    4: ("T", str, ihex),
    5: ("X", ihex, str),
    6: ("Z", ihex, str),
    7: ("Y", ihex, str),
}


tags = " ABCTXZY"

########################################################################
# Converting an tuple describing an unit vector to a sparse list 
########################################################################




def rand_str_to_scalar(p, value):
    """Map a character to a random integer modulo p

    For 'value we have:
    'u' means: return 1
    's' means: return 1 or -1 (modulo p)
    'n' means: return a random nonzerodivisor modulo p
    'r' means: return a random number modulo p

    The return value x satisfies 0 <= x < p.
    """
    if s == 'u': return 1
    if s == 's': return mod_rand_unit(p)
    if s == 'n': return mod_rand_invertible(p)
    if s == 'r': return randint(0, p-1)
    err = "Illegal scalar in tuple for MM vector"
    raise TypeError(err)



def gen_unit_std(p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    In the same way as function tuple_to_sparse(), this 
    function returns the flat numpy array 'a_indices'.

    It implements the standard case where the tag is followd
    by two indices i0, i1. The single entry of  numpy array 
    'a_indices' is coded as the integer

        (offset << 25) + (i0 << 14) + (i1 << 8) + scalar % p,

    where i0 and i1 run through the indices given by the input
    values i0 and i1. 'offset' is given by 'tag_offsets[tag]'
    and the boundaries for the indices i0 and i1 are given by 
    'tag_lengths[tag]'. 
    """   
    max0, max1 = tag_lengths[tag]
    if isinstance(i0, str):
        i0 = randint(0, max0 - 1)
    if isinstance(i1, str):
        i1 = randint(0, max1 - 1)
    if not (0 <= i0 < max0 and 0 <= i1 < max1):
        raise IndexError("Bad indices for MM vector tag %s" % tag)
    sc = scalar % p
    a = [tag_offsets[tag] + (i0 << 14) + (i1 << 8) + sc]
    return np.array(a, dtype = U32)
    

def gen_unit_A(p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tag 'A'.
    
    The input and the return value is as in gen_unit_std().
    This function is coded along the lines of function
    gen_unit_std(). It deals with the following peculiarity:

      -  if i0 < i1:  
             map (tag, i0, i1) to (tag, i1, i0)
    """
    if isinstance(i0, str):
        i0 = randint(0, 23)
    if isinstance(i1, str):
        i1 = randint(0, 23)
    i0, i1 = max(i0, i1), min(i0, i1)
    if not 0 <= i1 <= i0 < 24:
        raise IndexError("Bad indices for MM vector tag %s" % tag)
    sc = scalar % p
    a = [tag_offsets[tag] + (i0 << 14) + (i1 << 8) + sc]
    return np.array(a, dtype = U32)
    

def gen_unit_BC(p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tags 'B' and 'C'.
    
    The input and the return value is as in gen_unit_std().
    This function is coded along the lines of function
    gen_unit_A(). It deals with the following additional
    peculiarities:

      -  if i0 == i1:
             raise IndexError (This case is illegal)
      -  if i1 encodes a random index:
             ensure i1 != i0
    """
    if isinstance(i0, str):
        i0 = randint(0, 23)
    if isinstance(i1, str):
        i1 = randint(0, 23)
        while  i0 == i1:  i1 = randint(0, 23)
    i0, i1 = max(i0, i1), min(i0, i1)
    if not 0 <= i1 < i0 < 24:
        raise IndexError("Bad indices for MM vector tag %s" % tag)
    sc = scalar % p
    a = [tag_offsets[tag] + (i0 << 14) + (i1 << 8) + sc]
    return np.array(a, dtype = U32)


def gen_unit_D(p, scalar, tag, i0 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tag 'D'. ('D', i) is equivalent
    to ('A', i, i).
    """
    if isinstance(i0, str):
        i0 = randint(0, 23)
    if not 0 <= i0 < 24:
        raise IndexError("Bad indices for MM vector tag %s" % tag)
    a = [tag_offsets['A'] + i0 * 0x4100 + scalar % p]
    return np.array(a, dtype = U32)


def gen_unit_numeric(p, scalar, tag, i0 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    It deals with tag 'E' which interprets a single index as
    in the C version of the rep of the monster group.

    Therefore we use function
    
    mmgroup.structures.mm_aux_array_extern_to_sparse()
    """
    if isinstance(i0, str):
        i0 = randint(0, 196883)
    if not 0 <= i0 < 196884:
        raise IndexError("Bad indices for MM vector tag %s" % tag)
    a = np.array([i0], dtype = U32)
    try:
        mm_aux_index_extern_to_sparse(a, 1)
    except NameError:
        raise NotImplementedError(ERR_MM_LIN)   
    a[0] += scalar % p
    return a

def gen_unit_I(p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tag 'I'. ('I', i, j) is equivalent to

      ('A', i, i) +  ('A', j, j) - ('A', i, j) - 2 * ('B', i, j)

    for i != j.
    """   
    def entry(tag, i, j, factor = 1):
        sc = scalar * factor % p
        return tag_offsets[tag] + (i << 14) + (j << 8) + sc
    if isinstance(i0, str):
        i0 = randint(0, 23)
    if isinstance(i1, str):
        i1 = randint(0, 23)
        while  i0 == i1:  i1 = randint(0, 23)
    i0, i1 = max(i0, i1), min(i0, i1)
    if not 0 <= i1 < i0 < 24:
        raise IndexError("Bad indices for MM vector tag %s" % tag)
    a = [entry("A", i0, i0), entry("A", i1, i1),
            entry("A", i0, i1, -1),   entry("B", i0, i1, -2),   ]
    return np.array(a, dtype = U32)


tuple_to_sparse_functions = {
    'A': gen_unit_A,
    'B': gen_unit_BC,
    'C': gen_unit_BC,
    'T': gen_unit_std,
    'X': gen_unit_std,
    'Z': gen_unit_std,
    'Y': gen_unit_std,
    'D': gen_unit_D,
    'E': gen_unit_numeric,
    'I': gen_unit_I,
}



def tuple_to_sparse(p, *data):
    """Convert a tuple to a sparse representation of an MM vector 

    This function maps a tuple to the sparse representation of an
    MM vector, which is (usually) a unit vector. The tuple is 
    preceeded by the first argument p, which is a modulus. 
 
    The general format of such a tuple is 

       (scalar, tag, i0, i1,..)

    Here the 'tag' is a captial letter describing the unit vector,
    i0, i1, ... are indices which are specifying the unit vector,
    depending on the tag. The optional 'scalar' is an integer,
    which is reduced modulo p. It defaults to 1. The 'scalar' may
    also be a lower case letter; in that case a random scalar 
    is computed by function rand_str_to_scalar(p, scalar).

    The function returns the sparse representation of the vector
    to be created as a numpy array with dtype = np.uint32.

    If no data are given, the zero vector is created.
    """ 
    if len(data) == 0:
        return np.array([], dtype = U32)
    scalar = 1
    if isinstance(data[0], Integral):
        scalar, data = data[0] % p, data[1:]
    elif isinstance(data[0], str) and data[0] in "usnr":
        scalar, data = rand_str_to_scalar(p, data[0]), data[1:]
    try:
        tag, indices = data[0], data[1:]
        f = tuple_to_sparse_functions[tag]
    except:
        err = "Illegal or missing tag for creating MM vector"
        raise TypeError(err)
    return f(p, scalar, tag, *indices)




########################################################################
# Converting an index object to a sparse list 
########################################################################


def a_slice(index, max_index):
    """Map index to a tuple (shape, array)

    Here parameter 'index' is an one-dimensional index as used for
    indexing an axis of a numpy array. This means that index must be
    an integer, a slice object (obtained e.g. by evaluating colons
    ':' inside a square bracket), or an Iterable containing integers
    only.

    We assume 0 <= i < max_index for a valid index i, so all integers
    occuring as indices must satisfy this constraint.

    The function returns a pair (shape, indices).

    'indices' is a one-dimensional numpy array with dtype = np.uint32
    containing the requested indices. 
 
    'shape' is a tuple of length 0 if a the requested index is a 
    scalar. It is a tuple of length 1 containing the shape of the 
    array if the requested index is a one-dimensional numpy array. 
    
    """
    if isinstance(index, Integral):
        if -max_index <= index < max_index:
            return (), np.array([index % max_index], dtype = U32)
        else:
            raise IndexError("MM vector index out of range") 
    elif isinstance(index, slice):
        a = np.arange(*index.indices(max_index), dtype = U32)
        return (len(a),), a 
    elif isinstance(index, Iterable):
        if not all(isinstance(item, Integer) for item in index):
            raise TypeError("Bad index type for MM vector")
        if not all(0 <= item < max_index for item in index):
            raise IndexError("MM vector index out of range")
        a = np.array(index, dtype =  U32) 
        return (len(a),), a 
    else:
        raise TypeError("Bad index type for MM vector")
    



########################################################################
# Converting an tuple describing an index to a sparse list 
########################################################################


SLICE = slice(None)


def dbl_indices_to_sparse(tag, i0 = SLICE, i1 = SLICE):
    """Auxiliary function for function sparse_from_indices()

    In the same way as function sparse_from_indices(), this 
    function returns a pair 

          (shape, a_indices) .

    It implements the standard case where the tag is followed
    by two indices i0, i1. Each entry of the flat numpy array 
    'a_indices' is coded as the integer

        (offset << 25) + (i0 << 14) + (i1 << 8),

    where i0 and i1 run through the indices given by the input
    values i0 and i1. 'offset' is given by 'tag_offsets[tag]'
    and the boundaries for the indices i0 and j0 are given by 
    'tag_lengths[tag]'.
    """
    max0, max1 = tag_lengths[tag]
    if (isinstance(i0, Integral) and isinstance(i1, Integral)
            and -max0 <= i0 < max0 and -max1 <= i1 < max1): 
        c = tag_offsets[tag] + (i0 % max0 << 14) + (i1 % max1 << 8) 
        return (), U32(c).reshape(1)
    offset = tag_offsets[tag]
    d0, i0 = a_slice(i0, max0)
    d1, i1 = a_slice(i1, max1)
    i0 = np.expand_dims(i0, axis=1)
    i1 = np.expand_dims(i1, axis=0)
    a_out = np.array((i0 << 14) + (i1 << 8) + offset, dtype=U32)
    return d0 + d1, a_out.reshape(-1)


def D_indices_to_sparse(tag, index = SLICE):
    """Auxiliary function for function sparse_from_indices()

    This function deals with tag 'D'. ('D', i) is equivalent
    to ('A', i, i).
    """
    d0, index = a_slice(index, 24)
    a_out = 0x4100 * index + tag_offsets["A"]
    return d0, a_out.astype(U32)


def numeric_indices_to_sparse(tag, indices = SLICE):
    """Auxiliary function for function sparse_from_indices()

    It deals with tag 'E' which interpretes a single index as
    in the C version of the rep of the monster group.

    Therefore we use function
    
    mmgroup.structures.mm_aux_array_extern_to_sparse()
    """
    d0, a_indices = a_slice(indices, 196884)
    length = len(a_indices)
    if length:
        try:
            mm_aux_array_extern_to_sparse(a_indices, length)
        except NameError:
            raise NotImplementedError(ERR_MM_LIN)   
    return d0, a_indices


indices_to_sparse_functions = {
    'A': dbl_indices_to_sparse,
    'B': dbl_indices_to_sparse,
    'C': dbl_indices_to_sparse,
    'T': dbl_indices_to_sparse,
    'X': dbl_indices_to_sparse,
    'Z': dbl_indices_to_sparse,
    'Y': dbl_indices_to_sparse,
    'D': D_indices_to_sparse,
    'E': numeric_indices_to_sparse,
}


def sparse_from_indices(tag, *indices):
    """Convert index of MM vector to sparse representation.

    Although an MM vector is not a numpy ndarray, components 
    of an MM vector can be obtained in the same way as in the
    case of numpy arrays. Here the MM vector v is indexed as 
    v[tag, i0, i1,...], where 'tag' is a capital letter and 
    i0, i1, ... are numbers or slices or lists of integers as 
    used in numpy ndarrays.
    
    Here v[tag, i0, i1] returns an integer or a subset of
    the components of v as a numpy array. Coding  
    'v[tag, i0, i1] = A'  means that a number or an array-like 
    object A is converted to a scalar or a numpy array A1 of a 
    shape compatible to the shape of v[tag, i0, i1]. Then the 
    corresponding components of v are updated by those of A1. 

    This function converts the index of an MM vector
    given by (tag, *indices) to a pair

       (shape, a_indices).

    It returns that pair.

    Here 'shape' is a tuple describing the shape of a numpy
    array corresponding to the input (tag, *indices). So if
    the input describes a scalar, we have shape = ().

    'a_indices' is the flat numpy array containing the indices
    described by the input. We use the sparse form (with scalar
    equal to zero) for describing these indices. 

    """
    if isinstance(tag, str):
        return indices_to_sparse_functions[tag](tag, *indices)
    return numeric_indices_to_sparse('E', tag, *indices)


########################################################################
# Converting a sparse list entry to a standard form 
########################################################################

def purge_sparse_entry(n):
    """Return the standard form of the sparse index n.

    Return zero if the sparse index n is illegel.

    This is is a slow python function. The fast C functions purge
    sparse indices internally, whenever this is needed.
    """
    tag = n >> 25
    if 1 <= tag <= 3:   # Tags A, B, C
        # Change ('ABC', i0, i1) to ('ABC', i1, i0)  if  i0 < i1
        i0, i1 =  (n >> 14) & 0x7ff, (n >> 8) & 0x3f
        i0, i1 = max(i0, i1),  min(i0, i1) 
        n =  (n & 0xe0000ff) + (i0 << 14) + (i1 << 8)
        # max(i0, i1) >= 24 is illegal for tags 'ABC'
        # i0 == i1 is illegal for tags 'B' and 'C'
        return n if (i0 < 24) and (tag == 1 or i0 != i1) else 0
    elif tag == 4:      # Tag T
        # ('T', i0, i1) is illegal for i0 >= 759
        return n if ((n >> 14) & 0x7ff) < 759 else 0 
    elif 5 <= tag <= 7: # Tags X, Y, Z
        # ('XYZ', i0, i1) is illegal for i1 >= 24
        return n if (n & 0x3f00) < 0x1800 else 0
    else :              # Any other tags are illegal
        return 0

########################################################################
# Converting a sparse list to a list of tuples 
########################################################################

def sparse_to_tuples(a_indices):
    """Convert vector from sparse representation to list of tuples

    Here 'a_indices' is a numpy array representing a vector in
    sparse representation. The function returns a list of tuples

       (scalar, tag, i0, i1)

    Here (tag, i0, i1) represents a unit vector and the integer
    'scalar' is the coordinate of that unit vector.
    """
    return [(i & 0xff, tags[(i >>25) & 7], (i >> 14) & 0x7ff, 
         (i >> 8) & 0x3f) for i in a_indices if i & 0xe000000]


########################################################################
# Manipulating indices and data in a sparse list 
########################################################################


def sparse_add_data(p, shape, a_indices, data):
    """Augment indices (given in sparse representation) by data

    Function sparse_from_indices() returns indices in the form
    (shape, a_indices), where 'a_indices' is a flat numpy array
    containing the indices selected from an MM vector v. Here
    'shape' is the shape of that array.

    When updating these indices of v with data, we first augment 
    the array 'a_indices' with the arrray-like object 'data', 
    containing the data used for updating. Here the entries of 
    'data' are reduced modulo the characteristic 'p' of the 
    vector space.
    
    The array-like object 'data' is converted to a numpy array
    with the given 'shape'.
    """
    if isinstance(data, Number):
        a_indices += data % p
    else:
        if not isinstance(data, np.ndarray):
            data = np.array(data, dtype = np.int64)
        data = (data % p).astype(U32)
        if data.shape != shape:
            err = "Shape mismatch when updating mm vector"
            raise ValueError(err)
        a_indices += data.reshape(-1)
    return a_indices
    





     
def sparse_to_ndarray(p, shape, a_indices):
    """Convert vector in sparse representation to array

    This function returns the data contained in the vector
    'a_indices' in sparse representation to an numpy ndarray of 
    given 'shape'. The entries of the vector are reduced 
    modulo  p.   
    """
    a = np.reshape((a_indices & 0xff) % p, shape).astype(np.uint8)
    return a if len(shape) else int(a)


def sparse_to_str(p, a_indices):
    """Convert vector in sparse representation to a string

    This function returns the data contained in the vector
    'a_indices' in sparse representation as a string. The 
    entries of the vector are reduced modulo p.   
    """
    names = []
    for a in a_indices:
        tag, fmt0, fmt1 = str_tag[(a >> 25) & 7]
        i0, i1, v = (a >> 14) & 0x7ff, (a >> 8) & 0x3f, (a & 0xff) % p
        if tag and v:
            s = "-" if v == p - 1 else "+"
            if  1 < v < p - 1:  
                s += str(p)
            names.append("%s%s_%s_%s" % (s, tag, fmt0(i0), fmt1(i1)))            
    if len(names): 
        if names[0][:1] == "+": 
            names[0] = names[0][1:]  
        return "".join(names)
    return "0"

    