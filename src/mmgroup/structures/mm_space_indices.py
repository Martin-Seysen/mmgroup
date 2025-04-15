"""Components of a vector of the representation of the monster 

A vector of the representation of the monster contains 196884
entries which are integers modulo a small odd prime p. The index 
of such an entry is addressed as a tuple (tag, i0, i1).

Here 'tag' is a single capital letter specifying a group of 
indices and i0 and i1 are (usually) integers representing
elements of the Parker loop or the Golay cocode, depending 
on the tag. Details are specified elsewhere.

A tuple  ([scalar], tag, i0, i1) also specifies the multiple of
of scalar * u  of the unit vector u with index (tag, i0, i1).

If such a tuple is given by the user, it is first converted to
of 32-bit integer containing suitable bit fields for its
component scalar, tag, i0, i1. This conversion is done by the 
functions in this module. Indices i0, i1 may also be slices as 
used in numpy. They may also be elements of the Golay code, the 
Parker loop, or the Golay cocode, where appropiate. Furthermore,
an index may be "r" or "n", denoting an index selected at
random.

If a tuple contains one or two slices, it is converted to a
numpy array of type ``uint32``, which contains all entries of the
slice. The one-dimensional array unsigned of 32-bit integers
is passed to the to the fast C routines for reading and updating 
components of a vector.

"""

from numbers import Integral, Number
import warnings
from collections import defaultdict
from collections.abc import Iterable
from functools import partial
from random import randint
from importlib import import_module

import numpy as np


import mmgroup
from mmgroup import mat24 
from mmgroup.structures.abstract_rep_space import mod_rand_invertible
from mmgroup.structures.abstract_rep_space import mod_rand_unit

from mmgroup.structures.parse_atoms import ihex 
from mmgroup.structures.gcode import GCode, GcVector 
from mmgroup.structures.cocode import Cocode 
from mmgroup.structures.ploop import PLoop
from mmgroup.structures.suboctad import Octad

ERR_MM_LIN = "Linear access to MM vectors not supported"

try:
    from mmgroup.mm_op import mm_aux_array_extern_to_sparse
    from mmgroup.mm_op import mm_aux_index_leech2_to_sparse
except (ImportError, ModuleNotFoundError):
    mmgroup._warn(ERR_MM_LIN)
    def mm_aux_array_extern_to_sparse(*args, **kwds):
        raise NotImplementedError(ERR_MM_LIN)
    mm_aux_index_leech2_to_sparse =  mm_aux_array_extern_to_sparse


ERR_MM_SPARSE = "Sparse representation of MM vectors not supported"

try:
    from mmgroup.mm_op import mm_aux_mul_sparse
except (ImportError, ModuleNotFoundError):
    mmgroup._warn(ERR_MM_SPARSE)



########################################################################
# Collecting indices
########################################################################

INDEX_IS_SLICE = 1
INDEX_IS_RANDOM = 2

INDEX_ERROR_TYPE = "%sMM vector index of type '%s' is illegal for tag %s"
INDEX_ERROR_RAND = "Bad random description of MM vector index for tag %s"
INDEX_ERROR_RANGE = "%sMM vector index out if range for tag %s"
                     
POS_NAMES = defaultdict(str)
POS_NAMES.update({1: "First ", 2: "Second "})


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
# Converting numbers and Code elements to integers
########################################################################



def index_iterable(index, max_index, tag, pos):
    """Convert iterable object 'index' to an array of indices

    'index' must be an iterable object containing integers
    0 <= i < 'max_index' only. The function returns a numpy
    array of shape np.uint32 containing these indices in the 
    given order.

    'tag' and 'pos' are used in error messages only. They refer 
    to the tag to which the indes belong, an to the position of
    the indes (0 for the 1st, 1 for the 2nd index).  
    """
    if not all(isinstance(item, Integral) for item in index):
        err = "All MM vector indices in a list must be integers"
        raise TypeError(err)
    if not all(0 <= item < max_index for item in index):
        err = "%sMM vector indices out of range for tag %s"
        raise IndexError(err, POS_NAMES[pos], tag)
    return  INDEX_IS_SLICE, np.array(index, dtype =  U32), 0 


def a_slice(i_slice, max_index):
    """Convert slice object 'i_slice' to an array of indices

    The slice object 'i_slice' is interpreded as in a one-
    dimensional numpy array of length 'max_index'. The function
    returns the indices contained in that slice when applied
    to a one-dimensional numpy array.

    The function returns a numpy arrray of shape np.uint32.
    """
    return np.arange(*i_slice.indices(max_index), dtype = U32)

def index_I(tag, i, pos = 2):
    """Convert an index specifying a bit position

    0 <= 'i' < 24 is an index specifying a bit position in the 
    Golay code. If 'i' is of type Cocode or GcVector, corre-
    sponding to a bit vector of length one, then the position of 
    that bit is taken as the index. The function returns triple

      (index_type, a_indices, sign).

    Here 'index_type' is an integer interpreted as a bit field. 

    If index_type & INDEX_IS_SLICE is not zero then 'i'
    specifies a slice of indices. The numbers of these indices 
    returned in the numpy array 'a_indices'
    
    If index_type & INDEX_IS_SLICE is zero 'i' specifies a
    single index, and the numpy array 'a_indices' of length 1
    contains that indes.

    'sign' is always 0. It is for compatibility with functions 
    converting e.g. a Parker loop element (which comes with
    a natural sign) to a possible negated unit vector.
    
    'pos' ist used in error messages only. It denotes the 
    position of the index (0 for the 1st, 1 for the 2nd index).

    Both indices of tags A, B, and C are of this type and also 
    the second index of tags X, Y and Z. That tag must be given
    as Parameter 'tag'.
    """
    if isinstance(i, Integral):
        if 0  <= i < 24:
            return 0, i, 0
        raise IndexError(INDEX_ERROR_RANGE % (POS_NAMES[pos], tag))
    elif isinstance(i, str):
        if i in ["r", "n"]:
            return INDEX_IS_RANDOM, randint(0,23), 0
        raise IndexError(INDEX_ERROR_RAND % tag)
    elif isinstance(i, slice):     
        return INDEX_IS_SLICE, a_slice(i,24), 0
    elif isinstance(i, Cocode):
        if len(i) == 1: 
            return 0, i.syndrome().bit_list[0], 0  
        err = "MM vector index of class Cocode for tag %d must have weight 1"
        raise IndexError(err % tag)
    elif isinstance(i, GcVector):
        if len(i) == 1: 
            return 0, i.bit_list[0], 0  
        err = "MM vector index of class GcVector for tag %d must have weight 1"
        raise IndexError(err % tag)
    elif isinstance(i, Iterable):
        return index_iterable(i, 24, tag, pos)
    raise TypeError(INDEX_ERROR_TYPE % (POS_NAMES[pos], type(i), tag))


def index_T(tag, i):
    """Convert an index specifying an octad number

    0 <= 'i' < 759 is an index specifying the number of an octad
    in the Golay code. If 'i' is of type GCode or PLoop, corre-
    sponding  to a Golay code word or a Parker Loop element, that 
    element is taken, provided that it (or its complement) has
    length 8; otherwise we raise ValueError. 

    The function returns triple

      (index_type, a_indices, sign).

    Here index_type indicates whether the index has been selected
    ar random, or whether it refers to a slice, as in function 
    index_I().

    The index or the indices making up the slice are returned in
    the numpy array 'a_indices', as in function index_I(). In any
    case the numbers in that array are octad numbers.
    
    Note that a Parker loop element may refer to a negated unit 
    vector. In this case we return the number of the unit vector 
    and we return 'sign' = 1. It is illegal to specify slices 
    of Parker loop elements; so all entries of array 'a_indices'
    have the same sign.

    The first index of tag T is of this type. Parameter 'tag' is
    for compatibility with similar functions; ist must be 'T'.	 
    """
    if isinstance(i, Integral):
        if 0 <= i < 759:
            return 0, i, 0
        raise IndexError(INDEX_ERROR_RANGE % ("First ", tag))
    elif isinstance(i, str):
        if i in ["r", "n"]:
            return INDEX_IS_RANDOM, randint(0,758), 0
        raise IndexError(INDEX_ERROR_RAND % tag)
    elif isinstance(i, slice):     
        return INDEX_IS_SLICE, a_slice(i,759), 0
    elif isinstance(i, (GCode, PLoop)):
        oct = Octad(i) 
        return 0, oct.octad, oct.sign == -1
    elif isinstance(i, Iterable):
        return index_iterable(i, 759, tag, 1)
    raise TypeError(INDEX_ERROR_TYPE % ("First ", type(i), tag))

def index_suboctad(tag, i, index_type_T = 3, a_indices_T = None):
    """Convert an index specifying an suboctad number

    0 <= 'i' < 64 is an index specifying the number of a suboctad
    in the Golay cocode. Here a suboctad is a Golay cocode word
    which can be represented a subset of of an octad. For each
    octad there is a lexicographic numbering of its suboctad, and
    index 'i' refers to that lexicographic numbering. 

    This index type occurs as the second index of tag T only,
    so the first index of that tag gives us a reference octad. 
    If 'i' is of type Cocode or GcVector, corresponding to a 
    Golay cocode word, that cocode word is taken if it fits into 
    the reference octad; otherwise we raise ValueError. 

    The function returns triple

      (index_type, a_indices, sign).

    Here index_type indicates whether the index has been selected
    ar random, or whether it refers to a slice, as in function 
    index_I().

    The index or the indices making up the slice are returned in
    the numpy array 'a_indices', as in function index_I(). In any
    case the numbers in that array are suboctad numbers.

    Parameters ('index_type_T', 'a_indices_T') must be the values 
    (index_type, a_indices) returned by function index_T(), when
    called with the first index for tag T.
    
    The second index of tag T is of this type. Parameter 'tag' is
    for compatibility with similar functions; it must be 'T'.		 
    """
    if isinstance(i, Integral):
        if 0 <= i < 64:
            return 0, i, 0
        raise IndexError(INDEX_ERROR_RANGE % ("Second ", tag))
    elif isinstance(i, str):
        if i in ["r", "n"]:
            return INDEX_IS_RANDOM, randint(0,63), 0
        raise IndexError(INDEX_ERROR_RAND % tag)
    elif isinstance(i, slice):     
        return INDEX_IS_SLICE,  a_slice(i,64), 0
    elif isinstance(i, (GCode, Cocode, GcVector)):
        if index_type_T == 0:
            gcode = mat24.octad_to_gcode(a_indices_T) & 0xfff
            if isinstance(i, GCode):
                cocode = mat24.ploop_cap(gcode, i)
            else:
                cocode = Cocode(i).cocode
            sub = mat24.cocode_to_suboctad(cocode, gcode) & 0x3f
            return 0, sub, 0
        err = "Second vector index for tag T must be integer here"
        raise TypeError(err)
    elif isinstance(i, Iterable):
        return index_iterable(index, 64, tag, 2)
    raise TypeError(INDEX_ERROR_TYPE % ("Second ", type(i), tag))


def index_ploop(tag, i):
    """Convert an index specifying an Parker loop element

    0 <= 'i' < 2**11  is an index specifying the number of an
    element of the Parker loop. The Parker loop has 2**13 
    elements, with the first 2**11 elements being a transversal 
    of the quotient of the Parker Loop by its center. Flipping 
    bit 12 in the number of a Parker loop element always negates 
    the unit vector adressed by that element. Flipping bit 11 may
    either negate or keep that unit vetor, depending on the tag.
    In any case there are 2**11 positive and 2**11 negative unit 
    vectors to be distinguished, and an integer  0 <= 'i' < 2**11
    refers to a positive unit vector.

    The function returns the triple

      (index_type, a_indices, sign).

    Here index_type indicates whether the index has been selected
    ar random, or whether it refers to a slice, as in function 
    index_I().

    The index or the indices making up the slice are returned in
    the numpy array 'a_indices', as in function index_I(). In any
    case the numbers j in that array satisfy 0 <= j < 2**11.

    If 'i' is of type PLoop, corresponding  to a Parker loop 
    element, then the correct (positive or negative) unit vector
    is chosen, and we return 'sign' = 0 for a positve and 
    'sign' = 1 for a negative vector. Note that slices of Parker 
    loop  elements and integers >= 2**11 are illegal, so all 
    entries of the array 'a_indices' have the same sign.

    The first index of tags X, Y and Z is of this type. That tag 
    must be given as Parameter 'tag'.
    """

    if isinstance(i, Integral):
        if 0 <= i < 2048:
            return 0, i, 0
        raise IndexError(INDEX_ERROR_RANGE % ("Second ", tag))
    elif isinstance(i, str):
        if i in ["r", "n"]:
            return INDEX_IS_RANDOM, randint(0,2047), 0
        raise IndexError(INDEX_ERROR_RAND % tag)
    elif isinstance(i, slice):     
        return INDEX_IS_SLICE,  a_slice(i,2048), 0
    elif isinstance(i, (GCode, PLoop)):
        i = i.ord
        sign = (i >> 12) ^ (i >> 11) if tag == "Y" else (i >> 12)
        return 0, i & 0x7ff, sign & 1
    elif isinstance(i, Iterable):
        return index_iterable(index, 2048, tag, 2)
    raise TypeError(INDEX_ERROR_TYPE % ("Second ", type(i), tag))



def index_numeric(tag, i):
    """Convert an numeric index of a vector of the representation

    A vector in the monster group representation may also be 
    adressed by a numeric  index in the same way as a one-
    dimensional numpy array of length 196884.

    The function returns the triple

      (index_type, a_indices, sign).

    Here index_type indicates whether the index has been selected
    ar random, or whether it refers to a slice, as in function 
    index_I().

    The index or the indices making up the slice are returned in
    the numpy array 'a_indices', as in function index_I().

    'sign' is always zero. The 'tag' it for compatibility  with
    similar functions; it should always be 'E'.
    """
    if isinstance(i, Integral):
        if 0 <= i < 196884:
            return 0, i, 0
        raise IndexError(INDEX_ERROR_RANGE % ("First ", tag))
    elif isinstance(i, str):
        if i in ["r", "n"]:
            return INDEX_IS_RANDOM, randint(0,196883), 0
        raise IndexError(INDEX_ERROR_RAND % tag)
    elif isinstance(i, slice):     
        return INDEX_IS_SLICE, a_slice(i, 196884), 0
    elif isinstance(i, Iterable):
        return index_iterable(index, 196884, tag, 1)
    raise TypeError(INDEX_ERROR_TYPE % ("First ", type(i), tag))


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
    if value == 'u': return 1
    if value == 's': return mod_rand_unit(p)
    if value == 'n': return mod_rand_invertible(p)
    if value == 'r': return randint(0, p-1)
    err = "Illegal scalar in tuple for MM vector"
    raise TypeError(err)




    


INDEX_ERROR_SLICE = "Slices are illegal in MM vector index tuples"

def gen_unit_A(p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tag 'A'. So 'tag' must be 'A'.

    It implements the case for tag "A" where the tag is followd
    by two indices i0, i1. The function returns a numpy array 
    of shape np.uint32 with a single entry:

        (offset << 25) + (i0 << 14) + (i1 << 8) + scalar % p,

    where i0 and i1 are the indices given by the input
    values i0 and i1. 'offset' is given by 'tag_offsets[tag]'. 

    It deals with the following peculiarity:

      -  if i0 < i1:  
             map (tag, i0, i1) to (tag, i1, i0)
    """
    if i0 == 'r' and i1 == None:
        i1 ='r'
    index_type0, i0, _ = index_I(tag, i0, pos = 1)
    index_type1, i1, _ = index_I(tag, i1, pos = 2)
    if (index_type0 | index_type1) & INDEX_IS_SLICE:
        raise TypeError(INDEX_ERROR_SLICE)
    i0, i1 = max(i0, i1), min(i0, i1)
    sc = scalar % p
    a = [tag_offsets[tag] + (i0 << 14) + (i1 << 8) + sc]
    return np.array(a, dtype = U32)
    

def gen_unit_BC(p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tags 'B' and 'C'. So 'tag' must
    be 'B' or 'C'.  It implements the case for tags 'B' in
    the same way as function gen_unit_A() implements the
    case for tag 'A'. It deals with the following additional
    peculiarities:

      -  if i0 == i1:
             raise IndexError("This case is illegal")
      -  if i0 or i1 encodes a random index:
             ensure i1 != i0
    """
    if i0 == 'r' and i1 == None:
        i1 ='r'
    index_type0, j0, _ = index_I(tag, i0, pos = 1)
    index_type1, j1, _ = index_I(tag, i1, pos = 2)
    if (index_type0 | index_type1) & INDEX_IS_SLICE:
        raise TypeError(INDEX_ERROR_SLICE)
    is_rnd = (index_type0 | index_type1) & INDEX_IS_RANDOM
    while is_rnd and j0 == j1:
        _, j0, _ = index_I(tag, i0)
        _, j1, _ = index_I(tag, i1)
    if j0 == j1:
        err = "Vector indices 1 and 2 must be different for tag %s"
        raise ValueError(err % tag)
    j0, j1 = max(j0, j1), min(j0, j1)
    sc = scalar % p
    a = [tag_offsets[tag] + (j0 << 14) + (j1 << 8) + sc]
    return np.array(a, dtype = U32)


def gen_unit_D(p, scalar, tag, i0 = 'r', i1 = None):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tag 'D'. ('D', i) is equivalent
    to ('A', i, i).
    """
    index_type0, i0, _ = index_I(tag, i0, pos = 1)
    if (index_type0) & INDEX_IS_SLICE:
        raise TypeError(INDEX_ERROR_SLICE)
    a = [tag_offsets['A'] + i0 * 0x4100 + scalar % p]
    return np.array(a, dtype = U32)


def gen_unit_T(p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tag 'T'. It is implemented in
    the same fashion as function gen_unit_A(). Note that 
    the interpretation of parameter 'i1' depends on 'i0'.
    """
    if i0 == 'r' and i1 == None:
        i1 ='r'
    index_type0, i0, sign = index_T(tag, i0)
    index_type1, i1, _ = index_suboctad(tag, i1, index_type0, i0)
    if (index_type0 | index_type1) & INDEX_IS_SLICE:
        raise TypeError(INDEX_ERROR_SLICE)
    scalar = scalar * (-1)**sign
    a = [tag_offsets[tag] + i0 * 0x4000 + i1 * 0x100 + scalar % p]
    return np.array(a, dtype = U32)


def gen_unit_XYZ(p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tags 'X', 'Y' and 'Z'. It is 
    implemented in the same fashion as function gen_unit_A(). 
    Here parameter 'i0' and 'i1' are idependent.
    """
    if i0 == 'r' and i1 == None:
        i1 ='r'
    index_type0, i0, sign = index_ploop(tag, i0)
    index_type1, i1, _ = index_I(tag, i1)
    if (index_type0 | index_type1) & INDEX_IS_SLICE:
        raise TypeError(INDEX_ERROR_SLICE)
    scalar = scalar * (-1)**sign
    a = [tag_offsets[tag] + i0 * 0x4000 + i1 * 0x100 + scalar % p]
    return np.array(a, dtype = U32)


def gen_unit_numeric(p, scalar, tag, i0 = 'r', i1 = None):
    """Auxiliary function for function tuple_to_sparse()

    It deals with tag 'E'. So 'tag' should aways be 'E'. Tag 'E'
    means that 'i0' is interpreted as an index adressing a numpy 
    array of shape (196884,) in the same way as in the C version 
    of the representation of the monster group. 'i0' may be as
    slice in the same way as in a numpy array.
    
    Internally, a vector of the representation of the monster 
    group is not a standard numpy array. We use the C function
    
    mmgroup.structures.mm_aux_array_extern_to_sparse()

    for converting i0 to a standard index.
    """
    index_type0, i0, _ = index_numeric(tag, i0)
    if index_type0 & INDEX_IS_SLICE:
        raise TypeError(INDEX_ERROR_SLICE)
    a = np.array([i0], dtype = U32)
    mm_aux_array_extern_to_sparse(a, 1)
    a[0] += scalar % p
    return a
    

def gen_unit_IJ(sign, p, scalar, tag, i0 = 'r', i1 = 'r'):
    """Auxiliary function for function tuple_to_sparse()

    This function deals with tags 'I' and 'J'. Put f(sign, i, j) =

      ('A', i, i) + ('A', j, j) - ('A', i, j) - 2 * sign * ('B', i, j).

    Then ('I', i, j) and ('J', i, j) are equivalent to f(1, i, j) and
    f(-1, i, j), respecvtively, for i != j. The case i == j is illegal.
    """
    def A(i, j, factor = 1):
        sc = scalar * factor % p
        return int(tag_offsets["A"] + (i << 14) + (j << 8) + sc)
    b =  gen_unit_BC(p, -2 * sign * scalar, 'B', i0, i1)[0] 
    i0, i1 = (b >> 14) & 0x1f, (b >> 8) & 0x1f
    a = [A(i0, i0), A(i1, i1), A(i0, i1, -1), b]
    return np.array(a, dtype = U32)



def gen_unit_identity(p, scalar, *args):
    """Auxiliary function for function tuple_to_sparse()

    Generates the identity vector

    scalar * sum([('A', i, i) for i in range(24)])
    """
    start = tag_offsets["A"] + scalar % p
    d = (1 << 14) + (1 << 8)
    return np.array(range(start, start + 24 * d, d), dtype = U32)



def gen_vector_sparse(p, scalar, tag, data, p1=None):
    """Auxiliary function for function sparse_from_indices()

    It deals with tag 'S' where an array 'data' in sparse 
    representation (modulo 'p1') is multiplied with 'scalar'
    converted to a sparse array modulo 'p'.

    Therefore we use function ``mm_aux_mul_sparse``.
    """
    assert tag == 'S'
    if not p1:
        p1 = p
    a = np.array(data, dtype = np.uint32)
    if a.shape == (): 
        a = a.ravel()
    if p1 == p and scalar == 1:
        return a
    a1 = np.zeros(len(a), dtype = np.uint32)
    res =  mm_aux_mul_sparse(p1, a, len(a), scalar, p, a1)
    if res < 0:
        err = "Cannot reduce sparse MM vector modulo %d"
        raise ValueError(err % p)
    return a1[:res]


def gen_zero(*args):
    return np.zeros(0, dtype = np.uint32)

tuple_to_sparse_functions = {
    'A': gen_unit_A,
    'B': gen_unit_BC,
    'C': gen_unit_BC,
    'T': gen_unit_T,
    'X': gen_unit_XYZ,
    'Z': gen_unit_XYZ,
    'Y': gen_unit_XYZ,
    'D': gen_unit_D,
    'E': gen_unit_numeric,
    'I': partial(gen_unit_IJ, 1), 
    'J': partial(gen_unit_IJ, -1), 
    'U': gen_unit_identity,
    'S': gen_vector_sparse,
    '0': gen_zero,
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
# Converting an tuple describing an index to a sparse list 
########################################################################


SLICE = slice(None)


INDEX_ERROR_SLICE = "Cannot access a random index of an  MM vector"
PLOOP_SIGN_ERROR = "Negative Parker loop elements not allowed as indices" 

def indices_to_sparse(p, tag, type0, i0, type1, i1, sign = 0):
    """Auxiliary function for function sparse_from_indices()

    The function deals with adressing a vector v of the monster
    group representation  in the form v[tag, i0, i1]. Here
    tag, i0, and i1 refer to the index of a (possible negated) 
    component of v. i0 and i1 may also refer to slices of such
    components. 

    A function like index_ploop() or index_I() returns a triple
    (index_type, a_indices, sign) describing the information
    contained in a single index. Here 'index_type' is a bit field
    contaning the following information:

      - Does the index describe a single value or a slice?
      - Has the index been generated at random?  

    The array 'a_indices' is a numpy array containing all
    indices referred by a slice object. It contains one entry
    if the index refers to single component. If 'sign' & 1 == 1
    then the negated components have to be considered instead.

    Let (type0, i0, sign) be the output 
         (index_type, a_indices, sign)
    obtained by applying a function like index_ploop() to the
    first index and let  (type0, i0, dummy) be the output 
         (index_type, a_indices, sign)
    obtained by applying a function like index_I() to the
    second index. Depending on type0 and type1, the
    value v[tag, i0, i1] is either a scalar or an numpy of 
    dimension 1 or 2. 

    The function returns a pair (shape, a_out). Here shape is
    the expected shape of the array v[tag, i0, i1], with 
    shape = () if v[tag, i0, i1] is a scalar. 
    
    a_out is the flattened array contructed from arrays i0 and
    i1 containing the indices (j0, j1) for all j0 in i0 and
    j1 in sparse representation. 

    Remarks:
    
    - The sparse representation corresponding to v[tag, j0, j1] 
      is: 
           tag_offsets[tag]  + (i0 << 14) + (i1 << 8).

    - If 'sign' & 1 == 1 (meaning negated components) then we xor  
      the number p to that sparse representation. This makes sense 
      in case p = 2**k - 1. Then for any x with 0 <= x <= p we have 

           -x  =  x ^ p    (mod p).

      We support negated components in case p == 2**k-1  only.
      
    - If any index i0 or i1 has been generated at random, 
      we raise TypeErrror.

    - We use C-style, i.e. order = C, for reshaping numpy arrays.
    """
    if (type0 | type1) & INDEX_IS_RANDOM:
        raise TypeError(INDEX_ERROR_SLICE)
    offset = tag_offsets[tag] + (p & -sign)
    if (type0 & INDEX_IS_SLICE) == 0:
        if (type1 & INDEX_IS_SLICE) == 0:
            offset += (i1 << 8) + (i0 << 14)            
            return (), np.array([offset], dtype=U32)
        else:
            i1 = np.array(i1,  dtype=U32)
            offset += i0 << 14
            a_out = (i1 << 8) + offset
            return (len(i1),), a_out
    else:
        i0 = np.array(i0,  dtype=U32)
        if (type1 & INDEX_IS_SLICE) == 0:
            offset += i1 << 8
            a_out = (i0 << 14) + offset
            return (len(i0),), a_out
        else:
            i1 = np.array(i1,  dtype=U32)
            shape = (len(i0), len(i1))
            i0 = np.expand_dims(i0, axis=1)
            i1 = np.expand_dims(i1, axis=0)
            a_out = (i0 << 14) + (i1 << 8) + offset
            return shape, a_out.reshape(-1, order = 'C')
    

def ABC_indices_to_sparse(p, tag, i0 = SLICE, i1 = SLICE):
    shape0, a0, _ =  index_I(tag, i0, pos = 1)
    shape1, a1, _ =  index_I(tag, i1, pos = 1)
    return indices_to_sparse(p, tag, shape0, a0, shape1, a1)

def T_indices_to_sparse(p, tag, i0 = SLICE, i1 = SLICE):
    shape0, a0, sign =  index_T(tag, i0)
    shape1, a1, _ =  index_suboctad(tag, i1, shape0, a0)
    if sign and (p & (p+1)):
        raise NotImplementedError(PLOOP_SIGN_ERROR)
    return indices_to_sparse(p, tag, shape0, a0, shape1, a1, sign)

def XYZ_indices_to_sparse(p, tag, i0 = SLICE, i1 = SLICE):
    shape0, a0, sign =  index_ploop(tag, i0)
    shape1, a1, _ =  index_I(tag, i1)
    if sign and (p & (p+1)):
        raise NotImplementedError(PLOOP_SIGN_ERROR)
    return indices_to_sparse(p, tag, shape0, a0, shape1, a1, sign)




def D_index_to_sparse(p, tag, i0 = SLICE):
    shape0, a0, _ =  index_I(tag, i0, pos = 1)
    offset = tag_offsets["A"]
    if shape0 & INDEX_IS_RANDOM:
        raise TypeError(INDEX_ERROR_SLICE)
    if (shape0 and INDEX_IS_SLICE) == 0:
        offset += a0 * 0x4100
        return (), np.array([offset], dtype=U32)
    else:
        a_out = a0 * 0x4100 + offset
        return (len(a0),), a_out
        


def numeric_index_to_sparse(p, tag, i0 = SLICE):
    """Auxiliary function for function sparse_from_indices()

    It deals with tag 'E' which interpretes a single index as
    in the C version of the rep of the monster group.

    Therefore we use function
    
    mmgroup.structures.mm_aux_array_extern_to_sparse()
    """
    shape0, a0, _ =  index_numeric(tag, i0)
    if shape0 & INDEX_IS_RANDOM:
        raise TypeError(INDEX_ERROR_SLICE)
    if (shape0 and INDEX_IS_SLICE) == 0:
        shape = ()
        a_out = np.array([a0], dtype=U32)
    else:
        shape = (len(a0),)
        a_out = a0
    mm_aux_array_extern_to_sparse(a_out, len(a_out))
    return shape, a_out




indices_to_sparse_functions = {
    'A': ABC_indices_to_sparse,
    'B': ABC_indices_to_sparse,
    'C': ABC_indices_to_sparse,
    'T': T_indices_to_sparse,
    'X': XYZ_indices_to_sparse,
    'Z': XYZ_indices_to_sparse,
    'Y': XYZ_indices_to_sparse,
    'D': D_index_to_sparse,
    'E': numeric_index_to_sparse,
}



def xleech2_to_sparse(p, x):
    """Convert instance of class XLeech2 to to sparse representation
    """
    v = x.value
    vs = mm_aux_index_leech2_to_sparse(v & 0xffffff)
    if vs:
        vs += -(v >> 24) & p
        return (), np.array([vs], dtype=U32)
    E = "Index of type XLeech2 is not short in Leech lattice mod 2"
    raise ValueError(E)


def sparse_from_indices(p, tag, *indices):
    """Convert index of MM vector to sparse representation.

    Although an MM vector is not a numpy ndarray, components 
    of an MM vector can be obtained in the same way as in the
    case of numpy arrays. Here the MM vector v is indexed as 
    v[tag, i0, i1], where 'tag' is a capital letter and 
    i0, i1 are numbers or slices or lists of integers as 
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
        return indices_to_sparse_functions[tag](p, tag, *indices)
    from mmgroup.structures.xleech2 import XLeech2
    if isinstance(tag, XLeech2):
        return xleech2_to_sparse(p, tag)
    return numeric_index_to_sparse(p, 'E', tag, *indices)


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
         (i >> 8) & 0x3f) for i in map(int, a_indices) if i & 0xe000000]


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
        a_indices ^= data % p
    else:
        if not isinstance(data, np.ndarray):
            data = np.array(data, dtype = np.int64)
        data = (data % p).astype(U32)
        if data.shape != shape:
            err = "Shape mismatch when updating mm vector"
            raise ValueError(err)
        a_indices ^= data.reshape(-1)
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
    half = p >> 1
    for a in a_indices:
        tag, fmt0, fmt1 = str_tag[(a >> 25) & 7]
        i0, i1, v = (a >> 14) & 0x7ff, (a >> 8) & 0x3f, (a & 0xff) % p
        if tag and v:
            s, vabs = ("-", p - v) if v > half else ("+", v)
            s += (str(vabs) + '*') if vabs > 1 else ""
            names.append("%s%s_%s_%s" % (s, tag, fmt0(i0), fmt1(i1)))            
    if len(names): 
        if names[0][:1] == "+": 
            names[0] = names[0][1:]  
        return "".join(names)
    return "0"

    