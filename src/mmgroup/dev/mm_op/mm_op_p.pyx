
# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


cimport mm_op_p as op_p



def characteristics():
    """Return list of all supported moduli ``p`` 

    The function returns the list of all moduli ``p`` such that a
    representation of the Monster modulo ``p`` is supported.
    """ 
    cdef uint32_t i = 0
    cdef uint32_t j
    while op_p.MM_OP_P_TABLE[i]:
        i += 1
    return [op_p.MM_OP_P_TABLE[j] for j in range(i)] 


include "mm_basics.pxi"



######################################################################
### Python functions added for convenience
######################################################################




import numpy


PROTECT_OVERFLOW =  (17 << INT_BITS) // 19

def mmv_array(*args):
    """Create array with data type suitable for a Monster rep vector 

    This is a numpy array of shape given by the positional arguments
    and of type suitable for the monster group.

    For internal purposes only.
    """
    return numpy.zeros(args, dtype = numpy.uint64)

def mm_vector(p, *args):
    """Create zero vector or array suitable for a Monster rep vector 

    Such vectors are implemented as numpy arrays
    E.g.
    mm_vector(3)   # Create a monster rep vector in characteristic p=3
    mm_vector(3,5) # Create an array of 5 such vectors

    The most convenient way to use this function is e.g.:

    from op_p import MM
    mm3 = MM[3]        # Class containing monster functions for p=3
    v = mm3.vector()   # Create one monster rep vector v for p=3
    mm3.random_mmv(v)  # randomize the vector v
    w = mm3.vector(5)  # Create array of 5 such monster rep vectors 
    w[3] = v           # standard numpy operation
    mm3.get_mmv(v,3)   # get component 3 of v
    """
    assert (p & (p + 1) == 0) and 3 <= p < 256
    vsize = mm_aux_mmv_size(p)
    result = numpy.zeros(args + (vsize + 1,), dtype = numpy.uint64)
    copy = result.reshape(-1, vsize + 1)
    for i in range(copy.shape[0]):
         copy[i, vsize] = PROTECT_OVERFLOW
    return result






include "../pxd_files/mm_op_p.pxi"



