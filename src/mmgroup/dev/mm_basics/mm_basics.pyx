
# cython: language_level=3

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals




include "../pxd_files/mm_basics.pxi"




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

    from mm_op import MM
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

def mm_rng_make_seed(rng_no = 0, data = None):
    """Create a seed for the internal random generator

    Here rng_no is the number of the random generator, data is an
    array of bytes used for seeding it, see module mm_random.ske.
    """
    seed = numpy.zeros(266, dtype = numpy.uint8)
    cdef uint8_t[::1] seed_v_ = seed
    if data and len(data):
       a_data =  numpy.array(data, dtype =  numpy.uint8)
       length = len(a_data)
    else:
       length, a_data = 0, numpy.zeros(1, dtype =  numpy.uint8)
    cdef uint8_t[::1] a_data_v_ = a_data
    mm_rng_seed(seed_v_, rng_no, a_data_v_, len(a_data)) 
    del a_data
    return seed


def tag_dict_from_compressed(uint32_t p, d, a):
    cdef uint32_t[::1] pa = a
    cdef uint32_t len_a = len(a)
    cdef uint32_t tag, i, i0, i1, x
    for i in range(len_a):
        x = pa[i]
        tag = (x >> 25) - 1
        i0 = (x >> 14) & 0x7ff
        i1 = (x >> 8) & 0x3f
        x &= p
        try:
            d[("ABCTXZY"[tag], i0, i1)] = x
        except IndexError:
            raise IndexError("Bad compressed entry %s" % hex(pa[i]))        
        

def tag_dict_to_compressed(d, a):
    cdef uint32_t[::1] pa = a
    cdef uint32_t len_a = len(a)
    cdef uint32_t i = 0
    cdef uint32_t tag, i0, i1, x
    for  (t, i0, i1), x in d.items():
        tag = "0ABCTXZY".index(t)
        pa[i] = (tag << 25) + (i0 << 14) + (i1 << 8) + x
        i = i + 1
        if i > len_a:
            raise IndexError("Compressed array too short")
    return i

    
    





