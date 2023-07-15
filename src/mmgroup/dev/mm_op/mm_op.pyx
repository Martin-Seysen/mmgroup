
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


include "../pxd_files/mm_op_p.pxi"



