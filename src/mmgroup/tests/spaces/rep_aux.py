from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


"""Auxilary functions for group representations"""

import numpy as np 


def sgn(s):
    """return (-1)**(s)"""
    return 1 - ((s & 1) << 1)

def zero_vector(length, *data):
    """Return zero numpy vector of given length and type np.int32"""
    return np.zeros(length, dtype = np.int32)

 

#234567890123456789012345678901234567890123456789012345678901234567890

def pm_mat_from_function(f, l):
    """Create an  l times l matrix from a function f.

    The returned matrix m has entries m[i,j] = (-1)**f(i,j),  i, j 
    = 0,...,l-1. It is a numpy array of shape (l,l) and type np.int32.
    """
    a = np.zeros((l,l), dtype = np.int32)
    for i in range(l): 
        for j in range(l): 
            a[i,j] =  1 - 2 * (f(i,j) & 1)
    return a


def pm_diag_from_function(f, l):
    """Create an  l times l diagonal matrix from a function f.

    The returned matrix m has entries m[i,i] = (-1)**f(i),   
    i = 0,...,l-1, and zero entries for i != j.
    It is a numpy array of shape (l,l) and type np.int32.
    """
    a = np.zeros((l,l), dtype = np.int32)
    for i in range(l): 
        a[i,i] =  1 - 2 * (f(i) & 1)
    return a

