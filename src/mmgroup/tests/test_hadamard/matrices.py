"""Implement parity-adjusted Hadamard matrices

A parity-adjusted 2**n times 2**n Hadamard matrix M is defined as
follows:

Assmme that rows and columns are numbered with integers 0,...,2**n-1.
Then we put M[i,j] = -1**(parity(i & j) ^ (parity(i) & parity(j))), 
where operator '&' and '^' are defined on integers i, j as in C, and 
parity(i) returns the bit parity of i interpreted as a binary
number.

A standard Hadamard matrix H is defined by H[i,j] = -1**parity(i & j).

It turns out that all non-monomial blocks (except for one 3 times 3 
block) occuring in the representation of the Monster group are 
parity-adjusted 2**n Hadamard matrices of size 4, 6 or 64; or they 
are Kronecker products (i.e. tensor products) of two such matrices.

In this module we define the relevant parity-adjusted Hadamard 
matrices. Some of these matrices are taken from module
autogroup.NonMonomialOp_l.
"""

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import os
import sys
import collections
import numpy as np 


from mmgroup.tests.spaces.rep_aux import pm_mat_from_function
from mmgroup.tests.spaces.rep_aux import pm_diag_from_function
from mmgroup.tests.spaces.sparse_mm_space import NonMonomialOp_t
from mmgroup.tests.spaces.sparse_mm_space import NonMonomialOp_l

from mmgroup.bitfunctions import bitparity

########################################################################
# Create Hadamard matrices and some of their variation.
########################################################################


def identity_matrix(n):
    """Return the n times n identity matrix"""
    m = np.zeros((n,n), dtype=np.int32)
    for i in range(n): m[i,i] = 1
    return m

def hadamard_matrix(lg_n):
    """Return a 2**lg_n times 2**lg_n Hadamard matrix"""
    n = 1 << lg_n
    f = lambda i, j: bitparity(i & j) 
    return pm_mat_from_function(f, n)


def scal_prod_parity(i, j):
    """Return sign of entry i, j of parity-adjusted Hadamard matrix  

    Then such a matrix M has entries 

         M[i,j] = (-1) ** scal_prod_parity(i, j) .
    """
    return bitparity(i & j) ^ (bitparity(i) & bitparity(j))

def parity_hadamard_matrix(lg_n):
    """Return parity-adjusted Hadamard matrix of size 2**lg_n.

    The matrix is returned as a numpy array of shape (n, n) with
    n = 2**lg_n and of type np.int32. 
    """
    n = 1 << lg_n
    return pm_mat_from_function(scal_prod_parity, n)

# Next we collect the candidates for parity-adjusted Hadamard matrices

msym16 = NonMonomialOp_l.MSYM16
mdiag16 = NonMonomialOp_l.MDIAG16
mat_l_16 =  mdiag16 @ msym16 @ mdiag16 



