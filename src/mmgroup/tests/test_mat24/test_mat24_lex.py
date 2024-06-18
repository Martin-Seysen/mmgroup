"""Test lexical ordering of the Mathieu group M_24

The C function mat24_m24num_to_perm() in file mat24_functions.c 
converts a number 0 <= n < 244823040 to a permutation in the Mathieu 
group M_24. The C function mat24_perm_to_m24num() in the same file 
reverses that operation.

The numbering of the elements of M_24 provided by these functions is
lexicographic (based on the ordering of the permutations corresponding
to the elements of M_24). In this script we test if that ordering is
indeed lexicographic. Therefore we generate a large number of random
numbers  0 <= n < 244823040 and compute the permutations for these 
numbers. The we sort the list of the generated numbers and we check
that the list of permutations corresponding to that sorted list of 
numbers is also sorted.

Furthermore, function  mat24_m24num_to_perm() has two reference 
implementations in python and function mat24_perm_to_m24num()
has one reference implementation in python. For a small fraction 
of the test data we also check that the the (slow) reference 
implementations of the two C functions compute the same result as the
(fast) C functions.

""" 

from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals




from random import randint
import pytest

# Import the extension containing the C functions
from mmgroup import mat24
# Import two reference implementations for mat24_m24num_to_perm() 
from mmgroup.dev.mat24.mat24heptad import py_mat24_int_to_perm
from mmgroup.dev.mat24.mat24heptad import mat24_int_to_perm
# Import a reference implementation for mat24_perm_to_m24num()
from mmgroup.dev.mat24.mat24heptad import mat24_perm_to_int


MAT24_ORDER =  244823040 


def mat24lex_testcases(ntests = 1000):
    """Yield random test numbers 0 <= n < 244823040

    """
    # Yield a large number of random permutation numbers
    for  i in range(ntests):
        yield randint(0, MAT24_ORDER - 1)
    # Yield the lowest 10 permutation numbers
    yield from range(10)
    # Yield the highest 10 permutation numbers
    yield from range(MAT24_ORDER - 10, MAT24_ORDER)
    # Yield a random interval of 100 adjacent numbers
    r = randint(0, MAT24_ORDER - 100)
    for i in range(r, r+100):
        yield i


@pytest.mark.mat24
def test_mat24lex(ntests = 5000, verbose = 0):
    lindex = []; ldata = []
    for i, n in enumerate(mat24lex_testcases(ntests)):
        # Do  additional checks for the first few samples
        test_ref = i < 100
        # Compute permutation p with number n using the cython wrapper
        # mat24.m24num_to_perm() for function mat24_m24num_to_perm()
        p = mat24.m24num_to_perm(n)
        # For the first few random test data do:
        if test_ref:
            # Check result p against two slow reference implementations
            p_py_ref = py_mat24_int_to_perm(n)
            assert list(p) == list(p_py_ref), (n, p, p_py_ref)
            p_ref = mat24_int_to_perm(n)
            assert list(p) == list(p_ref), (n, p, p_ref)
        if verbose:
            print(i, n, p[:10])
        # Append pair (n i) to the list 'lindex' and put ldata[i] = p
        lindex.append((n, i))
           ### p[2] = randint(0,23) # This would destroy the order
        ldata.append(p)
        # Compute number n1 of permutation p using the cython wrapper
        # mat24.perm_to_m24num() for function mat24_perm_to_m24num()
        n1 = mat24.perm_to_m24num(p)
        # Then n == n1 must hold
        assert n1 == n, (hex(n),hex(n1), p)
        # For the first few random test data do:
        if test_ref:
            # Check result n1 against a slow reference implementation
            n1_ref = mat24_perm_to_int(p)
            assert n1_ref == n, (hex(n),hex(n1_ref), p)

    # Sort the index list 'lindex' by permutation numbers.
    lindex.sort()
    # Let n0 be the first permutation number in 'lindex' and
    # let p0 be the corresponding permutation.
    n0, index0 = lindex[0]
    p0 = ldata[index0]
    # For all subsequent permutation numbers n1 in the list 'lindex'
    # let p1 be the corresponding permutation.
    for n1, index1 in lindex[1:]:
        p1 = ldata[index1]
        if n1 > n0:
            # Check that permutation p1 is greater the the previous
            # permutation p0 (assuming that n1 is greater than the
            # previous number n0).
            assert p1 > p0
        elif n1 == n0:
            # Birthday paradox: n1 == n0 may (and will) happen
            assert p1 == p0
        else:
            raise ValueError("Something is going wrong here")
        # Update 'previous' number n0 and permutation p0
        n0, p0 = n1, p1


