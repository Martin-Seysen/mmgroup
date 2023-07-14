from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import os
import sys
import collections
import numpy as np
from random import randint 

import pytest





from mmgroup.dev.hadamard.hadamard_t import bitparity
from mmgroup.dev.mm_op.mm_op import MM_Op, INT_BITS

from mmgroup.tests.test_hadamard.matrices import hadamard_matrix
from mmgroup.tests.test_hadamard.matrices import msym16
from mmgroup.tests.test_hadamard.matrices import mat_l_16
from mmgroup.tests.test_hadamard.cythonize import SIZES 
from mmgroup.tests.test_hadamard.cythonize import PRIMES 
from mmgroup.tests.test_hadamard.cythonize import build, kill
from mmgroup.tests.spaces.sparse_mm_space import SparseMmV
from mmgroup.tests.groups.mgroup_n import MGroupNWord


NP_DTYPE = np.uint32 if INT_BITS == 32 else np.uint64


V = SparseMmV
G = MGroupNWord


################################################################
# Build external module with code to be tests
################################################################

test_hadamard = None

def build_code():
    """Build Cython code, import is as module test_hadamard

    This function builds the Cython module
    mmgroup.tests.test_hadamard.temp.test_hadamard
    an imports that module as test_hadamard.

    Since that build process takes some time, we do it on
    demand only.
    """
    global test_hadamard 
    if test_hadamard is None:
        build()
        print("importing  test_hadamard...")
        from mmgroup.tests.test_hadamard.temp import test_hadamard as t
        print("import done")
        test_hadamard = t

################################################################
# Auxiliary functions
################################################################

LOG_BIN = {4:2, 16:4, 64:6}

def pwr2modp(k, p):
    """Return 2**k mod p for any integer k"""
    if k < 0:
        assert p & 1
        return pow((p + 1) >> 1, -k, p)
    return pow(2, k, p)

################################################################
# Wrapper for functions in library test_hadamard.pyx
################################################################

def matrix_function(name, p, *data, **kwds):
    def pack_input_vector(data, dest, offset):
        for i in range(VECTOR_INTS):
            value = 0
            part = data[i * INT_FIELDS : (i+1) * INT_FIELDS]
            for j, x in enumerate(part):
                value |= (int(x) % p) << (j * FIELD_BITS)
            dest[offset + i] = value

    try:
        verbose = int(kwds["verbose"])
    except:
        verbose = 0
    basics = MM_Op(p = p)
    vlen = len(data[0])
    if name[-1:].isdigit():
        function_name = "mod%d_%s_test" % (p, name)
    else:
        function_name = "mod%d_%s%d_test" % (p, name, vlen)
    f =  test_hadamard.__dict__[function_name]
    FIELD_BITS = basics.FIELD_BITS
    INT_FIELDS = basics.INT_FIELDS
    VECTOR_INTS = max(1, vlen // INT_FIELDS)
    w_len = max(10, len(data)) * VECTOR_INTS
    w = np.zeros(w_len, dtype = NP_DTYPE)
    for i, d in enumerate(data):
        pack_input_vector(d, w, i * VECTOR_INTS)
    if verbose > 1:
        print("Input to function %s, p = %d, len = %d:" % (
           name, p, vlen))
        print([hex(x) for x in w[:VECTOR_INTS]])
    f(w)
    if verbose > 1:
        print("Output of function %s:" % name)
        print([hex(x) for x in w[:VECTOR_INTS]])
    result = np.zeros(vlen, dtype = np.int32)
    for i in range(vlen):
        index, sh = divmod(i, INT_FIELDS)
        o =  int(w[index]) >> (sh * FIELD_BITS)
        result[i] = ((o & p) % p)
    return result


################################################################
# Testing multiplication with Hadamard matrix
################################################################



def hadamard_function(p, vector, verbose = 0):
    return matrix_function("hadamard", p, vector, verbose=verbose)


def one_test_hadamard_function(p, vector, verbose = 0):
    if verbose:
        print("\nTest multiplication v * H with Hadamard matrix H, v =")
        print(vector)
    vector = np.array(vector, dtype = np.int32) % p
    obtained = hadamard_function(p, vector, verbose = verbose)
    k = LOG_BIN[len(vector)]
    h_matrix = hadamard_matrix(k) * pwr2modp(-k >> 1, p) % p
    expected = (vector @ h_matrix) % p
    ok = (obtained == expected).all()
    if verbose or not ok:
        if not ok:
            print("Multiplication v * H with Hadmard matrix H, v =")
            print(vector)
        print("Result obtained:")
        print(obtained)
        if not ok:
            print("Result expected:")
            print(expected)
            raise ValueError("Wrong result")


def hadamard_function_testcases(n_cases = 5):
    test_cases = [
        (3,  [1,0,0,0]),
        (3,  [0,1,0,0]),
        (3,  [0,0,1,0]),
        (3,  [0,0,0,1]),
    ]
    for t in test_cases:
        yield t
    for i in range(n_cases):
        for p in PRIMES:
            for l in SIZES:
                v = [randint(0, p - 1) for l in range(1 << l)]
                yield p, v 

@pytest.mark.compiler
@pytest.mark.slow
@pytest.mark.hadamard
def test_hadamard_function(n_cases = 5, verbose = 0):
    build_code()
    print("Test C code for multiplication with Hadamard matrix")
    for p, v in hadamard_function_testcases(n_cases):
        one_test_hadamard_function(p, v, verbose)
    print("Test passed")



################################################################
# Testing exchanging entries with odd parity
################################################################


def xch_parity_function(p, vector, verbose = 0):
    return matrix_function("xch_parity", p, vector, verbose=verbose)

def ref_xch_parity(vector):
    l = len(vector)
    res = np.zeros(l, dtype = np.int32)
    for i in range(l):
        j = l - i - 1 if bitparity(i) else i
        res[i] = vector[j]
    return res

def one_test_xch_parity_function(p, vector, verbose = 0):
    if verbose:
        print("\nTest exchanging entries of v with odd parity, v =")
        print(vector)
    vector = np.array(vector, dtype = np.int32) % p
    obtained = xch_parity_function(p, vector, verbose = verbose)
    k = LOG_BIN[len(vector)]
    expected = ref_xch_parity(vector) % p
    ok = (obtained == expected).all()
    if verbose or not ok:
        if not ok:
            s = "Exchanging entries of v (mod %d) with odd parity, v ="
            print(s % p)
            print(vector)
        print("Result obtained:")
        print(obtained)
        if not ok:
            print("Result expected:")
            print(expected)
            raise ValueError("Wrong result")

def xch_parity_function_testcases(n_cases = 5):
    test_p = [3, 7, 127]
    for p  in test_p:
        if not p in PRIMES:
            continue
        for i in range(64):
            lst = [0] * 64
            lst[i] = 1
            yield p, lst
    for i in range(n_cases):
        for p in PRIMES:
            v = [randint(0, p - 1) for l in range(64)]
            yield p, v 

@pytest.mark.compiler
@pytest.mark.slow
@pytest.mark.hadamard
def test_xch_parity_function(n_cases = 5, verbose = 0):
    build_code()
    print("Test C code for exchanging entries with odd parity")
    for p, v in xch_parity_function_testcases(n_cases):
        one_test_xch_parity_function(p, v, verbose)
    print("Test passed")



################################################################
# Testing operation t done by 64 times 64 matrices
################################################################



def op_t64_function(p, vector, exp, verbose = 0):
    exp = exp % 3
    if exp == 0: return vector
    x = [exp - 1]
    return matrix_function("op_t", p, vector, x, verbose=verbose)


    

def ref_op_t64(p, vector, exp):
    space = V(p)
    v = space()
    for i in range(64):
        v += int(vector[i]) * space('T',0, i)
    v *= G('t', exp)
    res = np.zeros(64, dtype = np.int32)
    for i in range(64):
        res[i] = v['T', 0, i]
    return res

def one_test_op_t64_function(p, vector, exp, verbose = 0):
    if verbose:
        print("\nTest op t64, exp = ", exp, ",  v =\n")
        print(vector)
    vector = np.array(vector, dtype = np.int32) % p
    obtained = op_t64_function(p, vector, exp, verbose = verbose)
    k = LOG_BIN[len(vector)]
    expected = ref_op_t64(p, vector, exp) 
    ok = (obtained == expected).all()
    if verbose or not ok:
        if not ok:
            s = "\nTest op t64 mod %d, exp = %d,  v ="
            print(s % (p, exp))
            print(vector)
        print("Result obtained:")
        print(obtained)
        if not ok:
            print("Result expected:")
            print(expected)
            raise ValueError("Wrong result")



def op_p64_function_testcases(n_cases = 5):
    for p  in PRIMES:
        for i in range(64):
            lst = [0] * 64
            lst[i] = 1
            yield p, lst, 1
            yield p, lst, 2
    for i in range(n_cases):
        for p in PRIMES:
            v = [randint(0, p - 1) for l in range(64)]
            yield p, v, 1 
            yield p, v, 2 


@pytest.mark.compiler
@pytest.mark.slow
@pytest.mark.hadamard
def test_op_p64_function(n_cases = 5, verbose = 0):
    build_code()
    print("Test operator t on tag 'T' done by 64 x 64 matrices")
    for p, v, exp1 in op_p64_function_testcases(n_cases):
        one_test_op_t64_function(p, v, exp1, verbose)
    print("Test passed")




################################################################
# Testing operation t done by 3 times 3 matrices
################################################################


def v24_ints(p):
    return MM_Op(p = p).V24_INTS

def int_fields(p):
    return MM_Op(p = p).INT_FIELDS


def op_t3_function(p, vector, exp, verbose = 0):
    if exp == 0: return vector
    assert len(vector) == 3 * int_fields(p)
    exp = exp % 3
    x = [exp - 1]
    return matrix_function("op_t3", p, vector, x, verbose=verbose)


def pr_a3(p, a):
    ld = len(str(p))
    fields =  int_fields(p)
    for i in range(0, 3*fields, fields):
        for j in range(fields):
            print("%*d" % (ld,a[i+j]), end = " ")
        print("")





def ref_op_t3(p, vector, exp):
    fields = int_fields(p)
    assert len(vector) == 3 * fields
    exp = exp % 3
    space = V(p)
    result = [None] * (3 * fields)
    for i in range(fields):
        v = space()
        for tag, j in [("A", 0), ("B", fields), ("C", 2 * fields)]:
            v += int(vector[i + j]) * space(tag, 1, 0)
        v *= G('t', exp)
        for tag, j in [("A", 0), ("B", fields), ("C", 2 * fields)]:
            result[i + j] = v[tag, 1, 0]
    return result



def one_test_op_t3_function(p, vector, exp, verbose = 0):
    if verbose:
        print("\nTest op t3, exp = ", exp, ",  v =\n")
        pr_a3(p, vector)
    vector = np.array(vector, dtype = np.int32) % p
    obtained = op_t3_function(p, vector, exp, verbose = verbose)
    expected = ref_op_t3(p, vector, exp) 
    ok = (obtained == expected).all()
    if verbose or not ok:
        if not ok:
            s = "\nTest op t3 mod %d, exp = %d,  v ="
            print(s % (p, exp))
            pr_a3(p,vector)
        print("Result obtained:")
        pr_a3(p, obtained)
        if not ok:
            print("Result expected:")
            pr_a3(p, expected)
            raise ValueError("Wrong result")


def op_p3_function_testcases(n_cases = 5):
    for p  in PRIMES:
        fields = int_fields(p)
        for i in range(3 * fields):
            lst = [0] * (3 * fields)
            lst[i] = 1
            yield p, lst, 1
            yield p, lst, 2
    for i in range(n_cases):
        for p in PRIMES:
            fields = int_fields(p)
            v = [randint(0, p - 1) for l in range(3 * fields)]
            yield p, v, 1 
            yield p, v, 2 


@pytest.mark.compiler
@pytest.mark.slow
@pytest.mark.hadamard
def test_op_p3_function(n_cases = 5, verbose = 0):
    build_code()
    print("Test operator t on tags 'ABC' done by 3 x 3 matrices")
    for p, v, exp1 in op_p3_function_testcases(n_cases):
        one_test_op_t3_function(p, v, exp1, verbose)
    print("Test passed")






################################################################
# test symmetric operation xi done by 64 times 64 matrices
################################################################

msym4a = np.array(
    [[1,1,1,1], [1,1,-1,-1], [1,-1,1,-1], [1,-1,-1,1]]
)

def op_xi_sym64_function(p, vector, verbose = 0):
    assert vector.shape == (16, 24)
    vector %= p    
    v1 = np.zeros( (16, 32), dtype = np.int32)
    v1[0:16, 0:24] = vector
    v1 = v1.reshape(16*32)
    v0 = np.zeros(8, dtype = np.int8)
    v2 = matrix_function("op_xi64", p, v1, [2], verbose=verbose)
    v2 = v2.reshape( (16, 32) )
    assert np.count_nonzero(v2[: , 24:]) == 0
    return v2[: , :24] % p


def ref_op_xi_sym64(p, vector):
    v = np.copy(vector)
    for i in range(16):
        for j in range(0, 24, 4):
            v[i, j:j+4] = v[i, j:j+4].dot(msym4a)
    for i in range(24):
        v[:,i] = mat_l_16.dot(v[:,i])
    q = pow((p + 1) >> 1, 3, p)
    return v * q % p
    


def one_test_op_xi_sym64_function(p, vector, verbose = 0):
    if verbose:
        print("\nTest op lsym64 mod %d,  v =\n" % p)
        print(vector)
    vector = np.array(vector, dtype = np.int32) % p
    obtained = op_xi_sym64_function(p, vector, verbose = verbose)
    expected = ref_op_xi_sym64(p, vector) 
    ok = (obtained == expected).all()
    if verbose or not ok:
        if not ok:
            print("\nTest op xi sym64 mod %d,  v =\n" % p)
            print(vector)
        print("Result obtained:")
        print(obtained)
        if not ok:
            print("Result expected:")
            print(expected)
            raise ValueError("Wrong result")


def op_xi_sym64_function_testcases(n_cases = 5):
    for p  in PRIMES:
        for i in range(16):
            for j in range(24):
                v = np.zeros((16,24), dtype = np.int32)
                v[i, j] = 1
                yield p, v
    for i in range(n_cases):
        for p in PRIMES:
            v = np.zeros((16,24), dtype = np.int32)
            for i in range(16):
                for j in range(24):   
                    v[i,j] = randint(0, p-1)
            yield p, v

@pytest.mark.compiler
@pytest.mark.slow
@pytest.mark.hadamard
def test_op_xi_sym64_function(n_cases = 5, verbose = 0):
    build_code()
    print("Test part of operator xi done by 64 x 64 matrices")
    for p, v in op_xi_sym64_function_testcases(n_cases):
        one_test_op_xi_sym64_function(p, v, verbose)
    print("Test passed")

################################################################
# test operation xi done by 64 times 64 matrices
################################################################


def op_xi64_function(p, vector, exp, verbose = 0):
    exp = exp % 3
    if exp == 0: return vector
    x = [exp - 1]
    assert vector.shape == (16, 24)
    vector %= p    
    v1 = np.zeros( (16, 32), dtype = np.int32)
    v1[0:16, 0:24] = vector
    v1 = v1.reshape(16 * 32)
    v2 = matrix_function("op_xi64", p, v1, x, verbose=verbose)
    v2 = v2.reshape( (16, 32) )
    assert np.count_nonzero(v2[: , 24:]) == 0
    return v2[: , :24] % p


def ref_op_xi64(p, vector, exp):
    exp = exp % 3
    if exp == 0: return vector
    space = V(p)
    v = space(0)
    for i in range(16):
        for j in range(24):
            v += int(vector[i,j]) * space('Y',i + 0x400, j)
    v *= G('l', exp)
    #print(v)
    v1 = np.zeros( (16, 24), dtype = np.int32)
    for i in range(16):
        for j in range(24):
            v1[i, j] = v['Y', i + 0x400, j]
    return v1
    

def one_test_op_xi64_function(p, vector, exp, verbose = 0):
    if verbose:
        print("\nTest op xi 64 mod %d, exp = %d,  v =\n" % (p, exp))
        print(vector)
    vector = np.array(vector, dtype = np.int32) % p
    obtained = op_xi64_function(p, vector, exp, verbose = verbose)
    expected = ref_op_xi64(p, vector, exp) 
    ok = (obtained == expected).all()
    if verbose or not ok:
        if not ok:
            print("\nTest op xi 64 mod %d, exp = %d, v =\n" % (p, exp))
            print(vector)
        print("Result obtained:")
        print(obtained)
        if not ok:
            print("Result expected:")
            print(expected)
            raise ValueError("Wrong result")


def op_xi64_function_testcases(n_cases = 5):
    if n_cases >= 100:
        for p  in PRIMES:
            for i in range(16):
                for j in range(24):
                    v = np.zeros((16,24), dtype = np.int32)
                    v[i, j] = 1
                    yield p, v, 1
                    yield p, v, 2
    for i in range(n_cases):
        for p in PRIMES:
            v = np.zeros((16,24), dtype = np.int32)
            for i in range(16):
                for j in range(24):   
                    v[i,j] = randint(0, p-1)
            yield p, v, 1
            yield p, v, 2

@pytest.mark.compiler
@pytest.mark.slow
@pytest.mark.hadamard
def test_op_xi64_function(n_cases = 5, verbose = 0):
    build_code()
    print("Test operator xi on tags 'YZ' done by 64 x 64 matrices")
    for p, v, exp in op_xi64_function_testcases(n_cases):
        one_test_op_xi64_function(p, v, exp, verbose)
    print("Test passed")






################################################################
# test operation xi done by 16 times 16 matrices
################################################################

def op_xi16_function(p, vector, exp, verbose = 0):
    exp = exp % 3
    if exp == 0: return vector
    x = [exp - 1]
    assert vector.shape == (4, 24)
    vector %= p    
    v1 = np.zeros( (4, 32), dtype = np.int32)
    v1[0:16, 0:24] = vector
    v1 = v1.reshape(4 * 32)
    v2 = matrix_function("op_xi16", p, v1, x, verbose=verbose)
    v2 = v2.reshape( (4, 32) )
    assert np.count_nonzero(v2[: , 24:]) == 0
    return v2[: , :24] % p



def ref_op_xi16(p, vector, exp):
    """Generate reference vector for 16 1 14 operation on tag 'A'"""
    assert (vector[:4,:4] == vector[:4,:4].T).all()
    # The previous assertion is a rather subtile symmetry condition
    # required for the python reference implementation to work.
    # Function symmetrize_test_matrix(vector) forces that symmetry.
    exp = exp % 3
    if exp == 0: return vector
    space = V(p)
    v = space()
    for i in range(4):
        for j in range(24):
            i0, j0 = max(i,j), min(i,j)
            if j >= 4 or i >= j:
                v += int(vector[i,j]) * space('A', i, j)
    v *= G('l', exp)
    #print(v)
    v1 = np.zeros( (4, 24), dtype = np.int32)
    for i in range(4):
        for j in range(24):
            i0, j0 = max(i,j), min(i,j)
            v1[i, j] = v['A', i0, j0]
    return v1


def one_test_op_xi16_function(p, vector, exp, verbose = 0):
    if verbose:
        print("\nTest op xi16 mod %d, exp = %d,  v =\n" % (p, exp))
        print(vector)
    vector = np.array(vector, dtype = np.int32) % p
    obtained = op_xi16_function(p, vector, exp, verbose = verbose)
    expected = ref_op_xi16(p, vector, exp) 
    ok = (obtained == expected).all()
    if verbose or not ok:
        if not ok:
            print("\nTest op xi 16 mod %d, exp = %d, v =\n" % (p, exp))
            print(vector)
        print("Result obtained:")
        print(obtained)
        if not ok:
            print("Result expected:")
            print(expected)
            raise ValueError("Wrong result")


def symmetrize_test_matrix(vector, i = 0):
   """Force symmetry of 'vector' required for function ref_op_xi16()"""
   i = i & ~ 3
   for j in range(4):
       for k in range(j):
           vector[i + j, k] = vector[i + k, j]


def op_xi16_function_testcases(n_cases = 5):
    if n_cases >= 0:
        for p in [7] + PRIMES:
            for i in range(4):
                for j in range(24):
                    v = np.zeros((4,24), dtype = np.int32)
                    v[i, j] = 1
                    symmetrize_test_matrix(v)
                    yield p, v, 1
                    yield p, v, 2
    for i in range(n_cases):
        for p in PRIMES:
            v = np.zeros((4,24), dtype = np.int32)
            for i in range(4):
                for j in range(24):   
                    v[i,j] = randint(0, p-1)
            symmetrize_test_matrix(v)
            yield p, v, 1
            yield p, v, 2


@pytest.mark.compiler
@pytest.mark.slow
@pytest.mark.hadamard
def test_op_xi16_function(n_cases = 5, verbose = 0):
    build_code()
    print("Test operator xi on tag 'A' done by 16 x 16 matrices")
    for p, v, exp in op_xi16_function_testcases(n_cases):
        one_test_op_xi16_function(p, v, exp, verbose)
    print("Test passed")



################################################################
# Clean up
################################################################


################################################################
# Main test program
################################################################



N_CASES = 20
if __name__ == "__main__":
    test_hadamard_function(N_CASES, verbose = 0)  
    test_xch_parity_function(N_CASES, verbose = 0)    
    test_op_p3_function(N_CASES, verbose = 0)
    test_op_p64_function(N_CASES, verbose = 0)
    test_op_xi16_function(N_CASES, verbose = 0)
    test_op_xi_sym64_function(N_CASES, verbose = 0)
    test_op_xi64_function(N_CASES, verbose = 0)
