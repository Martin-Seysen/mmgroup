from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from random import randint
import math

import pytest

from mmgroup.bitfunctions import bit_mat_random, bit_mat_inverse
from mmgroup.bitfunctions import bit_mat_mul, bit_mat_rank
from mmgroup.bitfunctions import lrange, lmap, bitlen, v2
from mmgroup.bitfunctions import bit_mat_orthogonal_complement
from mmgroup.bitfunctions import bit_mat_transpose
from mmgroup.bitfunctions import bit_mat_iszero


 
def inverse_testcases():
    testcases = [  [3,4,6,12], [2,3],
    ]
    NTESTS = 100
    for a in testcases:
         yield a
    for j in range(NTESTS):
        for n in lrange(1,13) + [24]: 
            a = bit_mat_random(n, n)
            yield a
    for j in range(200): yield  bit_mat_random(10, 10)


@pytest.mark.bitfunc
def test_inverse():
    print( "Testing matrix inversion ..." )
    MAXDIM_INVERSE_TEST = 20
    NSIGMA = 4.0
    elementary_testcases =  [  ([3,2], [3,2]) , ([2,4,1], [4,1,2]) ]
    for a, b in elementary_testcases:
        assert  bit_mat_inverse(a) == b
        if a != b: assert  bit_mat_inverse(b) == a


    n_cases = [0] *(MAXDIM_INVERSE_TEST+1) 
    n_ok = [0.0] *(MAXDIM_INVERSE_TEST+1) 
    for a in inverse_testcases():
        ok, n = False, len(a)
        if n > MAXDIM_INVERSE_TEST:
            continue
        n_cases[n] += 1
        try:
            inv = bit_mat_inverse(a)
            ok = True
            n_ok[n] += 1
        except:
            pass
        if ok:
            unit = [1 << i for i in range(n)]
            prod = bit_mat_mul(a, inv) 
            assert prod == unit, (lmap(hex,a), lmap(hex,inv), lmap(hex,prod))
    cases_displayed = [10]
    for n in range(1,MAXDIM_INVERSE_TEST):
        if n_cases[n] > 10:
            import math
            N, N_ok = n_cases[n], n_ok[n]
            p = 1.0
            for i in range(1,n+1): p *= 1.0 - 0.5**i
            mu =  N * p
            sigma = math.sqrt(N * p * (1.0 - p))
            if n in cases_displayed:
                print("Dim. %d, %d cases, ratio of invertibles: %.3f, expected: %.3f+-%.3f"  %(
                     n, N, N_ok/N, p, sigma/N))
            assert   abs(mu - N_ok) < NSIGMA *sigma, (n,p,mu,sigma,N,N_ok)
    print("Matrix inversion test passed")



def orthogonal_complement_testcases():
    testcases = [  ([6],5), ([3,4,6,12],None), ([2,3],None)
    ]
    for a, r in testcases:
         yield a, r
     
    for m in range(1,20):
        for n in range(1,20):
            a = bit_mat_random(m, n)
            yield a, n
   

@pytest.mark.bitfunc
def test_orthogonal_complement():
    for a, ncols in orthogonal_complement_testcases():
        ra = bit_mat_rank(a)
        ao = bit_mat_orthogonal_complement(a, ncols)
        prod = bit_mat_mul(a, bit_mat_transpose(ao))
        assert bit_mat_iszero(prod) , (a, ncols, ao)
        rao = bit_mat_rank(ao) 
        if ncols:
            assert ra + rao == ncols, (ra, rao, ncols, a, ao)
    print( "Test of othogonal complement computation passed" )


def v2_testcases():
    for i in lrange(32) + lrange(32,70,3):
         yield (1 << i)
         if i: yield (1 << i) - 1


def v2_old(x):
    v = 0
    assert x
    while x & 1 == 0:
        x, v = x >> 1, v + 1
    return v

@pytest.mark.bitfunc
def test_v2():
    for i in lrange(300)+[100000]:
        assert bitlen(1<<i) == i+1           
        assert bitlen((1<<i)-1) == i   
    for i in  v2_testcases():
         assert v2(i) == v2_old(i), (hex(i), v2(i), v2_old(i))
         pass


        
