from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals

from mmgroup.bitfunctions import bw24, lmap

import random
import numpy
from numpy import array, zeros, uint8, uint16, uint32

import pytest



from mmgroup.dev.mat24.mat24_ref import Mat24
from mmgroup.dev.mat24_xi.mat24_xi_ref import Mat24Xi
from mmgroup import mat24_xi  

        

#######################################################################
# Test functions
#######################################################################




#######################################################################
# test_mat24_xi_leech
#######################################################################


def xi_samples(ntests):
    data = [1,0x1000, 2, 0x100, 0x800001 ]
    for x in data:
         yield x  
    for i in range(ntests):
        yield random.randint(0, 0x1ffffff)



TEST_PAR = [(Mat24Xi, 200), (mat24_xi, 2000)]

@pytest.mark.mat24_xi
@pytest.mark.parametrize("mat24_xi_class, ntests", TEST_PAR )
def test_mat24_xi_leech(mat24_xi_class, ntests):
    xi = mat24_xi_class.xi_op_leech 
    mul =  mat24_xi_class.xi_mul_leech 

    xl_old =  [0,0,0]
    for x in xi_samples(ntests):
        xl = [x, xi(x,1), xi(x,2)]
        assert xl[2]  == xi(xl[1], 1) 
        assert xl[1]  == xi(xl[2], 2) 
        assert xi(xl[2], 1) == xi(xl[1], 2) == x 
        for i in range(1,3):
           assert mul(xl[i], xl_old[i]) == xi(mul(x, xl_old[0]), i)
        xl_old = xl


#######################################################################
# test_mat24_xi_short
#######################################################################


def is_short_vector_code(x):
    box, code = x >> 16, x & 0x7fff
    if box in [4,5]:
        return code & 0x1f < 24
    if box == 3:
        return code < 24576 # this is 384 * 64
    if box == 2:
        return code < 23040 # this is 360 * 64
    if box == 1:
        if code >=  2496: # this is 2 * 24 * 32 + 15 * 64 
             return False
        if code >= 1536:
             return True 
        j = code & 0x1f
        if j >= 24:
             return False
        i = (code >> 5) % 24
        return i != j
    return False

def normalize_short_vector_code(x):
    box, sign, code = x >> 16, x & 0x8000, x & 0x7fff
    if box == 1 and code < 1536 and (code & 31) < 24:
        b, i, j = code >= 768, (code >> 5) % 24, code & 31
        i, j = max(i,j), min(i,j)
        return 0x10000 + sign + 768 * b + 32* i + j 
    return x
    
 


def rand_short_vector(xmin = 0x10000, xmax = 0x60000):
    while 1:
        x = random.randint(xmin, xmax-1)
        if is_short_vector_code(x):
            return x
        


def short_samples(ntests):
    data = [ 0x10030,  0x10001, 0x35000, 0x3ddb0 ]
    for x in data:
        yield x
    for i in range(max(10, ntests // 50)):
        x = rand_short_vector(0x10000, 0x10000+1536)
        yield x
        yield x | 0x8000
    for i in range(ntests):
        yield rand_short_vector()
 

TEST_PAR = [(Mat24Xi, 200), (mat24_xi, 2000)]


          
@pytest.mark.mat24_xi
@pytest.mark.parametrize("mat24_xi_class, ntests", TEST_PAR )
def test_mat24_xi_short(mat24_xi_class, ntests):
    xi = mat24_xi_class.xi_op_leech 
    xi_short = mat24_xi_class.xi_op_short 
    to_leech = mat24_xi_class.xi_short_to_leech
    to_short = mat24_xi_class.xi_leech_to_short
    for x in short_samples(ntests):
        xl = to_leech(x)
        x_norm = normalize_short_vector_code(x)
        assert to_short(xl) == x_norm, lmap(hex, [x, xl, x_norm, to_short(xl)])
        for i in range(1,3): 
            continue           
            y = xi_short(x,i) 
            assert y, lmap(hex, [x, xl, xi(xl,i), y])
            ref = to_short(xi(xl,i))
            assert y == ref, lmap(hex, [x, xl, xi(xl,i), ref, y])



#######################################################################
# test_mat24_xi_ref
#######################################################################


@pytest.mark.mat24_xi
@pytest.mark.parametrize("mat24_xi_class, ref", [(mat24_xi, Mat24Xi)])
def test_mat24_xi_ref(mat24_xi_class, ref, ntests=200):       
    xi, xi_ref = mat24_xi_class.xi_op_leech, ref.xi_op_leech 
    mul, mul_ref =  mat24_xi_class.xi_mul_leech, ref.xi_mul_leech 

    x_old = 0
    for x in xi_samples(ntests):
        for i in range(0,3):
           res, ref_ = xi(x, i), xi_ref(x, i)
           assert res == ref_ , lmap(hex,[i, x,res, ref_ ])
        assert mul(x, x_old) == mul_ref(x, x_old)
        assert mat24_xi_class.xi_g_gray(x) == ref.xi_g_gray(x)
        assert mat24_xi_class.xi_w2_gray(x) == ref.xi_w2_gray(x)
        assert mat24_xi_class.xi_g_cocode(x) == ref.xi_g_cocode(x)
        assert mat24_xi_class.xi_w2_cocode(x) == ref.xi_w2_cocode(x)
        x_old = x
    xi_short = mat24_xi_class.xi_op_short 
    to_leech = mat24_xi_class.xi_short_to_leech
    to_short = mat24_xi_class.xi_leech_to_short
    xi_short_ref = ref.xi_op_short 
    to_leech_ref = ref.xi_short_to_leech
    for x in short_samples(ntests):
        xl = to_leech_ref(x)
        assert to_leech(x) == xl, lmap(hex,[x, xl,to_leech(x)])
        x_norm = normalize_short_vector_code(x)
        assert to_short(xl) == x_norm, lmap(hex,[xl, x, x_norm, to_short(xl)])
        for i in range(1,3):         
            ref = xi_short_ref(x,i)   
            assert xi_short(x,i) ==  ref, lmap(hex,[i, x, xi_short(x,i), ref])
            pass


