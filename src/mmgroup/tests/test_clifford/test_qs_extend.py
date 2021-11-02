from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix, qs_rand_matrix

from mmgroup.tests.test_clifford.test_qs_matrix import compare_complex


#####################################################################
# Test method rot_bits of class QStateMatrix
#####################################################################

rot_data = [
     (4, [0b000], 3, 2, 1),
     (4, [0b110_1001, 0b101_0111, 0b011_0100], 2, 1, 1),
]

def rot_testdata():
    for cols, data, rot, nrot, n0 in rot_data:
        yield QStateMatrix(0, cols, data), rot, nrot, n0
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                nrot = n0 = cols
                while n0 + nrot > cols:
                    n0, nrot = randint(0, cols), randint(0, cols)
                rot = randint(-nrot, nrot)
                m = qs_rand_matrix(0, cols, r)
                m.mul_scalar(randint(-8, 8), randint(0,7))  
                yield m, rot, nrot, n0

def rot_index(index, nqb, rot, nrot, n0):
    if  nrot <= 1: return index
    assert nrot + n0 <= nqb
    assert 0 <= index < 1 << nqb
    mask = ((1 << nrot) - 1) << n0
    i0, i1 = index & ~mask, index & mask
    rot %= nrot
    i1 = (i1 << rot) + (i1 >> (nrot - rot))
    return i0 + (i1 & mask)

       
@pytest.mark.qstate
def test_rot(verbose = 0):
    """Test function ``qstate12_rot_bits``. """
    for ntest, (m, rot, nrot, n0) in enumerate(rot_testdata()):
        if verbose:
            print("TEST rot bits %s" % (ntest+1))
            print(m)
            print("rot=", rot, ", nrot =", nrot, ", n0 =", n0)
        m1 = m.copy().rot_bits(rot, nrot, n0) 
        c, c1 = m.complex().ravel(), m1.complex().ravel() 
        c1_ref = np.copy(c)
        for i, x in enumerate(c):
            c1_ref[rot_index(i, m.ncols, rot, nrot, n0)] = x
        if verbose:
            print("output")
            print(m1)
        err = "Method rot has failed"
        compare_complex(c1_ref, c1, err)            
            
            
#####################################################################
# Test method rot_bits of class QStateMatrix
#####################################################################


xch_bits_data = [
     (4, [0b000], 2, 0b11),
     (4, [0b110_1001, 0b101_0111, 0b011_0100], 1, 0b101),
]

def xch_bits_testdata():
    for cols, data, sh, mask in xch_bits_data:
        yield QStateMatrix(0, cols, data), sh, mask
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                m = qs_rand_matrix(0, cols, r)
                m.mul_scalar(randint(-8, 8), randint(0,7))  
                if cols == 0:
                    yield m, 0, 0
                else:
                    sh = randint(0, cols - 1)
                    mask = randint(0, (1 << (cols - sh)) - 1)
                    mask &= ~(mask << sh)
                    yield m, sh, mask


def xch_bits_index(index, nqb, sh, mask):
    assert 0 <= index < 1 << nqb
    assert mask & (mask << sh) ==  0
    assert ((mask << sh) & (-1 << nqb)) == 0
    diff = mask & (index ^ (index >> sh))
    return index ^ diff ^ (diff << sh)


@pytest.mark.qstate
def test_xch_bits(verbose = 0):
    """Test function ``qstate12_xch_bits``. """
    for ntest, (m, sh, mask) in enumerate(xch_bits_testdata()):
        if verbose:
            print("TEST xch bits %s" % (ntest+1))
            print(m)
            print("sh=", sh, ", mask =", mask)
        m1 = m.copy().xch_bits(sh, mask) 
        c, c1 = m.complex().ravel(), m1.complex().ravel() 
        c1_ref = np.copy(c)
        for i, x in enumerate(c):
            c1_ref[xch_bits_index(i, m.ncols, sh, mask)] = x
        if verbose:
            print("output")
            print(m1)
        err = "Method xch_bits has failed"
        compare_complex(c1_ref, c1, err)            


#####################################################################
# Test methods extend and extend_zero
#####################################################################



extend_data = [
     (4, [0b000], 2, 2),
     (4, [0b110_1001, 0b101_0111, 0b011_0100], 1, 3),
]


def extend_testdata():
    for cols, data, j, nqb in extend_data:
        yield QStateMatrix(0, cols, data), j, nqb
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                m = qs_rand_matrix(0, cols, r)
                m.mul_scalar(randint(-8, 8), randint(0,7)) 
                j = randint(0, cols)
                nqb = randint(0, 5) 
                yield m, j, nqb



def extend_complex(c, j, nqb, zero = False):
    c1 = np.zeros(len(c) << nqb, dtype = complex)
    mask = (1 << j) - 1
    for i, x in enumerate(c):
        c1[(i & mask) + ((i & ~mask) << nqb)] = x
    if not zero:
        for i0 in range(0, len(c1), 1 << (nqb + j)):
            for i1 in range(1 << j):
                i = i0 + i1
                x = c1[i]
                c1[i : i + (1 << (nqb + j)) : 1 << j] = x
    return c1


@pytest.mark.qstate
def test_extend(verbose = 0):
    """Test function ``qstate12_extend``. """
    for ntest, (m, j, nqb) in enumerate(extend_testdata()):
        if verbose:
            print("TEST extend %s" % (ntest+1))
            print(m)
            print("j=", j, ", nqb =", nqb)
        m1 = m.copy().extend(j, nqb) 
        c, c1 = m.complex().ravel(), m1.complex().ravel() 
        c1_ref = extend_complex(c, j, nqb)
        if verbose:
            print("output after extend")
            print(m1)
        err = "Method extend has failed"
        compare_complex(c1_ref, c1, err)            

        m1 = m.copy().extend_zero(j, nqb) 
        c, c1 = m.complex().ravel(), m1.complex().ravel() 
        c1_ref = extend_complex(c, j, nqb, zero = True)
        if verbose:
            print("output after extend_zero")
            print(m1)
        err = "Method extend_zero has failed"
        compare_complex(c1_ref, c1, err)            


#####################################################################
# Test methods restrict, restrict_zero and sumup
#####################################################################


restrict_data = [
     (4, [0b0000], 2, 2),
     (4, [0b110_1001, 0b101_0111, 0b011_0100], 1, 3),
]


def restrict_testdata():
    for cols, data, j, nqb in restrict_data:
        yield QStateMatrix(0, cols, data), j, nqb
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                m = qs_rand_matrix(0, cols, r)
                m.mul_scalar(randint(-8, 8), randint(0,7)) 
                j = randint(0, cols)
                nqb = randint(0, cols - j) 
                yield m, j, nqb


def restrict_complex(c, j, nqb, zero = False):
    if not zero:
        mask = (1 << j) - 1
        c1 = np.zeros(len(c) >> nqb, dtype = complex) 
        for i in range(len(c1)):
            c1[i] = c[(i & mask) + ((i & ~mask) << nqb)]
    else:
        mask = ((1 << nqb) - 1) << j
        c1 = np.copy(c)
        for i, x in enumerate(c1):
            if i & mask: 
                c1[i] = 0
    return c1


def sumup_complex(c, j, nqb):
    c1 = c.reshape(-1, 1 << nqb, 1 << j)
    c1 = np.einsum("ijk->ik", c1)
    return c1.reshape(-1)

@pytest.mark.qstate
def test_restrict(verbose = 0):
    """Test function ``qstate12_restrict``. """
    for ntest, (m, j, nqb) in enumerate(restrict_testdata()):
        if verbose:
            print("TEST restrict %s" % (ntest+1))
            print(m)
            print("j=", j, ", nqb =", nqb)
        m1 = m.copy().restrict(j, nqb) 
        c, c1 = m.complex().ravel(), m1.complex().ravel() 
        c1_ref = restrict_complex(c, j, nqb)
        if verbose:
            print("output after restrict")
            print(m1)
        err = "Method restrict has failed"
        compare_complex(c1_ref, c1, err)            

        m1 = m.copy().restrict_zero(j, nqb) 
        c, c1 = m.complex().ravel(), m1.complex().ravel() 
        c1_ref = restrict_complex(c, j, nqb, zero = True)
        if verbose:
            print("output after restrict_zero")
            print(m1)
        err = "Method restrict_zero has failed"
        compare_complex(c1_ref, c1, err)    


        m1 = m.copy().sumup(j, nqb) 
        c, c1 = m.complex().ravel(), m1.complex().ravel() 
        c1_ref = sumup_complex(c, j, nqb)
        if verbose:
            print("output after sumup")
            print(m1)
        err = "Method sumup has failed"
        compare_complex(c1_ref, c1, err)    

        
