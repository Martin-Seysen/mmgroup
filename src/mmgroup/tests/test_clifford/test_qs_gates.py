from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix, qs_rand_matrix
from mmgroup.structures.qs_matrix import qs_pauli_matrix
from mmgroup.structures.qs_matrix import qs_ctrl_not_matrix
from mmgroup.structures.qs_matrix import qs_phi_matrix
from mmgroup.structures.qs_matrix import qs_ctrl_phi_matrix
from mmgroup.structures.qs_matrix import qs_hadamard_matrix


from mmgroup.tests.test_clifford.test_qs_matrix import compare_complex

from mmgroup.bitfunctions import bitparity

#####################################################################
# Test gate_not
#####################################################################


gate_not_data = [
     (4, [0b000], 0b11),
     (4, [0b110_1001, 0b101_0111, 0b011_0100], 0b11),
]

def gate_not_testdata():
    for cols, data, v in gate_not_data:
        yield QStateMatrix(0, cols, data), v
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                v = randint(0, (1 << cols) - 1)
                yield  qs_rand_matrix(0, cols, r), v

def complex_gate_not(c1, v):
    v &= len(c1) - 1
    c2 = np.copy(c1)
    for i, x in enumerate(c1):
        c2[i ^ v] = x
    return c2
        
@pytest.mark.qstate
def test_gate_not(verbose = 0):
    """Test function ``qstate12_gate_not``. """
    for ntest, (m, v) in enumerate(gate_not_testdata()):        
        if verbose:
            print("TEST gate not %s" % (ntest+1))
            print(m)
            print("m.shape=", m.shape, "v=", hex(v))
        m1 = m.copy().gate_not(v) 
        c, c1 = m.complex().ravel(), m1.complex().ravel()
        c1_ref = complex_gate_not(c, v)   
        if verbose:
            print("output")
            print(m1)
        err = "Method gate_not has failed"
        compare_complex(c1_ref, c1, err)            
        gm = qs_pauli_matrix(m.shape[1], v << m.shape[1]) 
        assert m @ gm == m1, [str(x) for x in (m, gm, m @ gm, m1)]         
            
            
#####################################################################
# Test gate_ctrl_not
#####################################################################


gate_ctrl_not_data = [
     (1, [0b0], 1, 0),
     (4, [0b110_1001, 0b101_0111, 0b011_0100], 0b11, 0b11),
]

def gate_ctrl_not_testdata():
    for cols, data, v1, v2 in gate_ctrl_not_data:
        yield QStateMatrix(0, cols, data), v1, v2
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                vc = randint(0, (1 << cols) - 1)
                v = randint(0, (1 << cols) - 1)
                while bitparity(vc & v):
                    v = randint(0, (1 << cols) - 1)
                m = qs_rand_matrix(0, cols, r)
                yield m, vc, v


def parity32(v):
    v ^= v >> 16
    v ^= v >> 8
    v ^= v >> 4
    v ^= v >> 2
    v ^= v >> 1
    return v & 1

def complex_gate_ctrl_not(c1, vc, v):
    vc &= len(c1) - 1
    v &= len(c1) - 1
    c2 = 0 * c1
    for i, x in enumerate(c1):
        c2[i ^ (v * parity32(vc & i))] += x
    return c2
  
@pytest.mark.qstate
def test_gate_ctrl_not(verbose = 0):
    """Test function ``qstate12_gate_ctrl_not``. """
    for ntest, (m, vc, v) in enumerate(gate_ctrl_not_testdata()):
        if verbose:
            print("TEST gate ctrl_not %s" % (ntest+1))
            print(m)
            print("vc=", hex(vc))
            print("v=", hex(v))
        m1 = m.copy().gate_ctrl_not(vc, v) 
        c, c1 = m.complex().ravel(), m1.complex().ravel()
        c1_ref = complex_gate_ctrl_not(c, vc, v)    
        if verbose:
            print("output")
            print(m1)
        err = "Method gate_ctrl_not has failed"
        compare_complex(c1_ref, c1, err) 
        gm = qs_ctrl_not_matrix(m.shape[1], vc, v)  
        assert m @ gm == m1, [str(x) for x in (m, m @ gm, m1)]         



#####################################################################
# Test gate_phi
#####################################################################



def gate_phi_testdata():
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                v = randint(0, (1 << cols) - 1)
                phi = randint(0,4); 
                m = qs_rand_matrix(0, cols, r)
                yield m, v, phi


phi_mult = [1, 1j, -1, -1j]

def complex_gate_phi(c1, v, phi):
    v &= len(c1) - 1
    f = phi_mult[phi & 3]
    c2 = np.copy(c1)
    for i, x in enumerate(c1):
        c2[i] = x * (f if parity32(i & v) else 1)
    return c2


@pytest.mark.qstate
def test_gate_phi(verbose = 0):
    """Test function ``qstate12_gate_phi``. """
    for ntest, (m, v, phi) in enumerate(gate_phi_testdata()):
        assert 0 <= v < 1 << m.shape[1]
        if verbose:
            print("TEST gate phi %s" % (ntest+1))
            print(m)
            print("v=", hex(v))
            print("phi=", phi & 3)
        m1 = m.copy().gate_phi(v, phi) 
        c, c1 = m.complex().ravel(), m1.complex().ravel()
        c1_ref = complex_gate_phi(c, v, phi)   
        if verbose:
            print("output")
            print(m1)
        err = "Method gate_phi has failed"
        compare_complex(c1_ref, c1, err)            
        gm = qs_phi_matrix(m.shape[1], v, phi)  
        assert m @ gm == m1, [str(x) for x in (m, m @ gm, m1)]
        if (phi & 3 == 2):         
            gm = qs_pauli_matrix(m.shape[1], v)  
            assert m @ gm == m1, [str(x) for x in (m, m @ gm, m1)]         


#####################################################################
# Test gate_ctrl_phi
#####################################################################

gate_ctrl_phi_data = [
     (1, [0b0], 1, 1),
     (2, [0b11], 0b10, 0b11),
]

def gate_ctrl_phi_testdata():
    for cols, data, v1, v2 in gate_ctrl_phi_data:
        yield QStateMatrix(0, cols, data), v1, v2
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                v1 = randint(0, (1 << cols) - 1)
                v2 = randint(0, (1 << cols) - 1)
                m = qs_rand_matrix(0, cols, r)
                yield m, v1, v1



def complex_gate_ctrl_phi(c1, v1, v2):
    v1 &= len(c1) - 1
    v2 &= len(c1) - 1
    c2 = np.copy(c1)
    for i, x in enumerate(c1):
        f = parity32(i & v1) * parity32(i & v2) 
        c2[i] = -x if f else x
    return c2


@pytest.mark.qstate
def test_gate_ctrl_phi(verbose = 0):
    """Test function ``qstate12_gate_ctrl_phi``. """
    for ntest, (m, v1, v2) in enumerate(gate_ctrl_phi_testdata()):
        if verbose:
            print("TEST gate ctrl phi %s" % (ntest+1))
            print(m)
            print("v1=", hex(v1))
            print("v2=", hex(v2))
        m1 = m.copy().gate_ctrl_phi(v1, v2) 
        c, c1 = m.complex().ravel(), m1.complex().ravel()
        c1_ref = complex_gate_ctrl_phi(c, v1, v2)   
        if verbose:
            print("output")
            print(m1)
        err = "Method gate_ctrl_phi has failed"
        compare_complex(c1_ref, c1, err)            
        gm = qs_ctrl_phi_matrix(m.shape[1], v1, v2)  
        assert m @ gm == m1, [str(x) for x in (m, m @ gm, m1)]         



#####################################################################
# Test gate_h
#####################################################################



gate_h_data = [
     (1, [0b0], 1),
     (2, [0b11], 0b11),
]

def gate_h_testdata():
    for cols, data, v in gate_h_data:
        yield QStateMatrix(0, cols, data), v
    for cols in list(range(6)):
        for r in range(cols+3):
            for i in range(3):
                v = randint(0, (1 << cols) - 1)
                m = qs_rand_matrix(0, cols, r)
                yield m, v



def complex_one_gate_h(c1, d):
    c2 = np.copy(c1)
    for i, x in enumerate(c1):
        f = -1 if bool(i & d) else 1
        c2[i] = f * c1[i] + c1[i ^ d]
    return c2


def complex_gate_h(c1, v):
    exp = 0.0
    d = 1
    while (d < len(c1)):
        if v & d:
           c1 = complex_one_gate_h(c1, d)
           exp -= 0.5
        d <<= 1
    return c1 * 2.0**exp


@pytest.mark.qstate
def test_gate_h(verbose = 0):
    """Test function ``qstate12_gate_ctrl_phi``. """
    for ntest, (m, v) in enumerate(gate_h_testdata()):
        if verbose:
            print("TEST gate ctrl phi %s" % (ntest+1))
            print(m)
            print("v=", hex(v))
        m1 = m.copy().gate_h(v) 
        c, c1 = m.complex().ravel(), m1.complex().ravel()
        c1_ref = complex_gate_h(c, v)   
        if verbose:
            print("output")
            print(m1)
        err = "Method gate_ctrl_phi has failed"
        compare_complex(c1_ref, c1, err)            
        gm = qs_hadamard_matrix(m.shape[1], v)  
        assert m @ gm == m1, [str(x) for x in (m, m @ gm, m1)]         




