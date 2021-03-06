from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix
from mmgroup.structures.qs_matrix import binary
from mmgroup.structures.qs_matrix import qs_pauli_matrix
from mmgroup.structures.qs_matrix import qs_unit_matrix
from mmgroup.structures.qs_matrix import pauli_vector_mul
from mmgroup.structures.qs_matrix import pauli_vector_exp

#####################################################################
# Create test matrices
#####################################################################


def rand_unitary_matrix(n, scalar = True):
    m = qs_unit_matrix(n)
    e = randint(-16,16) if scalar else 0
    m.mul_scalar(e, randint(0,7))
    if n == 0:
        return m
    mask = (1 << n) - 1
    for i in range(n+2):
        v = randint(0, mask-1)
        vc =  randint(0, mask) & ~v
        m.gate_ctrl_not(vc, v)
        v1 = randint(0, mask)
        v2 =  randint(0, mask)
        m.gate_ctrl_phi(v1, v2)
    m.gate_phi(randint(0, mask), 1)
    m.gate_phi(randint(0, mask), 2)
    m.gate_h(randint(0, mask))
    m.gate_not(randint(0, mask))
    for i in range(n+2):
        v = randint(0, mask-1)
        vc =  randint(0, mask) & ~v
        m.gate_ctrl_not(vc, v)
        v1 = randint(0, mask)
        v2 =  randint(0, mask)
        m.gate_ctrl_phi(v1, v2)
    return m.echelon()

def rand_pauli_vectors(n, length):
    mask = (1 << (2*n+2)) - 1
    return [randint(0, mask) for i in range(length)]

def create_conjugate_data():
    yield rand_unitary_matrix(0), [0,1,2,3]
    yield qs_unit_matrix(1), [0,1,2,3]
    
    for n in range(1, 13):
        for i in range(10):
            m = rand_unitary_matrix(n)
            v = rand_pauli_vectors(n, max(10, n))
            yield m, v



#####################################################################
# Test method pauli_conugate()
#####################################################################


    
@pytest.mark.qstate
def test_pauli_conjugate(verbose = 0):
    """Test the conjugation of Pauli matrix with unitary matrix"""
    for ntest, (m, v) in enumerate(create_conjugate_data()):
        n = m.shape[0]
        p = [qs_pauli_matrix(n, x) for x in v]
        pv = [x.pauli_vector() for x in p]
        ok = pv == v
        if verbose or not ok:
            mm = m.copy()
            print("\nTest %d: v = %s" % (ntest, binary(v[0],0,2*n+2)))
            print("m =", mm.reduce_matrix(), mm)
            print("Pauli matrix of v", p[0])
            if not ok:
                print("Input and recomputed Pauli vectors:")
                for i, x in enumerate(v):
                    bv =  binary(x,0,2*n+2)
                    bpv =  binary(pv[i],0,2*n+2)
                    print(bv, ",", bpv)
                if not ok:
                    err = "Error in recomputation of Pauli vector"
                    raise ValueError(err)
        mi = m.inv()
        assert m @ mi == qs_unit_matrix(n), (str(m), str(mi), str(m@mi))

        w = m.pauli_conjugate(v)
        w_ref = [ (m @ x @ mi).pauli_vector() for x in p]
        if verbose or w != w_ref:
            for i, w_i in enumerate(w):
                w_i_ref = w_ref[i]
                if w_i_ref != w_i:
                   s = "Pauli vector: %s, conjugated: %s, obtained: %s"
                   print (s % (hex(v[i]), hex(w_i_ref), hex(w_i)))
                   err = "Pauli vector conjugation error"
                   raise ValueError(err)
                elif verbose:
                   s = "Pauli vector: %s, conjugated: %s"
                   print (s % (hex(v[i]), hex(w_i)))
        

#####################################################################
# Create test matrices for Pauli group operation
#####################################################################


def create_pauli_vectors():
    yield 0, 2, 3
    yield 4, 0x1, 0x10
    yield 16, 0x100000001, 0x200010000
    for n in range(0,21):
        for i in range(4):
            v = rand_pauli_vectors(n, 2)
            yield n, v[0], v[1]

@pytest.mark.qstate
def test_pauli_multiplication(verbose = 0):
    for ntest, (n, v1, v2) in enumerate(create_pauli_vectors()):
        p1 = qs_pauli_matrix(n, v1)
        p2 = qs_pauli_matrix(n, v2)
        # Check product v1 * v2
        v3 = pauli_vector_mul(n, v1, v2)
        p3 = qs_pauli_matrix(n, v3)
        v3_ref = (p1 @ p2).pauli_vector()
        ok = v3 == v3_ref
        if verbose or not ok:
            print("Test %d, dim = %d" % (ntest, n))
            print("v1 = %s, v2 = %s, v1 * v2 = %s" % 
                (hex(v1), hex(v2), hex(v3_ref)) )
            if not ok:
                print("m1", p1)
                print("m2", p2)
                print("m1 * m2", p1 @ p2)
                print("pauli vector obtained: ", hex(v3))
                err = "Pauli vector multiplication failed"
                raise ValueError(err)
        # Check powers of v1
        p_e_list = [qs_unit_matrix(n), p1, p1 @ p1, p1.H]
        v_e_ref_list = [m.pauli_vector() for m in p_e_list]
        v_e_list = [pauli_vector_exp(n, v1, e) for e in range(4)]
        ok = v_e_list == v_e_ref_list
        if verbose or not ok:
            print("v1 = ", hex(v1))
            for e, v_ref in enumerate(v_e_ref_list):
                if ok:
                    print("v1 ** %d = %s" % (e, hex(v_ref)))
                else:
                    print("v1 ** %d = %s, obtained: %s, v1 = %s" % 
                         (e, hex(v_ref), hex(v_e_list[e]), v1))
            if not ok:
                err = "Pauli vector exponentiation failed"
                raise ValueError(err)
        



#####################################################################
# Test matrix exponentiation and order computation
#####################################################################


def create_exp_data():
    yield qs_unit_matrix(4), -3
    for n in range(1, 13):
        for i in range(10):
            m = rand_unitary_matrix(n)
            m.mul_scalar(randint(-2, 2))
            e = randint(-10000, 10000)
            yield m, e
   


def ref_power(m, e):
    if e > 1:
        m1 = ref_power(m, e >> 1)
        m1 = m1 @ m1
        return m1 @ m  if e & 1 else m1 
    if e == 1:
        return m
    if e == 0:
        return qs_unit_matrix(m.shape[0])
    if e < 0:
        mi = m.inv()
        assert m @ mi == qs_unit_matrix(m.shape[0])
        return  ref_power(mi, -e)
    


@pytest.mark.qstate
def test_matrix_power(verbose = 0):
    MAX_ORDER = (2**8-1)*(2**6-1)*2**10
    """Test the conjugation of Pauli matrix with unitary matrix"""
    for ntest, (m, e) in enumerate(create_exp_data()):
        me = m.power(e)
        me_ref = ref_power(m, e)
        ok = me == me_ref
        if verbose or not ok:
            mm = m.copy()
            print("\nTest %d: exponent = %d, " % (ntest, e))
            print("m =", m.reduce())
            print("Result of exponentiation:", me_ref)
            if not ok:
                print("Obtained:", me)
                raise ValueError("Matrix exponentiation failed")
        nqb = m.shape[0]
        m.mul_scalar(-m.lb_norm2())
        assert m.H == m.inv() 
        if nqb > 4:
            continue
        order = m.order(MAX_ORDER)
        if verbose:
            s = print("Scaled matrix m has order %d" % order)
        assert m.power(order) == qs_unit_matrix(m.shape[0])
                
