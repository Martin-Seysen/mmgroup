from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



from random import randint #, shuffle, sample
from functools import reduce
from operator import __or__

import numpy as np
import pytest

from mmgroup.structures.qs_matrix import QStateMatrix
from mmgroup.structures.qs_matrix import binary
from mmgroup.structures.qs_matrix import qstate12_reduce_matrix
from mmgroup.structures.qs_matrix import qstate_pauli_matrix
from mmgroup.structures.qs_matrix import qstate_unit_matrix

#####################################################################
# Create test matrices
#####################################################################


def rand_unitary_matrix(n, scalar = True):
    m = qstate_unit_matrix(n)
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
    yield qstate_unit_matrix(1), [0,1,2,3]
    
    for n in range(1, 13):
        for i in range(10):
            m = rand_unitary_matrix(n)
            v = rand_pauli_vectors(n, max(10, n))
            yield m, v



#####################################################################
# Test method qstate12_bit_matrix_t()
#####################################################################


    
@pytest.mark.qstate
def test_pauli_conjugate(verbose = 0):
    """Test the conjugation of Pauli matrix with unitary matrix"""
    for ntest, (m, v) in enumerate(create_conjugate_data()):
        n = m.shape[0]
        p = [qstate_pauli_matrix(n, x) for x in v]
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
        assert m @ mi == qstate_unit_matrix(n), (str(m), str(mi), str(m@mi))

        w = m.pauli_conjugate(v)
        w_ref = [ (mi @ x @ m).pauli_vector() for x in p]
        if verbose or w != w_ref:
            for i, w_i in enumerate(w):
                w_i_ref = w_ref[i]
                if w_i_ref != w_i:
                   s = "Pauli vector: %s, conjugated: %s, obtained: %s"
                   print (s % (hex(v[i]), hex(w_i_ref), hex(w_i)))
                   err = "Pauli vector conjgation error"
                   #raise ValueError(err)
                elif verbose:
                   s = "Pauli vector: %s, conjugated: %s"
                   print (s % (hex(v[i]), hex(w_i)))
        




