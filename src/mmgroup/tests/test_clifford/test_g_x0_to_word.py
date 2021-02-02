

from random import randint  
from functools import reduce
from operator import __or__, __xor__
import time

import numpy as np
import pytest

from mmgroup import MM, Xsp2_Co1
from mmgroup import mat24
from mmgroup.mat24 import MAT24_ORDER 
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.structures.xsp2_co1 import Xsp2_Co1_Word
from mmgroup.clifford12 import xsp2co1_mul_elem_word
from mmgroup.clifford12 import xsp2co1_elem_monomial_to_xsp
from mmgroup.clifford12 import xsp2co1_xspecial_vector
from mmgroup.clifford12 import xsp2co1_elem_to_word


#######################################################################
# Function monomial_to_word
#######################################################################


def print_monomial_4096(monomial, explain = False):
    print("Monomial operation of element (up to sign):") 
    if explain:
        print("Unit vector 1 << i is mapped (row 0) XOR (row i)")
    for i in range(13):
        bits = [(monomial[i] >> j) & 1 for j in range(11,-1,-1)]
        print("  " + "".join(map(str, bits)))

def monomial_to_word(elem, verbose = 0):
    assert isinstance(elem, Xsp2_Co1_Word)
    assert isinstance(verbose, int) , verbose   
    if verbose:
        print("Convert monomial element g of G_x0 to word. g is:")
        print(elem)
    a = []
    qs_i = elem.qs_t
    monomial = qs_i.monomial_row_op()
    if verbose:
        print_monomial_4096(monomial, explain = True)
    y = (monomial[12] & 0x7ff)
    mat24.matrix_from_mod_omega(monomial[1:])
    perm = mat24.autpl_to_perm(monomial[1:])
    perm = mat24.inv_perm(perm)
    pi = mat24.perm_to_m24num(perm)
    if pi:
        a.append(pi + 0xA0000000)
    if y:
        a.append(y + 0xC0000000)
    a = np.array(a, dtype = np.uint32)
    if verbose:
        print("Multiplier is:", MM.from_data(a))
        print_monomial_4096(monomial, explain = True)
    return a


def monomial_to_word_C(elem):
    a = np.zeros(2, dtype = np.uint32)
    elem_data = np.array(elem.data, dtype = np.uint64)
    len_a = xsp2co1_elem_monomial_to_xsp(elem_data, a)
    return a[:len_a]


#######################################################################
# Function elem_to_word
#######################################################################



def elem_to_word(elem, verbose = 1):
    assert isinstance(elem, Xsp2_Co1_Word)
    assert isinstance(verbose, int) , verbose   
    group = elem.group
    len_a = 0
    a0 = np.zeros(16, dtype = np.uint32)
    img_Omega = elem.xsp_conjugate(0x800000)

    if verbose:
        print("Convert element g of G_x0 to word. Element is:")
        print(elem)
        print("g conjugates 0x800000 to %s" % hex(img_Omega)) 

    len_a =  gen_leech2_reduce_type4(img_Omega, a0)
    a0 = list(a0[:len_a])
    elem_reduced = elem.copy().mul_data(a0)

    if verbose:
        print("Word w0 stabilizing Omega:")
        print(" ".join([hex(x) for x in a0]))
        print("After multiplication with w0, the element is:")
        print(elem_reduced)
    
    a1 = list(monomial_to_word(elem_reduced, verbose > 1))
    elem_reduced = elem_reduced.mul_data(a1)

    x = elem_reduced.as_xsp()
    x ^= mat24.ploop_theta(x >> 12)
    a = []
    if (x & 0x1fff000):
        a.append(0x30000000 + (x >> 12))
    if (x & 0xfff):
        a.append(0x10000000 + (x & 0xfff))
    a +=  [x ^ 0x80000000 for x in (a0 + a1)[::-1]]
    a = np.array(a, dtype = np.uint32)
    if verbose:
        print("Word is:\n " + str(MM.from_data(a)))
    return a
 

def elem_to_word_C(elem):
    a = np.zeros(10, dtype = np.uint32)
    elem_data = np.array(elem.data, dtype = np.uint64)
    len_a = xsp2co1_elem_to_word(elem_data, a, 0)
    return a[:len_a]


   
#######################################################################
# Yield test vectors as list of generator tuples
#######################################################################

from mmgroup.mat24 import MAT24_ORDER


MAX_TAG_ENTRY = {
   'd':0xfff, 'x':0x1fff, 'y':0x1fff, 'l':2, 'p': MAT24_ORDER-1
}
                 
                

def rand_G_x0_testword(t):
    """Yield element of G_{x0} as a list of pairs (tag, data) 

    ``t`` is a string of characters of generator tags, each character 
    must be in "dxypl"
    """
    w = []
    for tag in t:
         w.append((tag, randint(1, MAX_TAG_ENTRY[tag])))
    return w



def make_testwords(ntests = 100, monomial = True):
    test_words = [
       [],
       [('x', 0x800)],
       [('x', 0x1), ('y', 0x1), ('d', 0x800)],
       [('d', 0x800), ('y', 0x1)],
       [('d', 0x436), ('y', 0x12)],
       [('d', 0x837), ('y', 0x12)],
       [('x', 0x436), ('p', 0)],
       [('x', 0x436), ('p', 33445)],
       [('p',1234)],
    ]
    for w in test_words:
        yield w 
    tag_sets = ['d', 'dx', 'dxp', 'dxpy', 'dxpyl']
    for i in range(1,5):
         tag_sets.append('dxp' + 'lp' * i)
    for n in range(ntests):
        for t in tag_sets:
            if not monomial or not 'l' in t:
                 yield  rand_G_x0_testword(t)


#######################################################################
# Test function monomial_to_word
#######################################################################

@pytest.mark.qstate
def test_monomial_to_word(ntests = 10, verbose = 0):
    print("Test function test_monomial_to_word()")
    for i, w in enumerate(make_testwords()):
        if verbose:
            print("\nTest %d:" % (i+1))
            m = MM.word(*w)
            print("Testing word\n%s" % m)
        elem = Xsp2_Co1(*w)
        if verbose:
            print("This is element\n%s" % elem)
        data_i = monomial_to_word(elem, verbose > 1)
        elem_0 = elem.copy().mul_data(data_i)
        try:
            x = elem_0.as_xsp()
            ok = True
        except:
            raise
            ok = False
        if verbose or not ok:
            m_i = MM.from_data(data_i)
            print("Reduced inverse word\n%s" % str(m_i))
            if not ok:
                print("Product with inverse\n%s" % elem_0)
                err = "WTF"
                raise ValueError(err)
            else:
                print("Result is:", hex(x))
            print("Test %d ok" % (i+1))
        data_C = monomial_to_word_C(elem)
        assert (data_C == data_i).all(), (data_C, data_i)


#######################################################################
# Test function elem_to_word
#######################################################################

@pytest.mark.qstate
def test_elem_to_word(ntests = 50, verbose = 0):
    print("Test function test_elem_to_word()")
    for i, w in enumerate(make_testwords(monomial=False)):
        if verbose:
            print("\nTest %d:" % (i+1))
            m = MM.word(*w)
            print("Testing word\n%s" % m)
        elem = Xsp2_Co1(*w)
        if verbose:
            print("This is element\n%s" % elem)
        word = elem_to_word(elem, verbose > 1)
        if verbose:
            print("Reduced word\n%s" % MM.from_data(word))
        elem_1 = Xsp2_Co1.from_data(word)
        #elem_1 = Xsp2_Co1()
        #for w in word: elem_1 *= Xsp2_Co1.from_data([w])
        #for w in word: elem_1.mul_data([w])
        #elem_1.mul_data(word)
        if verbose:
            print("Recomputed element\n%s" % elem_1)
        ok = (elem == elem_1).all()
        if verbose or not ok:
            if not ok:
                print("Instead of to 1, the element reduces to:")
                print(elem * elem_1**(-1))
                err = "WTF"
                raise ValueError(err)
        word_C = elem_to_word_C(elem)
        assert (word_C == word).all(), (word_C, word)
            
            
#######################################################################
# benchmark
#######################################################################

def one_benchmark(ntests):
    tags = "dxyplplplplplp"
    samples = [Xsp2_Co1(*rand_G_x0_testword(tags)) for i in range(400)]
    img_omega = [s.xsp_conjugate(0x800000) for s in samples]
    a = np.zeros(10, dtype = np.uint32)
    t_start = time.process_time()
    for elem, img in zip(samples, img_omega):
        elem_data = elem._data
        for i in range(ntests):
            assert xsp2co1_elem_to_word(elem_data, a, img) >= 0
    t = time.process_time() - t_start
    return t / (ntests * len(samples))


def benchmark():
    t = 1.0e6 * one_benchmark(100)
    print("\nTime for function xsp2co1_elem_to_word: %.3f us" % t)


#######################################################################
# Main program 
#######################################################################

if __name__ == "__main__":
    #test_monomial_to_word()
    #test_elem_to_word()    
    benchmark()
    
    