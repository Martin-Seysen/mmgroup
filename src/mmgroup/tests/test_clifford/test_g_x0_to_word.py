

from random import randint  
from functools import reduce
from operator import __or__, __xor__
import time

import numpy as np
import pytest

from mmgroup import MM0
from mmgroup import mat24
from mmgroup.mat24 import MAT24_ORDER 
from mmgroup.generators import gen_leech2_reduce_type4
from mmgroup.structures.xsp2_co1 import Xsp2_Co1
from mmgroup.clifford12 import chk_qstate12
from mmgroup.clifford12 import xsp2co1_mul_elem_word
from mmgroup.clifford12 import xsp2co1_elem_monomial_to_xsp
from mmgroup.clifford12 import xsp2co1_xspecial_vector
from mmgroup.clifford12 import xsp2co1_elem_to_word
from mmgroup.clifford12 import xsp2co1_reduce_word
from mmgroup.clifford12 import xsp2co1_order_elem
from mmgroup.clifford12 import xsp2co1_half_order_word
from mmgroup.clifford12 import xsp2co1_order_word
from mmgroup.clifford12 import xsp2co1_power_elem



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
    assert isinstance(elem, Xsp2_Co1)
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
    assert isinstance(elem, Xsp2_Co1)
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
    len_a = chk_qstate12(xsp2co1_elem_to_word(elem_data, a))
    return a[:len_a]

def reduce_word_C(w):
    a = np.array(w, dtype = np.uint32)
    a1 = np.zeros(10, dtype = np.uint32)
    len_a1 = chk_qstate12(xsp2co1_reduce_word(a, len(a), a1))
    return a1[:len_a1]
    


   
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

@pytest.mark.xsp2co1
def test_monomial_to_word(ntests = 10, verbose = 0):
    print("Test function test_monomial_to_word()")
    for i, w in enumerate(make_testwords()):
        if verbose:
            print("\nTest %d:" % (i+1))
            m = MM0(w)
            print("Testing word\n%s" % m)
        elem = Xsp2_Co1(w)
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
            m_i = MM0("a", data_i)
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

@pytest.mark.xsp2co1
def test_elem_to_word(ntests = 50, verbose = 0):
    print("Test function test_elem_to_word()")
    for i, w in enumerate(make_testwords(monomial=False, ntests=ntests)):
        m = MM0(w)
        if verbose:
            print("\nTest %d:" % (i+1))
            print("Testing word\n%s" % m)
        elem = Xsp2_Co1(w)
        if i < 100:
            assert m == MM0(elem)
        if verbose:
            print("This is element\n%s" % elem)
        word = elem_to_word(elem, verbose > 1)
        if verbose:
            print("Reduced word\n%s" % MM0('a', word))
        elem_1 = Xsp2_Co1('a', word)
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
        word1_C = reduce_word_C(m.mmdata)
        assert (word1_C == word).all(), (word1_C, word)
            


#######################################################################
# Test computing the power of an element of G_x0
#######################################################################

def xsp2co1_ref_power(wx, e):
    """Foolproof exponentiation in G_x0"""
    assert isinstance(wx, Xsp2_Co1)
    if e > 1:
        h = xsp2co1_ref_power(wx, e >> 1)  
        return h * h * wx if e & 1 else h * h
    if e == 1:
        return wx
    if e == 0:
        return Xsp2_Co1()
    return (wx**(-1)) ** (-e)


def xsp2co1_fast_power(wx, e):
    """The safe exponentiation in G_x0 to be tested"""
    assert isinstance(wx, Xsp2_Co1)
    power = Xsp2_Co1()
    chk_qstate12(xsp2co1_power_elem(wx._data, e, power._data))
    return power

@pytest.mark.xsp2co1
def test_elem_power(ntests = 50, verbose = 0):
    print("Test function computation a of power in G_x0")
    for i, w in enumerate(make_testwords(monomial=False, ntests=ntests)):
        e = randint(-2**35, 2**35)
        #print("Test ", i, "e", e)
        wx = Xsp2_Co1(w)
        power =  xsp2co1_fast_power(wx, e)
        ref_power = xsp2co1_ref_power(wx, e) # wx ** e
        assert power == ref_power, e
        if i < 10:
            e1 = i - 5
            power =  xsp2co1_fast_power(wx, e1)
            assert power == xsp2co1_ref_power(wx, e1) 


#######################################################################
# Test computing the order of an element of G_x0
#######################################################################


def xsp2co1_ref_order(wx):
    assert isinstance(wx, Xsp2_Co1)
    o = wx.qs.order(120)
    if o & 1 == 0:
        o = o >> 1
    unit, pw = wx.group(), wx**o
    if pw == unit:
         return o
    for i in range(2):
        o, pw = 2*o, pw * pw
        if pw == unit:
            return o
    err = "Order of QStateMatrix object not found" 
    raise ValueError(err)

def xsp2co1_fast_order(wx, via_word = True):
    assert isinstance(wx, Xsp2_Co1)
    if not via_word:
        return chk_qstate12(xsp2co1_order_elem(wx._data))
    m = MM0(wx)
    return xsp2co1_order_word(m._data, m.length)

def xsp2co1_fast_half_order(wx):
    assert isinstance(wx, Xsp2_Co1)
    buf = np.zeros(10, dtype = np.uint32)
    m = MM0(wx)
    res = chk_qstate12(xsp2co1_half_order_word(m._data, m.length, buf))
    o, l = divmod(res, 256)
    assert 0 <= l <= 10
    out = Xsp2_Co1()
    chk_qstate12(xsp2co1_mul_elem_word(out._data, buf, l))
    return o, out


@pytest.mark.xsp2co1
def test_elem_order(ntests = 50, verbose = 0):
    neutral = Xsp2_Co1()
    print("Test function computation of order in G_x0")
    for i, w in enumerate(make_testwords(monomial=False, ntests=ntests)):
        wx = Xsp2_Co1(w)
        o_ref = xsp2co1_ref_order(wx)
        o = xsp2co1_fast_order(wx, via_word = i & 1)
        ok = o == o_ref
        if verbose or not ok:
            print("Test", i+1)
            print("g = ", MM(*w))
            print("order =", o)
            if not ok:
                print("expected:", o_ref)
                raise ValueError("Computation or order has failed")

        o1, invol = xsp2co1_fast_half_order(wx)
        assert o1 == o
        if (o & 1):
            assert invol == neutral
        else:
            assert xsp2co1_fast_power(wx, o >> 1) == invol
            assert invol * invol == neutral

#######################################################################
# benchmark
#######################################################################

def one_benchmark(ntests):
    tags = "dxyplplplplplp"
    samples = [Xsp2_Co1(*rand_G_x0_testword(tags)) for i in range(400)]
    a = np.zeros(10, dtype = np.uint32)
    t_start = time.process_time()
    for elem in samples:
        elem_data = elem._data
        for i in range(ntests):
            assert xsp2co1_elem_to_word(elem_data, a) >= 0
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
    
    