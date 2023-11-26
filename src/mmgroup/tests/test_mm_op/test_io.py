from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


from numbers import Integral
import numpy as np
from random import randint


import pytest

from mmgroup.mm_space import MMSpace, MMVector
from mmgroup.mm_space import characteristics
from mmgroup.mm_op import mm_aux_index_sparse_to_extern
from mmgroup.mm_op import mm_aux_index_sparse_to_intern
from mmgroup.mm_op import mm_aux_index_extern_to_sparse
from mmgroup.mm_op import mm_aux_index_intern_to_sparse
from mmgroup.mm_op import mm_aux_index_extern_to_intern
from mmgroup.mm_op import mm_aux_index_check_intern
from mmgroup.mm_op import mm_aux_mmv_extract_sparse
from mmgroup.mm_op import mm_aux_mmv_size

from mmgroup.tests.spaces.spaces import MMTestSpace
from mmgroup.tests.spaces.sparse_mm_space import SparseMmSpace
from mmgroup.tests.spaces.sparse_mm_space import SparseMmV
from mmgroup.tests.spaces.sparse_mm_space import SparseMmVector

from mmgroup.tests.chisquare import chisquare




#########################################################################
### Chi square test
#########################################################################

def chisquare_crosscor(p, b1, b2 = None, d = 1, m = 0):
    """chi-square test over cross correlation
    
    Let b1 and b2 be arrays of intgers (mod p) of the same length.
    the we put a[i] = b1[i] + b2[i + d] (mod p), with idex i + d 
    wrapping around. b2 may be 0, indicating an array of zeros.

    We return the chisquare test probability assumint that b1 and b2
    are equidistributed.

    In case m > 0 the entries a[i] are reduced modulo m. 
    """
    d %= len(b1)
    if b2 is None:
        a =  np.zeros(len(b1), dtype = np.uint32)
    else:
        assert len(b1) == len(b2)
        a =  np.array(np.concatenate((b2[d:], b2[:d])), dtype = np.uint32)
    a = (a + b1) % p
    f_exp = None
    if 0 < m < p:
        a %= m
        f = np.fromiter((x % m for x in range(p)), dtype = np.uint32)
        f_exp = np.bincount(f, minlength=m)
        f_exp = f_exp * len(a) / p
        p = m
    bins = np.bincount(a, minlength=p)
    chisq, prob = chisquare(bins, f_exp)
    return prob

def chisquare_crosscor_ok(p, b1, b2 = 0, d = 1, m = 0):
    for i in range(4):
        prob = chisquare_crosscor(p, b1, b2, d, m)
        if 0.01 < prob < 0.99: return True, prob
    return False, prob


#########################################################################
### Test random vectors and manipulation of them
#########################################################################

def do_test_sparse_rep(v):
    p = v.p
    space = MMVector
    sparse_space = SparseMmVector
    vsp = sparse_space(p, v)
    assert (v == space(p, vsp))
    assert (v == space(p, "S", v.as_sparse()))
    v2 = space(p, 'V', v.as_bytes())
    assert (v == v2)


def check_twin_ABC(a1, a2):
    assert 0 <= a1 <= 72*32, hex(a1)
    t1, a10 = divmod(a1, 24*32)
    t2, a20 = divmod(a2, 24*32)
    a1i, a1j = divmod(a10, 32)
    a2i, a2j = divmod(a20, 32)
    assert a1j < 24
    if a1i == a1j:
        assert t1 == 0, (hex(a1), hex(a2), hex(a1i), hex(a1j),  hex(t1))
        assert a2 == 0, (hex(a1), hex(a2), hex(a1i), hex(a1j), hex(a2))
    else:
        assert t1 == t2, (hex(a1), hex(a2), hex(a1i), hex(a1j), hex(t1), hex(t2))
        assert a1i == a2j, (hex(a1), hex(a2), hex(a1i), hex(a1j))
        assert a1j == a2i, (hex(a1), hex(a2), hex(a1i), hex(a1j))


def do_test_rep_conversion(v):
    p = v.p
    space = MMVector
    sparse_space = SparseMmVector
    ranges = [
        (0, 5, 1),
        (0, 50, 10),
        (randint(0,23), 852, randint(13, 17)),
        (randint(852,852+200), 196884, randint(390, 409)),
    ]
    for  start, stop, step in ranges:
        data_sparse =  v.as_bytes()[start:stop:step] 
        data = v[start:stop:step]
        assert len(data_sparse) == len(data), (len(data_sparse), len(data))
        assert (data_sparse == data).all(), (data_sparse, data)
        a = np.zeros(len(data), dtype = np.uint32)
        for i, index in enumerate(range(start,stop,step)):
            a[i] =  mm_aux_index_extern_to_sparse(index)             
            i_ext =  mm_aux_index_sparse_to_extern(a[i]) 
            assert i_ext == index, (hex(index), hex(a[i]), hex(i_ext))
            i_int1 =  mm_aux_index_extern_to_intern(index) 
            i_int2 =  mm_aux_index_sparse_to_intern(a[i]) 
            assert i_int1 == i_int1, (hex(index), hex(i_int1), hex(i_int2))
            i_sp1 =  mm_aux_index_intern_to_sparse(i_int1)
            assert i_sp1 == a[i], (hex(index), hex(i_sp1), hex(a[i]))
            twin = mm_aux_index_check_intern(i_int1)
            if i_int1 >= 72*32:
                assert twin == 0
            else:
                check_twin_ABC(i_int1, twin)
        mm_aux_mmv_extract_sparse(p, v.data, a, len(a))
        assert (a & 0xff == data).all(), (a & 0xff, data)



def check_direct_symmetric(v):
    for tag in "ABC":
        m = v[tag]
        assert (m == m.T).all()
        if tag in "BC":
            for i in range(24):
                assert m[i,i] == 0 
        assert v[tag, 3, 3] == m[3, 3]
        for k in range(4):
            i0, i1 = randint(0, 23), randint(0, 23)
            assert  v[tag, i0, i1] == m[i0, i1]
   
def do_test_random_io(p, verbose = False, to_sparse = 1):
    name_slices = [
        ("D", "D", 1),
        ("T", "T", 1),
        ("As",  slice(24,24+276), 1),
        ("Bs",  slice(24+276,24+2*276), 1),
        ("Cs",  slice(24+2*276,24+3*276), 1),
        ("ABCs",slice(0,24+3*276), 1), 
        ("X", "X", 1),
        ("Y", "Y", 1),
        ("Z", "Z", 1),
        ("A", "A", 0),
        ("B", "B", 0),
        ("C", "C", 0),
    ]
    chisqu_ok = {}
    for (name, _, do_chisqu_test) in name_slices: 
        chisqu_ok[name] = not do_chisqu_test
    for i in range(5):
        space = MMTestSpace(p)
        v = space(0)
        assert np.count_nonzero(v.data[:mm_aux_mmv_size(v.p)]) == 0
        v.check()
        v = space('R')
        v.check()
        check_direct_symmetric(v)
        if to_sparse and i == 0:
            do_test_sparse_rep(v)
        do_test_rep_conversion(v)
        
        for name, slice_, do_chisqu_test in name_slices:
            v1 = v[slice_]
            if name  in "DABCTXYZ":
                xx = v[name]
                assert (v1 % p == xx % p).all(), (p, name, slice_)
            if do_chisqu_test:
                v1 = v1.reshape(-1)
                prob1 =  chisquare_crosscor(p, v1,  m = 23)
                prob2 =  chisquare_crosscor(p, v1, v1, 7,  m = 11)
                if verbose:
                    print("p=%3d, %-8s: %.8f, %.8f" % (p,name,prob1,prob2))
            pr = [prob1, prob2]
            if 0.01 < min(pr) <= max(pr) < 0.99: 
                chisqu_ok[name] = 1 
    errors = ", ".join([ x for x in chisqu_ok if not chisqu_ok[x] ])
    if len(errors):
        print("\nChisqu test failed for p=%d and %s !!\n" % (p,errors))
        raise ValueError("Chisquare test failed for p=%d " % p)

@pytest.mark.mm_op
def test_all_random_io(verbose = 0):
    print("Chisquare test for randomization of an mm vector")
    for i in range(2):
        to_sparse = 0 if i else 1
        for p in characteristics():
            do_test_random_io(p, verbose, to_sparse)
    print("Test passed")



#########################################################################
### Test manipulation of symmetric tags A, B, C
#########################################################################


def do_test_sym_io(p, tag, verbose = 1):
    def array_data(shape,  data):
        if len(shape) == 0:
            return data 
        b1 = np.fromiter(data, dtype = np.int32).reshape(shape) 
        return b1         
    a_data = [
        [(slice(0,3), slice(0,3)), (3,3), 
                  (i*j+1 for i in range(3) for j in range(3))],
        [(slice(2,5), 17), (3,), 
                  (i for i in range(3))],
        [(slice(2,4), slice(2,4)), (2,2), 
                   (i*j+1 for i in range(2) for j in range(2))],
        [(slice(2,24), slice(2,24)), (22,22), 
                 (i*j+1 for i in range(22) for j in range(22))],
    ]               
    if verbose:
        print("p = %d, tag = %s" % (p,tag))
    space = MMTestSpace(p)
    b = np.zeros((24,24), dtype=np.int32)
    v = space(0)
    for slices, shape, data in a_data:
        b1 = array_data(shape, data)
        v[tag, slices[0], slices[1]] = b1
        b[slices[0], slices[1]] = b1 % p
        b1T = b1.T if len(shape) == 2 else b1 
        b[slices[1], slices[0]] = b1T % p
        
    if tag in "BC":
        for i in range(24): b[i,i] = 0
    vd = v[tag]
    assert (vd == vd.T).all()
    assert (b % p == vd).all(), (tag, b - vd, b, vd)

    for (sl0, sl1), shape, _ in a_data:
        eq = v[tag, sl0, sl1] == vd[sl0, sl1]
        if len(shape): eq = eq.all()
        assert eq

    abc = v[:852]
    v1 = space(0)
    v1[:852] = abc
    assert v == v1, (v, v1)


def do_test_sym_io_bad(p, tag, verbose = 1):
    def array_data(shape, value):
        if len(shape) == 0:
            return value
        return np.zeros(shape, dtype = np.int32) + value
    a_data = [
        [(slice(1,23), slice(2,24)), (22,22), 1],
        [(slice(0,11), slice(9,19)), (11,10), 2],
        [(1,23), (), 4],
        [(slice(7,11), 23), (4,),  8],
    ]
    if verbose:
         print("p = %d, tag = %s, bad data" % (p,tag))
    space = MMTestSpace(p)
    v = space(0)
    w = np.zeros( (24, 24), dtype = np.int32)
    for slices, shape, value in a_data:
        new_data = array_data(shape, value)
        v[tag, slices[0], slices[1]] = new_data
        w[slices[0], slices[1]] = new_data
        v.check()
        try:        
            new_data_t = new_data.T
        except:
            assert isinstance(new_data, Integral) or len(new_data.shape) < 2
            new_data_t = new_data
        w[slices[1], slices[0]] = new_data_t
        w = w % p
        if tag in "BC":
            for i in range(24): w[i,i] = 0
        assert (v[tag] == w).all(), (tag, v[tag] - w, v[tag], w)
           



@pytest.mark.mm_op
def test_sym_all(verbose = 0):
    print("Test io operations for tags ABC")
    for tag in "ABC":
        for p in characteristics():
            do_test_sym_io(p, tag, verbose = verbose)
            do_test_sym_io_bad(p, tag, verbose = verbose)
    print("Test passed")


#########################################################################
### Test manipulation of tags T, X, Y, Z
#########################################################################


def do_test_large_io(p, tag, verbose = 0):
    def array_data(shape,  data):
        if len(shape) == 0:
            return data 
        b1 = np.fromiter(data, dtype = np.int32).reshape(shape) 
        return b1         
    t_d = [
        (slice(0, 3, 1), slice(0, 5, 1)),
        (slice(18, 759, 20), slice(2,64, 3)),
        (slice(19,622, 70), 53),
        (52, slice(0,64,9)),
        (313,55), 
    ]
    x_d = [
        (slice(18, 2048, 79), slice(2,24, 7)),
        (slice(19, 1622, 70), 23),
        (1552, slice(0,24,9)),
        (1313,15), 
    ]
    lengths, a_data = ((759,64), t_d) if tag =='T' else ((2048,24), x_d)
    v = MMTestSpace(p)('R')
    for slices in a_data:
        prod, shape = 1, ()
        for l, s in zip(lengths, slices):
            if isinstance(s, slice):
                size_ = len(range(*s.indices(l)))
                prod *= size_
                shape = shape + (size_,)
        data = np.array(p * np.random.rand(prod), dtype = np.uint32)
        if len(shape) == 0:
            data = data[0]
        else:
            data = data.reshape(shape)
        v[tag, slices[0], slices[1]] = data
        eq = v[tag, slices[0], slices[1]] == data
        if len(shape): eq = eq.all()
        assert eq, (tag, v[tag, slices[0], slices[1]], data)
    v.check()

@pytest.mark.mm_op
def test_large_all(verbose = 0):
    print("Test io operations for tags TXYZ")
    for tag in "TXYZ":
        for p in characteristics():
            do_test_large_io(p, tag, verbose = verbose)
    print("Test passed")




