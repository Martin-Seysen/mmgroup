"""Test C function dealing with Leech lattice vectors type 2

"""

import pytest
import numpy as np

from mmgroup import Cocode, XLeech2, PLoop, Xsp2_Co1

from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_n_type_22
from mmgroup.generators import gen_leech2_find_v4_2xx
from mmgroup.generators import gen_leech2_reduce_2xx
from mmgroup.generators import gen_leech2_map_2xx
from mmgroup.generators import gen_leech2_u4_2xx
from mmgroup.generators import rand_get_seed

OMEGA = XLeech2(0x800000)
STD_V2 = XLeech2(Cocode([2,3])).ord & 0xffffff
STD_V233 = (XLeech2(Cocode([2])) * OMEGA).ord & 0xffffff
STD_V222 = XLeech2(Cocode([1, 2])).ord & 0xffffff



def G_x0_testdata(ntests):
    yield Xsp2_Co1(1)
    for i in range(ntests):
        yield Xsp2_Co1('r', 'G_x0')

@pytest.mark.gen_xi
def test_type22(verbose = 0):
    r"""Test generation of some type-2 vectors

    """
    data = set()
    for i in range(4600):
        w = gen_leech2_n_type_22(i)
        data.add(w)
        assert gen_leech2_type(w) == 2, (i, hex(w), gen_leech2_type(w))
        w2 = w ^ STD_V2
        assert gen_leech2_type(w2) == 2, (i, hex(w2), gen_leech2_type(w2))
    assert len(data) == 4600
    
    


@pytest.mark.gen_xi
def test_reduce_233(ntests = 1000, verbose = 0):
    tf_c = np.zeros(8, dtype = np.uint32)
    v100 = np.zeros(100, dtype = np.uint32)
    v100a = np.zeros(100, dtype = np.uint32)

    seed = rand_get_seed()
    for i, tf in enumerate(G_x0_testdata(ntests)):
        v2 = (XLeech2(STD_V2) * tf).ord & 0xffffff 
        v3 = (XLeech2(STD_V233) * tf).ord & 0xffffff
        v4 = gen_leech2_find_v4_2xx(v2, v3, seed)
        assert v4 >= 0
        vtype, v4 = divmod(v4, 0x1000000)
        assert vtype == 3
        assert gen_leech2_type(v4) == 4
        len_ = gen_leech2_reduce_2xx(v2, v3, v4, tf_c)
        assert len_ >= 0, hex(len_)
        tf = tf_c[:len_ & 0xff] 
        v2_tf = gen_leech2_op_word(v2, tf, len(tf)) & 0xffffff    
        v3_tf = gen_leech2_op_word(v3, tf, len(tf)) & 0xffffff
        if verbose:
            print("Test %d: v2 = %s, v3 = %s" % (i+1, hex(v2), hex(v3)))
        assert v2_tf == STD_V2   
        assert v3_tf == STD_V233
        if i < 50:
            gen_leech2_map_2xx(tf, len(tf), 3, v100)
            for j, v in enumerate(v100):
                assert gen_leech2_type(v) == 4
                assert gen_leech2_type(v ^ v2) == 2
                assert gen_leech2_type(v ^ v3) == 2
                assert gen_leech2_type(v ^ v2 ^ v3) == 2
            ret = gen_leech2_u4_2xx(v2, v3, seed, v100a)
            assert ret == 100
            assert  set(v100) == set(v100a)




@pytest.mark.gen_xi
def test_reduce_222(ntests = 1000, verbose = 0):
    tf_c = np.zeros(8, dtype = np.uint32)
    v891 = np.zeros(891, dtype = np.uint32)
    v891a = np.zeros(891, dtype = np.uint32)
    seed = rand_get_seed()
    for i, tf in enumerate(G_x0_testdata(ntests)):
        v2 = (XLeech2(STD_V2) * tf).ord & 0xffffff
        v3 = (XLeech2(STD_V222) * tf).ord & 0xffffff
        v4 = gen_leech2_find_v4_2xx(v2, v3, seed)
        assert v4 >= 0
        vtype, v4 = divmod(v4, 0x1000000)
        assert vtype == 2
        assert gen_leech2_type(v4) == 4
        len_ = gen_leech2_reduce_2xx(v2, v3, v4, tf_c)
        assert len_ >= 0, hex(len_)
        tf = tf_c[:len_ & 0xff]
        v2_tf = gen_leech2_op_word(v2, tf, len(tf)) & 0xffffff
        v3_tf = gen_leech2_op_word(v3, tf, len(tf)) & 0xffffff
        if verbose:
            print("Test %d: v2 = %s, v3 = %s" % (i+1, hex(v2), hex(v3)))
        assert v2_tf == STD_V2, (hex(v2_tf), hex(v3_tf), hex(STD_V222))
        assert v3_tf == STD_V222, (hex(v2_tf), hex(v3_tf), hex(STD_V222))
        if i < 50:
            gen_leech2_map_2xx(tf, len(tf), 2, v891)
            for j, v in enumerate(v891):
                assert gen_leech2_type(v) == 4, (j, hex(v))
                assert gen_leech2_type(v ^ v2) == 2, j
                assert gen_leech2_type(v ^ v3) == 2, j
                assert gen_leech2_type(v ^ v2 ^ v3) == 2, j
            ret = gen_leech2_u4_2xx(v2, v3, seed, v891a)
            assert ret == 891
            assert  set(v891) == set(v891a)

