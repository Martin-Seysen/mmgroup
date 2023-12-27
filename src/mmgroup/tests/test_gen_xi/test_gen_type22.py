"""Test C function dealing with Leech lattice vectors type 2

"""

import pytest
import numpy as np

from mmgroup import Cocode, XLeech2, PLoop, Xsp2_Co1

from mmgroup.generators import gen_leech2_type
from mmgroup.generators import gen_leech2_op_word
from mmgroup.generators import gen_leech2_n_type_22
from mmgroup.generators import gen_leech2_find_v4_233
from mmgroup.generators import gen_leech2_reduce_233
from mmgroup.generators import rand_get_seed

OMEGA = XLeech2(0x800000)
STD_V2 = XLeech2(Cocode([2,3])).ord & 0xffffff
STD_V3 = (XLeech2(Cocode([2])) * OMEGA).ord & 0xffffff


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
def test_reduce_233(verbose = 0):
    tf_c = np.zeros(8, dtype = np.uint32)
    seed = rand_get_seed()
    for i in range(1000):
        tf = Xsp2_Co1('r', 'G_x0')
        v2 = (XLeech2(STD_V2) * tf).ord & 0xffffff 
        v3 = (XLeech2(STD_V3) * tf).ord & 0xffffff 
        v4 = gen_leech2_find_v4_233(v2, v3, 5000, seed)
        assert v4 >= 0
        trials, v4 = divmod(v4, 0x1000000)
        if verbose:
            print("%d trials" % trials)
        assert gen_leech2_type(v4) == 4
        len_ = gen_leech2_reduce_233(v2, v3, v4, tf_c)
        assert len_ >= 0, hex(len_)
        tf = tf_c[:len_ & 0xff] 
        v2_tf = gen_leech2_op_word(v2, tf, len(tf)) & 0xffffff    
        v3_tf = gen_leech2_op_word(v3, tf, len(tf)) & 0xffffff
        if verbose:
            print("Testing v2 = %s, v3 = %s" % (hex(v2), hex(v3)))
        assert v2_tf == STD_V2   
        assert v3_tf == STD_V3 
        
    pass  

