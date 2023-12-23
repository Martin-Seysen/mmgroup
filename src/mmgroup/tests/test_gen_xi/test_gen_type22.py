"""Test C function dealing with Leech lattice vectors type 2

"""

import pytest

from mmgroup.generators import gen_leech2_n_type_22
from mmgroup.generators import gen_leech2_type

STD_V2 = 0x200


@pytest.mark.gen_xi
def test_type22(verbose = 0):
    r"""Test genration of some type-2 vectors

    """
    data = set()
    for i in range(4600):
        w = gen_leech2_n_type_22(i)
        data.add(w)
        assert gen_leech2_type(w) == 2, (i, hex(w), gen_leech2_type(w))
        w2 = w ^ STD_V2
        assert gen_leech2_type(w2) == 2, (i, hex(w2), gen_leech2_type(w2))
    assert len(data) == 4600
    
    
