from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


import pytest

from random import randint, shuffle, sample

from mmgroup import mat24
from mmgroup import GCode, Cocode, PLoop, PLoopZ
from mmgroup import Octad, SubOctad, GcVector 
from mmgroup import AutPL
from mmgroup import Parity
from mmgroup.mat24 import MAT24_ORDER

from mmgroup.mm_group import MMGroup
group = MMGroup()

#####################################################################
# Test embedding of Cocode and AutPL into MMGroup
#####################################################################

@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [100])
def test_autploop(n_cases):
    for i in range(n_cases):
        autp = AutPL(('d', 'r'), ('p', 'r'))
        g = group(autp)
        assert g ==  group(('d', autp.cocode), ('p',  autp.perm))
        assert g ==  group(('p',  autp))
        assert g ==  group.atom('p',  autp)
        coc, p_num = autp.cocode, autp.perm_num
        if coc and p_num:
            assert g.as_tuples() == autp.as_tuples() 
            assert g.as_tuples() == [('d', coc), ('p', p_num)] 
        Coc = Cocode(coc)
        h1 = range(9, 18)
        h2 = autp.perm[9:18]
        assert g ==  group(('d', Coc), ('p',  zip(h1, h2)))    
        assert g ==  group(('d', Coc), ('p',  dict(zip(h1, h2))))    


#####################################################################
# Test embedding of PLoop into MMGroup
#####################################################################

@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [100])
def test_xyz(n_cases):
    g1 = group()
    for i in range(n_cases):
        v = randint(0, 0x1fff)
        for tag in "xyz":
            g = group((tag, v)) 
            assert g ==  group((tag, PLoop(v))) 
            assert g ==  group.atom(tag, PLoop(v)) 
    group(('x', v))._mul(group(('y', v)), group(('z', v))) == g1
 

#####################################################################
# Test multiplication/inversion in groups AutPL and MMGroup
#####################################################################


@pytest.mark.ploop
@pytest.mark.parametrize("n_cases", [100])
def test_group_op(n_cases):
    for i in range(n_cases):
        a1 = AutPL(('d', 'r'), ('p', 'r'))
        c, p = randint(0, 0xfff), randint(0, MAT24_ORDER - 1)
        a2 = AutPL(('d', c), ('p', p))
        assert group(a1)._mul(group(a2)) == group(a1 * a2)
        assert group(a1**-1) == group(a1)**-1


