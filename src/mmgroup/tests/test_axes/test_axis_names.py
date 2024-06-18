import sys
from random import randint
from functools import partial
import numpy as np
from collections import defaultdict

import datetime
import time
import pytest

if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import MM0, MMSpace, MMV, characteristics, MMVectorCRT
from mmgroup import  XLeech2, PLoop, Cocode, Xsp2_Co1



def XLeech2_type2_elements():
    yield XLeech2(Cocode([1,2])) * MM0('r','N_x0')
    yield XLeech2(PLoop(list(range(8)))) * MM0('r','N_x0') 
    yield XLeech2(0, Cocode([1])) * MM0('r','N_x0')
    
def XLeech2_type2_elements_conjugation():
    for x in XLeech2_type2_elements(): 
        iclass, g = MM0(x).conjugate_involution_G_x0(x)
        assert iclass == '2A_x0'
        yield g**-1, x




@pytest.mark.axes
def test_axis_names(verbose =  0):
    g_0 = MM0('q', 'v+')
    assert g_0 == MM0('d', 'v+') == MM0('d', [2,3])
    assert g_0 * MM0('q', 'v-') == MM0('q', '-') == MM0('x', '-') == MM0('x', 0x1000)
    for p in characteristics() + [0]:
        V = MMV(p) if p else partial(MMVectorCRT, 15)
        std_axis = V('Axis', 'v+')
        if p == 0 and verbose:
             print(std_axis['A'])
             print(V("I", 3, 2)['A'])
        assert std_axis ==  V("I", 3, 2)
        opp_axis =  V('Axis', 'v-')
        assert opp_axis ==  V("I", 3, 2) * MM0('x', 0x200)
        for g, x in  XLeech2_type2_elements_conjugation():
            axis =  std_axis * g
            if verbose:
                print("g=", g)
                print("x=", x)
                if verbose > 1:
                    print("axis=", axis)
            assert axis ==  V('Axis', x)
            assert axis ==  V('Axis', x.ord)
            assert axis ==  V('Axis', g_0**g)
            

