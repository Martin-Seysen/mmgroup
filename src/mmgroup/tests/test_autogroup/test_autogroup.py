from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals



import sys
import warnings
from random import sample, randint
import re


from copy import deepcopy
from collections.abc import Sequence

import pytest



from mmgroup.tests.groups.free_group import FreeGroup











 
####################################################################
####################################################################
### Tests
####################################################################
####################################################################
 

@pytest.mark.auto_group
def test_autogroup():
    print("")
    g_bad = FreeGroup("abc", rewrite={"cb":"bc"} )       
    print(1, g_bad("cb"*30)) 
     
    print(2,  g_bad("bac") )

    g1 = FreeGroup("abc", commute = "ab", bad_commute="ab")
    w1 = g1("ab")
    EXP = 20
    print(3, w1**EXP)
    assert w1**EXP ==  g1("a")**EXP *  g1("b")**EXP 
    g2 = FreeGroup("abc", commute="ab")
    w2 = g2("ab")
    print(4, w2**30)
    assert w2**30 ==  g2("a")**30 *  g2("b")**30 
    w_long = w2**100
    assert w_long == 1 * w_long == w_long * 1
    assert w2 ** (-1)  ==  1 / w2
    assert w_long * w_long**(-1) == g2.neutral()

    g = FreeGroup("abc")
    w1 = g("A")
    w2 = g("a")
    print (5, w1, w2, w1*w2)
    assert w1 * w2 == g.neutral()

    g = FreeGroup("abc")
    w1 = g("AbC")
    w2 = g("cBBaac")
    print (6, w1, w2, w1*w2, w1**(-1)*w2, w1**(w2)) 
    assert  w2**(-1)*w1*w2 == w1**(w2)
    assert  w1*w2 == g("ABaac")
    






       