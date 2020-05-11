from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


from random import shuffle, randint
import pytest

from mmgroup.bitfunctions import lrange


from mmgroup.dev.mat24.benes_network import  benes_step
from mmgroup.dev.mat24.benes_network import  check_benes_step
from mmgroup.dev.mat24.benes_network import  build_benes24
from mmgroup.dev.mat24.benes_network import  benes_step_final




def do_benes_step(n, sh, p):
    res0, res1, q = benes_step(n, sh, p)
    check_benes_step(n, sh, p, res0, res1, q)



def exec_benes24(x, net):
    shifts = [1,2,4,8,16,8,4,2,1]
    for  (sh, msk) in zip(shifts, net):
        msk = (x ^ (x >> sh)) & msk
        x ^=  msk ^  (msk << sh) 
    return x
       



def benes24_testdata(ntests=1000):
    for i in range(10):
        p1 = [0,8,16]; shuffle(p1)
        p = [ [p1[j]] + lrange(8*j+1,8*j+8) for j in range(3) ]
        x = randint(0, 2**24-1)
        yield x, sum(p, [])

    for i in range(10):
        x = randint(0, 2**24-1)
        p0 = lrange(0,8); shuffle(p0)
        p1 = lrange(8,16); shuffle(p1)
        p2 = lrange(16,24); shuffle(p2)
        yield x , p0 + p1 + p2

    for i in range(ntests):
        x = randint(0, 2**24-1)
        p0 = lrange(0,24); shuffle(p0)
        yield x, p0
         




@pytest.mark.mat24
def test_benes24_network(ntests = 1000):
    print("24-bit permutation network test ...", end="")
    for x, p in benes24_testdata(ntests):
        net =  build_benes24(p[:])
        y = exec_benes24(x, net)
        for i in range(24):
            assert (x >> i) & 1 == (y >> p[i]) & 1, (i, hex(x), hex(y), p)  
    print(" passed")


@pytest.mark.mat24
def test_benes_steps():
    do_benes_step(4,0, [1,2,0,3])
    do_benes_step(4,1, [1,2,0,3])

    for t in range(10):
        for n in (4,8,16,32):
            for sh in range(5):
                if 1 << sh < n:
                    l = lrange(n)
                    shuffle(l)
                    do_benes_step(n, sh, l)

    for t in range(10):
        l = lrange(32)
        shuffle(l)
        for i in range(4):
            _1, _2, l = benes_step(32, i, l)
        benes_step_final(32, 4, l)
            
            
             