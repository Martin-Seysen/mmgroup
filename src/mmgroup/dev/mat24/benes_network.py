from __future__ import absolute_import, division, print_function
from __future__ import  unicode_literals


from random import shuffle, randint

from mmgroup.bitfunctions import lrange


def benes_step(n, sh, p):
    d = 1 << sh
    assert n == (n >> (sh + 1)) << (sh + 1)
    q = [None] * n
    for i in range(n):
        q[p[i]] = i                                               
    done = 0
    res0 = res1 = 0
    for i in range(n):
        j = i
        while not done & (1 << j):
            done |= 1 << j
            j = p[j]
            #next do: if (j & d): res1 |=  1 << (j & ~d)
            res1 |= ((j & d) >> sh) << (j & ~d)
            j = q[j ^ d]
            done |= 1 << j
            #next do: if (~j & d): res0 |=  1 << (j & ~d)
            res0 |= ((~j & d) >> sh) << (j & ~d)
            j = j ^ d    
    res2, res3 = res0 | (res0 << d), res1 | (res1 << d)
    for i in range(n):
        j = p[i ^ (((res2 >> i) & 1) << sh)] 
        q[i] = j ^ (((res3 >> j) & 1) << sh)
    return res0, res1, q


def benes_step_final(n, sh, p):
    d = 1 << sh
    res = 0
    assert n == (n >> (sh + 1)) << (sh + 1)
    for i in range(n):
        assert p[i] in [i, i^d]
        res |= ((i ^ p[i]) >> sh) << (i & ~d)
    return res
            


def benes_step_final3(n, sh, p):
    n3 = n//3
    assert n3 == 1 << sh
    res0 = res1 = res2 = 0
    for i in range(n3):
        t = p[i] >> sh
        t = 2 * t + ((p[i+n3] >> (sh + (t & 1))) & 1)
        t = (0x236407 >> (t << 2)) & 0xf
        res2 |=  (t & 1) << i
        res1 |=  ((t >> 1) & 1) << i
        res0 |=  ((t >> 2) & 1) << i
    return res0, res1, res2  



def check_benes_step(n, sh, p, res0, res1, q):
    d = 1 << sh
    def perm(n, res):
        p1 = [None] * n
        for i in range(n):
            if res & (1 << i):
                 p1[i], p1[i+d] = i+d, i
            elif i & d == 0:
                 p1[i], p1[i+d] = i, i+d
        return p1
    p0 = perm(n, res0)
    p1 = perm(n, res1)
    r = [None] * n
    for i in range(n):
        r[i] =  p1[p[p0[i]]]
    assert q == r, (q,r)
    for i in range(n):
        assert i & d == r[i] & d


def test_benes_step(n, sh, p):
    res0, res1, q = benes_step(n, sh, p)
    check_benes_step(n, sh, p, res0, res1, q)





def build_benes24(perm):
    """Compute modified Benes network for permutation of 24 entries

    The network consists of 9 layers. The returned array 'res' of 
    length 9 describes that network. In layer i, entry j is to be
    exchanged with entry  j + d[i], if bit j of the value res[i]
    is set. Here d[i] = 1,2,4,8,16,8,4,2,1 for i = 0,...,8.  
    """
    res = [None] * 9
    p = [None] * 24
    q = [None] * 24
    for i in range(24):
        p[i] = perm[i]       # copy permutation perm to p

    # The first and the last three layers are a standard Benes network. 
    # Do Benes netork looping algorithm steps for d = 1, 2, 4
    for sh in range(3):      
        d = 1 << sh
        for i in range(24):
            q[p[i]] = i      # q := inverse of p
        done = 0             # bit i marks that step i->p[i] is done
        res0 = 0             # initial looping transpositions
        res1 = 0             # final looping transpositions
        for i in range(24):  # Looping step for Benes network
            j = i            # j is a node not yet processed
            while not (done & (1 << j)): # while node j not done
                done |= 1 << j         # delare node j done
                j = p[j]               # j := permutation result p[j]
                # route node p[j] thru '0' part of inner Benes network
                # so we do: if (j & d): res1 |=  1 << (j & ~d)
                res1 |= ((j & d) >> sh) << (j & ~d)
                j = q[j ^ d]           # j = origin of buddy of p[j]
                done |= 1 << j         # declare that buddy done
                # route buddy thru '1' part of inner Benes network
                # so we do: if (~j & d): res0 |=  1 << (j & ~d)
                res0 |= ((~j & d) >> sh) << (j & ~d)
                j = j ^ d              # j = buddy of that origin
        res[sh] = res0       # save initial looping transposition
        res[8-sh] = res1     # save final looping transposition
        res0 |= res0 << d    # initial: exchange i with i^d if bit i set
        res1 |= res1 << d    # final: exchange i with i^d if bit i set
        for i in range(24):  # compute q = (initial) * p * (final)
            j = p[i ^ (((res0 >> i) & 1) << sh)] 
            q[i] = j ^ (((res1 >> j) & 1) << sh)
        for i in range(24):  # copy (initial) * p * (final) to p
            p[i] = q[i]
 
    # It remains to compute the 3 middle layers. They must compute
    # the permutation i -> p[i] with p[i] = i (mod 8). E.g. for i=0
    # we do the following transpositions, if (0, 8, 16) maps to:
    #
    #   ( 0,  8, 16):          (id)  *   (id)   *   (id)  // [1]
    #   ( 0, 16,  8):         (0,8)  *  (0,16)  *  (0,8)  // [0]     
    #   ( 8,  0, 16):         (0,8)  *   (id)   *   (id)  // [2]  
    #   ( 8, 16,  0):         (0,8)  *  (0,16)  *   (id)  // [3]
    #   (16,  0,  8):          (id)  *  (0,16)  *  (0,8)  // [4]  
    #   (16,  8,  0):          (id)  *  (0,16)  *   (id)  // [5]
    #          
    # For each permutation of (i, i+8, i+16) with compute a number j,
    # as indicated in square brackets above, from the bits 3 and 4 of
    # p[i] and p[i+8]. Then we use table look up for obtaining the
    # correct transpostions as in the list abvoe.   
    res0 = res1 = res2 = 0
    for i in range(8):
        j = p[i] >> 3
        j = 2 * j + ((p[i+8] >> (3 + (j & 1))) & 1)
        j = (0x236407 >> (j << 2)) & 0xf
        res2 |=  (j & 1) << i
        res1 |=  ((j >> 1) & 1) << i
        res0 |=  ((j >> 2) & 1) << i
    res[3], res[4], res[5] = res0, res1, res2
    return res










          