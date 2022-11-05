"""Generate random elements of subgroups of the Mathieu group ``M_24``


This is yet a fragment.

It implements the generation of random elements ofvarious subgroups
of the Mathieu group :math:`M_{24}`.

For mathematical background see **The mmgroup guide for develpers**,
section **Subgroups of the Mathieu group :math:`M_{24}`**.

"""

import sys
from random import randint
from collections import defaultdict
from collections.abc import Sequence
from numbers import Integral
import numpy as np
import pytest


if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import mat24
from mmgroup.mat24 import MAT24_ORDER 



"""
try:
    # Try importing the fast C function
    from mmgroup import mat24
    from mmgroup.mat24 import MAT24_ORDER 
except (ImportError, ModuleNotFoundError):
    # Import a pure python substitute for the mmgroup.mat24 extension
    # if the original extension has not been found. For background, 
    # see section 'Mock up C extension modules' in file
    # docs/source/conf.py.
    from mmgroup.dev.mat24.mat24_ref import  Mat24
    try:
        mat24 = Mat24    
        MAT24_ORDER = Mat24.MAT24_ORDER
        del Mat24
         
"""





RAND_2 = 1
RAND_o = 2
RAND_t = 4
RAND_s = 8
RAND_l = 16
RAND_3 = 32


RAND_USED = RAND_2 | RAND_o | RAND_t | RAND_s | RAND_l | RAND_3



def free_bit24(mask, r):
    # mask = bitmask, only a cleared bit may be returned
    # r = random number
    mask &= 0xffffff
    w = 24 - mat24.bw24(mask)
    r, k = divmod(r, w)
    bmask = mat24.spread_b24(1 << k, 0xffffff & ~mask)
    return mat24.lsbit24(bmask), mask | bmask, r


def remask(h2, new_mask):
    mask = 0xffffff & ~ new_mask
    for b in h2:
        mask |= 1 << b
    return mask

   
H1 = [3, 2, 1, 0, 5, 4, 8] 
 
def perm_from_suboctad(h2, mask, r): 
    # h2 = list of length <= 6 encoding a subset of an octad
    # mask = bitmask, only cleared bits may be used to complete an octad
    # r = random number
    mask_o = 0
    for b in h2:
        mask_o |= 1 << b
    mask = (mask | mask_o) & 0xffffff
    while len(h2) < 5:
        b, mask, r = free_bit24(mask, r)
        mask_o |= 1 << b
        h2.append(b)
    syn = mat24.syndrome(mask_o, 24)
    if len(h2) < 6:
        b, _, r = free_bit24(~syn, r)
        h2.append(b)
    b, _, r = free_bit24(mask_o | syn, r)
    h2.append(b)
    return mat24.perm_from_heptads(H1, h2)





def subgroup(mode, mode_subgroup, mode_of):
    if mode &  mode_subgroup == mode_subgroup:
        mode |= mode_of
    return mode


def complete_mat24_rand_mode(mode):
    mode = subgroup(mode, RAND_l,  RAND_o)
    mode = subgroup(mode, RAND_l | RAND_2,  RAND_t)
    mode = subgroup(mode, RAND_t | RAND_2,  RAND_l | RAND_o)
    mode = subgroup(mode, RAND_3 | RAND_t,  RAND_o | RAND_s)
    mode = subgroup(mode, RAND_3 | RAND_l,  RAND_o | RAND_s)
    mode = subgroup(mode, RAND_3 | RAND_t | RAND_l,  RAND_o | RAND_2)
    return mode

def complete_mat24_rand_mode_minimal(mode):
  for i in range(3):
    mode = subgroup(mode, RAND_l,  RAND_o)
    mode = subgroup(mode, RAND_l | RAND_2,  RAND_t)
    mode = subgroup(mode, RAND_t | RAND_2,  RAND_o | RAND_l)
    mode = subgroup(mode, RAND_3 | RAND_l,  RAND_s)
    mode = subgroup(mode, RAND_3 | RAND_t,  RAND_o | RAND_s)
    mode = subgroup(mode, RAND_3 | RAND_t | RAND_l, RAND_2)
    return mode


SUBGROUPS = [i for i in range(64) if complete_mat24_rand_mode(i) == i]
 

def complete_aff_trio(h1, h2, h3):
    return h1 ^ h2 ^ h3


AL = [0,1,2,3,4,5,7,6]

def complete_aff_line(h1, h2, h3):
    v = AL[h1] ^ AL[h2] ^ AL[h3]
    return AL[v]


def complete_perm(mode, r, verbose = 0):
    if mode & RAND_o:
        mask8 = 0xff
    elif mode & RAND_t:
        r, k = divmod(r, 3)
        mask8 = 0xff << (8 * k)
    else:
        mask8 = 0xffffff
    mask2 = mask3 =  mask4 = mask8
    if mode & (RAND_2 | RAND_3):
        if mode & RAND_s:
            mask4 = 0xf
        mask3 = 0xe if mode & RAND_3 else mask4
        mask2 = 0xc if mode & RAND_2 else mask3
    elif mode & RAND_s:
        mask = ~(mask8 & 0x111111) & 0xffffff
        b, _, r = free_bit24(mask, r)
        #print("mmm", hex(mask), b)
        mask2 = mask3 = mask4 = 15 << b
    if verbose > 2:
        print("complete_perm, mode = 0x%02x, r = %d" % (mode, r))
        print("  mask2-8: 0x%x 0x%x 0x%x 0x%x" % 
            (mask2, mask3, mask4, mask8)) 
    mask = ~mask2 

    b, mask, r = free_bit24(mask, r)
    h2 = [b]
    if (mode & (RAND_t | RAND_l) == (RAND_t | RAND_l)):
        h2 += [b ^ 1]
        mask = remask(h2, mask3)
        b, mask, r = free_bit24(mask, r) 
        h2 += [b, b ^ 1]
        mask = remask(h2, mask8)
        b, mask, r = free_bit24(mask, r)
        h2 += [b, b ^ 1]
    else:
        b, mask, r = free_bit24(mask, r)
        h2 += [b]
        mask = remask(h2, mask3)
        b, mask, r = free_bit24(mask, r)
        h2 += [b]
        mask = remask(h2, mask4)
        if (mode & (RAND_t | RAND_l)):
            f = complete_aff_trio if mode & RAND_t else complete_aff_line
            b = f(h2[0], h2[1], h2[2])
            h2 += [b]
            mask = remask(h2, mask8)
            b, mask, r = free_bit24(mask, r)
            h2 += [b]
            b = f(h2[0], h2[1], h2[4])
            h2 += [b]
        else:
            b, mask, r = free_bit24(mask, r)
            h2 += [b]
    if verbose > 2:
        print("  h2 =", h2)
    return perm_from_suboctad(h2, ~mask8, r)



                 
        
                
        
def mat24_perm_rand_local(mode, r, verbose = 0):
    mode = complete_mat24_rand_mode(mode)
    return complete_perm(mode, r, verbose)



def check_set(pi, start, diff):
    s = 0
    for i in range(start, 24, diff):
        for j in range(1, diff):
            s |= pi[i] ^ pi[i + j]
    return (s & -diff) == 0


def mat24_perm_in_local(pi):
    mode = 0
    s = (1 << pi[2]) | (1 << pi[3]) 
    if (s == 0xc):
        mode |= RAND_2
    s |= (1 << pi[1])
    if (s == 0xe):
        mode |= RAND_3
    s |= (1 << pi[0]) | (1 << pi[4])  | (1 << pi[5])
    s |= (1 << pi[6]) | (1 << pi[7]) 
    if (s == 0xff):
        mode |= RAND_o
    if check_set(pi, 8, 2):
        mode |= RAND_l
    if check_set(pi, 0, 4):
        mode |= RAND_s
    if check_set(pi, 0, 8):
        mode |= RAND_t
    return mode
    


class UIntUnionFind(Sequence):
    """Union-find algorithm for nonnegative integers less than n

    Here the given set is the set of singletons ``set(i)`` for
    integers ``0 <= i < n``, with ``n`` given by the constructor.

    ``self.union(i1, i2)`` unites the set containing elemnent ``i1``
    with the element containing ``i2``. 

    ``self[i]`` returns a representative of the set containg ``i``.
    
    """
    def __init__(self, n):
        """Initialize set with integers ``i``, with ``0 <= i < n``"""
        self.n = n
        self.a = np.array(range(n), dtype = np.uint32)
        self.weights = np.ones(n, dtype = np.uint32)
    
    def __getitem__(self, i):
        """Find and return the root of the set containing ``i``."""
        j = self.a[i]  # j = parent of i; fails if i is not in the set
        k = self.a[j]  # k = parent of j
        if j == k:     # optimize the case that j is already the root     
            return j
                    
        j, k = k, self.a[k]       # If not: do one more step and advance
        while j  != k:            # the path until the root k is found.
            j, k = k, self.a[k]
                    
        j = self.a[i]             # j is the parent of input i
        while i != j:             # compress the path from i to the root
            i, j, self.a[i] = j, self.a[j], k
        return k                  # return root

    def __len__(self):
        return n

    def union(self, i1, i2):
        """Unite the sets containing the integers ``i`` and ``j``"""
        k1, k2 = self[i1], self[i2]
        if k1 != k2:
            w1, w2 = self.weights[k1], self.weights[k2]
            if w1 > w2:
                self.a[k2] = k1
                self.weights[k1] += k2
            else:
                self.a[k1] = k2
                self.weights[k2] += k1

    def union_perm(self, perm):
        """Unite set containing ``i`` with set containing ``perm[i]``"""
        for j1, j2 in enumerate(perm):
            self.union(j1, j2)

    def bitmaps(self):
        """Return sets as a sorted list of bitmaps"""
        d = defaultdict(int)
        for i in range(self.n):
            d[self[i]] |= 1 << i
        return list(sorted(d.values()))
            

def check_transitivity(mode, transitivity_bitmaps):
    def check_subdivision(v, transitivity_bitmaps):
        for bm in transitivity_bitmaps:
            assert bm & v in [0, bm]
    super_mode = complete_mat24_rand_mode(mode)
    for mask, v in zip((RAND_2, RAND_3, RAND_o), (0xc, 0xe, 0xff)):
        if super_mode & mask:
            check_subdivision(v, transitivity_bitmaps)
             
        
        
def do_test_local(mode, n = 10, verbose = 0):
    super_mode = complete_mat24_rand_mode(mode)
    assert complete_mat24_rand_mode_minimal(mode) == super_mode 
    assert complete_mat24_rand_mode(super_mode) == super_mode
    generated = RAND_USED
    R_MAX = MAT24_ORDER - 1
    sets = UIntUnionFind(24)
    if verbose:
        print("Testing mode 0x%02x, super = 0x%02x" % (mode, super_mode))
    
    for i in range(n):
        r = randint(0, R_MAX)
        pi = mat24_perm_rand_local(mode, r, verbose)
        mat24.perm_check(pi)
        sets.union_perm(pi)
        in_group =  mat24_perm_in_local(pi)
        ok = (in_group & super_mode) == super_mode 
        generated &= in_group
        if verbose > 1 or not ok:
            print("pi =", pi)
            print("  mode = 0x%02x, super = 0x%02x, generates 0x%02x" %
                    (mode, super_mode, in_group))
            if not ok:
                err = "Permutation is not in expected subgroup"
                raise ValueError(err)
    ok = generated == super_mode
    transitivity_bitmaps = sets.bitmaps()
    check_transitivity(mode, transitivity_bitmaps)
    if verbose:
        print("Transitivity:", [hex(x) for x in transitivity_bitmaps])
    if verbose > 1 or not ok:
        print("Super = 0x%02x, generates 0x%02x" % (super_mode, generated))
        if not ok:
            err = "Permutations do not generate the expected subgroup"
            raise ValueError(err)
        print("")


@pytest.mark.mat24
def test_local():
    for mode in range(64):
        do_test_local(mode, n = 30, verbose = 0)



if __name__ == "__main__":
    test_local(verbose = 1)
    #print([hex(x) for x in SUBGROUPS])
    print(len(SUBGROUPS), "subgroups constructed")



            
        
                
         


         
 

             

