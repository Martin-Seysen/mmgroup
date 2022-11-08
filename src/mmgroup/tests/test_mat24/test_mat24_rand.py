"""Generate random elements of subgroups of the Mathieu group ``M_24``


This is yet a fragment.

It implements the generation of random elements ofv arious subgroups
of the Mathieu group :math:`M_{24}`. Here each subgroup is 
describend by an integer ``mode``, where each bit in integer
``mode`` corresponds to a certain subgroup of  :math:`M_{24}`
fixing a certain subset (or a set of subsets) of :math:`M_{24}`.
Setting several bits in ``mode`` means intersection of the 
corresponding subgroup.


For mathematical background see [1]


[1] **The mmgroup guide for deveolpers**, section 
    **Subgroups of the Mathieu group :math:`M_{24}`**.

"""

import sys
from random import randint
from collections import defaultdict
from collections.abc import Sequence
from operator import __and__, __or__
from numbers import Integral
from functools import reduce

import numpy as np
import pytest


if __name__ == "__main__":
    sys.path.append("../../../")


from mmgroup import mat24
from mmgroup.mat24 import MAT24_ORDER 


####################################################################
# Encoding the intersection of subgroups of :math:`M_{24}`
####################################################################


RAND_2 = 1
RAND_o = 2
RAND_t = 4
RAND_s = 8
RAND_l = 16
RAND_3 = 32

RAND_FLAGS = [RAND_2, RAND_o, RAND_t, RAND_s, RAND_l, RAND_3]
assert reduce(__and__, RAND_FLAGS) == 0


RAND_USED = reduce(__or__, RAND_FLAGS)


for letter, flag in zip("2otsl3", RAND_FLAGS):
    assert  mat24.MAT24_RAND[letter] == flag




def subgroup(mode, mode_subgroup, mode_of):
    if mode &  mode_subgroup == mode_subgroup:
        mode |= mode_of
    return mode


def py_mat24_complete_rand_mode(mode):
    """Complete an intersection of subgroups of :math:`M_{24}`

    Given an integer ``mode`` representing an intersection :math:`H`
    of subgroups of :math:`M_{24}`, the group :math:`H` may be 
    contained in more subgroups of :math:`M_{24}` encoded in the
    bits of ``mode``. This function sets all bits in ``mode``
    corresponding to groups containing  :math:`H`. Furthermore,
    the function clears all unused bits in parameter ``mode``.
    
    The function returns the modified parameter ``mode``.

    It computes precisely the implications shown in [1]. The main
    test in this file contains a computational proof that all valid 
    implications are consequences of the implications in [1].
    """

    mode = subgroup(mode, RAND_l,  RAND_o)
    mode = subgroup(mode, RAND_l | RAND_2,  RAND_t)
    mode = subgroup(mode, RAND_t | RAND_2,  RAND_l | RAND_o)
    mode = subgroup(mode, RAND_3 | RAND_t,  RAND_o | RAND_s)
    mode = subgroup(mode, RAND_3 | RAND_l,  RAND_o | RAND_s)
    mode = subgroup(mode, RAND_3 | RAND_t | RAND_l,  RAND_o | RAND_2)
    return mode

def py_mat24_complete_rand_mode_minimal(mode):
    """Alternative implementation of ``py_mat24_complete_rand_mode``

    This function is slower than ``py_mat24_complete_rand_mode``
    but closer to the statement of the implications in [1]. The 
    test below checks the equivalence of the two implementations-
    """
    for i in range(3):
        mode = subgroup(mode, RAND_l,  RAND_o)
        mode = subgroup(mode, RAND_l | RAND_2,  RAND_t)
        mode = subgroup(mode, RAND_t | RAND_2,  RAND_o | RAND_l)
        mode = subgroup(mode, RAND_3 | RAND_l,  RAND_s)
        mode = subgroup(mode, RAND_3 | RAND_t,  RAND_o | RAND_s)
        mode = subgroup(mode, RAND_3 | RAND_t | RAND_l, RAND_2)
    return mode


# List of all intersections ofsubgroups not ccontained in any
# other subgroup.  
SUBGROUPS = [
    i for i in range(64) if py_mat24_complete_rand_mode(i) == i
]
 

####################################################################
# Checking membership in subgroups of  :math:`M_{24}`
####################################################################




def check_in_set(pi, start, diff):
    """Check if permutation ``pi`` fixes a certain set.

    This set is the set 
    {{i,...,i + diff -1} for i = start, start + diff, start + 2*diff...}
    """ 
    s = 0
    for i in range(start, 24, diff):
        for j in range(1, diff):
            s |= pi[i] ^ pi[i + j]
    return (s & -diff) == 0




def py_mat24_perm_in_local(pi):
    """Compute some subgroups containing an element of :math:`M_{24}`

    Let :math:`p_1` be a permutation in :math:`M_{24}` given as an 
    array of 24 integers. The function computes a set of subgroups 
    of :math:`M_{24}` containing :math:`p_1`. These computations are 
    done for all subgroups corresponding to the flags defined in 
    ``enum mat24_rand_flags``. The function returns an integer 
    ``mode`` that is the combination of flags of type 
    ``enum mat24_rand_flags`` describing the
    subgroups of :math:`M_{24}` containing  :math:`H` .
    
    The function raises ValueError if :math:`p_1` is not
    in :math:`M_{24}`.
    """
    mat24.perm_check(pi)
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
    if check_in_set(pi, 8, 2):
        mode |= RAND_l
    if check_in_set(pi, 0, 4):
        mode |= RAND_s
    if check_in_set(pi, 0, 8):
        mode |= RAND_t
    return mode
    




####################################################################
# Generating a random element of an intersection of subgroups
####################################################################




def free_bit24(mask, r):
    """Select a 'random' cleared bit from ``mask``

    Here we consider the lowest 24 bits of parameter ``mask`` only.
    Let ``w`` be the number of such cleared bits and put
    ``r1 = r / w``, ``b = r % w``. Then the function computes the
    position ``k`` of the ``b``-th cleared bit in mask ``mask``.
    It returns the triple ``k, mask1, r1``, where ``mask1`` is
    obtained from mask ``mask`` by setting bit ``k`` in the mask.
    """
    # mask = bitmask, only a cleared bit may be returned
    # r = random number
    mask &= 0xffffff
    w = 24 - mat24.bw24(mask)
    r, k = divmod(r, w)
    bmask = mat24.spread_b24(1 << k, 0xffffff & ~mask)
    return mat24.lsbit24(bmask), mask | bmask, r


def remask(h2, new_mask):
    """Return amask computed from array ``h2`` and mask ``new_mask``


    The returned mask contains all bits in the array ``h2`` and also
    all bits not set in ``new_mask``. It contains bits at positions
    less than 24 only.
    """ 
    mask = 0xffffff & ~new_mask
    for b in h2:
        mask |= 1 << b
    return mask & 0xffffff;

   
H1 = [3, 2, 1, 0, 5, 4, 8] 
 
def perm_from_suboctad(h2, mask, r):
    """Complete image of subset of octad to an element of ``M_24``

    Let the array ``h2`` of length ``<= 6`` be an image of the
    prefix of the same length of the octad given by ``H1``.
    Then this function completes the mapping  from (a prefix of)
    ``H1`` to ``h2`` to a permutation in ``M_24`` and returns
    that permutation.

    If necessary,  ``h1`` is first completed to a subset of an 
    octad of size 6 of by selecting bits cleared in the mask 
    ``mask``. We use function ``free_bit24`` for this purpose; 
    and parameter ``r`` has the same meaning as in that function.
    """
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



                 
        
                
def py_mat24_perm_rand_local(mode, r, verbose = 0):
    """Generate random element of a subgroup :math:`M_{24}`

    The function generates an element of a subgroup :math:`H` of
    the Mathieu group  :math:`M_{24}`.  Here the parameter ``u_mode`` 
    is a combination of flags of type ``enum mat24_rand_flags`` 
    describing the group :math:`H` as an intersection  of subgroups 
    of :math:`M_{24}`. The generated permutation is stored in the 
    array ``p_out`` of length 24.
 
    Parameter ``u_rand`` is an integer describing the element of
    subgroup :math:`H` to be generated. Here ``u_rand`` is reduced
    modulo the order of  :math:`H`. In order to generate a uniform 
    random element of :math:`H`, the user should generate a uniform 
    random number ``0 <= u_rand < MAT24_ORDER``, where ``MAT24_ORDER`` 
    is the order of the group :math:`M_{24}`.
    """
    if (mode & RAND_USED) == 0:
        return mat24.m24num_to_perm(r % MAT24_ORDER)       
    mode = py_mat24_complete_rand_mode(mode)
    return complete_perm(mode, r, verbose)





####################################################################
# Checking classes of transitivity 
####################################################################



class UIntUnionFind(Sequence):
    """Union-find algorithm for nonnegative integers less than n

    Here the given set is the set of singletons ``set(i)`` for
    integers ``0 <= i < n``, with ``n`` given by the constructor.

    ``self.union(i1, i2)`` unites the set containing element ``i1``
    with the element containing ``i2``. 

    ``self[i]`` returns a representative of the set containing ``i``.
    """
    def __init__(self, n):
        """Initialize set with singletons ``{i}``, ``0 <= i < n``"""
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
                    
        j = self.a[i]             # Let j be the parent of input i
        while i != j:             # compress the path from i to the root
            i, j, self.a[i] = j, self.a[j], k
        return k                  # return root

    def __len__(self):
        return len(self.a)

    def union(self, i1, i2):
        """Unite the sets containing the integers ``i1`` and ``i2``"""
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
        for i in range(len(self.a)):
            d[self[i]] |= 1 << i
        return list(sorted(d.values()))
            

def check_transitivity(mode, transitivity_bitmaps):
    def check_subdivision(v, transitivity_bitmaps):
        for bm in transitivity_bitmaps:
            assert bm & v in [0, bm]
    super_mode = py_mat24_complete_rand_mode(mode)
    for mask, v in zip((RAND_2, RAND_3, RAND_o), (0xc, 0xe, 0xff)):
        if super_mode & mask:
            check_subdivision(v, transitivity_bitmaps)
             
 

####################################################################
# Testing crreation of random elements of a subgroup
####################################################################

       
        
def do_test_mat24_rand(mode, n, verbose = 0):
    """Test generation of random elements of subgroup given by ``mode``
 
    We generate ``n`` random elemnts ``g``  of the subgroup ``H`` of
    ``M_24`` using function  ``py_mat24_perm_rand_local``. Then we 
    test if ``g`` is actually in the required subgroup using function
    ``py_mat24_perm_in_local``.

    Given ``mode``, we compute a variable  ``super_mode`` =
    ``py_mat24_complete_rand_mode(mode)``. Then in variable 
    ``super_mode`` some (but not necessarily all) bits corresponding 
    to the subgroups ``K_i`` of ``M_{24}``  containing ``H`` are set.
    It has been shown in [1] that all inclusions of intersections of 
    subgroups  computed by function ``py_mat24_complete_rand_mode`` 
    are correct.

    One objective of this test function is to show that there are no
    further inclusions, i.e. that variable ``super_mode`` has maximum 
    possible bit weight. Therefore  we check that for each valid bit 
    cleared in ``super_mode``  there is an element ``g_i`` in ``H``,
    that is not contained in the group ``K_i`` corresponding to that
    bit. 

    So running this test for all possible combinations  of valid bits 
    of ``mode`` shows that there are no more inclusions than given
    in [1].
    
    Apart from this we make a few further tests. For subgroups of
    ``M_24`` fixing a subset of ``0,...,23`` we check that 
    transitivity of ``g`` is as expected. We also check that
    ``py_mat24_complete_rand_mode(mode)`` is equal to
    ``py_mat24_complete_rand_mode_minimal(mode)``. 

    We also check that the C functions ``mat24.perm_rand_local``, 
    ``mat24.perm_in_local``,  and ``mat24.complete_rand_mode`` are 
    equivalent to the  functions ``py_mat24_perm_rand_local``,
    ``py_mat24_perm_in_local``, and ``py_mat24_complete_rand_mode``,
    respecively.
    """
    super_mode = py_mat24_complete_rand_mode(mode)
    assert py_mat24_complete_rand_mode_minimal(mode) == super_mode 
    assert py_mat24_complete_rand_mode(super_mode) == super_mode
    assert mat24.complete_rand_mode(mode) == super_mode
    generated = RAND_USED
    R_MAX = MAT24_ORDER - 1
    sets = UIntUnionFind(24)
    if verbose:
        print("Testing mode 0x%02x, super = 0x%02x" % (mode, super_mode))
    
    for i in range(n):
        r = randint(0, R_MAX)
        pi = py_mat24_perm_rand_local(mode, r, verbose)
        pi_c = mat24.perm_rand_local(mode, r)
        assert pi_c == pi
        num_pi = mat24.perm_to_m24num(pi)
        assert mat24.m24num_rand_local(mode, r) == num_pi
        mat24.perm_check(pi)

        sets.union_perm(pi)
        in_group =  py_mat24_perm_in_local(pi)
        assert in_group == mat24.perm_in_local(pi)
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
def test_mat24_rand(verbose = 0):
    for mode in range(64):
        do_test_mat24_rand(mode, n = 40, verbose = verbose)

####################################################################
# The main test
####################################################################


if __name__ == "__main__":
    print("Starting test ...")
    test_mat24_rand(verbose = 0)
    #print([hex(x) for x in SUBGROUPS])
    print("passed")
    print(len(SUBGROUPS), "subgroups constructed")



            
        
                
         


         
 

             

