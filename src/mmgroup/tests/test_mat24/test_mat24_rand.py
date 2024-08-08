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
RAND_d = 64

RAND_FLAGS = [RAND_2, RAND_o, RAND_t, RAND_s, RAND_l, RAND_3, RAND_d]
assert reduce(__and__, RAND_FLAGS) == 0


RAND_USED = reduce(__or__, RAND_FLAGS)
assert RAND_USED == 127

for letter, flag in zip("2otsl3d", RAND_FLAGS):
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
    corresponding to groups containing  :math:`H`. 

    Corresponding C function:  mat24_complete_rand_mode
    """
    mode &= RAND_USED 
    mode_old = 0
    while (mode != mode_old):
        mode_old = mode;
        mode = subgroup(mode, RAND_d,  RAND_t); # [Iva09], Lemma 4.1.3
        mode = subgroup(mode, RAND_l,  RAND_o); # by definition
        mode = subgroup(mode, RAND_d | RAND_3, RAND_2); # by definition
        mode = subgroup(mode, RAND_l | RAND_t, RAND_d); 
        mode = subgroup(mode, RAND_l | RAND_2, RAND_o | RAND_t); 
        mode = subgroup(mode, RAND_t | RAND_2, RAND_o | RAND_l); 
        mode = subgroup(mode, RAND_l | RAND_3, RAND_o | RAND_s); 
        mode = subgroup(mode, RAND_t | RAND_3, RAND_o | RAND_s); 
        mode = subgroup(mode, RAND_o | RAND_t | RAND_d, RAND_l);
    return mode


# List of all intersections of subgroups not contained in any
# other subgroup.  
SUBGROUPS = [
    i for i in range(128) if py_mat24_complete_rand_mode(i) == i
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
    if check_in_set(pi, 8, 2) and mode & RAND_o:
        mode |= RAND_l
    if check_in_set(pi, 0, 4):
        mode |= RAND_s
    if check_in_set(pi, 0, 8):
        mode |= RAND_t
    if check_in_set(pi, 0, 2):
        mode |= RAND_d
    return mode
    




####################################################################
# Generating a random element of an intersection of subgroups
####################################################################



class RandPiS:
    """Equivalent to structure ``rand_pi_s`` in file ``mat24_random.c``
    """
    AL = [0,1,2,3,4,5,7,6]
    H1 = [3, 2, 1, 0, 5, 4, 8]
    def __init__(self, mode, rand):
        """Corresponding C function:  rand_pi_prepare
        """
        self.mode = py_mat24_complete_rand_mode(mode);
        self.rand = rand;
        self.bitmap = 0;
        self._debug = np.zeros(32, dtype = np.uint32)
        self._debug[0] = mode;
        self._debug[1] = self.rand;
        self._debug[2] = self.mode;
        self.h = [None] * 7

    def rand_int(self, n):
        """Corresponding C function:  rand_int
        """
        self.rand, i = divmod(self.rand, n)
        return i

    def add_mat24_point(self, index, mask):
        """Corresponding C function:  add_mat24_point
        """
        available = ~self.bitmap & mask & 0xffffff;
        w = mat24.bw24(available);
        if (w == 0):
            print(self.debug_info())
            print("index: %d, mask: 0x%06x, mode: %02x, h: %s" % (
                 index, mask, self.mode, self.h)) 
            ERR = "Could not add a point to a class RandPiS object"
            raise ValueError(ERR) 
        b = self.rand_int(w);
        bmask = mat24.spread_b24(1 << b, available);
        self.bitmap |= bmask;
        self.h[index] = b = mat24.lsbit24(bmask);
        self._debug[4 * index + 4] = self.bitmap;
        self._debug[4 * index + 5] = mask;
        self._debug[4 * index + 6] = b;

    @staticmethod
    def complete_aff_trio(h1, h2, h3):
        """Corresponding C function:  complete_aff_trio (roughly)
        """
        return h1 ^ h2 ^ h3

    @staticmethod
    def complete_aff_line(h1, h2, h3):
        """Corresponding C function:  complete_aff_line (roughly)
        """
        AL = RandPiS.AL
        return AL[AL[h1] ^ AL[h2] ^ AL[h3]]

    def find_img_0(self): 
        """Corresponding C function has same name 
        """
        fix = 0xffffff;
        if (self.mode & RAND_o): fix &= 0xff;
        if (self.mode & RAND_3): fix &= 0x0e;
        if (self.mode & RAND_2): fix &= 0x0c;
        return fix;

    def find_img_1(self): 
        """Corresponding C function has same name 
        """
        fix = 0xffffff;
        if (self.mode & RAND_t or self.mode & RAND_o):
            fix &= 0xff << (self.h[0] & 0xf8);
        self.mask_octad = fix;
        if (self.mode & RAND_s):
           fix &= 15 << (self.h[0] & 0xfc);
        self.mask_tetrad = fix;
        if (self.mode & RAND_3): fix &= 0x0e;
        if (self.mode & RAND_2): fix &= 0x0c;
        if (self.mode & RAND_d): fix &= 1 << (self.h[0] ^ 1);
        return fix;

    def find_img_2(self): 
        """Corresponding C function has same name 
        """
        fix = self.mask_tetrad;
        if (self.mode & RAND_3): fix &= 0x0e;
        return fix;

    def find_img_3(self): 
        """Corresponding C function has same name 
        """
        fix = self.mask_tetrad;
        h = self.h
        if (self.mode & RAND_d):
            fix &= 1 << (self.h[2] ^ 1);
        elif (self.mode & RAND_t):
            fix &= 1 << self.complete_aff_trio(h[0], h[1], h[2]);
        elif (self.mode & RAND_l):
            fix &= 1 << self.complete_aff_line(h[0], h[1], h[2]);
        return fix;

    def find_img_4(self): 
        """Corresponding C function has same name 
        """
        return self.mask_octad;

    def find_img_5(self): 
        """Corresponding C function has same name 
        """
        fix = self.syn = mat24.syndrome(self.bitmap, 0);
        h = self.h
        if (self.mode & RAND_d):
            fix &= 1 << (self.h[4] ^ 1);
        elif (self.mode & RAND_t):
            fix &= 1 << self.complete_aff_trio(h[0], h[1], h[4]);
        elif (self.mode & RAND_l):
            fix &= 1 << self.complete_aff_line(h[0], h[1], h[4]);
        return fix;

    def find_img_6(self): 
        """Corresponding C function has same name 
        """
        return ~(self.syn) & 0xffffff;

    def complete_perm(self):
        """Corresponding C function:  complete_perm
        """
        mask = self.find_img_0(); 
        self.add_mat24_point(0, mask);
        mask = self.find_img_1(); 
        self.add_mat24_point(1, mask);
        mask = self.find_img_2(); 
        self.add_mat24_point(2, mask);
        mask = self.find_img_3(); 
        self.add_mat24_point(3, mask);
        mask = self.find_img_4(); 
        self.add_mat24_point(4, mask);
        mask = self.find_img_5(); 
        self.add_mat24_point(5, mask);
        mask = self.find_img_6(); 
        self.add_mat24_point(6, mask);
        return mat24.perm_from_heptads(self.H1, self.h);


    @staticmethod
    def format_debug_info(data, as_str = True, source = None):
        if not as_str:
            return data
        s = "Debug info for function complete_perm%s:\n" % (
            " from " + str(source) if source else "")
        if len(data) == 0:
            return s + "(not available)\n"
        for i, x in enumerate(data[:32]):
            s += ("%08x%s" % (x, "\n" if i % 8 == 7 else " "))
        return s + "\n"      

    @staticmethod
    def debug_info_C(as_str = True):
        try:
            from mmgroup.mat24 import perm_rand_debug_info
            data = perm_rand_debug_info()
        except:
            data = []
        return RandPiS.format_debug_info(data, as_str, source = 'C')

    def debug_info(self, as_str = True):
        return self.format_debug_info(self._debug, as_str, "python")

                
def py_mat24_perm_rand_local(mode, u_rand, verbose = 0):
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
        return mat24.m24num_to_perm(u_rand % MAT24_ORDER) 
    p_i = RandPiS(mode, u_rand)
    perm = p_i.complete_perm()     
    if verbose:
         print("vvv", verbose)
         print(p_i.debug_info())
    return perm





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
# Testing creation of random elements of a subgroup
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
    transitivity of ``g`` is as expected.  

    We also check that the C functions ``mat24.perm_rand_local``, 
    ``mat24.perm_in_local``,  and ``mat24.complete_rand_mode`` are 
    equivalent to the  functions ``py_mat24_perm_rand_local``,
    ``py_mat24_perm_in_local``, and ``py_mat24_complete_rand_mode``,
    respecively.
    """
    super_mode = py_mat24_complete_rand_mode(mode)
    assert py_mat24_complete_rand_mode(super_mode) == super_mode
    assert mat24.complete_rand_mode(mode) == super_mode
    generated = RAND_USED
    R_MAX = MAT24_ORDER - 1
    sets = UIntUnionFind(24)
    if verbose:
        print("Testing mode 0x%02x, super = 0x%02x" % (mode, super_mode))
    
    pi_list = []
    for i in range(n):
        r = randint(0, R_MAX)
        pi = py_mat24_perm_rand_local(mode, r, verbose)
        try:
            pi_c = mat24.perm_rand_local(mode, r)
        except:
            print(RandPiS.debug_info_C())
            raise
        assert pi_c == pi
        num_pi = mat24.perm_to_m24num(pi)
        assert mat24.m24num_rand_local(mode, r) == num_pi
        mat24.perm_check(pi)

        sets.union_perm(pi)
        in_group =  py_mat24_perm_in_local(pi)
        in_group_C =  mat24.perm_in_local(pi)
        assert in_group == in_group_C, (hex(in_group), hex(in_group_C))
        ok = (in_group & super_mode) == super_mode 
        generated &= in_group
        if verbose > 1 or not ok:
            print("pi =", pi)
            print("  mode = 0x%02x, super = 0x%02x, generates 0x%02x" %
                    (mode, super_mode, in_group))
            if not ok:
                err = "Permutation is not in expected subgroup"
                raise ValueError(err)
        if i < 20:
            pi_list.append(pi) 
    ok = generated == super_mode
    transitivity_bitmaps = sets.bitmaps()
    check_transitivity(mode, transitivity_bitmaps)
    if verbose or not ok:
        print("Transitivity:", [hex(x) for x in transitivity_bitmaps])
    if verbose > 1 or not ok:
        if not ok:
            print(RandPiS.debug_info_C())
            print("  mode = 0x%02x, super = 0x%02x, generates 0x%02x" %
                    (mode, super_mode, generated))
            print('Some ermutations:')
            for i, pi in enumerate(pi_list):
                for j, x in enumerate(pi):
                    print("%2d%s" % (x, "\n" if j == 23 else " "), end = "")
        print("Super = 0x%02x, generates 0x%02x" % (super_mode, generated))
        if not ok:
            err = "Permutations do not generate the expected subgroup"
            raise ValueError(err)
        print("")


@pytest.mark.mat24
def test_mat24_rand(verbose = 0):
    for mode in range(0x80):
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



            
        
                
         


         
 

             

