"""Support for completing a permutation in the Mathieu group

A *heptad* is s subset of seven elements of the set on which that
Mathieu group acts. Under certain circumstances a mapping from one
heptad to another heptad can be completed to a unique element of
``mat_24``. The C function ``mat24_perm_from_heptads`` in file
``mat24_functions.c`` performs that task. It calls a simpler 
function ``mat24_perm_complete_heptad`` in that file which 
completes a mapping from a fixed heptad to variable heptad. 

This file contains python implementations of functions
``mat24_perm_complete_heptad`` and ``mat24_perm_from_heptads``.

The C function ``mat24_perm_to_m24num`` maps the elements of the 
Mathieu group ``Mat_24`` to the set if integers 
``0, ..., order(Mat_24)``. Function ``mat24_perm_complete_heptad`` 
is also used as a subroutine of function ``mat24_m24num_to_perm``
which compute the inverse of that mapping. This module also
contains python implementations of these functions.

Class ``HeptadCompleter`` is a table-providing class that is used
for code generation. It also exports the functionality of this
module.

Function ``hint_for_complete_heptads`` prints an explanation
of the implementation of function  ``mat24_perm_complete_heptad``,
including some precalulated data required for that implementation.

Function ``test_complete_octad`` test the correcteness of function
``mat24_perm_complete_heptad``. Here it suffices to check that
the identity permutation is generated correctly from the required
input values, and that the intersections of (possibly complemented)
octads and syndromes are indeed singletons if the ought to be
singletons.
"""


import sys
sys.path.append(r"C:\Data\projects\MonsterGit\src")

from mmgroup.dev.mat24.mat24tables import Mat24Tables as gc

#######################################################################
# Finding octads suitable for function complete_heptad()
#######################################################################

def lrange(*args):
    return list(range(*args))

def list_to_set24(l):
    """Convert a bit vector to an integer"""
    res = 0
    for x in l: res  ^= 1 << x
    return res & 0xffffff

def set24_to_list(v):
    """Convert an integer to a bit vector"""
    return [x for x in range(24) if v & (1 << x)]


_OCTAD_LIST = None
def octad_list():
    """Return the sorted list of all octads
    
    Each octad is given as a 24-bit integer describing a bit 
    vector in `vector`` representation. Octads are sorted 
    lexicographically, so lower bits of the integer are sorted
    with higher priority.        
    """
    global _OCTAD_LIST
    if _OCTAD_LIST is not None:
        return _OCTAD_LIST[:]
    l_oct = [set24_to_list(gc.octad_to_vect(o)) for o in range(759)]
    l_oct.sort()
    _OCTAD_LIST = [list_to_set24(l) for l in l_oct]
    return _OCTAD_LIST[:]

def superoctads(contained, not_contained, partial = [], weight = 0):
    """Return a specific list of octads
    
    The function returns a specific list of octads, where each 
    octad is given as a 24-bit integer describing a bit vector in 
    ``vector`` representation.
    
    An octad is included in the list if it satisfies the following
    requiremets.
     
      * All entries in the list ``contained ``must be contained 
        in the octad.
        
      * No entry in the list ``not_contained`` may be containe 
        in the octad.

      * Precisely ``weight`` entries of the list ``partial``
        must be in the octad.

    Each parameter ``contained `, ``not_contained``, or ``weight``
    may be a list of different integers < 24 or a 24-bit integer
    representing a bit vetor.    
    """
    if not isinstance(contained, int): 
        contained =  list_to_set24(contained)   
    if not isinstance(not_contained, int): 
        not_contained =  list_to_set24(not_contained)  
    if not isinstance(partial, int): 
        partial =  list_to_set24(partial)  
    return [x for x in octad_list() if (x & contained == contained)
        and (x & not_contained == 0)
        and (gc.bw24(x & partial) == weight) 
    ]

def l_superoctads(contained, not_contained, partial, weight):
    """Return a specific list of octads
    
    The function returns a spcific list of octads. Input parameters
    are as in function ``superoctads``, and also the same list of 
    octads is returns.
    
    But in contrast to function `superoctads`` each entry in that
    list of octad is represented as a sorted list of 8 integers
   
    THis function may be used for manually finding suitable octads 
    satisfying certain requiremnts.
    """
    l = superoctads(contained, not_contained, partial, weight)
    return [set24_to_list(v) for v in l]

def print_superoctads(contained, not_contained, partial, weight):
    """Print list of octads found by function ``superoctads``
    
    Input parameters are as in function ``superoctads``. Each octad
    is printed as a sorted list of 8 integers. Octads are printed
    innlexicographic order.
    """
    l = l_superoctads(contained, not_contained, partial, weight) 
    for x in l: print(x)
    if len(l) == 0: print("No octad found")

def hint_for_complete_heptads():
    """Print a hint for programming function ``complete_heptad``"""
    
    s = """\nHint for programming function ``complete_heptad()``    
    
Function ``complete_heptad()`` completes a partial 
permutation of 24 elements to a  permutation on 24 elements 
``0,..,24`` which is in the Mathieu group  ``Mat24``. Here the 
images of 7 elements of the set ``(0,..,24)`` must be given at 
the positions as indicated below. Then the constructed 
permutation is unique in the group `Mat24``.

Function ``complete_heptad()`` proceeds by completing pentads
(i.e. arbitrary subsets of size 5 of) to octads. If the images 
of the 5 elements of are known then the image of the octad is 
also known (up to permutations of the elements of the octad). 
By intersecting pentads and octads an their complements we 
will eventually split the set ``(0,..,24)`` into  singletons, 
so the the images of all elements of that set are known. 

We list suitable octads for performing the next intersection
step and also the singletons found after each step.
"""
    print(s)
    def add_to_found(found, new_found):
        print("\nSingletons found:", new_found, ". Useful octads:")
        return found + new_found
    found = [0,1,2,3,4,5,8]
    print("Given", found, ". Useful octads:")
    print_superoctads([9],lrange(14,24),found,5)
    found = add_to_found(found, [9])
    print_superoctads(lrange(8),[],found,6)
    print_superoctads([6],[7, 18, 19, 20, 21, 22, 23],found,5)
    found = add_to_found(found, [6,7])
    print_superoctads([10,12,14],[2,],found,5)
    found = add_to_found(found, [10,11,12,13,14])
    print_superoctads([15,16,18],[17,7],found,5)
    found = add_to_found(found, [16,17])
    print_superoctads([15,19,20],[17, 18],found,5)
    found = add_to_found(found, [15,18])
    print_superoctads([19,21,22],[15, 17, 20],found,5)
    found = add_to_found(found, [19,20])
    print_superoctads([21,23],[2,7,12,15,17,18,19,20,22],found,6)
    print("")




#######################################################################
# Auxiliary functions for function complete_heptad()
#######################################################################





def lsb24(x):
    """Return binary logarithm of ``x`` if ``x`` is a power of 2"""
    res = int(x).bit_length() - 1
    assert x == 1 << res, hex(x)
    return res

def odd_syn(v):
    """Return Golay cocode syndrom of a bit vector v

    The bit vector v must have odd parity.
    """
    coc = gc.vect_to_vintern(v)
    t = gc.syndrome_table[coc & 0x7ff ]
    return ( (1 << (t & 31)) ^ (1 << ((t >> 5) & 31)) ^
             (1 << ((t >> 10) & 31)) )


#######################################################################
# Function complete_heptad()
#######################################################################


""" MOG

    0      4      8   C=12   G=16   K=20
    1      5      9   D=13   H=17   L=21
    2      6   A=10   E=14   I=18   M=22
    3      7   B=11   F=15   J=19   N=23
"""

def mat24_complete_heptad(p_io):
    """Complete a permutation p given by p_io to an element of  Mat24.
   
    p must have length 24. Entries p[i], i = 0,1,2,3,4,5,8 must make
    up a valid umbral heptad, i.e. a heptad not contained in an octad.
    p[0],...,p[5] must be contained in an octad, p[8] must not be
    contained in that octad. The other entries of input p are ignored.
   
    It can be shown that such a permutation p can be completed to
    a unique element of Mat24.
   
    If a suitable permutation p an be found, the array p_io is
    completed so that it contains permutation p, and the function
    returns 0.
   
    Otherwise the function returns a nonzero value.
    """
    ## uint_fast32_t  err, s1, s5, s015, s3, s4, s8, s01234; 
    ## uint_fast32_t  s567, s67, s9AB, s9CD, s9, s6GH, s6;
    ## uint_fast32_t  sACE, sD, sFGI, sG, sFJK, sJLM, sALN;

    err = ((p_io[0] + 8) | (p_io[1] + 8) | (p_io[2] + 8) | 
      (p_io[3] + 8) | (p_io[4] + 8) | (p_io[5] + 8) | (p_io[8] + 8));
    err &= (0 - 0x20);
    s1 = 1 << p_io[1]; 
    s5 = 1 << p_io[5];
    s015 = (1 << p_io[0]) ^ s1 ^ s5;
    s3 = 1 << p_io[3]; 
    s4 = 1 << p_io[4]; 
    s8 = 1 << p_io[8];
    s01234 = s015 ^ s5 ^ (1 << p_io[2]) ^ s3 ^ s4;
    # if err == 0 then 0 <= s01234 < 0x1000000 has odd parity
    # octad = [0, 1, 2, 3, 4, 5, 6, 7]
    s567 = odd_syn(s01234);
    err |= (s01234) & s567;
    # if err == 0 then entries [0,1,2,3,4,6] are in an octad
    err |= (s01234 | s567) & s8;
    # if err == 0 then entry 8 is not in that octad
    err |= s5 ^ (s5 & s567);
    s67 = s567 & ~s5;
    # octad = [0, 1, 2, 3, 8, 9, 10, 11]
    s9AB = odd_syn(s01234 ^ s4 ^ s8);
    # octad = [0, 1, 4, 5, 8, 9, 12, 13]
    s9CD = odd_syn(s015 ^ s4 ^ s8);
    s9 = s9AB & s9CD;
    p_io[9] = lsb24(s9);
    # octad [1, 3, 5, 6, 8, 9, 16, 17]
    s6GH = odd_syn(s1 ^ s3 ^  s5 ^ s8 ^ s9);
    s6 = s67 &  s6GH;
    p_io[6] = lsb24(s6);
    p_io[7] = lsb24(s67 & ~s6GH);
    # still needed: 
    #   err, s1, s015, s3, s8, s01234, s9AB, s9CD, s9, s6GH, s6
    # octad [0, 2, 4, 6, 8, 10, 12, 14]
    sACE = odd_syn(s01234 ^ s1 ^ s3 ^ s6 ^ s8);
    p_io[10] = lsb24(s9AB & sACE);
    p_io[11] = lsb24(s9AB & ~sACE & ~s9);
    p_io[12] = lsb24(s9CD & sACE);
    sD = s9CD & ~sACE & ~s9;
    p_io[13] = lsb24(sD);
    p_io[14] = lsb24(sACE & ~s9AB & ~s9CD);
    # still needed: 
    #   err, s1, s015, s3, s8, s6GH, s6, sACE, sD
    # octad [0, 1, 5, 6, 13, 15, 16, 18]
    sFGI = odd_syn(s015 ^ s6 ^ sD);
    sG = s6GH & sFGI;
    p_io[16] = lsb24(sG);
    p_io[17] = lsb24(s6GH & ~s6 & ~sFGI);
    # octad [0, 1, 3, 5, 8, 15, 19, 20]
    sFJK = odd_syn(s015 ^ s3 ^ s8);
    p_io[15] = lsb24(sFGI & sFJK);
    p_io[18] = lsb24(sFGI & ~sG & ~sFJK);
    # octad [0, 3, 5, 6, 16, 19, 21, 22]
    sJLM = odd_syn(s015 ^ s1 ^ s3 ^ s6 ^ sG);
    p_io[19] = lsb24(sFJK & sJLM);
    p_io[20] = lsb24(sFJK & ~sFGI & ~sJLM);
    # octad [0, 1, 5, 6, 8, 10, 21, 23]
    sALN = odd_syn(s015 ^ s6 ^ s8); 
    p_io[21] = lsb24(sALN & sJLM);
    p_io[22] = lsb24(sJLM & ~sALN & ~sFJK);
    p_io[23] = lsb24(sALN & ~sACE & ~sJLM);
    return err;

#######################################################################
# Function based on function complete_heptad()
#######################################################################



STD_OCTAD = gc.vect_to_octad(0xff)

def mat24_int_to_perm(k):
    """Return the k-th permutation in the Mathieu group Mat24   

    Any integer 0 <= k < 244823040 is evaluated in mixed-radix with 
    bases  759, 8, 7, 6, 5, 4, 3, 16, with valence decreasing from 
    left to right. The first digit determines the number of the octad 
    which is the image  of the standard octad (0,...,7). The following 
    six digits determine an even permutation of that octad, up to a 
    permutation of the last two entries. The final digit determines 
    the image of element 8, i.e.  the first element not in the 
    standard octad. 

    The images of the remaining elements 6, 7 and 9,...,23  are 
    determined by calling function complete_heptad()
    """
    oct, k = divmod(k, 322560)
    if oct >= 759: return None
    oct -= 759 - STD_OCTAD
    oct += (oct >> 12) & 759 # give number 0 to standard octad
    #print("i2p", oct)
    oct = gc.octad_to_vect(oct);
    #print("i2p oct", hex(oct))
    p = [None]*24
    oct, j = 8 * oct, 0x8        
    for i in range(24):
        o = oct & 8
        p[(j >> o) & 0x1f] = i
        j += 1 << o
        oct >>= 1
    p[8] = p[8 + (k & 15)]
    #print("i2pfinal", k & 15)
    k >>= 4
    k *= (1 << 28) // 2520 + 1
    for i in range(6):
        k1 = i + (k >> 28)
        #print("i2p%d" % i, k >> 28)
        p[i], p[k1] = p[k1], p[i]  
        k = (k & 0xfffffff) * (7-i)
    mat24_complete_heptad(p)
    return p

def mat24_perm_to_int(p):
    """Convert permutation p in the Mathieu group Mat24 to an integer.

    This reverses member function int_to_perm(). The input permutation
    is not checked.
    """
    oct = sum(1 << x for x in p[:8])
    res = gc.vect_to_octad(oct) 
    #print("p2i oct", hex(oct))
    res -=  STD_OCTAD
    res += (res >> 12) & 759 
    #print("p2i", res)
    p1 = [24]*32
    oct, j = 8 * oct, 0x00        
    for i in range(24):
        o = oct & 8
        p1[i] = (j >> o) & 0x1f
        j += 1 << o
        oct >>= 1
    q, q_inv = [None]*8, [None]*8
    for i in range(8):
        j = p1[p[i] & 0x1f] & 7
        q[j] = i
        q_inv[i] = j
    for i in range(6):
        # exchange place i with place q_inv[i]
        j = q_inv[i]
        #q_inv[q[i]], q_inv[q[j]] = q_inv[q[j]], q_inv[q[i]]
        #q[i], q[j] = q[j], q[i]
        #assert q[:i] == q_inv[:i] == lrange(i)
        q_inv[q[i]] = q_inv[q[j]]
        q[j] = q[i]
        #print("p2i%d" % i, j-i)           
        res = res * (8 - i) + j - i
    #print("p2ifinal", p1[p[8] & 0x1f])  
    return 16 * res + p1[p[8] & 0x1f]
            
            
def mat24_perm_from_heptads(h1, h2):
    """Try to find a permutation p that maps heptad h1 to h2

    h1 and h2 must be lists of length 7 defining two umbral heptads,
    i.e. heptads not contained in an octad. If a permutation p in
    the Mathieu group Mat24 that maps h1 to h2 exists, it is unique. 

    Return permutation p if such a p exists an is unique,
    and return None otherwise.
    """
    # First find the special element of v h1 not contained in the octad
    v = 0
    for i in range(7):
        v ^= 1 << (h1[i] & 31)
    y = odd_syn(v)
    v =  lsb24(v & y)
    
    # Find position y of element v in h1
    y = 0
    for i in range(7): 
        y |= ((h1[i] != v) - 1) & i

    # Copy special element of h1 to position 8 of p1 and copy the other
    # elements of h1 to positions 0,...,6. Copy h2 similarly to p2
    p1 = h1[:7] + [None]*17
    p2 = h2[:7] + [None]*17
    p = [None] * 24
    p1[8] = p1[y]
    p1[y] = p1[6]
    p2[8] = p2[y]
    p2[y] = p2[6]

    # Complete p1 and p2 from heptad. Return error if any completion fails
    if mat24_complete_heptad(p1) | mat24_complete_heptad(p2):
        return None

    # If success, return p1**(-1) * p2
    for i in range(24):
        p[p1[i]] = p2[i]
    return p


#######################################################################
# Support for the code generator
#######################################################################


class HeptadCompleter(object):

    def __init__(*args):
        # Ignore arguments in constructor for compatibility
        # with older versions
        pass

    @staticmethod
    def complete_heptad(p_io):
        return  mat24_complete_heptad(p_io)

    @staticmethod
    def int_to_perm(k):
        return  mat24_int_to_perm(k)

    @staticmethod
    def perm_to_int(p):
        return  mat24_perm_to_int(p)

    @staticmethod
    def perm_from_heptads(h1, h2):
        return mat24_perm_from_heptads(h1, h2)

    @staticmethod
    def tables():
        """Return dictionary with the constant STD_OCTAD. """
        return {
            "MAT24_STD_OCTAD" : STD_OCTAD,
        }

    @staticmethod
    def directives():
        return { }


#######################################################################
# Testing this module
#######################################################################


def test_complete_octad():
    io = lrange(6) + [None,None,8] + [None]*15    
    complete_heptad(io)
    assert io == lrange(24)   
    print("The python function ``complete_heptad()`` is correct")

if __name__ == "__main__":
    hint_for_complete_heptads()
    test_complete_octad()
    
    