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
#######################################################################
# Finding octads suitable for function complete_heptad()
#######################################################################
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
#######################################################################
# Function complete_heptad()
#######################################################################
#######################################################################


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
    coc = gc.vect_to_vintern(int(v))
    t = int(gc.syndrome_table[coc & 0x7ff])
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
    for i, x in enumerate(p_io):
        if x is not None:
            p_io[i] = int(x)
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
#######################################################################
# Function based on function complete_heptad()
#######################################################################
#######################################################################


#######################################################################
# Compute the n-th element of the Mathieu group Mat24
#######################################################################


STD_OCTAD = gc.vect_to_octad(0xff)




def split_n_mat24(n):
    """Auxiliary function for function ``py_mat24_int_to_perm``

    This function converts the integer 0 <= n < 244823040 to
    mixed-radix notation with bases [16, 3, 20, 21, 22, 23, 24].
    It returns the list of the 8 corresponding digits, with the
    highest digit first.
    """
    n_list = []
    for d in [16, 3, 20, 21, 22, 23, 24]:
        n, r = divmod(int(n), d)
        n_list.append(r)
    n_list.reverse()
    return n_list

def py_mat24_int_to_perm(n):
    """Return the k-th permutation p in the Mathieu group Mat24   

    Here the elements of the Mathieu group Mat24 are numbererd
    in lexicographic order.

    Any integer 0 <= n < 244823040 is evaluated in mixed-radix with 
    bases  24, 23, 22, 21, 20, 3, 16, with valence decreasing from 
    left to right. The first five digits determine the images of
    the first five elements of the permutation. The next digit
    determines the image of entry 5. That entry has 3 possible 
    images, depending on the syndrome of the images of the first
    five elements. The final digit determines the image of element 8, 
    i.e.  the first element not in the standard octad. 

    The images of the remaining elements 6, 7 and 9,...,23  are 
    determined by calling function complete_heptad()

    This implementation of the function is easy to understand but
    not fast. Function ``mat24_int_to_perm`` is a faster version of
    this function.
    """
    n = int(n)
    assert 0 <= n < MAT24_ORDER
    # Compute the output permutation in ``p``
    p = [None] * 24
    # List of elements to be permuted by Mat24
    p_list = list(range(24))
    # Convert n to mixed-radix digit list ``n_list`` as described above.
    n_list = split_n_mat24(n)
    # Compute the images of the first five numbers 0,...,4 under p
    # from the the first five digits of the digit list.
    # After assigning an image, delete that image from ``n_list``,
    # in order to avoid repetitions.
    for i in range(5):
        index = n_list[i] 
        p[i] = p_list[index]
        del p_list[index]
    # Compute the syndrome these five images in the list ``syn``.
    # That syndrome has length 3.  
    bitmap = sum(1 << i for i in p[:5])
    cocode = gc.vect_to_cocode(bitmap)
    syn_tab = int(gc.syndrome_table[cocode & 0x7ff])
    syn = [(syn_tab >> i) & 31 for i in [0, 5, 10]] 
    # Select image of the number 5 from that syndrome,
    # using entry 5 f the digit list ``n_list``
    p[5] = syn[n_list[5]]
    # Let p_list be the list of number 0 <= i < 24 that are not images
    # any of the numbers of [0,1,2,3,4] and not in the syndrome.
    bitmap |= sum(1 << i for i in syn)
    p_list = [i for i in range(24) if ((1 << i) & bitmap) == 0]
    # Select image of the number 8 from that list,
    # using entry 6 f the digit list ``n_list``
    p[8] = p_list[n_list[6]]
    mat24_complete_heptad(p)
    return p


MAT24_ORDER =  244823040 
SH = 58
FACTOR = (1 << SH) // MAT24_ORDER + 1
FACTOR24 = 24 * FACTOR
DFIELD1 = 0x555555555555

def mat24_int_to_perm(n):
    """Return the k-th permutation in the Mathieu group Mat24   

    Here the elements of the Mathieu group Mat24 are numbererd
    in lexicographic order, starting with 0. 

    This is a fast version of function ``py_mat24_int_to_perm``.
    There is also a C version of this function. 
    """
    # See function ``py_mat24_int_to_perm`` for a simpler implementation
    # We document the tricks used to accelerate that implementation.
    assert 0 <= n < MAT24_ORDER
    # Compute the output permutation in ``p``
    p = [None] * 24

    # Compute ``n1`` such that ``k = n1 >> SH`` is the current digit in 
    # the mixed-radix representation of n as described in function    
    # ``py_mat24_int_to_perm``. For obtaining the next digit we subtract
    # k << SH from ``n1`` and then we  multiply ``n1`` by ``i``, 
    # for i = 23, 22, 21, 20, 3, 16 in that order. 
    n = int(n)
    n1 = FACTOR24 * n
    # Let k be the first digit of n. Store this in p[0] 
    p[0] = k = n1 >> SH  
    n1 -= k << SH      # process n1 for next digit
    # We also accumulate the images p[0],...,p[4] in a bitmap
    bitmap = 1 << k

    # For assigning the next images p[m], m = 1,2,3, we take the
    # m-th digit n_m from n and we put p[m] = n_m + d[m]. Here d[m]
    # is a difference that must be added to  n_m in order to avoid
    # repetitions in the images of p[m], 0 <= m < 5.  
    # Computation of d[m] is a bit tricky. We keep the array 
    # d[m], 0 <= m < 24 in the 64-bit integer d as a bit field, 
    # reserving two bits at positions 2*m, 2*m+1 for d[m]. After
    # selecting the image k = n_m, we put p[m] = k + d[k]. Then we
    # have to replace d[u] by d[u + 1] + 1 for u >= k. We don't 
    # change d[u] for u < j. We start with d[u] = 0 for all u. 
    # After putting p[0] = k, we put d[u] = 1 for u >= k.  
    d = DFIELD1 << (2*k)
    for i in [23, 22, 21]:
        # Obtain next digit k of n and adjust n1
        n1 = i * n1
        k = n1 >> SH
        n1 -= k << SH
        # put p[24-i] = j = k  - d[k]
        p[24-i] = j = k + ((d >> (2*k)) & 3)
        # Adjust bit field d: replace d[u] by d[u + 1] + 1 for u >= k.
        # Unfortunately, we would get an overflow in the last round. So
        # in the last round we replace replace d[u] by d[u + 1] instead,
        # and we remember the last value k in the variable ``last``. 
        mask = (1 << (2*k)) - 1
        if i > 21:
            d = (((d + DFIELD1) >> 2) & ~mask) + (d & mask) 
        else:
            d = ((d >> 2) & ~mask) + (d & mask) 
            last = k
        # Enter the image j into the bitmap
        bitmap |= 1 << j 

    # Assign the image p[4]. This works in the same way as assigning
    # the previous imgaes; but here we have to increment d[k] by one
    # in case  k >= last.
    n1 = 20 * n1
    k = n1 >> SH
    n1 -= k << SH
    p[4] = k + ((d >> (2*k)) & 3) + (k >= last)
    bitmap |= 1 << p[4] 

    # Now we have assigned the images p[m], 0 <= m < 5. An we have 
    # computed the bitmap of these images in ``bitmap``.  Next we 
    # compute the syndrome of these five images in ``syn``. Here
    # integer ``syn`` is to be interpreted as a list of three bit
    # fields of length 5. Each field contains an entry of the
    # syndrome. These entries are ordered. 
    cocode = gc.vect_to_cocode(bitmap)
    syn = int(gc.syndrome_table[cocode & 0x7ff])

    # Get next digit k from n. Assign the k-th entry of the
    # syndrome to p[5]
    k = (3 * n1) >> SH
    p[5] = (syn >>  (5 * k)) & 31 
   
    # Enter all entries of the syndrome into the bitmap
    bitmap |= (1 << (syn & 31)) | (1 << ((syn >> 5) & 31))
    bitmap |= (1 << ((syn >> 10) & 31))
    # Complement the bitmap. So the bitmap has 16 bits set.
    bitmap ^= 0xffffff

    # Compute ordered list of positions of the bits set in the bitmap    
    j = 0
    p1 = [None] * 24
    for i in range(24):
        p1[j] = i;              # write index i to output pos.
        j += (bitmap >> i) & 1; # increment output pos. if bitmap[i]=1

    # Let k (where k = n % 16) be the last digit of n. Assign the k-th
    # entry of that list to p[8].
    p[8] = p1[n & 15]

    # Now we have assigned p[m], m = 0,1,2,3,4,5,8; and we procceed
    # as in in function ``py_mat24_int_to_perm``.
    mat24_complete_heptad(p)
    return p








#######################################################################
# Compute the index of an element of the Mathieu group Mat24
#######################################################################




def mat24_perm_to_int(p1):
    """Convert permutation p1 in the Mathieu group Mat24 to an integer.

    This reverses member function mat24_int_to_perm(). The input 
    permutation is not checked.

    The function returns an undefined integer if permutation
    p1 is not in the Mathieu group.
    """
    p1 = [None if x is None else int(x) for x in p1]
    # We reverse to operation of function ``mat24_int_to_perm``.
    # We compute the digits of the result n in mixed-radix
    # representation as described in function ``mat24_int_to_perm``.
    # Digits are taken from the images p1[m], m = 0,1,2,3,4,5,8 in 
    # that order.
    # Get first digit k from p1[0]
    n = k = p1[0]
    # We also accumulate the images p1[0],...,p1[4] in a bitmap
    bitmap = 1 << k

    # For obtaining the next digits n_m from  p1[m], m = 1,2,3, we
    # put  n_m = p1[m] - d[m]. Here d[m] is a difference that must 
    # be subtracted n_m in order to ensure 0 <= n_m < 24-m.
    # Computation of d[m] is a bit tricky. We keep the array 
    # d[m], 0 <= m < 24 in the 64-bit integer d as a bit field, 
    # reserving two bits at positions 2*m, 2*m+1 for d[m]. After
    # selecting the image k = p1[m], we put n_m = k - d[k]. Then we
    # have to replace d[u] by d[u + 1] + 1 for u > k. We don't 
    # change d[u] for u < k. Since a permutation has no repetitions, 
    # we may ignore d[k]. We start with d[u] = 0 for all u. 
    # After putting n_0 = p1[0], we put d[u] = 1 for u >= n_0.  
    d = DFIELD1 << (2*k)
    for i in [23, 22, 21]:
        # Store current entry of p1 in k and adjust bitmap
        k = p1[24-i]
        bitmap |= 1 << k 
        # Set digit n_k = k - d[k] and compute next step for n
        n = i * n + k - ((d >> (2*k)) & 3)
        # Adjust bit field d: replace d[u] by d[u + 1] + 1 for u >= k.
        # Unfortunately, we would get an overflow in the last round. So
        # in the last round we replace replace d[u] by d[u + 1] instead,
        # and we remember the last value k in the variable ``last``.         
        if i > 21:
            d += DFIELD1 << (2*k) 
        else:
            last = k

    # Obtain digit n_4 from p[4]. This works in the same way as
    # assigning the previous imgaes; but here we have to increment
    # d[k] by one in case  k >= last.
    k = p1[4]
    bitmap |= 1 << k 
    n = 20 * n + k - ((d >> (2*k)) & 3) - (k >= last)

    # Now we have obtained the the digis n_m, 0 <= m < 5, from the
    # images p[m], 0 <= m < 5. An we have computed the bitmap of these 
    # images in ``bitmap``.  Next we compute the syndrome of these five
    # images in ``syn``. Here integer ``syn`` is to be interpreted as a 
    # list of three bit fields of length 5. Each field contains an 
    # entry of the syndrome. These entries are ordered. 
    # Also let syn1 = syn[1] be the middle entry of that list.
    cocode = gc.vect_to_cocode(bitmap)
    syn = int(gc.syndrome_table[cocode & 0x7ff])
    syn1 = (syn >> 5) & 31

    # Compute next digit n_5 = k from p[5]. Put k = 2 if p1[5] > syn1,
    # k = 1 if p1[5] == syn1, and  k = 0 if p1[5] < syn1.
    k = int(p1[5] > syn1) + int(p1[5] >= syn1)
    # Enter that digit into the accumulated value n
    n = 3 * n + k

    # Enter all entries of the syndrome into the bitmap
    bitmap |= (1 << (syn & 31)) | (1 << syn1)
    bitmap |= (1 << ((syn >> 10) & 31))

    # Delete all bits at positions >= p1[8] from the bitmap. The the
    # last digit of n ist the bit weight of that bit map.
    d = gc.bw24(((1 << p1[8]) - 1) & bitmap)
    # Enter the last digit into n
    n = 16 * n + p1[8] - d
    # Change n to zero if it is too large
    n & ((n >= MAT24_ORDER) - 1); 
    return n






#######################################################################
# Compute permitatio in the Mathieu group Mat24 from heptads
#######################################################################

           
            
def mat24_perm_from_heptads(h1, h2):
    """Try to find a permutation p that maps heptad h1 to h2

    h1 and h2 must be lists of length 7 defining two umbral heptads,
    i.e. heptads not contained in an octad. If a permutation p in
    the Mathieu group Mat24 that maps h1 to h2 exists, it is unique. 

    Return permutation p if such a p exists an is unique,
    and return None otherwise.
    """
    h1 = [None if x is None else int(x) for x in h1]
    h2 = [None if x is None else int(x) for x in h2]
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
#######################################################################
# Support for the code generator
#######################################################################
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
#######################################################################
# Testing this module
#######################################################################
#######################################################################


def test_complete_octad():
    io = lrange(6) + [None,None,8] + [None]*15    
    complete_heptad(io)
    assert io == lrange(24)   
    print("The python function ``complete_heptad()`` is correct")

if __name__ == "__main__":
    hint_for_complete_heptads()
    test_complete_octad()
    
    