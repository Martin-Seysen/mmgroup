r"""Compute classes of square roots of Q_x0 in G_x0.

Running this module creates a python file ``involution_samples.py``
containing samples of these classes including documentation.

Let :math:`Q_{x0}` be the  extraspecial normal subgroup of
structure :math:`2^{1+24}` of the group  math:`G_{x0}`. Then
:math:`G_{x0} / Q_{x0} \cong \mbox{Co}_1`,
and there is a natural homomrphism :math:`\phi` from  
:math:`Q_{x0}` to :math:`\Lambda / 2 \Lambda`,
i.e. to the Leech lattice modulo 2.  The kernel of
:math:`\phi` is :math:`Z(Q_{x0}) = \{x_{\pm 1}\}`. So the
group `\mbox{Co}_1` operates on :math:`Q_{x0} / Z(Q_{x0})`
by conjugation, an this operation is compatible with the
natural action of `\mbox{Co}_1` on :math:`\Lambda / 2 \Lambda`.

Ideally, we would like to separate all classes of square roots 
of :math:`Q_{x0}` in :math:`G_{x0}`, e.g. by computing the 
characters of the representations
:math:`299_x, 98280_x, 24_x \otimes 4096_x` of :math:`G_{x0}`.
But actually, we just separate the classes of involutions
in  :math:`G_{x0}`; and we generate a set of elements of
:math:`G_{x0}` runnning through all classes of square roots 
of :math:`Q_{x0}`.

According to [Wil13] there are five classes of involutions in
the group :math:`G_{x0}` that map to 2B involutions in the
monster. According to Table 2 in [Nor98] there are two classes 
of involutions in the group :math:`G_{x0}` that map to 
2A involutions in the monster. In this module we compute the
abovementioned characters of many involutions in :math:`G_{x0}`,
and we will see that these characters enable us to separate the 
seven classes of involutions in  :math:`G_{x0}` .

We store samples of different classes elements of
:math:`G_{x0}` that square up to elements of :math:`G_{x0}`,
(including the above characters of these samples) in file 
``involution_samples.py``. From this information we conclude 
that we can separate the seven classes of involutions 
in :math:`G_{x0}` by computing characters.

There are three classes of involutions in the factor group
group :math:`\mbox{Co}_1` of :math:`G_{x0}`, which are
labelled 2A, 2B, and 2C in the ATLAS :cite:`CCN+85`. 
Thus for any involution :math:`y \in G_{x0}` the coset
:math:`y Q_{x0}` must either be equal to :math:`Q_{x0}`
or lie in one of these classes of  :math:`\mbox{Co}_1`.

So it suffices to check the cosets :math:`g Q_{x0}` for
involutions, where :math:`g` runs through the neutral
element of :math:`G_{x0}`, and through one preimage in
:math:`G_{x0}` of the three classes of involutions
in :math:`\mbox{Co}_1`.

Running through all :math:`2^{25}` elements of a
coset :math:`g Q_{x0}` is still a formidable task. In
case :math:`g = 1` it suffices to inspect the
elements :math:`\{x_{\pm 1}\}`, and three other elements
in :math:`h \in Q_{x0}` such the the type of 
in :math:`\phi(h)` in :math:`\Lambda / 2 \Lambda` runs 
through the three possible nonzero values 2, 3, and 4. 

We deal with the preimages :math:`g` of involutions of
:math:`\mbox{Co}_1` as follows. 
For any :math:`h \in \mbox{Co}_1` let :math:`M(h)` be the 
matrix representing :math:`h` in the natural 24-dimensional
representation :math:`R` of :math:`\mbox{Co}_1` over 
:math:`\mathbb{F}_2` acting on :math:`\Lambda / 2 \Lambda`.
For :math:`h \in \mbox{Co}_1` let :math:`\imath(h)` be the
image of the matrix :math:`M(h) - 1`. For 
:math:`g \in G_{x0}` we put :math:`\imath(g) = \imath(h)`,
where :math:`h = g Q_{x0} \in \mbox{Co}_1`.

It is not difficult to show that for any :math:`g \in G_{x0}`
and :math:`v \in \imath(g))` there is an 
:math:`x \in \phi^{-1}(v) \subset Q_{x0}`,
such that  :math:`g x` is in the same class of :math:`G_{x0}` 
as :math:`g`. The proof of this statement relies on the
fact (mentioned above) that the operation of `\mbox{Co}_1` 
on :math:`Q_{x0} / Z(Q_{x0})` by conjugation is compatible 
with the action of `\mbox{Co}_1` on :math:`\Lambda / 2 \Lambda`.

Thus for finding representatives of all classes 
of :math:`G_{x0}` in the coset  :math:`g Q_{x0}` it suffices
to check :math:`g x` all :math:`x \in \phi^{-1}(L)`, where
:math:`L` is a maximal subspace of  
:math:`\Lambda / 2 \Lambda` with 
:math:`L \cap \imath(g) = \{0\}`.

If :math:`h \in \mbox{Co}_1` is of class  2A, 2B, or 2C then
:math:`\imath(h)` has dimension 8, 12, or 12, respectively.
So we have to check :math:`2^{17}, 2^{13}, 2^{13}` entries
in these three cases. Thus these checks are not too tedious; 
and they are done in this module.

We store some more information for all computed samples
(depending on the class of the sample only) in file
``involution_samples.py``. This information is used to
accelerate the computation of the character of the
representation  :math:`98280_x` in some cases, where that
computation would take a long time otherwise.
"""

import sys
import time
from math import floor
from random import randint, shuffle, sample
from collections import defaultdict
from time import sleep
import numpy as np

if __name__ == "__main__":
    sys.path.append("../../../")

from mmgroup import MM0, GCode, AutPL, PLoop, Cocode, Xsp2_Co1, XLeech2
from mmgroup.generators import gen_leech2_type
from mmgroup.clifford12 import xsp2co1_leech2_count_type2
from mmgroup.clifford12 import xsp2co1_traces_all


G = Xsp2_Co1


def std_hexads(*i):
    """Return a certain umbral hexad as a list h[i].

    Here ``0 <= i < 4 must hold``.
    For i = 0  return ``h[0] = [0, 4, 8, 12, 16, 20, 24]``. Otherwise 
    return ``h[i] = [4*j + 1 + (l[j + i]) mod 3], 0 <= i < 6``;
    where ``l[j] = [0, 1, 0, 1, 0, 1]``.
    
    If several parameters ``i1``, ``i2``, ... are given then the 
    corresponding lists  ``h[i1]``, ``[i2]``, are concatenated.

    Function ``check_std_hexads`` checks that the lists 
    ``std_hexads(0, i)`` are dodecads for 1 <= i <= 3.
    """
    if len(i) != 1:
        return sum([std_hexads(j) for j in i], [])
    i = i[0] & 3
    if i == 0:
        return list(range(0, 24, 4))
    l = [0, 1] * 3
    return [4*j + 1 + (l[j] + i) % 3 for j in range(6)]


def check_std_hexads(verbose = 0):
    """Check that some lists returned by ``std_hexads()`` are dodecads.

    Assertions in this function fail if this is not the case.
    See function ``std_hexads()`` for details.  
    """
    if verbose:
        print("Standard dodecads")
    for i in range(1, 4):
        l = std_hexads(0, i)
        assert set(GCode(l).bit_list) == set(l)
        if verbose:
            print(std_hexads(0,i))

check_std_hexads(verbose = 0)



def transversal_basis(g):
    r"""Compute basis of a transversal for :math:`g \in G_{x0}`

    Let ``g`` be an element of the group  :math:`G_{x0}`, such
    that ``Xsp2_Co1(g)`` constructs an instance of class
    `Xsp2_Co1``.

    Let :math:`\imath(g)` be the subspace of 
    :math:`\Lambda / 2 \Lambda` defined in the header of
    this module, and let :math:`L` be a complement of
    :math:`\imath(g)`, ad defined in ibid. Then the function 
    returns basis of :math:`\imath(g)` as a list of
    vectors in **Leech lattice encoding**.
    """
    a = Xsp2_Co1(g)._involution_invariants()[0]
    while len(a) and int(a[-1]) == 0:
        a = a[:-1] 
    while len(a) and  int(a[0]) & 0x4000000:
        a = a[1:]
    a = [int(x) & 0xffffff for x in a]
    bl = set(range(1,25)) - set([x.bit_length() for x in a])
    assert len(bl) == 24 - len(a)
    return [1 << (x - 1) for x in bl]

 

def vector_space(basis):
    r"""Map basis of vector space to vector space.

    Given a ``basis`` of a subspace of :math:`\Lambda / 2 \Lambda`
    in **Leech lattice encoding**, the function returns the list
    of all vectors spannned by that basis.
    """
    c = [0]
    for v in basis:
        c += [v ^ x for x in c]
    return c


def iter_Q_x0():
    r"""Iterator for all samples in :math:`Q_{x0}` to be checked

    This function performs the case ``g = 1`` in function
    ``iter_transversal```.
    """
    yield G()
    yield G("x", 0x1000)
    yield G("x", PLoop(range(8)))
    yield G("x", PLoop([0,4,8,12,16,20,2,7,10,]))
    yield G("x", 0x800)


def iter_transversal(g):
    r"""Iterator for all samples in :math:`g Q_{x0}` to be checked

    Let ``g`` be an element of the group :math:`G_{x0}` such
    that the coset :math:`g Q_{x0}` in  :math:`{Co}_1` has 
    order 1 or 2. Then the function iterates over the subset
    of :math:`g Q_{x0}` containing representatives of all 
    classes of :math:`G_{x0}` occuring in :math:`g Q_{x0}`,
    as described in the header of this module.
    """
    if g == G():
        yield from iter_Q_x0()
        return
    b = transversal_basis(g)
    for x in vector_space(b):
        yield g * G("q", x)    
        yield g * G("q", x) * G("x", 0x1000)   





# This is a representative of class 2A in :math:`\mbox{Co}_1`
y8 = G("y", PLoop(range(8)))

# This is a representative of class 2C in :math:`\mbox{Co}_1`
y12 =  G("y", PLoop(std_hexads(0, 1)))

# This is a the central involution in :math:`Q_{x0}`
neg =G('x', 0x1000)





def check_y12_involution_conjugates_to_its_negative():
    r"""Check involutions corresponding to class 2C 

    Let :math:`g \in G_{x0}` be such that :math:`g Q_{x0}`
    is of class 2C in :math:`\mbox{Co}_1`. We may e.g. put
    ``g = y12``.

    Then we check that all involutions in the coset
    :math:`g Q_{x0}` of :math:`G_{x0}` belong to the same 
    class of :math:`G_{x0}`, and that any such involution
    is in the same class as its nagative. (Here negation
    means multiplication with the central involution
    of :math:`G_{x0}`).
    """
    n_total = 0
    n_involution_transversals = 0
    n_good_tansversals = 0
    for g in iter_transversal(y12):
        n_total += 1
        #print(g.order())
        if g.order() == 2:
            n_involution_transversals += 1
            good = False
            #print(y12, g)
            g_neg = g * neg
            for i in range(1000000):
                c = G('x', 'r') * G('d', 'r') 
                g_conj = g ** c
                if g_conj == g_neg:
                    #print("YEEAH", i)
                    good = True
                    break
            if good: n_good_tansversals += 1
    assert  n_involution_transversals == n_good_tansversals  == 2
    assert n_total == 2**13, n_total


check_y12_involution_conjugates_to_its_negative()




def is_nice_permutation(pi):
    """Check that permutation centralizes the standard tetrad"""
    for i in range(0, 24, 4):
        k = pi[i] & -4
        if (k & -4 != i & -4): return False
        for j in range(i+1, i+4):
            if pi[j] & -4  != k: 
                return False
    return True


def get_perm(m):
    """Get atom with 'p' tag from monster elemment as permutation"""
    for tag, n in m.as_tuples():
        if tag == 'p':
             return AutPL(0, n).perm
    raise ValueError("Permutation in monster element not found")


def characters(g):
    a = np.zeros(4, dtype = np.int32)
    assert xsp2co1_traces_all(Xsp2_Co1(g)._data, a) == 0
    chi24, chisq24, chi4096, chi98280  = map(int, a[:4])
    chi299 = (chi24**2 + chisq24) // 2 - 1
    chi_M = chi299 + chi98280  + chi24 * chi4096
    return chi_M, chi299, chi24, chi4096, chi98280 
 



STD_FOURVOLUTION = None


def find_a_fourvolution(n_trials):
    """Try to find a fourvolution 

    See function ``find_fourvolution`` for the definition of a
    fourvolution. The function ties ``n_trials`` times.
    """
    for i in range(n_trials):
        e = y12 *  G('p', 'r') # G([("y",y12), ("p","r")])
        order = e.order()
        if  order % 4 == 0:
            v = (e ** (order//4))
            if v * v == neg:
                chi = characters(v)
                if chi[0] == -13 and is_nice_permutation(get_perm(v)):
                    return v
    return None


def find_fourvolution(verbose = 0):
    r"""Return a certain 'fourvolution' in G_x0 as an element of MM

    Here a 'fourvolution' is an element v of G_x0 auch that the
    the image v Q_x0 of v in Co_1 is in class 2B of Co_1. Then 
    any preimage of v Q_x0 in Co_0 squares to -1. The returned 
    element of v has the following properties:

       vQ    is in   y12 * pi * Q,

    where y12 is an element y_d of N_x0 with, where d is a Golay 
    code word that intersects with each tetrad of the standard 
    sextet in a set of two elements.
    Here pi is a permutation that centralizes the standard sextet.
    The standard sextet contains the tetrad [0,1,2,3].
    """
    global STD_FOURVOLUTION 
    if STD_FOURVOLUTION is not None:
        return STD_FOURVOLUTION  
    print("Searching for a 'fourvolution' in G_x0...")
    N_TRIALS, N_SINGLE = 1000000, 1000
    try:
        from multiprocessing import Pool, cpu_count
        n_cpus = max(1, min(32, cpu_count()))
        n_list = [N_SINGLE] * n_cpus
        for i in range(N_TRIALS // N_SINGLE + 1): 
            with Pool(processes = n_cpus) as pool:
                results = pool.map(find_a_fourvolution, n_list)
            pool.join()
            results = [x for x in results if x is not None]
            if len(results):
                v = results[0]
                #print("Found %s after %d rounds", v, i)
                break
    except:
        print("Parallel computation of fourvolution failed!")
        v = find_a_fourvolution(N_TRIALS)
        if v is None:
            raise ValueError("No suitable fourvolution found")
    print("found" , v)
    if verbose:
         print("y part")
         print(GCode(v.mmdata[0]).bit_list)
         print("permutation part")
         print(AutPL(0, v.mmdata[2] & 0xfffffff).perm)
    STD_FOURVOLUTION  = v
    return v


#find_fourvolution()


def invariant_count_type2(iv):
    v0  =  (int(iv[0]) >> 24) & 7 if len(iv) else 0
    while len(iv) and int(iv[-1]) == 0:
        iv = iv[:-1] 
    if len(iv) == 1:
        return gen_leech2_type(iv[0])
    if (v0 & 4 == 4):
         return xsp2co1_leech2_count_type2(np.copy(iv), len(iv))
    if len(iv) == 12:
        v1  =  (int(iv[1]) >> 24) & 7 
        if v0 == 2 and v1 == 1:
            cp = np.copy(iv[1:])
            cp[0] = 0
            return xsp2co1_leech2_count_type2(cp, 11)
    return 0
    

def invariant_type(g):
    iv, _1, ortho = Xsp2_Co1(g)._involution_invariants()
    iv_len = 12
    while iv_len and iv[iv_len -1] == 0:
        iv_len = iv_len - 1
    v0, v1 =  (int(iv[0]) >> 24) & 7, (int(iv[1]) >> 24) & 7
    ct2 = invariant_count_type2(iv)
    return iv_len, v0, v1, gen_leech2_type(ortho), ct2


def display_involution_invariants(g):
    iv, v1, v0 = Xsp2_Co1(g)._involution_invariants()
    print("Involution invariants of", MM(g), ":")
    for i in range(12):
        v = int(iv[i])
        if i and v == 0:
            break
        print("0x%06x 0x%01x 0x%06x" % 
            ((v >> 32) & 0xffffff, (v >> 24) & 0xf, v & 0xffffff)
        )
    print("Orthogonal: 0x%06x, 0x%06x\n" % (v1, v0))



#display_involution_invariants(MM("M0<y_0fh*x_1003h*d_5h>"))
#1/0 


   

def show_characters(g0, class_name):
    output = ["# Characters for Co_1 class %s" % class_name]
    if g0 == G():
        transversal = iter_Q_x0()
    else:
        transversal = iter_transversal(g0)
    all_characters = []
    for nn, g in enumerate(transversal):
        o = g.order(), XLeech2(g**2).type
        o_chi =  o, characters(g), invariant_type(g)
        if not o_chi in all_characters:
            all_characters.append(o_chi)
            chi = o_chi[1] 
            x = chi[0] - chi[1] - chi[2] * chi[3]
            output.append(str([list(o_chi), str(MM0(g))]) + ",")
        #if nn & 0xffff == 0: print(".")
    return "\n".join(output) + "\n"



DOC = r"""Classes of square roots of Q_x0 in G_x0.

The list INVOLUTION_SAMPLES shows invariants and representatives g of
all classes in the subgroup G_x0 of the monster that square up to an 
element of Q_x0. Let A be the matrix corresponding to the operation 
of g on the Leech lattice mod 2. Then g operates trivially on the
preimage of \ker (A - 1) in Q_x0 up to sign. Here the negated
elements of  Q_x0 are given by a linear form s on \ker (A - 1).
The linear form s is also defined in the subspace \im (A - 1) of
\ker (A - 1).

There is a nonsingular symmetric bilinear form B on \im (A - 1)
called the Wall parametrization. The orthogonal complement s^\perp
of s (with regard to B) is a vector in \im (A - 1).

Let (\ker (A - 1))^+ be the subspace of the kernel where s vanishes.
Let (\im (A - 1))^+ be the orthogonal complement of (\ker (A - 1))^+,
and let (\im (A - 1))^- be the affine space 
(\im (A - 1))^+  \setminus \im (A - 1).


For each class of preimages in G_x0 of an involution in Co_1
we store a list of invariants of that class.

After the list, the character of the class in the rep 98280_x 
of G_x0 follows.

After that, a representative of g of the class is printed.


The list of invariants contains the following three tuples,
and an appended character value:

Tuple 1:
   - Order of the class
   - Type Leech lattice mod 2 vector corresponding to squared class

Tuple2:
   - Character in the 196883-dimensional rep of the monster
   - Character in the 299-dimensional 299_x rep of Co_1
   - Character in the 24-dimensional rep 24_x of Co_0
   - Character in the 4096-dimensional rep 4096_x
   - Character in the 98280-dimensional monomial rep 98280_x

   Representations 24_x and 4096_x are defined up to sign only, 
   but the tensor product of these two representations is well 
   defined. We normalize a pair of characters (chi_24, chi_4096)
   so that the first onzero entry of that pair is positive.

Tuple 3:
   - Minimum of 12 and the dimension of  (\im (A - 1))^+
   - Row 0, column bits 26, 25, 24  of invariant matrix ``invar``
   - Row 1, column bits 26, 25, 24  of invariant matrix ``invar``
   - Type of s^\perp in Leech lattice mod 2
   - Related to the number of type-2 vectors in \im (A - 1)
     For class 1A  in Co_1
         Type of vector in (\im (A - 1))^- if present, otherwise 0
     For class 2A  in Co_1
         Number of type-2 vectors in (\im (A - 1))^-
     For class 2C  in Co_1
         Number of type-2 vectors in (\im (A - 1))^+
         This number is not calculated and set to zero if
         (\im (A - 1))^- contains no vectors of type 2 or 4.
     For all other classes:
         Not calculated, and set to zero
     
   Column bit 26 in row 0 is one iff (\im (A - 1))^+ is strictly 
   greater than (\im (A - 1)); this is calculated for classes 1A 
   and 2A  in Co_1 only.

   Column bit 25 in row 0 is 1 for class 2C in Co_1 and 0 for 
   all other classes in Co_1. It is 1 iff \im (A - 1) contains
   any type-3 vectors.

   Column bits 24 are set to one in the following cases only:

   Class 1A in Co_1:  
      In row 0 iff g corresponds to a type-3 vector.

   Class 2A in Co_1:
      In row 0 iff s is 0 and (\im (A - 1))^- contains a type-3 vector. 
      Then all vectors in (\im (A - 1))^- are of type 3.

      In row 1 iff s is not 0.

   Class 2C in Co_1:  
      In row 0 iff s is equal to the linear form x -> type(x) (mod 2).

      In row 1 iff s is not 0 and bit 24 in row 0 is not set to 1.

   Class 2B in Co_1:
      In row 0 iff s is not 0. 


"""


def get_character_data():
    t_start = time.time()
    g4 = find_fourvolution()
    char_samples = [
            (G(), "1A"), (y8, "2A"), (y12, "2C"), (g4, "2B")
    ]
    print("Searching for samples of involutions in subgroup G_x0")
    try:
        from multiprocessing import Pool
        with Pool(processes = 4) as pool:
            data = pool.starmap(show_characters, char_samples)
        pool.join()
    except:
        print("Parallel computation of data has failed!")
        data = [show_characters(g, name) for g, name in char_samples]
    t = time.time() - t_start
    print("Involution data computed in %.2f seconds" % t)
    return data

def print_invariants(file = None):
    data = get_character_data()
    is_py = not file is None
    if file is None:
        file = sys.stdout
    do_open =  isinstance(file, str)
    if do_open:
        file = open(file, "wt")
    if is_py:
        s = "This file has been created automatically, do not change!!!"
        print("# " + s + "\n", file = file)
        print('r"""', end = "", file = file)
    print(DOC, file = file)
    if is_py:
        print('"""\n\nINVOLUTION_SAMPLES = [',file = file)
    for record in data:
        print(record,  file = file)
    if is_py:
        print(']', file = file)
    if do_open:
        file.close()
        sleep(0.01)





if __name__ == "__main__":
    print_invariants(file = "involution_samples.py")


