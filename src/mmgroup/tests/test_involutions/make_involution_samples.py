import sys
import time
from math import floor
from random import randint, shuffle, sample
from collections import defaultdict
import numpy as np

sys.path.append(r"C:\Data\projects\MonsterGit\src")

from mmgroup import MM, AutPL, PLoop, Cocode, Xsp2_Co1
from mmgroup.generators import gen_leech2_type
from mmgroup.clifford12 import xsp2co1_involution_invariants
from mmgroup.clifford12 import xsp2co1_orthogonal_complement


G = MM

def transversal_space(g):
    a = Xsp2_Co1(g)._involution_invariants()[0]
    while len(a) and int(a[-1]) == 0:
        a = a[:-1] 
    while len(a) and  int(a[0]) & 0x7000000:
        a = a[1:]
    a = [int(x) & 0xffffff for x in a]
    bl = set(range(1,25)) - set([x.bit_length() for x in a])
    assert len(bl) == 24 - len(a)
    return [1 << (x - 1) for x in bl]
 


def iter_Q_x0():
    yield G()
    yield G(("x", 0x1000))
    yield G(("x", PLoop(range(8))))
    yield G(("x", PLoop([0,4,8,12,16,20,2,7,10,])))
    yield G(("x", 0x800))

def iter_transversal(g):
    if g == G():
        yield from iter_Q_x0()
        return
    b = transversal_space(g)
    yield g * G(("q", 0x1000000))
    c = [0]
    for v in b:
        c += [v ^ x for x in c]
    for x in c:
        yield g * G(("q", x))     




y8 = G(("y", PLoop(range(8))))

y12 =  G(("y", PLoop([0,4,8,12,16,20,2,7,10,])))


#from fourv import find_fourvolution
#g4 = find_fourvolution()
#Possible value for g4:
g4 = G("M<y_3afh*d_480h*p_41443>")
assert g4.chi_G_x0()[0] == -13 and g4.order() == 4




def invariant_count_type2(iv):
    v0  =  (int(iv[0]) >> 24) & 7 if len(iv) else 0
    while len(iv) and int(iv[-1]) == 0:
        iv = iv[:-1] 
    if (v0 & 4 != 4):
        return 0
    data = [int(iv[0]) & 0xffffff]
    for i in range(1, len(iv)):
        if (int(iv[i]) >> 24) & 1 == 0:
            x = int(iv[i]) & 0xffffff
            data += [x ^ y for y in data]
    s = sum(((gen_leech2_type(x) & 0xf0) == 0x20 for x in data))
    return s
    

def invariant_type(g):
    iv, _1, ortho = Xsp2_Co1(g)._involution_invariants()
    iv_len = 12
    while iv_len and iv[iv_len -1] == 0:
        iv_len = iv_len - 1
    v0, v1 =  (int(iv[0]) >> 24) & 7, (int(iv[1]) >> 24) & 7
    ct2 = invariant_count_type2(iv)
    return iv_len, v0, v1, gen_leech2_type(ortho) >> 4, ct2


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



#display_involution_invariants(MM("M<y_0fh*x_1003h*d_5h>"))
#1/0 


   

def show_characters(g0, file):
    characters = []
    for nn, g in enumerate(iter_transversal(g0)):
        o = (g.order(), (g**2).type_Q_x0())
        o_chi =  o, g.chi_G_x0(), invariant_type(g)
        if not o_chi in characters:
            characters.append(o_chi)
            chi = o_chi[1] 
            x = chi[0] - chi[1] - chi[2] * chi[3]
            print([list(o_chi), x, str(g)], ",", file = file)
            if o_chi[2][:3] == (4,1,4):
                #display_involution_invariants(g)
                pass
        #if nn & 0xffff == 0: print(".")



DOC = """Classes of square roots of Q_x0 in G_x0.

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


The list of invariant contains the following three tuples:

Tuple 1:
   - Order of the class
   - Type Leech lattice mod 2 vector corresponding to squared class

Tuple2:
   - Character in the 196883-dimensional rep of the monster
   - Character in the 299-dimensional 299_x rep of Co_1
   - Character in the 24-dimensional rep 24-x of Co_0
   - Character in the 4096-dimensional rep 4096_x

   The last two reps are defined up to sign only, but the tensor
   product of these two reps is well defined as a rep of G_x0.

Tuple 3:
   - Minumum of 12 and the dimension of  (\im (A - 1))^+
   - Row 0, column bits 26, 25, 24  of invariant matrix ``invar``
   - Row 1, column bits 26, 25, 24  of invariant matrix ``invar``
   - Type of s^\perp in Leech lattice mod 2
   - Number of type-2 vectors in (\im (A - 1))^-; this is
     calculated for classes 1A and 2A  in Co_1 only.
     
   Column bit 26 in row 0 is one iff (\im (A - 1))^- is strictly 
   greater than (\im (A - 1)); this is calculated for classes 1A 
   and 2A  in Co_1 only.

   Column bit 25 in row 0 is 1 for class 2B in Co_1 and 0 for 
   all other classes in Co_1. It is 1 iff \im (A - 1) contains
   type-3 vectors.

   Column bits 24 are set to one in the following cases only:

   Class 1A in Co_1:  
      In row 0 iff g corresponds to a type-3 vector.

   Class 2A in Co_1:
      In row 0 iff s is 0 and (\im (A - 1))^- contains a type-3 vector. 
      Then all vectors in (\im (A - 1))^- are of type 3.

      In row 1 iff s is not 0.

   Class 2B in Co_1:  
      In row 0 iff s is equal to the linear form x -> type(x) (mod 2).

      In row 1 iff s is not 0 and bit 24 in row 0 is not set to 1.

   Class 2C in Co_1:
      In row 0 iff s is not 0. 


"""



def print_invariants(file = None):
    is_py = not file is None
    if file is None:
        file = sys.stdout
    do_open =  isinstance(file, str)
    if do_open:
        file = open(file, "wt")
    if is_py:
        s = "This file has been created automatically, do not change!!!"
        print("# " + s + "\n", file = file)
        print('"""', end = "", file = file)
    print(DOC, file = file)
    if is_py:
        print('"""\n\nINVOLUTION_SAMPLES = [',file = file)
    print("#Characters for Co_1 class 1A", file = file)
    show_characters(G(), file = file)
    print("\n#Characters for Co_1 class 2A", file = file)
    show_characters(y8, file = file)
    print("\n#Characters for Co_1 class 2B", file = file)
    show_characters(y12, file = file)
    print("\n#Characters for Co_1 class 2C", file = file)
    show_characters(g4, file = file)
    if is_py:
        print(']', file = file)
    if do_open:
        file.close()




"""Output:
Possible orders, characters, and involution invariants of x*Q_x0

Characters for Q_x0, x has Co_1 class 1A
[1, (196883, 299, 24,  4096), (0, 0, 0,   0)], 98280, M<1> 
[2, (   275, 299, 24, -4096), (0, 0, 0,   0)], 98280, M<x_1000h> 
[2, (  4371, 299, 24,     0), (4, 0, 0,   1)], 4072, M<x_80fh> 
[4, (   275, 299, 24,     0), (5, 0, 0,   0)], -24, M<x_0ae0h> 
[2, (   275, 299, 24,     0), (4, 0, 0,   0)], -24, M<x_800h> 

Characters for y_octad, x has  Co_1 class 2A
Transversal dimension =  16
[2, (4371,  43, 8,  256), (0, 0, 0,  0)], 2280,  M<y_80fh> 
[4, ( 275,  43, 8,    0), (4, 1, 4, 16)],  232,  M<y_80fh*d_1h> 
[2, ( 275,  43, 8,    0), (4, 0, 0, 16)],  232,  M<y_80fh*d_3h> 
[4, (  51,  43, 8,    0), (4, 1, 2,  2)],    8,  M<y_80fh*d_801h> 
[4, (  19,  43, 8,    0), (4, 1, 4,  0)],  -24,  M<y_80fh*d_803h> 
[4, (  19,  43, 8,    0), (5, 0, 0,  0)],  -24,  M<y_0fh*x_1003h*d_5h> 
[2, ( 275,  43, 8, -256), (0, 0, 0,  0)], 2280,  M<y_100fh> 

Characters for y_dodecad, Co_1 class 2B
Transversal dimension =  13
[4, ( 275,  11, 0,   64), (2, 0, 0, 0)], 264, M<y_0ae0h> 
[4, (  51,  11, 0,    0), (2, 1, 2, 0)],  40, M<y_2e0h*x_1400h> 
[4, (  19,  11, 0,    0), (2, 1, 4, 0)],   8, M<y_2e0h*x_1600h*d_501h> 
[4, ( -13,  11, 0,    0), (2, 1, 4, 0)], -24, M<y_2e0h*x_1700h*d_701h> 
[8, (  11,  11, 0,    0), (2, 1, 3, 0)],   0, M<y_0ae0h*d_800h> 
[2, ( 275,  11, 0,    0), (3, 0, 4, 0)], 264, M<y_0ae0h*d_200h> 
 
Characters for Co_1 class 2C
Transversal dimension =  12
[4, ( -13, -13, 0,   64), (0, 0, 0, 0)], 0, M<y_3afh*d_480h*p_41443>]
[4, ( -13, -13, 0,    0), (1, 0, 4, 0)], 0, M<y_3afh*x_0ch*d_2dh*p_41443>

"""


if __name__ == "__main__":
    print_invariants(file = "involution_samples.py")


