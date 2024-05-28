# This file has been created automatically, do not change!!!

r"""Classes of square roots of Q_x0 in G_x0.

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

INVOLUTION_SAMPLES = [
# Characters for Co_1 class 1A
[[(1, 0), (196883, 299, 24, 4096, 98280), (0, 0, 0, 0, 0)], 'M0<1>'],
[[(2, 0), (275, 299, 24, -4096, 98280), (0, 0, 0, 0, 0)], 'M0<x_1000h>'],
[[(2, 0), (4371, 299, 24, 0, 4072), (1, 4, 0, 0, 2)], 'M0<x_80fh>'],
[[(4, 0), (275, 299, 24, 0, -24), (1, 5, 0, 0, 3)], 'M0<x_0ae0h>'],
[[(2, 0), (275, 299, 24, 0, -24), (1, 4, 0, 0, 4)], 'M0<x_800h>'],

# Characters for Co_1 class 2A
[[(2, 0), (4371, 43, 8, 256, 2280), (8, 0, 0, 0, 0)], 'M0<y_80fh>'],
[[(2, 0), (275, 43, 8, -256, 2280), (8, 0, 0, 0, 0)], 'M0<y_0fh>'],
[[(4, 4), (275, 43, 8, 0, 232), (9, 4, 1, 4, 16)], 'M0<y_80fh*d_1h>'],
[[(2, 0), (275, 43, 8, 0, 232), (9, 4, 0, 0, 16)], 'M0<y_80fh*d_3h>'],
[[(4, 2), (51, 43, 8, 0, 8), (9, 4, 1, 2, 2)], 'M0<y_80fh*d_801h>'],
[[(4, 4), (19, 43, 8, 0, -24), (9, 4, 1, 4, 0)], 'M0<y_80fh*d_803h>'],
[[(4, 0), (19, 43, 8, 0, -24), (9, 5, 0, 0, 0)], 'M0<y_0fh*x_1003h*d_406h>'],

# Characters for Co_1 class 2C
[[(4, 4), (275, 11, 0, 64, 264), (12, 2, 0, 0, 0)], 'M0<y_0ae0h>'],
[[(2, 0), (275, 11, 0, 0, 264), (12, 3, 0, 4, 0)], 'M0<y_0ae0h*d_20h>'],
[[(8, 3), (11, 11, 0, 0, 0), (12, 2, 1, 3, 132)], 'M0<y_0ae0h*d_800h>'],
[[(4, 4), (19, 11, 0, 0, 8), (12, 2, 1, 4, 136)], 'M0<y_2e0h*x_1001h*d_0eh>'],
[[(4, 4), (-13, 11, 0, 0, -24), (12, 2, 1, 4, 120)], 'M0<y_2e0h*x_1011h*d_21h>'],
[[(4, 2), (51, 11, 0, 0, 40), (12, 2, 1, 2, 152)], 'M0<y_2e0h*x_1020h*d_0fh>'],

# Characters for Co_1 class 2B
[[(4, 0), (-13, -13, 0, 64, 0), (12, 0, 0, 0, 0)], 'M0<y_9d0h*d_700h*p_21289010>'],
[[(4, 4), (-13, -13, 0, 0, 0), (12, 1, 0, 4, 0)], 'M0<y_9d0h*d_791h*p_21289010>'],

]
