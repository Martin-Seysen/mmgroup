r"""
Computation in the real Clifford group
--------------------------------------

Backround from quantum computing
................................

Let :math:`V = \mathbb{F}_2^n` be a Boolean vector space and 
:math:`\mathbb{T}` be the unit circle in 
the set :math:`\mathbb{C}` of the complex numbers. Then 
:math:`\mathbb{F}_2^n` is an additive and :math:`\mathbb{T}` 
is a multiplicative Abelian group. For any mapping 
:math:`f: V \rightarrow \mathbb{T}` define 
the mapping :math:`\Delta f: V \times V ->  \mathbb{T}` by

.. math::

   \Delta f(x,y) = f(x+y) \cdot f(x)^{-1} \cdot f(y)^{-1} * f(0) .

A mapping :math:`g: V \rightarrow \mathbb{T}`  is bilinear if 

.. math::

   g(x_1 + x_2, x_3) = g(x_1, x_3) \cdot g(x_2, x_3) ,  \\
   g(x_1, x_2 + x_3) = g(x_1, x_3) \cdot g(x_1, x_3) .

Then we obviously have :math:`g(x_1,x_2) = \pm 1` and 
:math:`g(0,x_2) = g(x_1,0) = 1`.
A function :math:`f: V \rightarrow \mathbb{T}` is called a 
*quadratic* if :math:`\Delta f` is 
bilinear. For a quadratic function :math:`f`  we also require 
the functional values  :math:`f(x)` to be 8-th roots of unity.


Let :math:`Q: F_2^m \rightarrow  \mathbb{T}` be  a quadratic function. 
Let :math:`A: F_2^m \rightarrow  F_2^n` be an affine mapping. 
Let :math:`e` be an integer. Then we define a function 
:math:`f = f(e, A, Q): F_2^n \rightarrow \mathbb{C}` by

.. math::
   :label: quadratic_mapping

   f(x) = 2^{e/2} \cdot \sum_{y \in F_2^m: A(y) = x} Q(y) .  

We call a function :math:`f`  statisfying  :eq:`quadratic_mapping`
a *quadratic mapping*.


By definition, quadratic mappings are invariant under multiplication.
A function :math:`f: F_2 \rightarrow \mathbb{C}` hs a natural
interprestion as a vector in :math:`\mathbb{C}^2`. Similarly, a
function :math:`f: F_2^n \rightarrow \mathbb{C}` has a natural
interpreation as a tensor in the tensor product 
:math:`(\mathbb{C}^2)^{\otimes n}`. With this interpretation,
tensor products and tensor contractions of quadratic function
are quadratic functions by definition.

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

In the theory of quantum computation, a state vector representing the 
state of n qubits can be written as a vector in CC**(2**n), where the 
2**n unit vectors of CC**n are labeled by the elements of F_2**n. Such 
a state vector v can be considered as a function  F_2**n -> CC mapping 
the labels of the unit vectors to the corresponding coordinates of the 
vector v. We call this function the coordinate function of the 
vector v. Coordinate functions of a complex 2**n \times 2**n matrix 
(describing a transition from one n-qubit state to another) are 
defined similarly. In the theory of quantum computing the vector space
CC**(2**n) is considered as an an n-fold tensor product of CC*2.


A vector (or a matrix) is called a quadratic state vector (or a 
quadratic state matrix) if its coordinate function is a quadratic 
state function. It is not difficult to show that tensor products and 
tensor contractions of quadratic state vectors are quadratic state 
vectors, which can easily be computed in practice. Since matrix 
multiplication can be obtained from tensor products and contraction, 
the matrix product of two quadratic state matrices is also a 
quadratic state matrix. 

The following Theorem is less obvious:

Theorem 1:
For every quadratic state mapping g: F_2**n -> CC there is a 
function f(e, A, Q) satisfying (1) with g = f(e, A, Q)  such that
the affine mapping  A: F_2**m -> F_2**n is injective.
Given g as any quadratic state mapping, quantities e, A and 
Q: F_2**m -> T  with  g = f(e, A, Q) and m <= n can be computed 
in time polynomial in the input size. 

Sketch proof
For the proof of Theorem 1 we have to calculate a sum 
   S(y) =  sum  Q(x) * (-1)**B(x,y) ,   for x in F_2**m.
Here Q: F_2**m -> T is a quadratic mapping,  B is a symmetric bilinear 
form on F_2**m \times F_2**n. The sum S has to be calculated as a 
function of y for all y in F_2**n. Such sums can be computed as 
generalized Gauss sums by completing the square. With not too much
difficulty we can show that there is an affine subspace W of F_2**n 
such that the restiction of S is a W quadratic mapping (up to a scalar 
factor 2**(e1/2) * w**e2, with w an 8th root of unity), and that 
S(y) = 0 for y not in W.  We omit the details here.
q.e.d.



The complex Clifford group C_n is a group which is defined in [AaGo08] 
and [NeRS00]. It has a unitary representation in CC*(2**n).

Theorem 2
The unitary complex quadratic state matrices in CC*(2**n)
form a representation of the Clifford group C_n.

Proof
It is easy to see that all generators of C_n in [AaGo08] or [NeRS00] 
are unitary complex quadratic state matrices. The group C'_n of such 
matrices is closed under multplication. It is obviously closed under 
computing the inverse, which is the conjugate transpose for a unitary 
matrix. Thus C_n is a subgroup of C'_n. By Theorem 2 the group C_n is 
finite. By [NeRS00], Theorem 6.5, all finite supergroups of C_n in the 
unitary group are generated by C_n and a multiple of the unit matrix 
by a root of unity. Comparing the scalar multiples of the unit matrix 
in C_n and C'_n yields C_n = C'_n.  
q.e.d. 


Call a vector in CC*(2**n) a stablizer state vector if it has norm 1
and if it is also a quadratic state vector. By Theorem 1 the product
of a stabilizer state vector with a unitary quadratic state matrix is 
a stabilizer state vector. Thus Theorem 2 implies that the Clifford 
group C_n stabilizes the stablizer state vectors. It is not difficult 
to see that our stablizer state vectors coincide with the stabilizer 
states in [AaGo08].

Theorems 1 and 2 give a direct way to describe quadratic state 
vectors and matrices, and hence also Clifford group elements and 
stabilizer states in space quadratic in the number n of qubits. 
Operations on these objects such as matrix multiplication and 
computing scalar products can be done in time polynomial in n.


Our implementation of the real Clifford group
---------------------------------------------

After this motivation we feel free to use the language of quantum
computing. But we restrict this implementation to the case of rather 
small real quadratic state functions. More specifically we require
m + n < 64 in (1) and that all intermediate binary vectors fit in
a 64-bit integer. In practice the restricts us to 12 qubits, which
is sufficient to similate the subgroup 2**{1+24}.Co_1 of the monster.

Then the quadratic mapping Q in (1) can be considered as a real 
quadratic form F_2**m -> F_2, with  the embedding F_2 -> T given by 
k -> (-1)**k. Note that Q can be realized as an alternating bilinear  
form B(Q) with  B(Q)(x,y) = Q(x+y) - Q(x) - Q(y) + Q(0), augmented by 
a diagonal vector containing Q(e_k) for all unit vectors e_k of 
F_2**m. For a non-unit vector y the value Q(y) can be computed using 
Q(x1 + x2) =  Q(x1) + Q(x2) + B(Q)(x1, x2).
 
For the full set of complex quadratic functions we would have to 
consider mappings Q: ZZ_4, with the group of integers modulo 4, and 
the embedding F_2 -> ZZ_4 given by x -> 2x (as an embedding of 
additive groups), and the embedding ZZ_4 -> T given by x -> i**x,
where i = sqrt(-1). Then B(Q) is still a symmetric bilinear form and 
the diagonal vector has entries in ZZ4.   

We represent A as an (1 + m) times n matrix. All vector and matrix 
indices start with 0 as in the C language. Matrix A is considered to
have a constant invisible column (with index, say, -1) which has entry 
1 in row and 0 elsewhere. This column corresponds to an ancilla qubit 
with constant value 1. We represent B as an (m + 1) \times (m + 1) 
matrix. Evaluating the quadratic function Q at vector x, for x in 
F_2**m is done by computing (1,x)**t * (matrix B) * (1,x), where '*' is 
matrix multiplication, '**t' is transposition and (1,x) is the vector 
with first coordinate 1 and the other coordinates taken from vector x.    

Since we always evaluate Q at a vector (1,x) in F_2**(1+m), we need
not store the functional values Q(e_k) of Q at the unit vectors. We
have Q((1,e_k)) = Q((0,e_k)) + Q((1,0)) + B((0,e_k),(1,0)); so 
flipping both bits Q((0,e_k)) and B((0,e_k),(1,0)) does not change 
Q(1,e_k). This means that by tacitly assuming Q(e_k,0) = 0, it
suffices to store the bit B((0,e_k),(1,0)).  

In this module we support state vectors, operation of gates on state
vectors, and tensor contraction of state vectors. Here we alway contract 
over the first n qubits of a state vector. So we support tensor 
oparetions on the qubit level only. Tensor products and matrices over 
spaces of vector states are out of the scope of this module. We also 
support writing down stat vectors modulo smaall odd numbers, since we 
need  this feature for calculations in the monster group. 


Some naming conventions for C functions
---------------------------------------

The C functions in this module do operations on the bit matrices A 
and Q used to define a real state vector.

The state vector is desribed by a structure containing (at least) 
the following components

typedef struct {
    uint32_t nrows;    // Number m of rows of bit matrix A
    uint32_t ncols;    // Number n of columns of bit matrix A and Q
    int32_t  factor;   // A scaling factor, see below
    uint64_t *data;    // Pointer to the data bits of A and Q
} qbstate12_type;

Then it is implied that matrix Q has m rows and m columns. A and Q are 
simply concatenated to a matrix M. So M[i, j] refers to A[i, j] for 
j < ncols and to Q[i, j - ncols] for j >= ncols. The scaling factor is
an integer f which is interpreted as the constant factor
   (-1)**(f & 1) * sqrt(2)**(f >> 1) .
For complex state vectors we would have to store powers of sqrt(2)
multiplied by an 8th root of unity.

Data bits of M are stored in a strided array M' of unsigned integers 
of type uint64_t. (see e.g. the Python numpy package for documentation 
of a strided array). So the bit M[i,j] is stored in bit j % 64 of the 
array M'[i, j  / 64], and the position of the integer M'[i, k] is at 
byte  (i * stride[1] + k * stride[0]) relative to the pointer pointer 
to the data bits in the structure. Here bit j of an uint64_t is the
bit with valence 2**j. 

The current implementation requires nrows + ncols <= 64,  so that 
stride[1] is fixed to 1 and we need not care about stride[0]. With 
this restriction we may process transformation matrices for up to 12 
qubits, which is sufficient for implementing the subgroup 
2_+**{1+24}.Co_1 of the monster. In a future version the strides 
should be more flexible, so that more then 12 qubits may be processed. 
Also the full complex Clifford group for a given number of qubits 
could be supported in a future version. 

C functions supporting this module are prefixed with 'qbstate12_'.
Unless otherwise stated, these functions return an int32_t, where
a nonegative value is interpreted as success, and a negative value 
is intepreted as failure. Depending on the function, a nonnegative 
return value may e.g. mean an index for a matrix A, M or Q.  

Most functions take a pointer to a structure of type qbstate12_type
with name pqs, pqs1, pqs2,... . Other typical parameter names are

nqb           Number of qubits, i.e. of columns of matrix A.
nrows         Number of rows of matrix A, Q, and M.
i, i1,...     Index of a row of A, Q or M, starting with 0.
j, j1,..      index of a column of A, with a column of matrix A
              corrsesponding to a qubit, starting with 0.
              If appropriate, an index  j >= ncols refers to
              column (j - ncols) of matrix Q
pv, pv1,...   pointer to a row or column vector of matrix A, Q, or M,
              of type uint64*


Referneces
See file references.txt.


Side remark
-----------

Our theory of quadratic mappings from F_2**n to T can be generalized
to polynomial mappings of higher degree. Then state mappings of
higher degree can be defined similarly. We claim that mappings
of degree d >= 2 can be used to describe quantum computing circuits
with Clifford gates an phase gates of phase exp(2 \pi i * 2**(-d)).
Already in the cubic case d = 3 the group of unitary cubic state
matrices is no longer discrete. In contrary, it unleashes 
(asymtocially) the full power of quantum computing. While matrix
products can still be stored in polynomial space for n qubits and 
degree d, an analogue for Theorem 1 is no longer valid even for 
cubic state mappings.
************************************************************************/

"""



r"""

/*************************************************************************
Updating the quadratic form Q of a state vector

Recall that we implement the quadratic form Q as symmetric matrix, where
the entries of Q are integers modulo 4. The off-diagonal entries of 
matrix Q are known modulo 2 only, but we always have Q[i1,i2] = Q[i2,i1] 
(modulo 4). The elementry operation in function qbstate12_row_op()
means adding row i1 to row i2 of Q and also adding column i1 to 
column i2 of Q. For the off-diagonal elements this is equivalent to a 
row and column operation on a binary matrix. The story is a bit more
complicated for the diagonal elements. Here we have:

     Q[i2, i2] <-  Q[i2, i2] + Q[i1, i1] + 2 * Q[i1, i2]  (mod 4).

For any  0 <= i,j < nrows let  let Q[i,i] = 2*qh_i + q_i_i,  and
Q[i,j] = q_i_j,  for bits qh_i, and q_i_j equal to 0 or 1. Then the  
transformation of Q[i2, i2] given above may be written as:

    qh_i2    <-  qh_i2 ^ qh_i1 ^ q_i1_i2 ^ (q_i1_i1 & q_i2_i2);
    q_i2_i2  <-  q_i2_i2 ^ q_i1_i1;
    q_i2_j   <-  q_i2_j ^ q_i1_j,        for j != i2;
    q_i_j    <-  q_i_j,                  for i != i2, j != i2;
    qh_i_i   <-  qh_i_i,                 for i != j2;

where bit operations '^' and '&' are as in the C language. This means 
that a matrix with entries q_i_j, 0 <= i,j < nrows,  transforms as a 
bilinear form over the binary field. Using the bit operations on a
computer, such a transformation can easily be implemented; and we are 
left with just the computation of the additional bit qh_i2.

In the case of the real Clifford group (which we have implemented here),
there is a further simplification, namely, that all bits q_i_i are zero.

Note that row 0 of our symmetric matrix Q corresponds to a ancilla
qubit with fixed value 1. Mathematically this means that we have to
compute Q(v) = v**t * Q * v for vectors v with component v[0] = 1 only. 
We have seen that flipping the two bits qh_i and q_0_i (with q_i_0 
identical to q_0_i) simultaneously does not change the value Q(v) for 
any such vector v. So for i2 > 0 we need not store the bit qh_i2, and 
we flip bit q_0_i2 whenever the bit qh_i2 is to be flipped. 

We remark that flipping bit qh_0 corresponds to a global phase change 
by an angle of pi (or to multiplication of the unitary matrix 
represented by A and Q with -1). This done by updating component 
'factor' of the structure referring to A and Q. Similarly, changing
q_0_0 means multiplying the global phase by sqrt(-1) or -sqrt(-1). 

This explains the mathematical theory behind the row/column operations
on the quadratic form Q in function qbstate12_row_op() and also in 
function qbstate12_pivot(). The implementations of these two functions 
are still  a bit tricky. Note that in these two functions the relevant 
bits q_i_j of matrix Q are given as follows:

     q_i_j = (m[i] >> (pqs->ncols + j)) & 1,  where m = pqs->data .

For matrix A, with A[i,j] given by A[i,j] = (m[i] >> j) & 1, the 
function qbstate12_row_op() simply adds row i1 to row i2. Here 
addition means '^', i.e. bitwise addition modulo 2.
*************************************************************************/
"""






r"""
/*************************************************************************
Summing up the kernel of the transformation matrix A

Assume that the higest row A[n] of A is zero. Assume further that there
is an i > 0 with Q[i,n] != 0, and that i is the highest such n. Then we 
first add row n to all rows of Q where Q[i,n]= 0. This means pivoting 
column n of Q with row i of Q, and also doing the appropriate row 
operations on A. After pivoting Q loos like this:

      x ... x  ... 0  (0)
        ... x  ... 0          Here we have drawn rows and columns 0,
      x xxx x  xxx 1  (i)     i and n of Q completely, and x denotes
        ... x  ... 0          a bit with an unknown value
      0 000 1  000 x  (n)

Depending on the symmetric matrix after pivoting column n, we have
the following cases:

Case 1: Q[n,n] != 0.

In this case (which does not occur in the real Clifford group) we have 
Q[n,n] = 1, Q[j,n] = 0 for j < n, and hence
Q((v,0)) = sqrt(-1) * Q((v,1)) for any v in F_2**(n). Thus
Q((v,0)) + Q((v,1)) = (1 + sqrt(-1)) * Q((v,0)). So we may delete row n
of A and Q, and column n of Q, and multiply the state vector with  
1 + sqrt(-1). This operation does not change the state vector.

Case 2: For some i with 0 < i < n we have  Q[i,n] != 0.

In this case we pivot column i of Q in the same way as column n has
previously been pivodetd, so that we obtain:

      x ... x  ... 0  (0)
        ... 0  ... 0          Since we might have changed row and column
      x 000 0  000 1  (i)     i, we cannot control bit Q[0,i], due to 
        ... 0  ... 0          the effect of the operation in the comment
      0 000 1  000 0  (n)     after function qbstate12_row_op(). 

Case 2a:  There is a 0 < i < n, such that row A[i] is zero

Then we have to sum over the 4 possible combinations of bit i and bit 
n of vector v. This amounts to computing the sum of (-1)**(i*n + i*x)
for all i, n = 0,1. Clearly, this sum is equal to 2 for any x = 0,1. 
So we may delete rows i and n of A and Q, and also  the corresponding
columns of Q, and multiply the state vector by 2, without changing the 
state vector.

Case 2b:  There is no 0 < i < n, such that row A[i] is zero

Thn there is an i with  Q[i,n] = 0 and row i of A is nonzero.
In this case we have Q((v,0)) = Q((v,1)) if bit i of vector v in
F_2**(n) is zero and  Q((v,0)) = -Q((v,1)) is that bit is equal to one. 
So these values sum up to two, if bit i of v is zero, and they cancel
otherwise. So we may delete row i, since the sum is 0 whenever bit i 
of v is set, without changing the state vector. After deleting that 
row, we have   Q((v,0)) = Q((v,1)), so that we may also delete row Q 
and multiply the state vector by 2, Thus, computationally, Case 2b is 
treated exactly in the same way as Case 2a.  

Case 3: There is no i with 0 < i < n such that Q[i,n] != 0.

Then matrix Q looks as follows:

       x  xxx y  (0)     Here we have drawn rows and columns 0 and n
       x  ... 0          completely, and the operation to be done
       y  000 0  (n)     depends on the bit y = Q[0,n] = Q[n,0].
 
Case 3a: After pivoting column n, the bit Q[0,n] is zero

Then Q((v,0)) = Q((v,1)) for all v. Thus we may delete row Q 
and multiply the state vector by 2.

Case 3b: After pivoting column n, the bit Q[0,n] is one

Then Q((v,0)) = -Q((v,1)) for all v, where bit 0 of v is equal to one.
Bit 0 of v corresponds to an ancilla qubit with fixed value 1. So the 
terms Q((v,0)) and Q((v,1)) cancel in all relevant cases. In this case 
the state vector must be set to zero.

*************************************************************************/

"""

r"""
/*************************************************************************
Applying a Z or CZ gate to a state

Applying an CZ gate to the qubits j1 and j2 of a state qs means to negate 
all entries of the state vector where both qubits, j1 and j2 are equal to
one. Ignoring the global scalar factor in (1), our state qs is given as
a sum of classical states by 

   qs = sum_{x in  F_2**n}  f(x) * |x>,  where

   f(x) =  sum_{y in F_2**m: A(y) = x} Q(y) ,   

for the quadratic state mapping  f: F_2**n -> CC  given by (1). Here (1)
is as in the comment at the beginning of the module. So we have
to negate all terms of that sum where A(y) is in K, with K the
linear subset of entries x in F_2**n, such that both qubits j1 and j2 
of x are set to one.  

Let A[j1] and A[j2] be the j1-th and the j2-th column vector of submatrix 
A of the state. A(y) is in K if (A[j1],y) = (A[j2],y) = 1, with (.,.) the 
scalar product in  F_2**m. So we have

   (A[j1],y) = (A[j2],y) = 1  <==>   (A[j1],y) * (A[j2],y) = 1
                              <==>
       sum_{0 <= k1,k2 < m}  A[k1,j1] * A[k2,j2] * y[k1] *y[k2] = 1,

with y1 = y1[0],...,y1[m-1], y2 = y2[0],...,y2[m-1], and A[k,j] the
entry of A in row k and column j.  The sum in the last line is a 
quadratic form Q1 on F_2**m. Adding that quadratic form to the quadratic 
form Q on F_2**m (modulo 2) achieves the required phase changes for 
the given CZ operation. So we have to replace Q by Q + Q1.

Therefore we put 

   Q[k1,k2] = Q[k1,k2] + A[k1,j1] * A[k2,j2] + A[k1,j2] * A[k2,j1]. 

This achieves the required addition Q = Q + Q1 for all off-diagonal
elements of Q with k1 != k2. For the diagonal elements of Q these two
additions modulo 2 cancel. The diagonal entries of Q have to be
modified as follows:

    Q[k1,k1] = Q[k1,k1] + 2 * A[k1,j1] * A[k1,j2]   (mod 4)

According to our general philosophy discussed at the beginning of
the module, we may replace the last operation by

    Q[k1,0] =  Q[0,k1] = Q[0,k1] +  A[k1,j1] * A[k1,j1] (mod 2)

for k1 > 0. The operation  Q[0,0] = Q[0,0] + 2 * A[0,j1] * A[0,j2]
is equivalent to negating the global scalar factor of the matrix,
if A[0,j1] * A[0,j2] = 1.

An unconditional Z gate is applied to j1 as follows:
We simply set the column vector A[j2] to (1,0,...,0). This is justified
by the fact that row 0 of matrices A and Q corresponds to an ancilla
qubit with fixed value 1. 

*************************************************************************/
"""




r"""
/*************************************************************************
Applying a Hadamard operation to qubit j of a state.

We first pivot column j of submatrix A of the state so that, after 
pivoting,  column j of A has at most one nonzero entry in row i. After
pivoting, A(y) is of shape |u,1,u>  iff such an i exists and y[i] = 1. 
Otherwise it is of shape  |u,0,u>. Then we replace matrices A and Q by 
new matrices A' and Q' as follows.

For obtaining A', we first add a new zero row to A and to Q, and we 
copy all entries for the old rows of A' from A. Of course, we also add 
a zero column to the quadratic form Q'. Let ix be the index of the new 
row of A'. Then we put A'[i,j] = 0,  A'[ix,j] = 1, and
Q'[i,ix] = Q'[ix,i] = A[i,j], without changing any other entry of A'.



Proof of correctness
--------------------

As in the previous comment, our state qs is given as a sum of classical 
states by 

   qs = sum_{x in  F_2**n}  f(x) * |x>,  where

   f(x) =  sum_{y in F_2**m: A(y) = x} Q(y) ,           (H_gate 1)

up to a scalar factor. We write  |u,b,u>  for a state where  qubit j
of that state is equal to the bit b, and the other bits are equal to 
a (usually irrelevant) bit vector u. Then our Hadamard gate maps

   
        |u,0,u>  to  (|u,0,u> + |u,1,u>) / sqrt(2) 
  and   |u,1,u>  to  (|u,0,u> - |u,1,u>) / sqrt(2)  .   (H_gate 2)

Let H(s) be the result of this mapping applied to the state s.
In the following description we ignore the factor 1/sqrt(2). 

Case 1: A[i,j] = 1  for exactly one i > 0
-----------------------------------------

For a column vector y in F_2^m let (y,c) the vector in F_2^(m+1)
with entry c in the last row ix, and the same entries as in y in 
the other rows. Then A(y) = |u,b,u>  implies A'(y,0) = |u,0,u>, for 
b = 0, 1, and for all y in F_2**m. This yields the first term 
|u,0,u>  of H(|u,b,u>) in (H_gate 2), for all terms x = |u,b,u> 
occuring in the sum given by (H_gate 1).

It remains to enter the second term (-1)**b * |u,1,u> of H(|u,b,u>) 
in  (H_gate 2) into the sum given by (H_gate 1). By definition of
A' and Q' we have A'[ix,j] = 1, Q'[i,ix] = Q'[ix,i] = 1; and all
other entries in row ix are zero. Thus  A(y) = |u,b,u>  implies 
A'(y,1) = |u,1,u>, and we have  Q'(y,1) = Q(y) * (-1)**b. 

So the term of the state vector obtained from (y,1) is equal to 

   Q(y) * (-1)**b * |u,1,u> ,  

as required by (H_gate 1) and (H_gate 2).

Case  2: A[i,j] = 0  for all i > 0
----------------------------------

Let b = A[0,j]. In case b = 1 the proof of correctness is the same as 
in Case 1, replacing i by 0. In case b = 0, we have A(y) = |u,0,b>.
Then A'(y,0) = |u,0,b>, A'(y,1) = |u,1,b>, and Q'(y,0) = Q'(y,1) = Q(y)
for all y in F_2**m, as required.

*************************************************************************/
"""
 



r"""

/*************************************************************************
Equating states for tensor contraction in function qbstate12_prep_mul()

Function qbstate12_prep_mul() takes two states qs1 and qs2 and equates 
the first nqb columns, so that they will satify conditions (QS1) and
(QS2) on return of the function.


We focus on the operations on the submatrices A1 and A2 of the states 
qs1 and qs2, which are row operations. These row operations are 
performed by the functions qbstate12_pivot(), qbstate12_xch_rows() and
qbstate12_del_rows(). These functions also process the quadratic forms
Q1 and Q2 contained in the states qs1 and qs2 appropriately, so that 
we need not bother about Q1 and Q2.

In the main loop, a counter 'col_pos' counts from 0 to 'nqs' - 1.
There is also a variable 'row_pos', initalized to one, which markes
the rows of A1 and A2 which have already been processed. 

After the loop iteration with counter = 'col_pos', we claim that the
columns  0, ..., 'col_pos' - 1 of A1 and A2 are equal, and that in all 
rows with index >= 'row_pos', the entries in these columns are zero. 
We also claim that the submatrices of A1 and A2 that consist of rows 
1,...,'row_pos' - 1 and of columns 0,...,'col_pos' - 1 have full rank. 

Entries in rows < 'row_pos' have already been processed. Such rows may 
have nonzero entries and they may no longer be changed independently 
without violating our claim. However, in later iterations of 
the loop we may do row operations on A1 and A2 simultaneously or
delete the same rows of A1 and A2.

In the loop iteration with counter = 'col_pos' we first compute the
index i1 of the highest nonzero bit in column 'col_pos' of A1 and also 
the index i2 of the highest nonzero bit in column 'col_pos' of A2.
The we try to pivot matrices A1 and A2 with these two bits in order
to equate the columns 'col_pos' of the two matrices. The details of
that pivoting process are a bit hairy, since we must maintain the
assertions in our claim for any of the bits in rows with index less
than 'row_pos' of matrices A1 and A2. We consider 4 main cases:

Case 1:  i1 >= 'row_pos',   i2 >= 'row_pos'. 
------
In this case we simply pivot matrix A1 with row i1 and matrix A2 with
row i2, both in column 'col_pos'. After pivoting, all entries of A1 
and A2 in column 'col_pos' are zero, except for A1[i1, col_pos] and
A1[i2, col_pos]. Then we exchange row i1 of A1 with row 'row_pos' and
also row i2 of A2 with row 'row_pos'. Afterwards the columns 
'col_pos' of the matrices A1 and A2 are equal. We increment 'row_pos'
by one to indicate that row 'row_pos' has nonzero entries.

Note that A1[i1, j] = A2[i2, j] = 0 for all j < 'col_pos', so that 
the columns j of of A1 and A2 are not changed for  j < 'col_pos'.
Thus our claim hold for the bits in all columns j < 'col_pos'. 
After pivoting and exchanging rows of A1 and A2, columns 'col_pos'
of A1 and A2 are equal.

After incrementing 'row_pos' and 'col_pos' and exchanging rows, 
matrices A1 and A2 still have full rank, since then we have
A1[row_pos-1, col_pos-1] = A2[row_pos-1, col_pos-1] = 1. 



Case 2:  i1 >= 'row_pos',   i2 < 'row_pos'. 
----
In this case we may pivot matrix A1 with row i1 in column 'col_pos', 
but we may not pivot matrix A2 with row i2 in that column, since 
that pivoting operation may violate our claim. So we pivot A1 with 
row i1 in such a way that  A1[i, col_pos] = A2[i, col_pos] for all i, 
except for i = i1. Then we delete row i1 from matrix A1.
Finally we increment 'col_pos' by one.

This operation does not change columns j of A1 and A2 for 
j < 'col_pos', and, after pivoting, the columns 'col_pos' of matrices
A1 and A are equal. So our claim still holds after this operation. 

Deleting row i1 from matrix A1 is justified as follows:

Note that then we have A1[i, col_pos] = A2[i, col_pos] = 0 for 
i > 'row_pos'. 

Let A1' and A2' be the submatrices of A1 and A2 that consist of the
first 'nqb' columns of A1 and A2, respectively. 

Fact: 
Row i1 of A1' (after pivoting) is not a linear combination of the
other rows of A1' and the rows of A2'.

Proof of the fact:
Let A1*, A2* be the submatrices of A1' and A2' that that consist of
columns 0,...,'col_pos' - 1. By our claim the lowest 'row_pos' - 1
rows of  A1* and  A2* are linear independent and equal, and the
other rows of  A1* and  A2* are zero. Thus any linear combination
of rows of A1' and  A2' that sums up to row i1 of A1 necessarily
consists of the same set of rows  0,...,'row_pos' - 1 of matrices A1' 
and A2', and of an arbitrary set of other rows of A1' and A2'.   
Considering column 'col_pos' of such a linear combination, we see
that the only nonzero terms in that column can come from rows with 
index >= 'row_pos', since the terms coming from row i, i < 'row_pos',
of A1' and A2' cancel.  After pivoting we have 
A1[i, col_pos] = A2[i, col_pos] = 0 for i > 'row_pos' and  
A2[row_pos, col_pos] = 0. Thus the only nonzero term in the set  
{A1[i,'col_pos'],  A2[i,'col_pos'], i >='row_pos'} is A1[i1, col_pos]. 
Hence row i1 of A1' cannot be a linear combination of the other rows 
of A1' or A2'. This proves the fact.


For the required tensor contraction we only have to consider those 
linear combinations of rows of A1' which are equal to any linear 
combination of rows of A2'. Assume that A1' contains a row i1 which 
cannot be expressed as a linear combination of the other rows of A1' 
and the rows of A2'. Then a linear combination of rows of A1' 
containing row i1 cannot be equal to any linear combination of rows 
of A2'. So we may delete row i1 from A1 in that case.


Case 3:  i1 < 'row_pos',   i2 >= 'row_pos'. 
------

This is similar to Case 2, exchanging the role of A1 and A2. 


Case 4:  i1 < 'row_pos',   i2 < 'row_pos'. 
-------
Here we have to consider the following subcases

Case 4.1: columns 'col_pos' of A1 and A2 are equal

In this case we may simply increment 'col_pos' without any other
action; and our claim about the submatrices consisting of rows 
1,...,'row_pos' - 1 and of columns 0,...,'col_pos' - 1 of A1 and A2
remains valid.

Case 4.2: columns 'col_pos' of A1 and A2 are equal except for row 0

In this case we delete all rows from both, A1 and A2, and stop.
This means that we set both states, qs1 and qs2 to zero. 

Since row 0 of A1 or A2 correspond to an ancilla qubit with fixed
value 1, we only have to consider linear combinations of rows of 
A1 or A2 containing row 0. In this case such linear combinations
of rows of A1 and of A2 always differ in column 0. Thus the tensor
contraction of the states qs1 and qs2 is zero in this case. 

Case 4.3: A1[i, col_pos] != A2[i, col_pos] for some i > 0

Let i be the idnex of the highest row with 
A1[i, col_pos] != A2[i, col_pos], and let v be the sum of the column
vectors of A1 and A2 in column 'col_pos' (modulo 2). 

Then we add row i of A1 to all rows j of A1 where v[j] = 1. Similarly,
we add row i of A2 to all rows j of A2 where v[j] = 1. Then we delete
row i in both A1 and also from A2. We also increment 'col_pos'

Clearly, this operation preservs the assertions in our claim. After
these row operations and incrementing 'col_pos', the submatrices of 
A1 and A2 consisting of rows 0, ..., 'row_pos' - 1 and columns
0, ..., 'col_pos' - 1 are equal, with the exception 
A1[i, col_pos] !=  A2[i, col_pos]. A simlar argument as in Case 2
shows that we may delete row i from A1 and A2 without changing
the tensor contraction of the two states qs1 and qs2.
 

Remarks

In Cases 2 and 3 we delete a row from A1 or A2 by exchanging it with
the last row of the matrix and dercementing the number of rows of
the matrix. This operation requires no postprocessing. In case 4 we
delete the same row i from A1 and A2 by marking it for deletion.
Therfore we set bit i in the vector 'deleted'. We also zero the
entries of row i in A1 and in A2 to avoid useless pivoting. At the
end we use function qbstate12_del_rows() to delete the marked rows
from A1 and A2. 
 
*************************************************************************/

"""



