r"""
Computation in the Clifford group 
---------------------------------

In this section we present effective algorithms for computing
in the complex Clifford group of structure 
:math:`\frac{1}{2} ( 2^{1+2n}_+ \times Z_8 ). \mbox{Sp}_{2n}(2)`
defined in  :cite:`NRS01`. Here a factor :math:`\frac{1}{2}` 
means identification of central subgroups of order :math:`2` 
of the factors of a direct product as in :cite:`Atlas`.

The complex Clifford group  :math:`\mathcal{X}_{n}` has a
unitary complex representation of dimension :math:`2^n` 
which has been studied in the theory of quantum computation,
see e.g. :cite:`AG04`. For calculations in the monster group
it would be sufficient to deal with the subgroup 
:math:`\mathcal{C}_{n}` of :math:`\mathcal{X}_{n}` consisting
of the real matrices of that representation. However, the
extra effort for extending the real to the complex case 
is marginal, and using our algorithms for the complex 
Clifford group might be useful in the theory of quantum 
computation.

Usually, in quantum physics it suffices to calculate in a
unitary group up to a scalar multiple of the unit matrix. 
Therefore we cannot simply cut an paste existing algorithms
for calculating in :math:`\mathcal{X}_{n}` used in quantum
physics. 
We keep the following exposition independent of the theory
of quantum computation, but we do not hide the fact that
the main ideas are strongly influenced by that theory.

We present an algorithm that performs matrix multiplication
in  :math:`\mathcal{X}_{n}` in :math:`O(n)^3` bit operation.
Therefore we store an element of  :math:`\mathcal{X}_{n}` 
in :math:`O(n)^2` bits.  In case  :math:`n=12`
this is considerably more efficient than dealing with
complex :math:`4096 \times 4096` matrices.

We will introduce certain complex-valued functions defined
on a Boolean vector space :math:`\mathbb{F}_2^n` which we 
will call *quadratic mappings*. Such a function has a natural
interpretation as a vector in :math:`\mathbb{C}^{2^n}`.
It will turn out that quadratic mappings are closed under 
tensor products and under tensor contraction. Since matrix
multiplication is just a special case of tensor contraction,
we may consider quadratic mappings on     
:math:`\mathbb{F}_2^n \times \mathbb{F}_2^n` as a monoid
of :math:`2^n \times 2^n` matrices closed under matrix 
multiplication. It turns out that the unitary matrices in
this monoid are just a representation of the Clifford group  
:math:`\mathcal{X}_{n}`.



Quadratic mappings
..................

Let :math:`V = \mathbb{F}_2^n` be a Boolean vector space and 
:math:`\mathbb{T}` be the unit circle in 
the set :math:`\mathbb{C}` of the complex numbers. Then 
:math:`\mathbb{F}_2^n` is an additive and :math:`\mathbb{T}` 
is a multiplicative Abelian group. For any mapping 
:math:`f: V \rightarrow \mathbb{T}` define the mapping 
:math:`\Delta f: V \times V \rightarrow  \mathbb{T}` by

.. math::

   \Delta f(x,y) = f(x+y) \cdot f(x)^{-1} \cdot f(y)^{-1} \cdot f(0) .

A mapping :math:`g: V \rightarrow \mathbb{T}`  is bilinear if 

.. math::

   g(x_1 + x_2, x_3) = g(x_1, x_3) \cdot g(x_2, x_3) ,  \\  
   g(x_1, x_2 + x_3) = g(x_1, x_3) \cdot g(x_1, x_3) .

Then we obviously have :math:`g(x_1,x_2) = \pm 1` and 
:math:`g(0,x_2) = g(x_1,0) = 1`. Thus there is a unique symmetric
bilinear form :math:`\beta(g): V \times V \rightarrow \mathbb{F}_2`
with :math:`g(x_1, x_2) = (-1)^{\beta(g)(x_1, x_2)}`.
A function :math:`q: V \rightarrow \mathbb{T}` is called 
*quadratic* if  :math:`q(0) = 1` and  :math:`\Delta q` is bilinear.
For a quadratic function :math:`q` all functional values of
:math:`q` are fourth roots of unity.   
We write :math:`\beta(q)` for :math:`\beta(\Delta  q)`. Put

.. math::
   R_8 = \{ 2^{e/2} \cdot w \mid e \in \mathbb{Z}, 
   w \in \mathbb{C}, w^8 = 1\} \cup \{0\} \; .

Let :math:`e \in R_8`, let :math:`q: \mathbb{F}_2^m \rightarrow  
\mathbb{T}` be a quadratic function, and let
:math:`a: \mathbb{F}_2^m \rightarrow \mathbb{F}_2^n` 
be an affine mapping. Then we define a function  
:math:`f = f(e, a, q): \mathbb{F}_2^n \rightarrow \mathbb{C}` 
by 

.. math::
   :label: quadratic_mapping

   f(e,a,q)(x) = e \cdot
   \sum_{\{y \in \mathbb{F}_2^m \mid a(y) = x\}} q(y) .  

We call a function :math:`f(e,a,q)` satisfying  
:eq:`quadratic_mapping` a *quadratic mapping*. We call the triple  
:math:`(e, a, q)` a *representation* of the quadratic mapping  
:math:`f(e, a, q)`. Occasionally we also consider quadratic
mappings :math:`f(e, a, q)` where :math:`q` is a scalar multiple 
of a  quadratic function with :math:`q(0) \in R_8`. We sometimes
abbreviate :math:`f(e, a, q)(x)` to :math:`f(x)` if the meaning
of :math:`e, a, q` is clear from the context.

By definition, quadratic mappings are invariant under multiplication.
A function :math:`f: \mathbb{F}_2 \rightarrow \mathbb{C}` has a 
natural interpretation as a vector in :math:`\mathbb{C}^2`. 
Similarly, a function 
:math:`g: \mathbb{F}_2^n \rightarrow \mathbb{C}` 
has a natural interpretation as a tensor in the tensor product 
:math:`(\mathbb{C}^2)^{\otimes n}`. If :math:`\mathbb{C}^2` has
basis :math:`(b_0, b_1)` then the tensor corresponding to
:math:`g` has coordinate :math:`g(i_1,\ldots,i_n)` at the basis
vector :math:`b_{i_1} \otimes \ldots \otimes b_{i_n}` of 
:math:`(\mathbb{C}^2)^{\otimes n}`. We call :math:`g` the
*coordinate function* of the tensor.

A vector (or a tensor) in a space over :math:`\mathbb{F}_2^n` is 
called a *quadratic state vector* (or a quadratic state tensor) if 
its coordinate function is a quadratic mapping. Obviously tensor 
products and tensor contractions of quadratic state vectors are 
quadratic state vectors, which can easily be computed in practice. 
A matrix is a special kind of a tensor, and matrix multiplication 
is a special case of the contraction of a tensor product. Thus the 
matrix product of two quadratic state matrices is also a quadratic 
state matrix. 


For any affine mapping 
:math:`a: \mathbb{F}_2^m \rightarrow \mathbb{F}_2^n`
there is a unique linear mapping :math:`l(a)` that differs from 
:math:`a` by the constant :math:`a(0)`. We write :math:`\ker a`
for :math:`\ker l(a)`. We call a representation :math:`(e, a, q)`
of a quadratic mapping *injective* if  :math:`\ker a =  0`.
We have:



Lemma 1:

Every quadratic mapping has an injective representation.

Proof

Let 
:math:`e \in R_8, a : \mathbb{F}_2^m \rightarrow \mathbb{F}_2^n`, 
and :math:`q` a quadratic function on  :math:`\mathbb{F}_2^m` such 
that :math:`g = f(e, a, q)` and :math:`m` is minimal. 
Assume :math:`\ker a \neq 0`. Then we construct a tuple 
:math:`(e', a', q'`) with :math:`g = f(e', a', q')`, such that the
domain of  :math:`a'` and of  :math:`q'` is  a proper affine 
subspace of the domain of :math:`a`.


Let :math:`v \in \ker A \setminus  \{0\}` and :math:`W` 
be a subspace of :math:`\mathbb{F}_2^m` with   
:math:`\left<v \right> \oplus W = \mathbb{F}_2^m`. Then for 
:math:`w \in W` we have 
:math:`q(w+v) = q(v) \cdot (-1)^{h(w)} \cdot r`, 
with :math:`h` a linear form on :math:`W` given by
:math:`h(w) = \beta(q)(v, w)` and :math:`r = q(v)/q(0)`
a fourth root of unity. Since :math:`v \in \ker A`, we have

.. math::

   f(e,a,q)(x) = e \cdot
   \sum_{\{y \in W \mid a(y) = x\}} q(y) + q(y + v)  
   = e \cdot
   \sum_{\{y \in W \mid a(y) = x\}} q(y) t(h(y)) \; ,
   
where :math:`t(z) = 1 + r \cdot (-1)^z` for 
:math:`z \in \mathbb{F}_2`. Let :math:`W'` be the 
support of :math:`t \circ h`, i.e.
:math:`W' = \{y \in W \mid t(h(y)) \neq 0\}`. In case 
:math:`W' = \emptyset` we have :math:`f(e,a,q) = 0`.
Otherwise it suffices to show that :math:`W'` is an affine 
subspace of :math:`W`, and that the restriction of
:math:`t \circ h` to :math:`W'` is a quadratic function up to
a scalar multiple in  :math:`R_8`.
 
Since :math:`r` is a fourth root ot unity, :math:`t(0)` and
:math:`t(1)` are in :math:`R_8`. If :math:`h` is zero then 
:math:`t \circ h` is a constant function on :math:`W`. 

If :math:`r` is an imaginary fourth root of unity then
:math:`t(0) \neq 0` and :math:`t(1)/t(0) = \pm \sqrt{-1}`. 
Thus :math:`t` is a quadratic function on :math:`\mathbb{F}_2` 
up to the factor :math:`t(0)`. Since :math:`h` is linear, 
:math:`t \circ h` is also a quadratic function on :math:`W`. 

If :math:`r = \pm 1, h \neq 0`, then the function
:math:`t \circ h` is constant on :math:`W'`, and we have
:math:`W' = \ker h` if :math:`r = 1` and 
:math:`W' = W \setminus  \ker h`  if :math:`r = -1`.

q.e.d.



The complex Clifford group :math:`\mathcal{X}_n` is a group which 
is defined in  :cite:`AG04` and :cite:`NRS01`. It has a unitary 
representation in  :math:`(\mathbb{C}^2)^{\otimes n}`.

Lemma 2

The unitary complex quadratic state matrices representing
endomorphisms of the space :math:`\mathbb{C}^{2^n}` make up a 
representation of the complex Clifford group  
:math:`\mathcal{X}_n` .

Sketch Proof

It is easy to see that all generators of :math:`\mathcal{X}_n` 
in  :cite:`NRS01` are unitary complex quadratic 
state matrices. The group  :math:`\mathcal{X}'_n` of such 
matrices is closed under multiplication. It is obviously closed under 
computing the inverse, which is the conjugate transpose for a unitary 
matrix. Thus :math:`\mathcal{X}_n`  is a subgroup of 
:math:`\mathcal{X}'_n` . By Lemma 1 the group :math:`\mathcal{X}'_n` 
is  finite. By Theorem 6.5 in :cite:`NRS01`  all finite supergroups 
of :math:`\mathcal{X}_n` in the unitary group are generated by 
:math:`\mathcal{X}_n` and a scalar multiple of the unit matrix 
by a root of unity. Comparing the scalar multiples of the unit matrix 
in :math:`\mathcal{X}_n` and :math:`\mathcal{X}'_n` 
yields :math:`\mathcal{X}_n` = :math:`\mathcal{X}'_n`.  

q.e.d. 


Background from the theory of quantum computing
...............................................


In the theory of quantum computation a state vector representing the 
state of :math:`n` qubits can be written as a vector in 
:math:`(\mathbb{C}^2)^{\otimes n}`, where the :math:`2^n` basis vectors 
of that space are labelled by the elements of  :math:`\mathbb{F}_2^n`. 
Here the :math:`n` qubits correspond to the :math:`n` factors
:math:`\mathbb{F}_2` of   :math:`\mathbb{F}_2^n`. 
In  :cite:`AG04` the unit vectors which are also  quadratic state 
vectors are called  *stabilizer states*.  By Lemma 1 the product
of a stabilizer state  with a unitary quadratic state matrix is 
a stabilizer state. Thus Lemma 2 implies that the Clifford 
group :math:`\mathcal{X}_n` stabilizes the stabilizer state vectors. 
In the theory of quantum computation this fact is known as the
Gottesman-Knill theorem. 

In :cite:`AG04` there is a fast algorithm for calculating in the
Clifford group   :math:`\mathcal{X}_n`. As usual in quantum theory,
this algorithm ignores scalar factors in the matrix representation
of :math:`\mathcal{X}_n`. This means that we have to create our
own algorithm for computing in  :math:`\mathcal{X}_n`. 

We remark that in the  graphical ZX-calculus (which is used for 
describing linear maps between qubits) is an appropriate setup
for 'explaining' the operations defined in the following 
sections, see https://en.wikipedia.org/wiki/ZX-calculus .




Implementation of quadratic mappings
....................................

Let :math:`g = f(e, a, q): \mathbb{F}_2^n \rightarrow \mathbb{C}` 
with :math:`e \in R_8`, 
:math:`a : \mathbb{F}_2^m \rightarrow \mathbb{F}_2^n` an affine 
mapping, and :math:`q: \mathbb{F}_2^m \rightarrow  \mathbb{T}`. 
We implement the quadratic mapping :math:`g` as a triple 
:math:`(e, A, Q)`. 

Here :math:`A` is an :math:`(m + 1) \times n` bit matrix 
representing :math:`a`, and :math:`Q` is a symmetric 
:math:`(m + 1) \times (m + 1)` bit matrix representing :math:`q`. 
All vector and matrix indices start with :math:`0` as usual in 
the C language. For 
:math:`x = (x_1,\ldots, x_m) \in \mathbb{F}_2^m`, 
:math:`x_0 = 1` and bit matrices :math:`A, Q` as above we put:

.. math::
    a(x) =(x_0,\ldots, x_m)  \cdot A  \quad ,  \quad
    q(x) = \exp \left(\pi \sqrt{-1} /2 \cdot \
    \sum_{j,k = 0}^m Q_{j,k} x_j x_k  \right) \; .

We write :math:`\sqrt{-1}` for the imaginary unit,
and we use the letters :math:`i, j, \ldots` as indices.
    
We stress that all bits :math:`Q_{j,k}, x_j, x_k` in that sum
must be interpreted as integers equal to  :math:`0` or :math:`1` 
(modulo :math:`4`). But since :math:`Q` is a symmetric matrix, it 
suffices to know the off-diagonal entries of :math:`Q` modulo 
:math:`2` if we assume  :math:`Q_{j,k} = Q_{k,j}`. 
Due to the condition :math:`q(0) = 1` we only consider matrices 
:math:`Q` with :math:`Q_{0,0} = 0`. It is easy to check that we
can encode any quadratic function
:math:`q: \mathbb{F}_2^m \rightarrow  \mathbb{T}` as a symmetric
:math:`(m + 1) \times (m + 1)` bit matrix :math:`Q`
(with :math:`Q_{0,0} = 0`) uniquely as above. Also, we can encode
any affine mapping
:math:`a : \mathbb{F}_2^m \rightarrow \mathbb{F}_2^n` as an 
:math:`(m + 1) \times n` bit matrix :math:`A` (not uniquely)
as above.

Given :math:`e` and matrices :math:`A, Q` as above, we define 
:math:`f(e,A,Q) = f(e,a,q)`, where :math:`a` and :math:`q` 
are as in the last equation.

We encode a number :math:`e \in R_8 \setminus \{0\}` as an
integer :math:`e'` such that 
:math:`e = \exp(e' \pi \sqrt{-1} / 4) 
\cdot 2^{\lfloor e'/16 \rfloor / 2}`,
with :math:`\lfloor x \rfloor` the greatest integer :math:`\leq x`.
We encode the constant quadratic mapping :math:`0` as a matrix
:math:`A` with zero rows.


Let :math:`(b_1, \ldots, b_m)` be the standard basis of 
:math:`\mathbb{F}_2^m`.
For :math:`1 \leq i, j \leq m, i \neq j`, define the linear 
transformation 
:math:`T_{i,j} : \mathbb{F}_2^m \rightarrow \mathbb{F}_2^m` 
by:

.. math::
   T_{i,j}(b_j) = b_j + b_i \, , \,   T_{i,j}(b_k) = b_k
   \quad \mbox{for} \quad k \neq j \; . 

We also define the (affine) translation 
:math:`T_{0,j} : \mathbb{F}_2^m \rightarrow \mathbb{F}_2^m` by
:math:`T_{0,j}(x) = x + b_j` for all :math:`x \in \mathbb{F}_2^m`. 

Let :math:`q'` be the quadratic function :math:`q \circ T_{i,j}` 
given by :math:`x \mapsto q(T_{i,j}(x))`. Let :math:`Q'` be the
symmetric bit matrix representing :math:`q'`.
Then for :math:`i, j > 0, i \neq j` we have

.. math::
   \begin{array}{rcll}
   Q'_{i,0}  = Q'_{0,i}  & = & Q_{i,0} + Q_{j,0} + Q_{i,j} 
   + Q_{i,i} \cdot Q_{j,j}    \, ,   \\      
   Q'_{i,i} &  = & Q'_{i,i} + Q_{j,j}   \, ,   \\  
   Q'_{i,k} = Q'_{k,i} & = &  Q_{i,k} + Q_{j,k}  &
   \mbox{for} \quad  k > 0, k\neq i  \, ,   \\ 
   Q'_{k,l} & = & Q_{k,l} &
   \mbox{for all other pairs} \quad (k, l) \; .
   \end{array}

The quadratic function 
:math:`q \circ T_{0,j}, \;  j > 0` is equal to :math:`e \cdot q'`
with :math:`e = \exp(\pi \sqrt{-1} \cdot(Q_{0,i} + Q_{i,i}/2))`.
Here  :math:`q'` is the quadratic function represented by the bit
matrix  :math:`Q'` with


.. math::
   \begin{array}{rcll}
   Q'_{0,0} & = & 0   \, ,   \\   
   Q'_{0,k} = Q'_{k,0} &  = &  Q_{0,k} + Q_{j,k}  &
   \mbox{for} \quad  k > 0 \, ,   \\ 
   Q'_{k,l} & = & Q_{k,l} &
   \mbox{for all other pairs} \quad (k, l) \; .
   \end{array}

For an affine mapping
:math:`a: \mathbb{F}_2^m \rightarrow \mathbb{F}_2^n` the mapping 
:math:`a' = a \circ T_{i,j}` with :math:`i \geq 0, j > 0, j \neq i`
is represented by the matrix :math:`A'` given by:

.. math::
   \begin{array}{rcll}
   A'(i, l) & = &  A(i, l) + A(j, l) \, ,   \\ 
   A'(k, l) & = &  A(k, l)  &  \mbox{for} \quad  k \neq j  \; .
   \end{array}


This means that adding row and column :math:`j >0` to row and 
and column :math:`i \geq 0` of a matrix  :math:`Q` representing 
a quadratic function changes that matrix to a matrix :math:`Q'`
representing :math:`Q \circ T_{i,j}`, up to a scalar factor 
(which is a fourth root of unity) and some corrections
required for the entries 
:math:`Q'_{0,j}, Q'_{j,0}` and :math:`Q'_{0,0}` of
:math:`Q'`. Similarly, adding row :math:`j` to row :math:`i` 
of a matrix  :math:`A` representing an affine mapping from 
:math:`\mathbb{F}_2^m` to :math:`\mathbb{F}_2^n` changes
matrix  :math:`A` to a matrix representing the affine mapping
:math:`A \circ T_{i,j}`.
 

We obviously have 
:math:`f(e, A, Q) = f(e, A \circ T_{i,j}, Q \circ T_{i,j})`.
So we may perform a row operation on matrix :math:`A`, and 
also a row and a column operation on matrix `Q` without
changing `f(e, A, Q)`, (up to a few corrections in line and
column :math:`0` of :math:`Q`). 
  
If  :math:`f(e, A, Q)` is in *echelon form*, (as defined in the
following section), then it remains in echelon form after 
applying an operation  :math:`T_{i,j}, i < j`.



Reducing the representation of a quadratic mapping
..................................................

We call a representation  :math:`(e, A, Q)`  of a quadratic 
mapping *reduced* if the following conditions hold:

  * :math:`A` has no nonzero rows.

  * Let :math:`A_1` be the matrix obtained from :math:`A` by 
    prepending the column vector :math:`(1,0,\ldots,0)^\top`. 
    Then the leading coefficient of a row of :math:`A_1`
    is always strictly  to the right of the leading coefficient 
    of the row above it.
    
  * Each column containing a leading coefficient of a row
    has  zeros in all its other entries.  
  
Here the *leading coefficient* of a row of matrix :math:`A` is 
the first nonzero entry in that row. 

Obviously, a reduced representation :math:`(e, A, Q)` also
injective. It is easy to see that any quadratic mapping has a
unique reduced representation. 


In this section we present an algorithm for converting a 
representation of a quadratic function to a reduced 
representation. Function ``qstate12_reduce`` in module
``qstate12.c`` implements this algorithm.


A matrix is in *row echelon form* if 

  * All rows consisting of only zeros are at the bottom.
  
  * The leading coefficient of a nonzero row is always strictly 
    to the right of the leading coefficient of the row above it.
    
A matrix  :math:`A` is in *reduced row echelon* form if  it is 
in echelon form, and each column containing a leading coefficient 
of a row has  zeros in all its other entries.  See
https://en.wikipedia.org/wiki/Row_echelon_form .

In order to reduce a representation :math:`(e, A, Q)` of a
quadratic mapping :math:`f(e, A, Q)` we apply several 
operations on the components :math:`e, A, Q` that do not change
:math:`f(e, A, Q)`. Let :math:`T_{i,j}` be the operation on 
:math:`e, A, Q` defined in the last section. Let 
:math:`X_{i,j}, i, j > 0` be the operation on :math:`e, A, Q` 
that exchanges rows :math:`i` with row :math:`j` of :math:`A` 
and :math:`Q`, and also column :math:`i` with column :math:`j` 
of :math:`Q`. Neither of these two operations changes 
:math:`f(e, A, Q)`. 

Given :math:`A`, let :math:`A_1 = A_1(A)` be obtained from 
:math:`A` as above. By applying a sequence of operations 
:math:`T_{i,j}`, :math:`X_{i,j}` we may convert :math:`A_1` 
to reduced echelon form. If :math:`A_1` is in reduced
echelon form and contains no zero rows at the bottom then the 
representation :math:`(e, A, Q)` is reduced. Otherwise we use 
the following algorithm repeatedly to remove the zero rows 
from the bottom of :math:`A`.

Let :math:`i` be the index of the last row of :math:`A` and assume
:math:`A_{i,j}=0`  for all :math:`j`.

Let :math:`i'` be the highest index with  :math:`Q_{i,i'}=1`.
If such an :math:`i'` exists then we add row :math:`i'` of 
:math:`A` to all  rows :math:`k` where   :math:`Q_{k,i}=1` holds
and we adjust :math:`Q`. Afterwards we have  :math:`Q_{i',i}=1`
for at most one index  :math:`i'`.



Case 1:  :math:`Q_{i',i}=0` for all :math:`i'`

The we may delete the last row of :math:`A` , adjust :math:`Q`,
and double  :math:`e` without changing :math:`f(e, A, Q)`.


Case 2:  :math:`Q_{0,i}=1` 

Then :math:`f(e, A, Q)` is the constant function :math:`0`.
  
  
Case 3.  :math:`Q_{i',i}=1, 0 < i' < i`

Then we add row :math:`i` of  :math:`A` to all  rows 
:math:`k \notin \{i,i'\}` where   :math:`Q_{k,i'}=1` holds
and we adjust :math:`Q`. Afterwards we have
:math:`Q_{k,i'} = Q_{i',k} = Q_{k,i} = Q_{i,k} = 0` for all
:math:`k \notin \{i,i'\}`,  :math:`Q_{i,i} = 0` , and 
:math:`Q_{i,i'} = Q_{i',i} = 1`. So we may delete rows
:math:`i`  and :math:`i'`  of :math:`A`, adjust  :math:`Q`,
and double :math:`e` without changing :math:`f(e, A, Q)`.


Case 4:  :math:`Q_{i,i}=1`

Then :math:`f(e, A, Q)` is not changed if we delete the last row 
:math:`i` of :math:`A`, adjust :math:`Q`, and multiply :math:`e` 
by :math:`1 + \sqrt{-1}`.


Remark

The algorithm sketched in this subsection yields an alternative
proof of Lemma 1.


Extending a quadratic mapping
.............................


Let :math:`g: \mathbb{F}_2^{n} \rightarrow \mathbb{C}` be a 
quadratic mapping with :math:`g = f(e,A,Q)`, where  :math:`A`
is an :math:`(m+1) \times n`  and :math:`Q` is an 
:math:`(m+1) \times (m+1)` bit matrix. Define 
:math:`g': \mathbb{F}_2^{n+1}\rightarrow \mathbb{C}` by
 
.. math::
   g'(x_0,\ldots, x_{n-1}, x_n) = g(x_0,\ldots, x_{n-1}) \; .
   
Thus  :math:`g'(x)` does not depend on the last bit of :math:`x`.
   
Then we have :math:`g = f(e,A',Q')` for matrices :math:`A',Q'`
defined as follows. :math:`A'` is obtained form :math:`A'`  by
appending a zero row at the bottom and a zero column at the
right and changing  the rightmost lowest entry  
:math:`A'_{m+1,n}` to 1.  :math:`Q'` is obtained form :math:`Q'`  
by appending a zero row at the bottom and a zero column at the
right.

Of course, we may insert a zero column at any position :math:`j` 
of matrix :math:`A` instead at the end, and perform the other 
modifications of :math:`A` math :math:`Q` as above. Let
:math:`A^{j}` and math :math:`Q^{j}` be the modified matrices
:math:`A` math :math:`Q`. Then we have

.. math::
   f(e,A^{j},Q^{j})(x_0,\ldots, x_{j-1}, x_j, \ldots, x_n) = 
   f(e,A,Q)(x_0,\ldots, x_{j-1}, x_{j+1}, \ldots, x_n) \; .

Function ``qstate12_extend`` in module ``qstate12.c`` 
implements the extension of a quadratic mapping.

Restricting a quadratic mapping
...............................

Let :math:`g: \mathbb{F}_2^{n} \rightarrow \mathbb{C}` with
:math:`g = f(e,A,Q)` be as in the last section. Define 
:math:`\hat{g}^{(j)}: \mathbb{F}_2^{n}\rightarrow \mathbb{C}` 
by
 
.. math::
   \hat{g}^{(j)}(x_0,\ldots, x_j, \ldots, x_{n-1}) = 
   \left\{
   \begin{array}{ll}
   g(x_0,\ldots, x_j, \ldots, x_{n-1}) & 
   \quad \mbox{if} \quad x_j = 0    \\ 
   0   & \quad \mbox{if} \quad  \mbox  x_j = 1
   \end{array}
   \right.

Function ``qstate12_restrict_zero`` in module ``qstate12.c`` 
implements the restriction of a quadratic mapping.

We can compute matrices :math:`\hat{A}^{j}` and math
:math:`\hat{A}^{j}` with  
:math:`\hat{g}^{(j)} = f(\hat{e}^{j},\hat{A}^{j},\hat{Q}^{j})` 
as follows.

Case 1:  :math:`A_{i,0}=0` for all :math:`i`

The we put 
:math:`(\hat{e}^{(j)},\hat{A}^{(j)},\hat{Q}^{(j)}) = (e, A, Q)`.


Case 2: :math:`A_{i,0}=1`, :math:`A_{i,j}=0` for :math:`i>0`

Then :math:`\hat{g}^{(j)}` is the constant function :math:`0`
and we may put :math:`\hat{e}^{(j)} = 0`.


Case 3: There is an :math:`i > 0` with :math:`A_{i,j} = 1`

Then we add row :math:`i` of matrix :math:`A` to all rows
:math:`k` of :math:`A` where :math:`A_{k,j} = 1`, and we 
adjust :math:`Q` and :math:`e`  so that :math:`g` is not changed. 
Let :math:`(e', A', Q')` be the adjusted triple :math:`(e, A, Q)`.
Then :math:`f(e', A', Q') = g(e, A, Q)` and column :math:`j` of
:math:`A` has precisely one nonzero entry in row :math:`i`.  

We obtain :math:`\hat{A}^{(j)}` from :math:`A` by deleting 
row :math:`i`, and  :math:`\hat{Q}^{(j)}` from :math:`Q` by 
deleting row  and column :math:`i`. We put  
:math:`\hat{e}^{(j)} = e`.

Remark

We obtain the restriction of :math:`g` from 
:math:`\mathbb{F}_2^n` to 
:math:`\mathbb{F}_2^{} \times \{0\} \times \mathbb{F}_2^{n-1-j}`
as :math:`g( \hat{e}^{(j)}, A', \hat{Q}^{(j)})` where :math:`A'` 
is obtained from  :math:`\hat{A}^{(j)}` by deleting column 
:math:`j`. This special case of a restriction is implemented in 
function ``qstate12_restrict_zero`` in module ``qstate12.c``. 



Products and tensor products of quadratic mappings
..................................................

TODO: Check documentation from this point on!!!

In this section we present an algorithm for multiplying
quadratic mappings. Later we will use this algorithm for
tensor contraction of quadratic state matrices an also for 
matrix multiplication of quadratic state matrices.

For :math:`\lambda = 1, 2` let  
:math:`g^{(\lambda)} : \mathbb{F}_2^{n_\lambda}
\rightarrow   \mathbb{C}` be quadratic mappings.
for :math:`j \leq \min(n_1, n_2)` we define a mapping 
:math:`(g^{(1)} \odot g^{(2)})_j : \mathbb{F}_2^{n_1+n_2-j}
\rightarrow  \mathbb{C}` by

.. math::
    (g^{(1)} \odot g^{(2)})_j(x, x_1, x_2) = 
    g^{(1)}(x, x_1) \cdot   g^{(2)}(x, x_2) \; ,  \quad 
    x \in \mathbb{F}_2^{j} \; , \quad
    x_\lambda \in \mathbb{F}_2^{x_\lambda-j} \; .

Obviously, :math:`(g^{(1)} \odot g^{(2)})_n` is a quadratic 
mapping. Condsidering the corresponding quadratic state 
vectors we see that
:math:`(g^{(1)} \odot g^{(2)})_0` is just the tensor product
:math:`g^{(1)} \otimes g^{(2)}` of  
:math:`g^{(1)}` and :math:`g^{(2)}`. 




For actually computing :math:`(g^{(1)} \odot g^{(2)})_j` we
may extend the mapping 
:math:`g^{(1)}: \mathbb{F}_2^{j} \times \mathbb{F}_2^{n_1-j}
\rightarrow  \mathbb{C}`
to a mapping
:math:`g^{(1')}: \mathbb{F}_2^{j} \times \mathbb{F}_2^{n_1-j}
\times \mathbb{F}_2^{n_2-j} \rightarrow  \mathbb{C}`, with
:math:`g^{(1')}`  not depending on the last factor
:math:`\mathbb{F}_2^{n_2-j}`, as described in one of the
last sections. Similarly, we may extend :math:`g^{(2)}` to
a mapping 
:math:`g^{(2')}: \mathbb{F}_2^{j} \times \mathbb{F}_2^{n_1-j}
\times \mathbb{F}_2^{n_2-j} \rightarrow  \mathbb{C}`, with
:math:`g^{(2')}`  not depending on the factor
:math:`\mathbb{F}_2^{n_2-j}` in the middle.

Then we simply have 
:math:`(g^{(1)} \odot g^{(2)})_j = g^{(1')} \cdot g^{(2')}`.


In the next section we give an effective algorithm for computing 
the product :math:`h_1 \cdot h_2` of two quadratic mappings 
:math:`h_1 ,  h_2:  \mathbb{F}_2^k \rightarrow \mathbb{C}`.


The appropriate setup for 'explaining' the operation
:math:`(. \odot .)_j` is the ZX-calculus. Since we only need
this operation for fast matrix multiplication, we do not go
into details here.



Multiplication of quadratic mappings
....................................


For  :math:`\lambda = 1,2` let
:math:`g^{(\lambda)} : \mathbb{F}_2^{n_\lambda}` be 
quadratic mappings and :math:`j \leq \min(n_1, n_2)`
as in the last section. Let
:math:`g^{(\lambda)} = 
f\big(e^{(\lambda)}, A^{(\lambda)}, Q^{(\lambda)} \big)` 
be a reduced representation of :math:`g^{(\lambda)}`.
We assume that :math:`A^{(\lambda)}` is in reduced
echelon form. 


For any  :math:`j \leq \min(n_1, n_2)` there are quadratic
mappings :math:`g^{(\lambda,j)} = 
f\big(e^{(\lambda,j)}, A^{(\lambda,j)}, Q^{(\lambda,j)} \big)` 
with 
:math:`(g^{(1,j)} \odot g^{(2,j)})_j` =
:math:`(g^{(1)} \odot g^{(2)})_j`, where the first 
:math:`j` columns of :math:`A^{(1,j)}` and :math:`A^{(2,j)}`
are equal, and both, :math:`A^{(1,j)}` and :math:`A^{(2,j)}`
are in reduced echelon form.
We put  :math:`A^{(\lambda,0)} = A^{(\lambda)}`.

Next we give an algorithm for computing
:math:`\big(e^{(\lambda,j)}, A^{(\lambda,j)}, Q^{(\lambda,j)}\big)`
from
:math:`\big(e^{(\lambda,j-1)}, A^{(\lambda,j-1)}, 
Q^{(\lambda,j-1)}\big)`.

Asumme that 
:math:`e^{(\lambda,j-1)}, A^{(\lambda,j-1)}, 
Q^{(\lambda,j-1)}` satisfy the conditions above.


Case 1:  Both, :math:`A^{(1,j-1)}` and :math:`A^{(2,j-1)}`,
have a row with leading coefficient in column  :math:`j`.

Since :math:`A^{(1,j-1)}` and :math:`A^{(2,j-1)}` are in reduced
echelon form and the first :math:`j-1` columns of these two
matrices are equal, we conclude that the first  :math:`j` 
columns of  :math:`A^{(1,j)}` and :math:`A^{(2,j)}` are equal.

So we may put 
:math:`\big(e^{(\lambda,j)}, A^{(\lambda,j)}, Q^{(\lambda,j)}\big) =
\big( e^{(\lambda,j-1)}, A^{(\lambda,j-1)}, Q^{(\lambda,j-1)}\big)` 
for :math:`\lambda  = 1,2`.

Case 2: Only :math:`A^{(1,j-1)}` has a row with leading 
coefficient in column  :math:`j`.

Assume that this row of :math:`A^{(1,j-1)}` has index :math:`i`. 
We add row :math:`i` to all rows :math:`k` of :math:`A^{(1,j-1)}` 
where :math:`A^{(1,j-1)}_{k,j} \neq  A^{(2,j-1)}_{k,j}`. 
We also adjust  :math:`Q^{(1,j-1)}` and  :math:`e^{(1,j-1)}`
so that :math:`g^{(1,j-1)}` is not changed. Write 
:math:`A^{(1,j-1)'}` for the matrix obtained from
:math:`A^{(1,j-1)}` by performing these row operations.
Let  :math:`Q^{(1,j-1)'}` and  :math:`e^{(1,j-1)'}` be the
adjusted quantities such that
:math:`f\big(e^{(1,j-1)'}, A^{(1,j-1)'}, Q^{(1,j-1)'}\big)`
= :math:`f\big(e^{(1,j-1)}, A^{(1,j-1)}, Q^{(1,j-1)}\big)`.

We obtain  :math:`A^{(1,j)}` from :math:`A^{(1,j-1)'}` by deleting 
row :math:`i` of  :math:`A^{(1,j-1)'}`. We obtain :math:`Q^{(1,j)}` 
from :math:`Q^{(1,j-1)'}` by deleting row and column :math:`i`.  
We put :math:`e^{(1,j)} = e^{(1,j-1)'}`. We put
:math:`\left(e^{(2,j)},A^{(2,j)},Q^{(2,j)}\right)` =  
:math:`\left(e^{(2,j-1)},A^{(2,j-1)},Q^{(2,j-1)}\right)`.

By construction, matrix :math:`A^{(1,j)}`  is in reduced echelon
from and the first :math:`j` columns of  :math:`A^{(1,j)}` and 
:math:`A^{(2,j)}` are equal.

It remains to show that deleting row :math:`i` of matrix 
:math:`A^{(1,j-1)'}` does not change 
:math:`(g^{(1)} \odot g^{(2)})_j`. 

Let :math:`D^{(\lambda)}` be the submatrix of 
:math:`A^{(\lambda,j-1)'}` that consists of the first :math:`j`
columns of :math:`A^{(1,j-1)'}`. For computing 
:math:`(g^{(1)} \odot g^{(2)})_j` we only have to consider
rows of matrix :math:`D^{(1)}` that are linear combinations of 
rows of matrix :math:`D^{(2)}`, excluding row :math:`0` of both 
matrices. By construction, :math:`D^{(1)}` and :math:`D^{(2)}` 
are in reduced echelon form, row :math:`i` of :math:`D^{(1)}` has 
its leading coefficient in column :math:`j`, and in 
:math:`D^{(2)}` there is no row with leading coefficient in 
column :math:`j`. Thus row  :math:`i` of :math:`D^{(1)}`  is
not a linear combination of the rows of  :math:`D^{(2)}`,
ignoring row :math:`0` of :math:`D^{(2)}`.



Case 3: Only :math:`A^{(2,j-1)}` has a row with leading 
coefficient in column  :math:`j`.

This case is symmetric to case 2, exchanging the role of
:math:`A^{(1,j-1)}` and :math:`A^{(2,j-1)}`.

Case 4:
Neither :math:`A^{(1,j-1)}` nor  :math:`A^{(2,j-1)}` has a 
row with leading coefficient in column  :math:`j`.

Case 4.1: Columns :math:`j` of :math:`A^{(1,j-1)}` and  
:math:`A^{(2,j-1)}`  are equal.

Then we may proceed as in case 1.

Case 4.2: Column :math:`j`of :math:`A^{(1,j-1)}` and  
:math:`A^{(2,j-1)}` are equal except in row :math:`0`.

Assume :math:`A^{(\lambda,j-1)}_{0,j} =  \lambda + y \pmod{2}`. 
Then 
:math:`g^{(\lambda,j-1)}(x_0, \ldots, x_{j-1}, x_j, \ldots) = 0` 
if  :math:`x_j \neq \lambda + y \pmod{2}`. Hence
:math:`(g^{(1)} \odot g^{(2)})_j` is zero for all arguments
and we may put :math:`e^{(1,j)}  = e^{(2,j)}  = 0`.



Case 4.3: There is an :math:`i>0` with 
:math:`A^{(1,j-1)}_{i,j}  \neq A^{(2,j-1)}_{i,j}` 


Let :math:`i` be the index of the highest row with
:math:`A^{(1,j-1)}_{i,j}  \neq A^{(2,j-1)}_{i,j}` .

For :math:`\lambda = 1,2` we add row :math:`i` to all rows 
:math:`k` of :math:`A^{(\lambda,j-1)}` where 
:math:`A^{(1,j-1)}_{k,j} \neq  A^{(2,j-1)}_{k,j}`. 
We also adjust  :math:`Q^{(\lambda,j-1)}` and  
:math:`e^{(\lambda,j-1)}` so that :math:`g^{(\lambda,j-1)}` is 
not changed. Write :math:`A^{(\lambda,j-1)'}` for the matrix 
obtained from :math:`A^{(\lambda,j-1)}` by performing these 
row operations. Let  :math:`Q^{(\lambda,j-1)'}` and  
:math:`e^{(\lambda,j-1)'}` be the adjusted quantities such that
:math:`f\big(e^{(\lambda,j-1)'}, A^{(\lambda,j-1)'}, 
Q^{(\lambda,j-1)'}, \big)`
= :math:`f\big(e^{(\lambda,j-1)}, A^{(\lambda,j-1)}, 
Q^{(\lambda,j-1)} \big)`.

As in Case 2 we obtain  
:math:`A^{(\lambda,j)}` from :math:`A^{(\lambda,j-1)'}` 
by deleting row :math:`i` of  :math:`A^{(\lambda,j-1)'}`. We obtain 
:math:`Q^{(\lambda,j)}` from :math:`Q^{(\lambda,j-1)'}` by deleting 
row and column :math:`i`.  We put 
:math:`e^{(\lambda,j)} = e^{(\lambda,j-1)'}`.

A similar argument as in case 2 shows that matrices 
:math:`A^{(1,j)}` and  :math:`A^{(2,j)}` are as required and 
that deleting row :math:`i` from 
:math:`A^{(1,j-1')}` and  :math:`A^{(2,j-1')}` does not change
:math:`(g^{(1)} \odot g^{(2)})_j`.  


Remark

In case :math:`n_1 = n_2 = n` we may compute the product of
:math:`g^{(1)}` and :math:`g^{(2)}`  as follows:

.. math::
   g^{(1)} \cdot g^{(2)}  = (g^{(1)} \odot g^{(2)})_n
   = f\big( e^{(1,n)} \cdot e^{(2,n)}, A^{(1,n)}, 
   Q^{(1,n)} \odot    Q^{(2,n)} \big) \; . 
   
If the symmetric :math:`(m+1) \times (m+1)` bit matrices 
:math:`Q^{(\lambda)}, \lambda = 1,2`, represent the quadratic
functions 
:math:`q^{(\lambda)} : \mathbb{F}_2^m \rightarrow\mathbb{T}`,
then :math:`Q^{(1)} \odot  Q^{(2)}` is the symmtric
:math:`(m+1) \times (m+1)` bit matrix representing the
quadratic function :math:`q^{(1)} \cdot  q^{(2)}`.
The entries  :math:`Q^{(1 \odot 2)}_{i,j}` of 
:math:`Q^{(1)} \odot  Q^{(2)}` are given as follows:

.. math::
    \begin{array}{ll}
    Q^{(1 \odot 2)}_{i,j} = Q^{(1)}_{i,j} + Q^{(2)}_{i,j}
    & \quad \mbox{for} \quad i, j > 0 \; , \\ 
    Q^{(1 \odot 2)}_{i,0} =  Q^{(1 \odot 2)}_{0,i}
    = Q^{(1)}_{0,i} + Q^{(2)}_{0,i} + 
    Q^{(1)}_{i,i} \cdot  Q^{(2)}_{i,i} 
    & \quad \mbox{for} \quad i > 0 \; ,  \\ 
    Q^{(1 \odot 2)}_{0,0} =  0 \; ,
    \end{array}
  
   
   
The corrections in row and column :math:`0` of 
:math:`Q^{(1 \odot 2)}` are necessary, since the diagonal 
entries of :math:`Q^{(1)}`and :math:`Q^{(2)}` are to be 
interpreted modulo :math:`4`.
   
So our algorithm allows us to multiply guadratic mappings 
effectively. 

When we want to compute the expression
:math:`(g^{(1)} \odot g^{(2)})_j` defined in the last section,
then in the algorithm above there may be colums of matrices
:math:`A^{(\lambda)}, \lambda = 1, 2` that araise from 
extending the functions :math:`g^{(\lambda)}`. Here the extended
function :math:`g^{(\lambda)}` does not depend on the bit 
corresponding to such a column of matrix :math:`A^{(\lambda)}`. 
Expointing that independence leads to a considerable simplifiction 
of the aglorithm given above. We omit the details here.


Tensor contraction
..................

Let :math:`g^{(\lambda)}, \lambda = 1, 2` be as in the
last section. We also may consider
:math:`g^{(\lambda)}` as a tensor in the space
:math:`(\mathbb{C}^2)^{ \otimes j} \otimes 
(\mathbb{C}^2)^{\otimes n_\lambda-j}`. 
Then 
:math:`g^{(1)} \otimes g^{(2)} \in (\mathbb{C}^2)^{\otimes j}
\otimes  (\mathbb{C}^2)^{\otimes n_1 - j}
\otimes  (\mathbb{C}^2)^{\otimes j}
\otimes (\mathbb{C}^2)^{\otimes n_2 -j}`.
For the contraction :math:`(g^{(1)} \otimes g^{(2)})_n`
of :math:`g^{(1)} \otimes g^{(2)}` over
the two spaces :math:`(\mathbb{C}^2)^{\otimes j}` we have

.. math::
    (g^{(1)} \otimes g^{(2)})_j (x_1, x_2) =
    \sum_{x \in   \mathbb{F}_2^{j}}
    (g^{(1)} \odot g^{(2)})_n(x, x_1, x_2) \; , \quad 
    x_\lambda \in \mathbb{F}_2^{n_\lambda - j} \; . 

    
In the last section we have presented an effective algorithm for
computing :math:`h = (g^{(1)} \odot g^{(2)})_n` as a quadratic 
mapping :math:`h = f(e, A, Q)`.  Computing the sum in the 
expression given above corresponds to computing  
:math:`h' = f(e, A', Q)`, where :math:`A'` is obtiend from
:math:`A` by ropping the :math:`n` leftmost columns.  This means 
that we have to reduce the representation :math:`(e, A', Q)` of 
:math:`h'`.


Applying a (controlled) not gate to a quadratic mapping
.......................................................

In the theory of quantum computing we may apply so-called
*gates* to a quadratic state vector in 
:math:`(\mathbb{C}^2)^{\otimes n}`. For our puposes a gate
is a linear operation on  
:math:`(\mathbb{C}^2)^{\otimes n} = (\mathbb{C}^2)^{\otimes k}
\otimes  (\mathbb{C}^2)^{\otimes n-k}`
which may be written as a tensor product of a unitary
:math:`2^k \times 2^k` matrix `G` and a
:math:`2^{n-k} \times 2^{n-k}` identitiy matrix for a small
number :math:`1 \leq k \leq 2`. Here we may permute the factors
:math:`\mathbb{C}^2` of :math:`(\mathbb{C}^2)^{\otimes n}` 
arbitrarily before the decompostion into a tensor product as 
above.

A *not* gate operating in qubit :math:`j` maps a  state 
:math:`g` to a state :math:`g'` with
:math:`g'(x) = g(x + e_j)`, where
:math:`e_j = (0,\ldots,0,1,0,\ldots,0)` and  the component
:math:`1` is at position `j`. 
A *not* gate operating in qubit :math:`j` is implemented
for :math:`g = f(e,A,Q)` by flipping the bit :math:`A_{0,j}`.
The C function ``state12_gate_not`` implements a sequence of 
not gates.

A *controlled not* gate is a gate that negates a target qubit
:math:`j \neq j` controlled by a qubit :math:`j'`. Such a 
gate maps a state :math:`g` to a state :math:`g'` with
:math:`g'(x) = g(x + \langle e_{j'},x \rangle \cdot e_j)`,
where :math:`\langle .,. \rangle` is the scalar product 
of bit vectors. Such a gate is implemented for 
:math:`g = f(e,A,Q)` by adding column :math:`j'` of 
:math:`A` to  column :math:`j` of :math:`A`. The C function 
``state12_gate_ctrl_not`` implements a generalization of a
controlled not gate.


In the following two sections we discuss more types of gates.
Altogether, these gates generate the Clifford group 
:math:`\mathcal{X}_n` on :math:`n` qubits.

Applying a (controlled) phase gate to a quadratic mapping
.........................................................


Applying a phase :math:`\phi` gate  to qubit :math:`j` of
a state :math:`g= f(e,A,Q)` changes the state :math:`g` to 
a state :math:`g'` with 

.. math::
    g'(x_0,\ldots,x_j,\ldots,x_{n-1}) =
    \exp(\phi  x_j \sqrt{-1}) \cdot
    g(x_0,\ldots,x_j,\ldots,x_{n-1}) \; .
    
We consider only phase :math:`\phi` which are multiples of
:math:`\pi/2`. For an :math:`(m+1) \times n` matrix 
:math:`A` let :math:`A_j` be the :math:`j`-th column of
matrix :math:`A`. Let  :math:`A_{-1}` be the column vector 
:math:`(1,0,\ldots,0)^\top` with :math:`m+1` entries.
Then a phase :math:`\pi` gate on qubit :math:`j` maps
:math:`f(e,A,Q)` to 

.. math::
    f\left( (-1)^{A_{0,j}} \cdot e, A, Q \odot A_{-1} A_j ^\top
    \odot    A_j A_{-1}^\top \right) \; .
    
Here we consider :math:`A_j, A_{-1}` as a :math:`(m+1) \times 1`
matrices, so that the matrix product :math:`A_j A_{-1}^\top` 
is an :math:`(m+1) \times (m+1)` matrix. Operator :math:`\odot`
is as in section *Multiplication of quadratic mappings*.

A phase :math:`\pi/2` gate on qubit :math:`j` maps
:math:`f(e,A,Q)` to 
 
.. math::
    f\left( \sqrt{-1}^{A_{0,j}} \cdot e, A, Q \odot A_j A_j ^\top
    \right) \; .


Applying a controlled phase :math:`\pi` gate  to qubits
:math:`j` and :math:`j'` of a state :math:`g= f(e,A,Q)` changes 
the state :math:`g` to  a state :math:`g'` with 

.. math::
    g'(\ldots,x_j,\ldots,x_{j'},\ldots ) =
    (-1)^{x_j x_{j'}} \cdot
    g(\ldots,x_j,\ldots,x_{j'},\ldots) \; .
    
    
A conrtolled phase :math:`\pi` gate on qubit :math:`j` and
:math:`j'` maps  :math:`f(e,A,Q)` to     


.. math::
    f\left( (-1)^{A_{0,j} \cdot A_{0,j'}} \cdot e, A, 
    Q \odot A_j A_{j'}^\top \odot A_{j'} A_j^\top  \right) \; .


We leave the proofs of these statements to the reader.

Applying a Hadamard gate to a quadratic mapping
...............................................

A Hadamard gate at qubit :math:`j` is a a mapping that changes
a quadratic mapping :math:`g` to another quadratic mapping 
:math:`1/\sqrt{2} \cdot g'` with

.. math::
    g'(x_0,\ldots,x_{j-1},x_j,x_{j+1},\ldots,x_{n-1}) =
    g(x_0,\ldots,x_{j-1},0,x_{j+1},\ldots,x_{n-1}) + (-1)^{x_j} \cdot 
    g(x_0,\ldots,x_{j-1},1,x_{j+1},\ldots,x_{n-1}) \; .
    
    
We implement the application of a Hadamard gate on qubit :math:`j`
to a quadatic mapping :math:`g` represented as :math:`(e, A, Q)` 
as follows.

We append a zero row at :math:`A` and also a zero row and a
zero column at :math:`Q`. Let :math:`i` be the index of the 
appended row and column. Then we put 
:math:`Q_{i,k} = Q_{k,i} = A_{k,j}`, :math:`A_{k,j} = 0`
for all :math:`k \neq i`, and  :math:`A_{i,j} = 1`. Let 
:math:`A', Q'` be the modified matrices :math:`A', Q'`. 
Then :math:`g' = f(e, A', Q')`.

The correctenss of this algorithm can be seen as follows.
W.l.o.g we assume that :math:`j` is the last index :math:`n-1`.
Let :math:`x = (x_0,\ldots, x_j) \in \mathbb{F}_2^n`,
:math:`y  = (y_0,\ldots, y_m) \in   \mathbb{F}_2^{m+1}`
with :math:`y_0 = 1`,  and assume :math:`y \cdot A = x`
for the matrix product :math:`y \cdot A`. 
Then 

.. math::
   (y, b) \cdot A' = (x_0,\ldots,x_{n-1},b) \quad 
   \mbox{for} \quad b \in \mathbb{F}_2 \; .

Let :math:`q, q'` be the quadratic mappings given by :math:`Q, Q'`.
Then 

.. math::
    q'(y, b) = (-1)^{b \cdot \langle y, A_j \rangle} \cdot q(y)
    = (-1)^{b \cdot x_j} \cdot q(y) \; ,
    
where :math:`A_j` is the :math:`j`-th column of :math:`A`. 
Thus

.. math::
   f(e, A',Q')(x_0, \ldots, x_{j-1}, x_j) =
   g(x_0, \ldots, x_{j-1}, 0) + 
   (-1)^{x_j} \cdot g(x_0, \ldots, x_{j-1}, 1) \; .


The C function ``state12_gate_h`` implements a sequence of 
Hadamard gates.

C functions dealing with quadratic state vectors
................................................


The C functions in this module do operations on quadratic state
vectors given by triples :math:`(e, A, Q)` as defined above.
Here component :math:`e` encodes the number  
:math:`\exp(e \pi \sqrt{-1} / 4) \cdot 
2^{\lfloor e/16 \rfloor / 2}`, and
:math:`A` is an  :math:`(1+m) \times n` bit matrix.
:math:`Q` is a symmetric :math:`(1+m) \times (1+m)` bit matrix 
representing an symmetric bilinear form. We always have
:math:`Q_{0,0}=0`. Put :math:`m'=1+m`. Matrices :math:`A` and 
:math:`Q` are concatenated to an :math:`m' \times (n+m')` matrix
:math:`M` with :math:`M_{i,j} = A_{i,j}` for :math:`j < n` and
:math:`M_{i,j} = Q_{i-n,j}` for :math:`j \geq n`. Matrix
:math:`M` is encoded in a one-dimensional array of unsigned
64-bit integers. Here bit :math:`j` of entry :math:`i` 
corresponds to :math:`M_{i,j}`, with bit :math:`0` the least
significant bit.

A quadratic state vector is described by a structure containing 
the following components:

.. code-block:: c

  typedef struct {
    uint32_t maxrows; // No of entries allocated to component data
    uint32_t nrows;   // No m' = m + 1 of rows of bit matrix A
    uint32_t ncols;   // No n of columns of bit matrices A and Q
    int32_t  factor;  // A number e encoding a scaling factor
    uint64_t *data;   // Pointer to the data bits of matrix M
    uint32_t shape1;  // (Will be discussed in the next session)
  } qbstate12_type;

Recall that a quadratic state vector :math:`v` of type
``qbstate12_type`` with component ``ncols = n`` models a complex 
vector in a vector space  :math:`V` of dimension :math:`2^n`, and 
that the basis of ``V`` is labelled by the elements of the Boolean
vector space :math:`\mathbb{F}_2^n`. In C and python programs
we represent the element :math:`(x_0, \ldots, x_{n-1})` of
:math:`\mathbb{F}_2^n` by the integer 
:math:`\sum_{0 \leq i < n} 2^i \cdot x_i`. This leads to a natural
representation of ``v`` as a one-dimensional complex array of
length :math:`2^n`, starting with index ``0``.

The zero state is encoded as a matrix with :math:`m'=0` rows.
We do not update the entries :math:`Q_{i,0}`, so the 
corresponding bits in compoment ``data`` of the structure
are garbage. One may use the C function ``qstate12_check`` to
set these bits to their proper values.


The current implementation requires ``n + m <= 63``.  
This can easily be generalized to larger Clifford 
groups by reserving an array of several integers for each row 
of matrix :math:`M`. Here we also leave the details to the reader.

C functions supporting this module are prefixed with ``qbstate12_``.
Unless otherwise stated, these functions return an ``int32_t``, 
where a nonegative value is interpreted as success, and a negative 
value is intepreted as failure. Depending on the function, a 
nonnegative return value may e.g. mean an index for a matrix
:math:`A`, :math:`M`, or :math:`Q`.

Typical names for parameters of functions in this module are:

   ================== ================================================
   ``pqs, pqs1, ...`` Pointer to structure of type ``qbstate12_type``
   ``nqb``            Number of qubits, i.e. of columns of matrix 
                      :math:`A`.
   ``nrows``          Number of rows of matrix :math:`A`, :math:`M`, 
                      and :math:`Q`.
   ``i, i1, ...``     Index of a row of matrix :math:`A`, :math:`M`,  
                      or and :math:`Q`, starting with 0.
   ``j, j1, ...``     Index of a column of matrix :math:`A`, with a 
                      column of :math:`A`, corrsesponding to a qubit, 
                      starting with ``j = 0``.
                      If appropriate, an index  ``j >= ncols`` refers 
                      to column ``(j - ncols)`` of matrix  :math:`Q`.
   ``pv, pv1,...``    Pointer to a row or column vector of matrix 
                      :math:`A`, :math:`M`, or  :math:`Q`.
   ================== ================================================


Quadratic state matrices
........................

A *quadratic state matrix* :math:`S` of shape :math:`(n_0, n_1)` 
is an element of the tensor product :math:`V_0 \otimes V_1` 
with the basis vectors of :math:`V_k` being indexed by 
:math:`\mathbb{F}_2^{n_k}, k = 0, 1`, such that the
coordinate function of :math:`S` is a quadratic mapping.
That coordinate function is a function
:math:`\mathbb{F}_2^{n_0} \times \mathbb{F}_2^{n_1}
\rightarrow \mathbb{C}`.

Thus :math:`S` corresponds to a complex 
:math:`2^{n_0} \times 2^{n_1}` matrix. We may implement
:math:`S` as a quadratic state vector of :math:`n_0 + n_1` 
qubits, augmented by an information about its shape
:math:`(n_0, n_1)` . We let the  :math:`n_0` qubits with
high indices  correspond to the rows of  :math:`S`; and
we let the  :math:`n_1` qubits with low indices  correspond 
to the columns  of  :math:`S`. For implementing a quadratic
state matrix in C, we store the part :math:`n_1` of its shape 
in component ``qstate1`` of a structure ``qs`` of type
``qstate12_type`` and we compute the part  :math:`n_0`  as 
``qs.ncols - qs.qstate1``.

In python we implement a *quadratic state matrix* as a instance of 
class ``QStateMatrix`` in module ``mmgroup.structures.qs_matrix``.
A matrix of shape :math:`(0,n)` corresponds to a row vector of 
dimension :math:`2^{n}` and a matrix of shape :math:`(n,0)` 
corresponds to a column vector of dimension :math:`2^{n}`.
We have seen above that the unitary quadratic state 
matrices of shape :math:`(n,n)` form the Clifford group 
:math:`\mathcal{X}_{n}`. We have also
discussed fast algorithms for multiplication and inversion
of such matrices. So class ``QStateMatrix`` supports fast
computation in Clifford group :math:`\mathcal{X}_{n}`.
Our implementation requires :math:`n \leq 12`, which is
sufficient for computing in the subgroup
:math:`2^{1+24}.\mbox{Co}_1` of the monster group.

Function ``qstate12_matmul`` in file ``qmatrix12.c`` multiplies
two quadratic state matrices.


Reducing a quadratic state matrix
.................................

We define a special *reduced matrix representation* :math:`(e, A, Q)` 
of a quadratic state matrix :math:`S` of shape :math:`(n_0, n_1)` 
that differs slightly from the reduced representation of a 
quadratic state vector. That representation satisfies the 
following conditions:

 1. We require that the representation :math:`(e, A, Q)` is in 
    (not necessarly reduced) echelon form as described in section
    *Reducing the representation of a quadratic mapping*.
    
 2. Let :math:`A'` be the bit matrix obtained from :math:`A` by
    removing the first  :math:`n_1` columns from  :math:`A`.
    We require that  a permutation of the rows of matrix
    :math:`A'`, excluding row :math:`0`, is in (not necessarily 
    reduced) echelon form.
    
    Note that the first :math:`n_1` columns of the bit matrix
    :math:`A` correspond to the :math:`2^{n_1}` columns of 
    the complex matrix :math:`S`; and the remaining  
    :math:`n_0` columns of :math:`A` correspond to the
    :math:`2^{n_0}` rows of :math:`S`.
        
 3. Let :math:`I_0` be the set of the rows of the bit matrix 
    :math:`A` such that all bits in columns :math:`\geq n_1` 
    of that row are zero. Let :math:`I_1` be the set of the 
    rows of :math:`A` such that all bits in columns 
    :math:`< n_1` of that row are zero. We exclude row 
    :math:`0` from :math:`I_0` and from :math:`I_1`. 
    
    If the bit matrix :math:`A` has :math:`m'` rows then
    the bit matrix :math:`Q` is a symmetric 
    :math:`m' \times m'` bit matrix. Let :math:`Q'` be the 
    submatrix of :math:`Q` that consists of all rows with
    index in :math:`I_0` and of all columns with index in
    :math:`I_1` .

    We require that submatrix :math:`Q'` of :math:`Q`  has at 
    most one nonzero entry in each row and in each column. 

    Note that  :math:`I_0 \cap I_1 = \emptyset`, since
    :math:`A` is in echelon form.

We may obtain a *reduced matrix representation* :math:`(e, A, Q)` 
of a quadratic state matrix :math:`S` of shape :math:`(n_0, n_1)` 
as follows:

  * Starting from a reduced representation :math:`(e_0, A_0, Q_0)` 
    of a quadratic state vector :math:`S` we may obtain a 
    representation :math:`(e_1, A_1, Q_1)` of the matrix  :math:`S` 
    satifying condition *(2.)* by applying a seqence of  
    transformations :math:`T_{i,j}, i < j`. Then 
    :math:`A_1` is also in echelon form. :math:`T_{i,j}` 
    is defined in section *Implementation of quadratic mappings*. 

  * We apply a sequence of transformations :math:`T_{i,j}, i < j`,
    with :math:`i, j \in I_0` or :math:`i, j \in I_1` to 
    :math:`(e_1, A_1, Q_1)`. These transformations preserve the 
    echelon form of :math:`A_1` and also property *(2.)*. Since 
    :math:`I_0 \cap I_1 = \emptyset`, these transformations act as
    row and column operations on the submatrix :math:`S'`
    of  :math:`S`. So we may achieve property *(3.)* by a sequence
    of such transformations, thus obtaining a suitable
    representation :math:`(e, A, Q)` of  :math:`S`.



Using a *reduced matrix representation* :math:`(e, A, Q)` of a
quadratic state matrix :math:`S` has a variety of advantages.

 * We can easily compute the rank of :math:`S` as follows:

   For :math:`(e, A, Q)` let :math:`I_0, I_1, S'` be defined as 
   above. Let  :math:`I_2` be the subset :math:`I_0` containing 
   alls rows of :math:`A` such that the corrsponding row of 
   bit matrix  :math:`Q'` is zero. Then the binary logarithm of
   the rank of :math:`S` is equal to the number of rows of
   matrix :math:`A` with are neither in  :math:`I_1` nor in 
   :math:`I_2`. Here we have to exclude row :math:`0` of :math:`A`.  

   We omit the proof of this fact since we do not need it
   for our purposes.

 * That representation can be used for decomposing :math:`S`
   into a product :math:`M_1 \cdot H \cdot M_2` of quadratic 
   state matrices, where :math:`M_1, M_2` are monomial and 
   :math:`H` is a Hadamard-like matrix. By a Hadamard-like
   matrix we mean a matrix obtained from the unit matrix by
   applying Hadamard gates to the row vectors.

   If :math:`S` has shape :math:`(n,n)` then such a decomposition
   reduces the complexity of multiplying :math:`S` with an
   arbitrary complex vector form :math:`O(4^n)` to
   :math:`O(n \cdot 2^n)`.

   Essentially, we have used such a decomposition for 
   multiplying the non-monimial part of generator :math:`\xi` 
   of the monster :math:`\mathbb{M}` with a vector of our 
   representation of :math:`\mathbb{M}`. Since this special
   case is discussed elsewhere, we do not go into details here.

 * In the next section we will introduce the Pauli group 
   :math:`\mathcal{P}_{n}`, which is an important normal 
   subgroup of :math:`\mathcal{X}_{n}`. Using the reduced 
   matrix representation of an element :math:`S` of 
   :math:`\mathcal{X}_{n}` we can conugate any element
   of :math:`\mathcal{P}_{n}` with :math:`S` in  
   :math:`O(n^2)` bit operations.
   
   With quadratic state matrices, a general matrix 
   multiplication in the Clifford group :math:`\mathcal{X}_{n}` 
   costs :math:`O(n^3)` bit operations.

   

Function ``qstate12_reduce_matrix`` in file ``qmatrix12.c``
converts any representation of a quadratic state matrix to a
reduced matrix representation.

The Pauli group
...............

The *Pauli group* :math:`\mathcal{P}_{n}` of :math:`n` qubits
is the normal subgroup of the Clifford group 
:math:`\mathcal{X}_{n}` generated by the not gates, the phase 
:math:`\pi` gates in :math:`\mathcal{X}_{n}`, and by the scalar 
multiples of the unit matrix by a fourth root of unity. It has 
structure :math:`\frac{1}{2}(2_+^{1+2n} \times Z_4)`, exponent
math:`4` and order :math:`2^{2n+2}`.
 
We represent an element of :math:`\mathcal{P}_{n}` as a product
of :math:`2n+2` generators. Each generator may have exponent 
:math:`0` or :math:`1`. The sequence of these exponents are   
stored as a bit vector as follows:

  * Bit :math:`2n+1` corresponds to multiplication with 
    the scalar :math:`\sqrt{-1}`.

  * Bit :math:`2n` corresponds to  multiplication with
    the scalar :math:`-1`.

  * Bit :math:`n+i` with :math:`0 \leq i < n` corresponds to 
    a not gate applied to qubit :math:`i`.

  * Bit :math:`i` with :math:`0 \leq i < n` corresponds to 
    a phase :math:`\pi` gate applied to qubit :math:`i`.


Factors are ordered by bit positions, with the most significant 
bit position occuring first. In the C language we represent bit 
vectors a integers as usual.

All generators commute and have order :math:`2` except for the
following cases:

  * A phase :math:`\pi` gate anticommutes with a not gate
    applied to the same qubit, i.e their commutator is the
    scalar  :math:`-1`.

  * Of course, the scalar :math:`\sqrt{-1}`
    squares to the scalar :math:`-1`.


Functions ``qstate12_pauli_vector_mul`` and 
``qstate12_pauli_vector_exp`` in module ``qs_matrix.c`` 
perform multiplication and exponentiation in the Pauli group
:math:`\mathcal{P}_{n}`.


Given an element :math:`p` of the Pauli group and a
*reduced matrix representation* :math:`(e, A, Q)` of a
unitary quadratic state matrix :math:`S`, we can quickly 
compute the conjugate  :math:`S^{-1} \cdot p \cdot S` as 
follows:

  * Left multiply :math:`S` with :math:`p` by applying the 
    appropriate gates to :math:`S`. This affects row 
    :math:`0` of the bit matrix :math:`A` and  row and 
    column :math:`0` of the symmetic bit  matrix :math:`Q` 
    only.

  * Restore the original values  :math:`A[0,j]` or all
    :math:`j \geq n` and the original values 
    :math:`Q[0,k] = Q[k,0]` for all :math:`k > n` by applying 
    appropriate transformations :math:`T_{0,j}` to the 
    modified representation :math:`(e, A, Q)`. 
    This does not change the value :math:`p \cdot S` of 
    the complex matrix computed in the previous step.
    :math:`T_{0,j}` is defined in section
    *Implementation of quadratic mappings*.

  * We may restore the the remaining original values of row 
    and column :math:`0` of :math:`Q` by applying phase 
    :math:`\pi` gates to :math:`S`. These gate operations 
    correspond to a right multiplication with a element 
    :math:`p_1` of the Pauli group.

  * We may restore the the remaining original values of row 
    :math:`0` of :math:`A` by applying not gates to
    :math:`S`. These gate operations correspond to a right
    multiplication with a element :math:`p_2` of the Pauli
    group.

  * Up to a known scalar factor we have obtained an equation
    :math:`S = p  \cdot S \cdot p_2 \cdot p_1`, with 
    :math:`p_1, p_2` in the  Pauli group. With this equation
    the requested conjugation is an easy computation in the 
    Pauli group.

Function ``qstate12_pauli_conjugate`` in module 
``qs_matrix.c`` performs this conjugation.


"""



import sys

class QSstate_tables:
    directives = {}
    def __init__(self):
        self.tables = {
            "QSTATE_DOC":  sys.modules[__name__],
        }



