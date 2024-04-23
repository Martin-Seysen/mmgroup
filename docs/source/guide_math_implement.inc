

.. _implement-autpl-label:

Implementing  Automorphisms of the Parker loop
----------------------------------------------
 
A standard automorphism :math:`\pi` of the Parker loop 
:math:`\mathcal{P}` is implemented as an array of ``12`` integers 
of type ``uint32_t``. The lower ``13`` bits of the :math:`i`-th entry 
of that array contain the image of the :math:`i`-th basis element 
:math:`(b_i,0)` of the Parker loop, where :math:`b_i` is the 
:math:`i`-th basis vector of the Golay code :math:`\mathcal{C}`. These 
images describe the automorphism :math:`\pi` uniquely. Each element of
:math:`\mathcal{P}` has a unique representation as a tuple
:math:`d = (d, \lambda)`, :math:`d \in \mathcal{C}`, 
:math:`\lambda \in \mathbb{F}_2`.


A key task is to compute the image :math:`(d, \lambda) \cdot \pi` of 
an element  :math:`(d, \lambda)` of the Parker loop under the 
automorphism :math:`\pi`.
By Lemma 4.1 in :cite:`Seysen20` there is a quadratic form
:math:`q_\pi` with associated bilinear form :math:`\theta_\pi`
on the Golay code satisfying

.. math::
   \theta_\pi  \, = \, \theta^\pi + \theta \;, \quad
   \mbox{where} \quad \theta^\pi(d,e) = \theta(d^\pi,e^\pi) \;, \\
   (d, \lambda)^\pi = (d^{\pi}, \lambda + q_\pi(d)) \; ,	 


for any :math:`(d, \lambda) \in  \mathcal{P}`. So the image
:math:`(d, \lambda) \cdot \pi` can be computed if :math:`q_\pi(d)`
can be computed. A functional value of :math:`q_\pi(d)` can easily be 
computed from the associated bilinear from :math:`\theta_\pi`` if 
the values :math:`q_\pi(b_i)` are known for all basis vectors 
:math:`b_i`. The values :math:`q_\pi(b_i)` can be computed 
from the sign bits of the images if the basis vectors.

So it suffices to compute the bilinear form :math:`\theta_\pi` as a  
:math:`12 \times 12` bit matrix. The first term :math:`\theta` of 
:math:`\theta_\pi` does not depend on :math:`\pi` and can be stored 
as a constant. The second term  :math:`\theta^\pi` of that matrix 
is computed as follows:

Let :math:`c_i` = :math:`(b_i)^\pi`, and let :math:`C` be the matrix 
with entries :math:`c_{i,j}` such that :math:`c_i` = 
:math:`\sum_j c_{i,j} b_j`. So :math:`C` is just the matrix of the 
images of the basis vectors, ignoring the signs of the Parker loop.

Row  :math:`i`, column :math:`j` of matrix :math:`\theta^\pi` is equal 
to  :math:`\theta(c_i,c_j)`.  Since :math:`\theta(c_i,c_j)` is linear
in its second argument, we have 

.. math::
   \theta^\pi(c_i, c_j) \, = \, \sum_k c_{j,k}  \theta(c_i,b_k) \; ,\\
   \theta^\pi(c_i)  \, = \,  \theta(c_i) \cdot C^\top \; .

Here  :math:`\theta^\pi(c_i)` is just the :math:`i`-th row of the bit 
matrix :math:`\theta^\pi`. We also store the values  :math:`\theta(d)` 
for all :math:`d \in \mathcal{C}` in a table. Thus the bilinear form 
:math:`\theta_\pi` can essentially be computed as a product of two
:math:`12 \times 12` bit matrices. 
 
We store the lower triangular submatrix of the bit matrix 
:math:`\theta_\pi` in bits ``13, ..., 24`` of the array representing
the automorphism :math:`\pi` of the Parker loop.


.. _implement-gen-mm:


Implementing generators of the Monster group
--------------------------------------------


The operation :math:`y_{f} \cdot x_{e} \cdot x_\epsilon`
........................................................

The operation  :math:`g = y_{f} \cdot x_{e} \cdot x_\epsilon` on
:math:`\rho_p` is coded in function ``mm_op<p>_xy`` in the
automatically generated C file ``mm<p>_op_xy.c`` for the modulus ``p``.
The C function performs that operation for arbitrary
:math:`e, f \in \mathcal{P}` and :math:`\epsilon \in \mathcal{C}^*` 
in a single pass.
That operation is monomial and the formula for it can easily be
deduced from the information in :cite:`Seysen20`.

For  :math:`v \in \rho_p` 
the C program calculates the components in :math:`V_A, V_B, V_C \ldots` 
of the result  :math:`v \cdot g` in their natural order. For 
implementing such a program it is useful to have a formula for the 
operation of 
:math:`g^{-1} =  x_\epsilon \cdot x_{\bar{e}} \cdot y_{\bar{f}}`.
In the sequel we state the operation of  :math:`g^{-1}` on
:math:`\rho_p`.


In this subsection we use the operator ':math:`\oplus`' for the 
*programmer's multiplication* 
in the Parker loop :math:`\mathcal{P}`, which we define by:

.. math::
  (\tilde{d_1}, \lambda_1 ) \oplus  (\tilde{d_2}, \lambda_2) =
  (\tilde{d_2} + \tilde{d_2} , \,
  \lambda_1 + \lambda_2) \; , \qquad
  \tilde{d_1}, \tilde{d_2}
  \in \mathcal{C}, \, \lambda_1,\lambda_2 \in \mathbb{F}_2 \; .


Then ':math:`\oplus`' is a simple XOR operation on a computer; and
the standard product :math:`d_1 d_2` in the Parker loop is given by
:math:`d_1 d_2 = (-1)^{\theta(d_1, d_2)} d_1 \oplus d_2` for 
:math:`d_1, d_2 \in \mathcal{P}`. 
Note that ':math:`\oplus`'  depends on the selected 
cocycle :math:`\theta`. We may store 
the cocycles :math:`\theta(d)` in a table, with one 12-bit entry 
representing the element
:math:`\theta(d Z(\mathcal{P}))`
of the Golay cocode for each of the 2048 cosets :math:`d Z(\mathcal{P})` of 
:math:`Z(\mathcal{P})  = \{\pm 1, \pm \Omega\}` in :math:`\mathcal{P}`.
So multiplication in :math:`\mathcal{P}` is easy.


Put  
:math:`d^{[0]} = d^+, d^{[1]} = d^-`, and
:math:`X^+_{d,\delta} = X^{\vphantom{+}}_{\Omega d,\delta}` if
:math:`x^{\vphantom{+}}_{\Omega d,\delta}` is short,  
:math:`X^+_{d,\delta} = X^{\vphantom{+}}_{d,\delta}`
otherwise. Then

.. math::
  X_{d,i}^+ &  \stackrel{g^{-1}}{\longrightarrow} 
   (-1)^{s_{{X}\vphantom{X^x}}} \cdot  X_{d\oplus f,i}^+ \, ,
  \quad \mbox{with} \; \\
  {s_X} & = P(f) + P(ef) + (|\epsilon| + 1) P(d) + P(def) + \left< e, i\right> + 
  \big< d, i^{|\epsilon|} \epsilon  A(e,f) \theta(f) \big>  \; ; \\
  (d^{[\tau]} \otimes_1 i) &  \stackrel{g^{-1}}{\longrightarrow}   
   (-1)^{s_{{Y\!Z}\vphantom{X^x}}} \cdot 
  d^{[\sigma]}_{d \oplus  e \oplus  f^{\sigma + 1}}  \otimes_1 i \, ,
  \quad \mbox{with} \; \;  \sigma = \tau +  |\epsilon|  \pmod{2} \; , \\  
  {s_{Y\!Z}} & =  (\sigma+1)\theta(f,e) + \sigma P(f) + P(de) + P(def)
  +  \left<f, i\right>  +
  \left< d, \epsilon \theta(e) \theta(f)^{\sigma+1} \right>    \; ; \\
  X^+_{d \cdot \delta} &  \stackrel{g^{-1}}{\longrightarrow}  
     (-1)^{s_{{T}\vphantom{X^x}}} \cdot  X^+_{d \cdot \delta \delta'} \, , 
  \quad  \mbox{for} \; |d| = 8, \; \; \delta \mbox{ even} \, ,
  \quad \mbox{with}  \; \delta' = A(d,f) \, , \; \;   \mbox{and} \\
  s_T &=  P(e) + P(de) +  \left<  d, \epsilon \right> + 
    \left<ef,\delta\right> + |\delta||\epsilon|/2  \; ; \\
  X_{\Omega^m \cdot ij}  &  \stackrel{g^{-1}}{\longrightarrow}
    (-1)^{ m |\epsilon| + \left< ef, ij \right> }
  \cdot   X_{\Omega^n \cdot ij}  \; ,
  \quad \mbox{with} \; \;  n = m + \left<f, ij \right> \, ;  \\  
  (ij)_1    &  \stackrel{g^{-1}}{\longrightarrow} 
    (-1)^{\left< f, ij \right>} (ij)_1  \; .

The sign bits :math:`s_X, s_{YZ},` and :math:`s_T` can also be
computed as follows:

.. math::

  s_X & =
  |\epsilon| \left<d, i \right>
  + \left< e, i \right>
  + |\epsilon| P(d)
  + \left<ef, \theta(d)  \right>
  + \left<d, \epsilon \theta(e) \right> + P(f)
  \, ;  \\
  s_{YZ} & =
  \left<f,  i \right>
  + \left<f, \theta(d) \right>
  + \left<d, \epsilon \theta(ef) \theta(f)^{\sigma + 1} \right> \\
  & + \left<f,  \theta(e) \right>
  + \sigma \left<e, \theta(f) \right>
  + (\sigma + 1) P(f)
  \, ;  \\
  s_{T} & =
  |\delta||\epsilon|/2
  + \left<ef,\delta\right>
  + \left<e, \theta(d)  \right>
  + \left<d,  \epsilon \theta(e)  \right>
  \; .



The operation of the triality elements  :math:`\tau^e`
.......................................................

The operation  :math:`g = \tau^e` on :math:`\rho_p` is coded in
function ``mm_op<p>_t`` in the automatically generated C file 
``mm<p>_op_t.c`` for the modulus ``p``. Since  :math:`\tau` has
order 3, it suffices to present the formulas for :math:`e = \pm 1`.
Note that :math:`\tau` fixes  :math:`(ii)_1`. From the discussion
in :cite:`Seysen20`, Section 8, we obtain for :math:`i \neq j`:

.. math::
   \begin{array}{ccccccc}
  (ij)_1  
  & \stackrel{\tau}{\longrightarrow} &
  X_{ij} -  X^+_{ij}
  & \stackrel{\tau}{\longrightarrow} &
  X_{ij} +  X^+_{ij} 
  & \stackrel{\tau}{\longrightarrow} & 
  (ij)_1  \, , \\
  2 X_{ij}  
  & \stackrel{\tau}{\longrightarrow} &
  (ij)_1 + X_{ij} +  X^+_{ij}
  & \stackrel{\tau}{\longrightarrow} &
  (ij)_1 + X_{ij} -  X^+_{ij} 
  & \stackrel{\tau}{\longrightarrow} & 
  2 X_{ij} \, , \\
  2 X^+_{ij}   
  & \stackrel{\tau}{\longrightarrow} &
  (ij)_1 - X_{ij} -  X^+_{ij}
  & \stackrel{\tau}{\longrightarrow} &
  -(ij)_1 + X_{ij} -  X^+_{ij} 
  & \stackrel{\tau}{\longrightarrow} & 
  2 X^+_{ij} \, , \\
  X^+_{d \cdot i}
  & \stackrel{\tau}{\longrightarrow} &
  (-1)^{\left<d, i \right>} (d^- \otimes_1 i)
  & \stackrel{\tau}{\longrightarrow} &
  (-1)^{P(d)} (d^+ \otimes_1 i)
  & \stackrel{\tau}{\longrightarrow} &
  X^+_{d \cdot i} \, , \\
  X^+_{d,\delta}
  & \stackrel{\tau}{\longrightarrow} &
  \textstyle\frac{1}{8} 
  \sum_\epsilon \textstyle{(-)_{\delta,\epsilon}} X^+_{d,\delta}
  & \stackrel{\tau}{\longrightarrow} &
  \textstyle\frac{1}{8} 
  \sum_\epsilon \textstyle{(-)_{\epsilon,\delta}} X^+_{d,\delta}
  & \stackrel{\tau}{\longrightarrow} &
  X^+_{d,\delta} \, .
  \end{array}

The last formula describes the operation of :math:`\tau` on
:math:`X^+_{d,\delta}`, where :math:`d` is a preimage of an octad
in the Parker loop, and :math:`\delta` is an even element of the
Golay cocode represented by a subset of that octad. The summation
runs over all even subsets :math:`\epsilon` of that octad; and
:math:`(-)_{\delta,\epsilon}` means
:math:`(-1)^{\left|\delta \cap \epsilon \right| + \left|\delta \right| /2}`
as in :cite:`Seysen20`. We will briefly comment on the implementation
of the last formula.

For any octad :math:`d` as above we consider the row vector
:math:`X^+_{d}` with 64 entries :math:`X^+_{d, \delta}`, where
:math:`\delta` runs over the 64 even subsets of the octad corresponding
to  :math:`d`. Then operation of :math:`\tau` on  :math:`X^+_{d}` can 
be written as
   
.. math::
     X^+_{d} 
     \stackrel{\tau}{\longrightarrow} 
     \textstyle\frac{1}{8} X^+_{d} \cdot D \cdot H
     \stackrel{\tau}{\longrightarrow} 
     \textstyle\frac{1}{8}  X^+_{d} \cdot H \cdot D
     \stackrel{\tau}{\longrightarrow} 
     X^+_{d} \; .

Here matrix :math:`D` is a diagonal matrix with an entry
:math:`(-1)^{\left|\delta \right| /2}` in the row and column
corrsponding to :math:`\delta`. Matrix :math:`H` is a Hadamard matrix
with an entry :math:`(-1)^{\left|\delta \cap \epsilon \right|}` in
the row and column corresponding to :math:`\delta` and
:math:`\epsilon`, respectively.
    
When implementing the muliplication with matrices :math:`D` and
:math:`H` we have to be aware of the numbering of the entries 
:math:`\delta` used in Section :ref:`mmrep-label`. Using that
numbering, it turns out that matrix :math:`D` has diagonal entries
:math:`(-1)^{w(i(i+1)/2)}`, with :math:`i` running from 0 to 63,
and :math:`w(j)` the bit weight of the binary representation
of a nonnegative integer :math:`j`. Matrix :math:`H` is a
(commuting) product of a permutation matrix :math:`P`, and a
standard :math:`64 \times 64` Hadamard matrix :math:`H_{64}`
according to Sylvester's construction.
Matrix :math:`P` fixes entry :math:`i` if  :math:`w(i)` is
even and exchanges entry :math:`i` with entry :math:`63 - i`
otherwise.  

Here multiplication with matrix :math:`H_{64}` is the most
time-consuming operation. We can speed this up by using the
decomposition:

.. math::
    H_{64} = H_2 \otimes H_2 \otimes H_2 \otimes H_2
    \otimes H_2 \otimes H_2 \, , \quad \mbox{with} \;
    H_2 = \left(
    \begin{array}{cc}
    1 & 1 \\
    1 & -1
    \end{array}
    \right) \, .


Monomial operation of the generators  :math:`\xi^e`
.......................................................

Generator :math:`\xi` is an element of the subgroup :math:`G_{x0}`
of the monster. The basis vectors with tags ``B, C, T, X`` (as
explained below) of the representation  :math:`\rho` of the
Monster span a subspace :math:`98280_x` of :math:`\rho`. The basis
vectors with the tags listed above (together with their opposite
vectors) correspond to the short elements of the normal subgroup
:math:`Q_{x0}` of :math:`G_{x0}` of structure :math:`2^{1+24}`.
The group :math:`G_{x0}` operates on these basis vectors in the
same way as it operations on :math:`Q_{x0}` by conjugation. Thus
:math:`G_{x0}` (and hence also :math:`\xi^{\pm 1}`) operates
monomially on the basis vectors of :math:`98280_x`.  

We have :math:`Q_{x0} / \{\pm 1\} \cong \Lambda / 2 \Lambda`. An
element  :math:`x` of :math:`Q_{x0}` is short if it corresponds to
a vector in the Leech lattice :math:`\Lambda` of norm 4. The 
operation of :math:`\xi^{\pm 1}` on  :math:`Q_{x0}` (and hence also
of :math:`98280_x`) is given by Lemma 9.5 in :cite:`Seysen20`.

In  section :ref:`mmrep-label` we use the following names and tags 
for the basis vectors of :math:`98280_x`.


    ================ === =================================== ========
    Name             Tag Entries                             Remarks
    ================ === =================================== ========
    :math:`X^+_{ij}`  B  i, j;  0 <= j < i < 24              (1)
    :math:`X_{ij}`    C  i, j;  0 <= j < i < 24              (1)
    :math:`X_{o,s}`   T  o, s;  0 <= o < 759, 0 <= s < 64    (2)
    :math:`X_{d,j}`   X  d, j;  0 <= d < 2**11, 0 <= j < 24  (1),(3)
    ================ === =================================== ========

Remarks

(1)  i and j, 0 <= i,j < 24  refer to basis vectors of the
     Boolean vector space in which the Golay code is defined.
(2)  o is one of 759 octads, s is one of 64 even subsets of
     octad d (modulo its complement in d), as described in
     section :ref:`octads_label`.
(3)  d with 0 <= d < 2048 refers to the Parker loop element d. 
     

We group these basis vectrors into five boxes (labelled 1, ..., 5)
with each box containing at most 24576 entries.
Element :math:`\xi` permutes these boxes as follows:

.. math::
  
   \mbox{Box1} \rightarrow \mbox{Box1} \, ,  \; 
   \mbox{Box2} \rightarrow \mbox{Box2} \, ,  \; 
   \mbox{Box3} \rightarrow \mbox{Box4} \rightarrow    
   \mbox{Box5} \rightarrow \mbox{Box3} \, . 


The mapping from the basis vectors to entries in boxes is:

   ===================== ===     =======================
   Basis vector          Box     Entry
   ===================== ===     =======================
   B[i,j]                 1         0 +  32 * i + j
   C[i,j]                 1       768 +  32 * i + j
   T[o,s], o < 15         1      1536 +  64 * o + s             
   T[o,s], 15 <= o < 375  2        64 * (o - 15) + s           
   T[o,s], o >=  375      3        64 * (o - 375) + s
   X[d,j], d <  1024      4        32 *  d + j
   X[d,j], d >= 1024      5        32 * (d - 1024) + j
   ===================== ===     =======================

This subdivision looks weird, but is has quite a few advantages:

  * The lower index (j, or s) has stride 1, and the stride of the
    higher index (i, o, or d) is the lowest possible power of two. So 
    accessing an entry is easy. In the C code for a representation
    derived from :math:`\rho` the entries of a vector of such a 
    representation are strided in the same way. So the tables 
    computed under these assumptions can be used in the C code
    for the operator :math:`\xi`.

  * Boxes are permuted as above, so the operation of
    :math:`\xi^{\pm 1}` can be stored in a set ot 5 tables
    encoding mappings between boxes.

A similar, but finer subdivision of the whole space :math:`\rho`
(and not only of the subspace :math:`98280_x`) is given in 
:cite:`Iva09`, section 3.4.

We further subdivide a box into *clusters*, *rows*, and *columns*.
Here a row always contains 32 columns; with either all or just the
first 24 columns being used. A cluster contains a variable nubmer of
rows, and a box contains a variable numer of clusters, as indicated
in the  following table. We also assign names to the boxes.


 === ==== ======= ======== ====== ============
 Box Name Size    Clusters Rows   Columns used
 === ==== ======= ======== ====== ============
  1   BC     2486        1     78           32
  2   T0    23040       45     16           32 
  3   T1    24576       64     12           32 
  4   X0    24576       64     16           24 
  5   X1    24576       64     16           24 
 === ==== ======= ======== ====== ============

If a box has :math:`c` clusters, :math:`m` rows, and :math:`n` used
columns, then the triple :math:`(c, m, n)` is called the *shape*
of the box. 

It turns out that :math:`\xi` always maps the :math:`i`-th cluster
of a box to the :math:`i`-th cluster of a (possibly different) box.
This greatly simplifies the implementation of the mapping. So it
suffices to specify the mapping as a collection of permutations inside
a cluster. Since a cluster has size at most 2486, such a permutation
can be stored in an array of 16-bit integers, with indices relative
to the beginning of the cluster.

An optimized implementation of a mapping :math:`\xi^{\pm 1}` between
two boxes depends on the shape of these two boxes. We want to use
the same program code for both mappings, using different parameters
in both cases. Therefore we implement five procedures, labelled
with a *stage* counting form 0 to 4. These procedures map the boxes
as follows: 


===== ============ ================= ==============================
Stage :math:`\xi:` :math:`\xi^{-1}:` Mapping of shapes:
===== ============ ================= ==============================
    0   BC -> BC       BC -> BC       (1, 78, 32)  ->  (1, 78, 32)
    1   T0 -> T0       T0 -> T0      (45, 16, 32)  -> (45, 16, 32)
    2   T1 -> X0       T1 -> X1      (64, 12, 32)  -> (64, 16, 24)
    3   X0 -> X1       X1 -> X0      (64, 16, 24)  -> (64, 16, 24)
    4   X1 -> T1       X0 -> T1      (64, 16, 24)  -> (64, 12, 32)
===== ============ ================= ==============================

The functions in module ``gen_xi_functions.c`` deal with the 
operation of the generator :math:`\xi` of the Monster on
:math:`Q_{x0}` and :math:`98280_x`. They are also used for 
computing tables describing the monomial operation of 
:math:`\xi^{\pm 1}` on :math:`98280_x`.
Module ``mmgroup.dev.generators.gen_xi_ref`` is a pure python 
substitute for this set of C functions; but calculating the tables
with python takes a rather long time.

Class ``MM_TablesXi`` in module ``mmgroup.dev.mm_basics.mm_tables_xi``
contains 3-dimensional arrays that can be used for computing the
monomial operation of :math:`\xi^{\pm 1}`. In each array the first
index (running from 0 to 4) indicates the stage of the computation.
as shown in the previuos table.
 
The following table decribes the entries of the arrays 
PERM_TABLES, SIGN_TABLES, OFFSETS, and SHAPES in that class.

.. table:: Arrays encoding the monomial operation  :math:`\xi^{\pm 1}`
  :widths: 20, 20, 60

  =========== ============= ===============================================================
  Array name  Shape         Meaning
  =========== ============= ===============================================================
  PERM_TABLES (5, 2, c*m*n) Entry (s, e, h*m*n + i*n + j) is the index of the preimage of 
                            output entry (h,i,j) in stage s under :math:`\xi^{e+1}`; 
                            with [c, m, n] = SHAPES[s, 1].

  SIGN_TABLES (5, 2, c*m)   Bit j of entry (s, e, h*m + i) is the sign of the preimage of 
                            output entry (h,i,j) in stage s under :math:`\xi^{e+1}`; 
                            with [c, m, n] = SHAPES[s, 1].

  OFFSETS     (5, 2, 2)     Entries (s, e, 0) and (s, e, 1) are the offests of the source
                            box and the destination box when computing :math:`\xi^{e+1}`
                            in stage s.

  SHAPES      (5, 2, 3)     Entries (s, 0) and (s, 1) are triples encoding the shape of 
                            the source and the destiniation box at stage s. 
  =========== ============= ===============================================================


Here we abbreviate the phrase "the entry in cluster h, row i, column j,
of the output vector computed in stage s" to
"output entry (h,i,j) in stage s". Entries in array PERM_TABLES are
indices relative to the beginning of the current cluster.  


The following code block multiplies a vector ``v_in`` in the
representation :math:`\rho` of the Monster with :math:`\xi^e` for
:math:`e=1, 2` and stores the monomial part of the result in the
vector ``v_out``. It uses the arrays described in the table above. 
The function in that code block is tested in module
``mmgroup.tests.test_mm_op.test_table_xi``.

.. code-block:: python

    from mmgroup.dev.mm_basics.mm_tables_xi import MM_TablesXi
    tbl = MM_TablesXi()
    PERM_TABLES = tbl.PERM_TABLES
    SIGN_TABLES = tbl.SIGN_TABLES
    OFFSETS = tbl.OFFSETS
    SHAPES = tbl.SHAPES

    def map_xi(v_in, e, v_out):
        r"""Compute v_out = v_in * \xi**e, for e = 1, 2.

        We compute the monomial part of v_out only.
        v_in is a array of integers corresponding to a vector in the
        representation \rho of the Monster. Vector v_out is
        an empty vector of integers of the same size.
        """
        for stage in range(5):
            box_in = v_in[OFFSETS[stage][e-1][0]:]
            box_out = v_out[OFFSETS[stage][e-1][1]:]
            shape_in = SHAPES[stage][0]
            shape_out = SHAPES[stage][1]
            cluster_in_size = shape_in[1] * 32
            cluster_out_size = shape_out[1] * 32
            cluster_perm_size = shape_out[1] * shape_out[2]
            cluster_sign_size = shape_out[1]
            perm_table = PERM_TABLES[stage][e-1]
            sign_table = SIGN_TABLES[stage][e-1]
            for cluster in range(shape_out[0]):
                cluster_in = box_in[cluster * cluster_in_size:]
                cluster_out = box_out[cluster * cluster_out_size:]
                cluster_perm = perm_table[cluster * cluster_perm_size:]
                cluster_sign = sign_table[cluster * cluster_sign_size:]
                for i in range(shape_out[1]):
                    for j in range(shape_out[2]):
                        x = cluster_in[cluster_perm[shape_out[2]*i + j]]
                        x = (-1)**(cluster_sign[i] >> j) * x
                        cluster_out[32 * i + j] = x


Non-monomial operation of the generators  :math:`\xi^e`
.........................................................


This subsection is yet to be written!



Conjugation of  :math:`\tilde{x}_d x_\delta` with :math:`y_e`
.............................................................
  

Sometimes we have to conjugate an element :math:`\tilde{x}_d x_\delta`, 
where :math:`\tilde{x}_d = x_d x_{\theta(d)}`,  with :math:`y_e`. We have 

.. math::
     y_e^{-1} \tilde{x}_d x_\delta y_e & = 
     x_{-1}^\alpha x_\Omega^\beta \tilde{x}_d \tilde{x}_e^{|\delta|} 
     x_\delta  x_\epsilon \, , \\
     \alpha & = \theta(d,e) + \langle e,\delta\rangle^{1 + |\delta|}
     + \mbox{sign}(e) \, ,\\
     \beta &= \theta(e,d) + \langle e,\delta\rangle
     + {P(e)}^{|\delta|} \, , \\
     \epsilon &= A(d,e) +  {\theta(e)}^{|\delta|} \, .



