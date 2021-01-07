

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


Implementing generators of the Monster group
--------------------------------------------


The operation :math:`y_{f} \cdot x_{e} \cdot x_\epsilon`
........................................................

The operation  :math:`g = y_{f} \cdot x_{e} \cdot x_\epsilon` on
:math:`\rho_p` is coded in function ``mm_op<p>_xy`` in the
automatically generated C file ``mm_op<p>_xy.c`` for the modulus ``p``.
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



Computations in the Leech lattice modulo 2
------------------------------------------

In this section we describe the operation of the Conway group
:math:`Co_1` on the type-4 vectors in :math:`\Lambda / 2 \Lambda`,
where :math:`\Lambda` is the Leech lattice.
It is well known that such a type-4 vector corresponds to a 
coordinate frame in the Leech lattice in standard coordinates,
see e.g. :cite:`CS99`, :cite:`Iva99`. 

This operation is important for computing in the subgroup
:math:`G_{x0}` of structure :math:`2^{1+24}.\mbox{Co}_1` of the 
monster, as described in :cite:`Con85` and :cite:`Seysen20`.
Let :math:`Q_{x0}` be the normal subgroup of :math:`G_{x0}`
of structure :math:`2^{1+24}`.
The group :math:`\mbox{Co}_1 = G_{x0}/Q_{x0}` is the automorphism
group of :math:`\Lambda / 2 \Lambda`.  We assume that 
:math:`\mbox{Co}_1` and also :math:`G_{x0}` operate on  
:math:`\Lambda / 2 \Lambda` by right multiplcation. Let 
:math:`\Omega` \in :math:`\Lambda / 2 \Lambda` be the 
type-4 vector corresponding to the the standard coordinate frame 
in :math:`\Lambda`.
The stabilizer of :math:`\Omega` is a subgroup :math:`N_{x0}` of  
:math:`G_{x0}` of structure :math:`2^{1+24}.2^{11}.M_{24}`. Thus 
the set of type-4 vectors in  :math:`\Lambda / 2 \Lambda`
corresponds the set of right cosets
:math:`\mbox{Co}_1 \setminus 2^{11}.M_{24}`. We have 
:math:`G_{x0} \setminus N_{x0} \cong \mbox{Co}_1 \setminus 2^{11}.M_{24}`.
So identifying the right coset :math:`N_{x0} h` for a
:math:`h \in G_{x0}` reduces to the computation of the frame
:math:`\Omega h`.

In our construction of the monster, the element :math:`h` may 
given as an arbitrary word in the generators of :math:`G_{x0}`
(modulo :math:`Q_{x0}`). So it is important to find a short word 
:math:`g` in the the generators of :math:`G_{x0}` with 
:math:`g \in N_{x0} h`. This can be achieved by applying a sequences
:math:`g' = g'_1,\ldots,g'_k` of automorphisms of
:math:`\Lambda / 2 \Lambda` to :math:`\Omega h` such that  
:math:`\Omega h g' = \Omega`, and each :math:`g'_i` corresponds to 
a generator of  :math:`G_{x0}` (modulo :math:`Q_{x0}`).


We assume that the reader is familiar with the
Conway group and the geometry of the Leech lattice as described in
:cite:`Iva99`, section 4.1 - 4.7.



Orbits of the group :math:`N_{x0}` in the Leech lattice mod 2
.............................................................

In this subsection w  describe the orbits of :math:`N_{x0}` on
:math:`\Lambda / 2 \Lambda`. The type of a vector in 
:math:`v \in  \Lambda` is the halved scalar product 
:math:`\frac{1}{2} \langle v, v \rangle`. The type of a vector 
in :math:`\Lambda / 2 \Lambda` is the type of its shortest
representative in :math:`\Lambda`. Each vector in  
:math:`\Lambda / 2 \Lambda` has type 0, 2, 3 or 4; and the 
group :math:`G_{x0}` is transtive on the vectors of any of these
types.

The orbits of the groups :math:`N_{x0}` on the vectors of
type 2,3, and 4 on :math:`\Lambda` have been described in 
:cite:`Iva99`,  Lemma 4.4.1. Since :math:`N_{x0}` acts monomially
on the  the lattice :math:`\sqrt{8} \Lambda` given in the
standard Euclidean basis, an orbit of :math:`N_{x0}` on
:math:`\Lambda` can be described by the *shape* of its vectors in 
:math:`\sqrt{8} \Lambda`. Here the shape of a vector is the 
multiset of the absolute values of the coordinates of the vector. 
E.g. a vector of shape :math:`(3^5 1^{19})` has 5 coordinates 
with absolute value 3 and 19 coordinates with absolute value 1.

A vector of type 2 or 3 in :math:`\Lambda / 2 \Lambda` has a unique
representative of the same type in :math:`\Lambda` (up to sign); so 
its shape is uniquely defined.  A vector of type 4 in 
:math:`\Lambda / 2 \Lambda` has :math:`2 \cdot 24` representatives 
of type 4 in :math:`\Lambda` which are orthogonal except when equal 
or opposite.   

The table at Lemma 4.4.1. in :cite:`Iva99` assigns a name and a 
shape to each orbit of  :math:`N_{x0}` on the vectors of type
2, 3, and 4 in  :math:`\Lambda`. The table at  Lemma 4.6.1. 
in :cite:`Iva99` assigns a name and a shape to each orbit of  
:math:`N_{x0}` on the vectors of type 4 in  
:math:`\Lambda / 2 \Lambda`. We repreduce this information
for the orbits of :math:`N_{x0}` on :math:`\Lambda / 2 \Lambda`
in the following table. Here we also assign a subtype
(which is an 2-digit number) to each orbit.  Here the first 
digit specifies the type of the orbit and the second digit is
used to distinguish between orbits of the same type. 

In :cite:`Con85`, section 7 and `Sey19`, section section 6.1,
an element of :math:`\Lambda / 2 \Lambda` is also described
as a product :math:`x_d \cdot x_\delta \in Q_{x0} / \{\pm 1\}`,  
with :math:`x_d \in \mathcal{C}`, 
:math:`x_\delta \in \mathcal{C}^*`. Here with :math:`\mathcal{C}`
is the Golay code and :math:`\mathcal{C}^*` is the Golay cocode.
The :math:`\Lambda / 2 \Lambda`-orbit of :math:`x_d \cdot x_\delta`
can also be computed from :math:`d` and :math:`\delta`
as indicated in the table:


.. math::
    \begin{array}{|c|c|c|c|c|c|c|}
    \hline 
    \mbox{Subype} & \mbox{Name} & \mbox{Shape}  & |d| & |\delta| &
    \langle d , \delta\rangle & \mbox{Remark} \\
    \hline
    00 & 
        & (0^{24})  & 
    0 & 0  & 0 \\     
    \hline
    20 & 
    \Lambda_2^4 & (4^2 0^{22}) & 
    0, 24  &  2  & 0 &    \\     
    \hline
    21 & 
    \Lambda_2^3 & (3 \,  1^{23}) & 
    \mbox{any} & 1 & |d| / 4 &  \\     
    \hline
    22 & 
    \Lambda_2^2 & (2^8 0^{16}) & 
    8, 16  &  \mbox{even} & 0 &  1. \\     
    \hline
    31 & 
    \Lambda_3^5 & (5 \, 1^{23}) & 
    \mbox{any} & 1 & |d| / 4 + 1 &  \\     
    \hline
    33 & 
    \Lambda_3^3 & (3^3 1^{21}) & 
    \mbox{any} & 3 & |d| / 4 &  \\     
    \hline
    34 & 
    \Lambda_3^4 & (4 \, 2^{8} 0^{15}) & 
    8, 16 & \mbox{even}  &  1 &  \\    
    \hline
    36 & 
    \Lambda_3^2 & (2^{12} 0^{12}) & 
    12 & \mbox{even}  &  0 &  \\    
    \hline
    40 & 
    \bar{\Lambda}_4^{4a} & (4^4 0^{20}) & 
    0, 24 & 4 &  0 &  \\    
    \hline 
    42 & 
    \bar{\Lambda}_4^{6} & (6 \, 2^7 0^{16}), (2^{16} 0^8) & 
    8, 16 & \mbox{even} &   0 & 2.\\    
    \hline
    43 & 
    \bar{\Lambda}_4^{5} & (5 \, 3^2 1^{21}), (3^{5} 1^{16}) & 
    \mbox{any} & 3 & |d| / 4 + 1 &  \\    
    \hline
    44 & 
    \bar{\Lambda}_4^{4b} & (4^2 2^8 0^{14}), (2^{16} 0^8) & 
    8, 16 & \mbox{even} &   0 & 3.\\    
    \hline 
    46 & 
    \bar{\Lambda}_4^{4c} & (4 \, 2^{12} 0^{11}) & 
    12 & \mbox{even} &   1 &  \\    
    \hline 
    48 & 
    \bar{\Lambda}_4^{8} & (8 \,  0^{23}) & 
    24 & 0 &  0 &  \\    
    \hline 
    \end{array}


Remarks

   1. :math:`|\delta|/2 = 1 + |d|/8 \pmod{2}`, 
      :math:`\delta \subset d \Omega^{1 + |d|/8}`
      for a suitable representative :math:`\delta` of the cocode.

   2. :math:`|\delta|/2 = |d|/8 \pmod{2}`, 
      :math:`\delta \subset d \Omega^{1 + |d|/8}`
      for a suitable representative :math:`\delta` of the cocode.

   3. None of the conditions stated in Remarks 1 and 2 hold.
  

Here column *Subtype* lists the numbers that we use for describing 
the orbit. Columns  *Name* and *Shape* list the names and the shapes
of the orbits as given in :cite:`Iva99`, Lemma 4.1.1  and 4.6.1.
Columns :math:`|d|` and  :math:`|\delta|` list conditions on the
weigth of a Golay code word :math:`d` and of (a shortest 
representative of) the Golay cocode word :math:`\delta`, 
respectively. Column :math:`\langle d, \delta \rangle` lists 
conditions on the scalar product of  :math:`|d|` and  
:math:`|\delta|`. All this information can easily be derived from  
:cite:`Iva99` and  :cite:`Con85` (or :cite:`Seysen20`).

The table provides enough information for effectively computing
the subtype of an element :math:`x_d \cdot x_\delta`. Function
``gen_leech2_type`` in file ``gen_leech.c`` compute that subtype.


Operation of the group :math:`G_{x0}` on the Leech lattice
..........................................................


Yet to be documented!

From type 46 to type 44
.......................

Yet to be documented!


Here the part :math:`x_d` of the vector in the Leech lattice mod 2 
corresponds to a dodecad. We map that dodecad to the dodecad 
:math:`(0,1,2,3,8,9,14,15,18,19,22,23)` by applying a permutation
:math:`\pi \in 2^{11}.M_{24}`. 

