

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



.. _computation-leech2:


Computations in the Leech lattice modulo 2
------------------------------------------


In this section we describe the operation of the Conway group
:math:`\mbox{Co}_1` on :math:`\Lambda / 2 \Lambda`,
where :math:`\Lambda` is the Leech lattice. The description of
these orbits is taken from  :cite:`Seysen22`, Appendix A.

This operation is important for computing in the subgroup
:math:`G_{x0}` of structure :math:`2^{1+24}.\mbox{Co}_1` of the 
monster, as described in :cite:`Con85` and :cite:`Seysen20`.
Let :math:`Q_{x0}` be the normal subgroup of :math:`G_{x0}`
of structure :math:`2^{1+24}`.
The group :math:`\mbox{Co}_1 = G_{x0}/Q_{x0}` is the automorphism
group of :math:`\Lambda / 2 \Lambda`.  We assume that 
:math:`\mbox{Co}_1` and also :math:`G_{x0}` operate on  
:math:`\Lambda / 2 \Lambda` by right multiplication. Let 
:math:`\Omega` \in :math:`\Lambda / 2 \Lambda` be the 
type-4 vector corresponding to the the standard coordinate frame 
in :math:`\Lambda`.
The stabilizer of :math:`\Omega` is a subgroup :math:`N_{x0}` of  
:math:`G_{x0}` of structure :math:`2^{1+24}.2^{11}.M_{24}`. Thus 
the set of type-4 vectors in  :math:`\Lambda / 2 \Lambda`
corresponds the set of right cosets
:math:`2^{11}.M_{24} \backslash \mbox{Co}_1`. We have 
:math:`N_{x0} \backslash G_{x0} \cong 2^{11}.M_{24} 
\backslash \mbox{Co}_1`.
So identifying the right coset :math:`N_{x0} h` for a
:math:`h \in G_{x0}` reduces to the computation of the frame
:math:`\Omega h`.

In our construction of the monster the element :math:`h` may be
given as an arbitrary word in the generators of :math:`G_{x0}`
(modulo :math:`Q_{x0}`). So it is important to find a short word 
:math:`g` in the the generators of :math:`G_{x0}` with 
:math:`g \in N_{x0} h`. This can be achieved by applying a sequences
:math:`g' = g'_1,\ldots,g'_k` of automorphisms of
:math:`\Lambda / 2 \Lambda` to :math:`\Omega h` such that  
:math:`\Omega h g' = \Omega`, and each :math:`g'_i` corresponds to 
a generator of  :math:`G_{x0}` (modulo :math:`Q_{x0}`).


We assume that the reader is familiar with the Conway group 
:math:`\mbox{Co}_1` and the geometry of the Leech lattice as 
described in :cite:`Iva99`, section 4.1 - 4.7.



Orbits of the group :math:`N_{x0}` in the Leech lattice mod 2
.............................................................

In this subsection we describe the orbits of :math:`N_{x0}` on
:math:`\Lambda / 2 \Lambda`. The type of a vector in 
:math:`v \in  \Lambda` is the halved scalar product 
:math:`\frac{1}{2} \langle v, v \rangle`. The type of a vector 
in :math:`\Lambda / 2 \Lambda` is the type of its shortest
representative in :math:`\Lambda`. Each vector in  
:math:`\Lambda / 2 \Lambda` has type 0, 2, 3 or 4; and the 
group :math:`\mbox{Co}_1` is transitive on the vectors of any of 
these types.

The orbits of the groups :math:`N_{x0}` on the vectors of
type 2,3, and 4 on :math:`\Lambda` have been described in 
:cite:`Iva99`,  Lemma 4.4.1. :math:`N_{x0}` acts monomially
on the  the lattice :math:`\sqrt{8} \Lambda`, which has integers
coordinates in the standard Euclidean basis. Thus an orbit of 
:math:`N_{x0}` on :math:`\Lambda / 2 \Lambda` can be described by 
the *shapes* of the shortest vectors in the corresponding orbit of
:math:`\sqrt{8} \Lambda`. Here the shape of a vector is the 
multiset of the absolute values of the coordinates of the vector. 
E.g. a vector of shape :math:`(3^5 1^{19})` has 5 coordinates 
with absolute value 3 and 19 coordinates with absolute value 1.

A vector of type 2 or 3 in :math:`\Lambda / 2 \Lambda` has two
opposite representatives of the same type in :math:`\Lambda`; so 
its shape is uniquely defined.  A vector of type 4 in 
:math:`\Lambda / 2 \Lambda` has :math:`2 \cdot 24` representatives 
of type 4 in :math:`\Lambda` which are orthogonal except when equal 
or opposite. It is well known that a type-4 vector in 
:math:`\Lambda / 2 \Lambda` corresponds to a 
coordinate frame in the Leech lattice in standard coordinates,
see e.g. :cite:`CS99`, :cite:`Iva99`. 

  
The table at Lemma 4.4.1. in :cite:`Iva99` assigns a name and a 
shape to each orbit of  :math:`N_{x0}` on the vectors of type
2, 3, and 4 in  :math:`\Lambda`. The table at  Lemma 4.6.1. 
in :cite:`Iva99` assigns a name and one or more shapes to each 
orbit of :math:`N_{x0}` on the vectors of type 4 in  
:math:`\Lambda / 2 \Lambda`. We reproduce this information
for the orbits of :math:`N_{x0}` on :math:`\Lambda / 2 \Lambda`
in the following table. Here we also assign a subtype (which is 
a 2-digit number) to each orbit. The first digit of the subtype
specifies the type of the orbit and the second digit is used to 
distinguish between orbits of the same type.

Let :math:`\mathcal{C}` be the Golay code and :math:`\mathcal{C}^*` 
is the Golay cocode. For 
:math:`d \in \mathcal{C}, \delta \in \mathcal{C}^*` let 
:math:`\lambda_d, \lambda_\delta \in \Lambda / 2 \Lambda` be
defined as in :cite:`Seysen20`, Theorem 6.1. Conway's definition
of :math:`\lambda_d, \lambda_\delta` in :cite:`Con85` differs
slightly from that definition. Each 
:math:`v \in  \Lambda / 2 \Lambda` has a unique decomposition 
:math:`v = \lambda_d + \lambda_\delta, d \in \mathcal{C}, 
\delta \in \mathcal{C}^*`. The subtype of a vector 
:math:`\lambda_d + \lambda_\delta \in \Lambda / 2\Lambda` can
be computed from :math:`d` and :math:`\delta`,
as indicated in the following table:


.. math::
    \begin{array}{|c|c|c|c|c|c|c|}
    \hline 
    \mbox{Subtype} & \mbox{Name} & \mbox{Shape}  & |d| & |\delta| &
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
    \bar{\Lambda}_4^{5} & (5 \, 3^2 1^{21}), (3^{5} 1^{19}) & 
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
      :math:`\delta \subset d \Omega^{1 + |d|/8}` for a suitable 
      representative :math:`\delta` of the cocode element.

   2. :math:`|\delta|/2 = |d|/8 \pmod{2}`, 
      :math:`\delta \subset d \Omega^{1 + |d|/8}` for a suitable 
      representative :math:`\delta` of the cocode element.

   3. None of the conditions stated in Remarks 1 and 2 hold.
  

Here column *Subtype* lists the two-digit number that we use for 
describing the orbit. Columns  *Name* and *Shape* list the names and 
the shapes of the orbits as given in :cite:`Iva99`, Lemma 4.1.1 and 
4.6.1. Columns :math:`|d|` and  :math:`|\delta|` list conditions 
on the weight of a Golay code word :math:`d` and of (a shortest 
representative of) the Golay cocode element :math:`\delta`, 
respectively. Column :math:`\langle d, \delta \rangle` lists 
conditions on the scalar product of  :math:`|d|` and  
:math:`|\delta|`. All this information can easily be derived from  
:cite:`Iva99` and  :cite:`Con85` (or :cite:`Seysen20`).

The table provides enough information for effectively computing
the subtype of an element :math:`\lambda_d + \lambda_\delta` from
:math:`d` and :math:`\delta`. Function ``gen_leech2_type`` in file 
``gen_leech.c`` computes that subtype. It returns e.g. the subtype
46 as the hexadecimal number 0x46.

The following table may be helpful for memorizing the second digit
of a subtype of a vector in :math:`\Lambda / 2 \Lambda`.

.. math::
    \begin{array}{|c|l|}
    \hline 
    \mbox{Subtype} & \mbox{Description}  \\
    \hline 
    x0 & \mbox{Contains a vector $\lambda_d+\lambda_\delta$ with $d=0$;
    $\delta$ even} \\ 
    \hline 
    x1 & \mbox{Vectors $\lambda_d+\lambda_\delta$ with $|\delta|=1$} \\ 
    \hline 
    x2 & \mbox{Vectors $\lambda_d+\lambda_\delta$ with $|d| = 8, 16$; 
    and  $\delta \subset d \Omega^{1 + |d|/8}$, $|\delta|$ even}\\ 
    \hline 
    x3 & \mbox{Vector $\lambda_d+\lambda_\delta$ with $|\delta|=3$} \\ 
    \hline 
    x4 & \mbox{Vectors $\lambda_d+\lambda_\delta$ with $|d| = 8, 16$;
    not of subtype $x2$, $|\delta|$ even} \\ 
    \hline 
    x6 & \mbox{Vectors $\lambda_d+\lambda_\delta$ with $|d| = 12$; 
    $|\delta|$ even} \\ 
    \hline 
    x8 & \mbox{The vector $\lambda_d+\lambda_\delta$ with $|d| = 24$; 
    $\delta = 0$} \\ 
    \hline 
    \end{array}



Operation of the group :math:`G_{x0}` on the Leech lattice
..........................................................

The generators :math:`x_d, x_\delta, y_d, x_\pi, \xi`, with
:math:`d \in \mathcal{C}, \delta \in \mathcal{C}^*, \pi \in
\mbox{Aut}_{St} \mathcal{P}` operate on the subgroup 
:math:`G_{x0}` of the monster and also on 
:math:`\Lambda / 2 \Lambda`, as described in :cite:`Seysen20`. 
Here :math:`x_d, x_\delta` operate trivially on 
:math:`\Lambda / 2 \Lambda`. :math:`y_d` and :math:`x_\pi` 
act as permutations and sign changes on the coordinates of 
:math:`\Lambda`, respectively; they do not change the subtype. 

From the *Leech graph* in :cite:`Iva99`, section 4.7, we see
how the generators :math:`\xi` and :math:`\xi^2` may change
the subtype of a type-4 vector in  :math:`\Lambda / 2 \Lambda`.
The following figure shows a subgraph of that graph. It shows
the transitions of subtypes of type-4 vectors (due to transformations
by  :math:`\xi` or :math:`\xi^2`) that we actually use.


 
.. tikz:: Transitions of subtypes of type-4 vectors in the Leech lattice mod 2

   nodes={draw, % General options for all nodes
   line width=1pt,
   anchor=center, 
   text centered,
   rounded corners,
   minimum width=1.5cm, minimum height=8mm
   }, 
   \matrix (m) [matrix of  nodes,row sep=2.0ex,column sep=2em,
   minimum width=2em, nodes={rounded corners=6pt, draw}]
   {   
   43 &    &  42   &     \\
   &      &       &   40 &  48  \\
   46 	&      &  44 \\
   };         	
   \path[->] 
   (m-1-1) edge (m-1-3) 
   (m-1-1) edge (m-3-3) 
   (m-3-1) edge (m-3-3) 
   (m-1-3) edge (m-2-4) 
   (m-3-3) edge (m-2-4) 
   (m-2-4) edge (m-2-5) 
   ;







We will provide some more information about the operation of 
:math:`\xi^k` on :math:`\Lambda`. 

We say that a vector :math:`w \in \mathbb{Z}^n` is of shape
:math:`(m^\alpha 0^{n-\alpha} \bmod 2m)` if :math:`w` has
:math:`\alpha` coordinates equal to :math:`m` :math:`\pmod{2m}`
and :math:`n-\alpha` coordinates equal to 
:math:`0` :math:`\pmod{2m}`.
For a vector :math:`v = (v_0,\ldots, v_{23}) \in \mathbb{R}^{24}`
define 
:math:`(v_{4i}, v_{4i+1}, v_{4i+2}, v_{4i+3}) \in \mathbb{R}^4` 
to be the :math:`i`-th column of :math:`v`. This definition is
related to the Miracle Octad Generator (MOG) used for the
description of the Golay code and the Leech lattice, see
:cite:`CS99`, chapter 11. The following lemma has been shown
in :cite:`Seysen22`.

Lemma 

Let :math:`v \in \sqrt{8} \Lambda`, and let 
:math:`v^{(k)} = v \cdot {\xi^k}`.
Let :math:`w` be a column of :math:`v`, and let :math:`w^{(k)}`
be the corresponding column of :math:`v^{(k)}`. Then
:math:`|w^{(k)}| = |w|`. If :math:`w` has shape
:math:`(m^4 \bmod{2m})` then there is a unique :math:`k \in \{1,2\}`
such that :math:`w^{(k)}` has shape :math:`(0^4 \bmod{2m})`.
If :math:`w` has shape :math:`(m^2 0^2 \bmod{2m})` then 
:math:`w^{(k)}` has shape :math:`(m^2 0^2 \bmod{2m})` for all
:math:`k \in \mathbb{Z}`.


Computing a transversal of :math:`G_{x0} \backslash N_{x0}`
............................................................

The cosets in  :math:`G_{x0} \backslash N_{x0}` correspond to the
type-4 vectors in the Leech lattice modulo 2. Thus for computing 
a transversal of :math:`G_{x0} \backslash N_{x0}` we have to compute
an element of  :math:`G_{x0}` that maps the standard frame
:math:`\Omega` in :math:`\Lambda/ 2 \Lambda` to a vector  :math:`v`
for a given vector :math:`v` of type 4 in :math:`\Lambda/ 2 \Lambda`.

In the next few subsections we show how to map a vector of subtype
4x to a vector of subtype 48 by a sequence of generators
:math:`x_\pi, \xi^k` of length at most 6. The proofs are a bit 
involved, but the C function ``gen_leech2_reduce_type4_vector`` in 
file ``gen_leech.c`` constructs that sequence very fast. Then the
inverse of the element computed by function 
``gen_leech2_reduce_type4_vector`` is the requested element of
the transversal.



If we select a permutation in the Mathieu group :math:`M_{24}` then
we actually have to select an element of the group
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` of standard automorphisms
of the Parker loop described in :ref:`aut_ploop_label`. Whenever we have 
to select a permutation mapping a certain set to another set, we want
to keep the inverse of the selected element of 
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` as simple as possible.
More specifically, for that inverse we select the standard 
representative in :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` of the 
least feasible permutation in lexicographic order. We use the
function ``mat24_perm_from_map`` in file ``mat24_functions.c``
for finding a suitable permutation in :math:`M_{24}`.




Incorporating the baby monster into the computation of a transversal
....................................................................

The capability of computing in a subgroup :math:`H^+` of structure 
:math:`2.B` of the monster :math:`\mathbb{M}` (where :math:`B` is
the baby monster) is crucial for the success of our project. 
Although not strictly necessary,
we want to select elements of :math:`H^+` in our transversal of 
:math:`G_{x0} \backslash N_{x0}` whenever possible. 
In this subsection we explain how to do so. This is a bit technical 
and may be skipped at first reading. 

Let :math:`\beta = {\small \frac{1}{\sqrt{8}}}(e_2 - e_3) \in \Lambda`,
where  :math:`e_i` is the  :math:`i`-th unit vector in 
:math:`\mathbb{R}^{24}`. Then :math:`\beta` is of type 2 in the
Leech lattice :math:`\Lambda`, and the centralizer  :math:`H^+`
of :math:`x_\beta` in :math:`\mathbb{M}` has structure  :math:`2.B`, 
and the centralizer :math:`H` of :math:`x_\beta` in :math:`G_{x0}` 
has structure :math:`2^{1+22}.\mbox{Co}_2`. Our goal is find a 
representative of a coset in  :math:`N_{x0} \backslash G_{x0}` 
in :math:`H` (centralizing the element :math:`x_\beta` of 
:math:`\mathbb{M}`) whenever possible. 

The coset in :math:`N_{x0} \backslash G_{x0}` corresponding to 
a vector :math:`v \in \Lambda / 2 \Lambda` of type 4 has a nonempty 
intersection with :math:`H` if and only if :math:`v + \beta` is of 
type 2 and orthogonal to :math:`\beta` in the real Leech lattice 
:math:`\Lambda`, see  :cite:`Seysen22` for details. 

If this is the case then we compute an element of 
:math:`H \subset G_{x0}` 
that maps :math:`v + \beta` to the type-2 vector 
:math:`\Omega + \beta` and fixes :math:`\beta`. 

Note that the generator :math:`\xi` fixes :math:`\beta`. A permutation
:math:`\pi \in M_{24}` fixes  :math:`\beta` if and only if it fixes
the set :math:`\{2,3\}` of points.



From subtype 46 to subtype 44
.............................

Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of 
subtype 46. Then :math:`d` is a dodecad. Assume that :math:`d`
contains one column of the MOG. Then up to signs and permutation
of the columns of the MOG, and up to permutation of the entries of
a column, there is a vector :math:`v \in \sqrt{8}\Lambda` with 
:math:`v/\sqrt{8} = \lambda_r \pmod{2\Lambda}` 
that has MOG coordinates

   
.. math::
    \begin{array}{|c|c|c|c|c|c|}
    \hline 
    2 & 0 & 2 & 0 & 0 & 0 \\
    \hline 
    2 & 0 & 2 & 0 & 0 & 0  \\ 
    \hline 
    2 & 0 & 4 & 2 & 2 & 2 \\ 
    \hline 
    2 & 0 & 0 & 2 & 2 & 2  \\ 
    \hline 
    \end{array} \quad .
    
By the Lemma given above  there is a 
:math:`k \in \{1,2\}` such that column 0 of the the vector 
:math:`w = v \cdot {\xi^k}` has shape :math:`(4 \, 0^3)`. The other
columns of :math:`w` have the same shape as the corresponding 
columns of :math:`v`. Thus  :math:`w` has shape 
:math:`(4 \, 2^8 0^{15})`. So using the table given above we see
that  :math:`w` has type 44. 

Note that both, the dodecad :math:`d` and its complement, 
contain exactly one column of the MOG. The union of these two
columns is a grey even octad :math:`o`.

With a little extra effort we can show 
:math:`k = 2 - \langle o, \delta \rangle`, where
:math:`\langle ., . \rangle: \mathcal{C} \times \mathcal{C}^*
\rightarrow \{0,1\}` is the scalar product.

If dodecad  :math:`d` does not contain a column of the MOG
then we select a permutation that maps the first four entries
of  dodecad  :math:`d` to the points  :math:`0,1,2,3`, i.e.
to the first column of the MOG. Then we proceed as above.






From subtype 43 to subtype 42 or 44
.....................................


Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of subtype 43.
Then :math:`|\delta| = 3`. The three points of :math:`\delta` can 
lie in one, two, or three different columns of the MOG.

If all entries of  :math:`\delta` lie in the same column of the MOG
then up to signs and permutation of the columns of the MOG, and up 
to permutation of the entries of a column, there is a vector 
:math:`v \in \sqrt{8}\Lambda` with 
:math:`v/\sqrt{8} = \lambda_r \pmod{2\Lambda}` 
that has MOG coordinates

.. math::
    \begin{array}{|c|c|c|c|c|c|}
    \hline 
    5 & 1 & 1 & 1 & 1 & 1 \\
    \hline 
    3 & 1 & 1 & 1 & 1 & 1  \\ 
    \hline 
    3 & 1 & 1 & 1 & 1 & 1 \\ 
    \hline 
    1 & 1 & 1 & 1 & 1 & 1  \\ 
    \hline 
    \end{array}  \quad .

By the Lemma given above  there is a :math:`k \in \{1,2\}` such that 
one column of :math:`w = v \cdot {\xi^k}` (and hence all columns) 
have even entries. Column 0 of :math:`w` has squared norm 
:math:`44 = 6^2 + 2^2 + 2^2 + 0^2`. That decomposition of
:math:`44` into a sum of four even squares is unique, so 
column 0 has shape :math:`(6 \, 2^2 0)`. The other columns
have shape  :math:`(2 \, 0^3)`. So :math:`w` is of
shape  :math:`(6 \, 2^7 0^{15})` and hence of subtype 42.
 
With a little extra effort we can show 
:math:`k = 2 - \langle d, \omega \rangle`, where :math:`\omega`
is the standard tetrad represented by :math:`\{0,1,2,3\}`, and 
:math:`\langle ., . \rangle` as above. 

If the entries of  :math:`\delta` lie in two different columns of 
the MOG then up to signs and permutations as above,  the 
corresponding vector in :math:`v \in \sqrt{8}\Lambda` has MOG 
coordinates


.. math::
    \begin{array}{|c|c|c|c|c|c|}
    \hline 
    3 & 5 & 1 & 1 & 1 & 1 \\
    \hline 
    3 & 1 & 1 & 1 & 1 & 1  \\ 
    \hline 
    1 & 1 & 1 & 1 & 1 & 1 \\ 
    \hline 
    1 & 1 & 1 & 1 & 1 & 1  \\ 
    \hline 
    \end{array}  \quad .


By a similar argument as in the first case there is a 
:math:`k \in \{1,2\}` such that the first two columns of 
:math:`w = v \cdot {\xi^k}` have shape  :math:`(4 \, 2\, 0^2)`
and  :math:`(4 \, 2^3)`, and that the other columns of :math:`w`
have shape :math:`(2 \, 0^3)`. Thus :math:`w` is of subtype 44.

If the entries of  :math:`\delta` lie in three different columns 
of  the MOG then we apply a permutation in :math:`M_{24}` that maps
:math:`\delta` to  :math:`\{1,2,3\}`, and proceed as in the first
case.


From subtype 44 to subtype 40
.............................


Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of subtype 44. 
Then :math:`d` is an octad or a complement of
an octad. We call that octad :math:`o`.
If the cocode word :math:`\delta` is a duad, let :math:`c_0, c_1` 
be the two points in that duad. Otherwise, let :math:`c` by any 
tetrad intersecting the octad :math:`o` in two points, so that
:math:`c` is equal to the sextet :math:`\delta` modulo the cocode
:math:`\mathcal{C}^*`. Let :math:`c_0, c_1` be the two points in  
:math:`c` which are not in the octad  :math:`o`. 

Assume first that :math:`o` is a grey even octad, i.e. a union of
two columns of the MOG, and that the points  :math:`c_0, c_1`
are in the same column of the MOG.

Then up to signs and permutation of the columns of the MOG, and up 
to permutation of the entries of a column, there is a vector 
:math:`v \in \sqrt{8}\Lambda` with 
:math:`v/\sqrt{8} = \lambda_r \pmod{2\Lambda}` 
that has MOG coordinates

.. math::
    \begin{array}{|c|c|c|c|c|c|}
    \hline 
    2 & 2 & 4 & 0 & 0 & 0 \\
    \hline 
    2 & 2 & 4 & 0 & 0 & 0  \\ 
    \hline 
    2 & 2 & 0 & 0 & 0 & 0 \\ 
    \hline 
    2 & 2 & 0 & 0 & 0 & 0  \\ 
    \hline 
    \end{array}

Similar to the previous cases, we can show that there is a unique
:math:`k \in \{1,2\}` such that :math:`w = v \cdot {\xi^k}` has 
shape :math:`(4^4 0^{20})`. So :math:`w` is of subtype 40.

One can show :math:`k = 2 - \sigma`, where
:math:`\sigma = |c \cap w | \pmod{2}`, :math:`\sigma \in \{0,1\}` ,
and :math:`w` is any of the two columns of the MOG contained in
the octad  :math:`o`. 


If the octad :math:`o` and the two points  :math:`c_0, c_1` do not 
satisfy the condition given above then we first construct a 
permutation that maps  :math:`c_0, c_1` to  :math:`2, 3`, and four 
points of  :math:`o` to :math:`4, 5, 6, 7` as follows.

We compute the syndrome of the set  :math:`c_0, c_1, o_0, o_1, o_2`.
That syndrome intersects with  :math:`o` in exactly one point
:math:`o_j`. The we construct a permutation in :math:`M_{24}` that 
maps the hexad :math:`(c_0, c_1, o_0, o_1, o_2, o_j)` to the hexad
:math:`(2, 3, 4, 5, 6, 7)`. Such a permutation exists, since both
hexads are subsets of an octad.


From subtype 42 to subtype 40
.............................

Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of subtype 42. 
Then :math:`d` is an octad or a complement of an octad. We call that
octad :math:`o` as in the case of type 42. If that octad is grey,
we map vector :math:`\lambda_r` to a vector of type 40 in the same 
way as in that case.

Otherwise we first map the lowest four points of octad  :math:`o` to 
the set :math:`0, 1, 2, 3`, using a permutation in  :math:`M_{24}`.


From subtype 40 to subtype 48
.............................


Let :math:`\lambda_r` be of subtype 40. Then
:math:`\lambda_r = \alpha \lambda_\Omega + \lambda_\delta`,
:math:`\alpha = 0, 1`, for some 
:math:`\delta \in \mathcal{C}^*`, :math:`|\delta|=4`.
If :math:`\delta` is equal to the standard tetrad  :math:`\omega` 
represented by :math:`\{0,1,2,3\}` then there is a unique power of
:math:`\xi` that maps :math:`\lambda_r` to :math:`\lambda_\Omega`, 
which is of type 48.

Otherwise let :math:`c` be the tetrad containing the point :math:`0`
that corresponds to :math:`\delta`. The we first apply a permutation
in  :math:`M_{24}` that maps  :math:`c` to  :math:`\{0,1,2,3\}`.


From subtype 21 to subtype 22, fixing :math:`\beta`
...................................................

Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of subtype 21,
and orthogonal to :math:`\beta` in the real Leech lattice. 
Then :math:`|\delta| = 1`.

We  proceed as in the case of subtype 43. Then there is a
:math:`k \in \{1,2\}` such that one column of 
:math:`w = v \cdot {\xi^k}` (and hence all columns) have even 
entries. It is easy to see that :math:`w` is of subtype 22.


From subtype 22 to subtype 20, fixing :math:`\beta`
...................................................

Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of subtype 22,
and orthogonal to :math:`\beta` in the real Leech lattice.  
Then :math:`d` is an octad or a complement of an octad. We call that 
octad :math:`o` as in the cases of subtype 42 or 44. For subtype 22
the cocode word :math:`\delta` can be considered as a subset of
octad :math:`o`.

The set :math:`\{2,3\}` either contained in or disjoint to octad
:math:`o`. 

If octad :math:`o` is grey and even then we may directly apply a power
of :math:`\xi` that transforms :math:`\lambda_r` to a vector of
subtype 22, similar to the operation in the case of subtype 42.
Otherwise we proceed as follows.

If the set :math:`\{2,3\}` is disjoint from :math:`o` we proceed
exactly as in the case of subtype 44, replacing :math:`c_0` and 
:math:`c_1` by  :math:`2` and :math:`3`. This operation fixes the
points :math:`2` and :math:`3`.

Otherwise we first fix the points :math:`2` and :math:`3` and map 
the lowest two points of :math:`o \setminus \{2,3\}` to :math:`0` 
and :math:`1`, using a permutation in  :math:`M_{24}`.


From subtype 20 to the vector :math:`\beta + \Omega`, fixing :math:`\beta`
..........................................................................

Let :math:`\lambda_r` be of subtype 20, and orthogonal to 
:math:`\beta` in the real Leech lattice.  Then
:math:`\lambda_r = \alpha \lambda_\Omega + \lambda_\delta`,
:math:`\alpha = 0, 1`, for some 
:math:`\delta \in \mathcal{C}^*`, :math:`|\delta|=2`. In case
:math:`\delta =  \{2,3\}` we must have :math:`\alpha = 1`, i.e.
:math:`\lambda_r = \beta + \Omega`, and we are done. Otherwise both 
points in :math:`\delta` are different from :math:`2` and :math:`3`. 

In case :math:`\delta = \{0, 1\}` the transformation :math:`\xi^k`
with :math:`k = 2- \alpha` maps :math:`\lambda_r` to :math:`\beta`.

Otherwise we first apply a permutation in :math:`M_{24}` that maps
:math:`\delta` to :math:`\{0, 1\}`, and fixes :math:`2` and :math:`3`.




Orbits of the group :math:`N_{x0}` in the  group :math:`Q_{x0}`
---------------------------------------------------------------

Function ``gen_leech2_reduce_n`` in file ``gen_leech_reduce_n.c`` 
computes a standard representative of the orbit of a nonzero vector
in the extraspecial subgroup :math:`Q_{x0}` of :math:`G_{x0}` under
the action of the group :math:`N_{x0}`. Here :math:`N_{x0}` acts on
the group :math:`Q_{x0}` of structure  :math:`2_+^{1+24}` by 
conjugation.

The quotient of :math:`Q_{x0}` by its center  :math:`\{x_{\pm 1}\}` 
is isomophic to  :math:`\Lambda / 2 \Lambda`, i.e. to the Leech 
lattice mod 2. So we may define the subtype of any 
:math:`v \in Q_{x0}` as the subtype of the element :math:`v x_{\pm 1}`
of  :math:`\Lambda / 2 \Lambda`. 

With one exception the orbit of a nonzero vector :math:`v \in Q_{x0}` 
under  :math:`N_{x0}` depends on the subtype of :math:`v` only. There are
exactly two elements :math:`x_{\pm 1}` of  :math:`Q_{x0}` of subtype
00, which are of course in different :math:`N_{x0}` orbits.

For each nonzero orbit of :math:`N_{x0}` we choose a 
representative  :math:`x_d x_\delta`, 
:math:`d \in \mathcal{P}, \delta \in \mathcal{C}^*` (depending
on the subtype of the orbit) according to the following table:


.. math::
    \begin{array}{|c|c|c|}
    \hline 
    \mbox{Subtype} &  x_d  & \delta \\
    \hline 
    00   & x_{-1}         & \{\}  \\
    20   & x_1            & \{2,3\} \\
    21   & x_1            & \{0\} \\
    22   & x_o            & \{\} \\
    31   & x_\Omega       & \{0\} \\
    33   & x_\Omega       & \{1,2,3\} \\
    34   & x_o            & \{0,8\} \\
    36   & x_D            & \{\} \\
    40   & x_1            & \{0,1,2,3\} \\
    42   & x_\Omega x_o   & \{\} \\
    43   & x_1            & \{1,2,3\} \\
    44   & x_o            & \{8, 9\} \\
    46   & x_D            & \{0, 12\} \\
    48   & x_\Omega       & \{\} \\
    \hline 
    \end{array}

:math:`x_o` is the (positive) element of the Parker loop
corresponding to the stadard octad :math:`\{0,1,2,3,4,5,6,7\}`.

:math:`x_D` is the (positive) element of the Parker loop
corresponding to the stadard dodecad 
:math:`\{0,4,8, 13,14,15, 17,18,19, 21,22,23\}`.

This table can be computed by function ``gen_leech2_reduce_n_rep`` 
in file ``gen_leechc_reduce_n.c``.  

.. _check_equality_monster:

Checking equality of monster elements and membership in :math:`2^{1+24}.\mbox{Co}_1`
-------------------------------------------------------------------------------------

Checking equality of monster elements
.....................................

For checking equality of two elements of the monster we may use 
two nonzero vectors :math:`w_{71}` and :math:`w_{94}` in the 
198883-dimensional representation of the monster that are fixed by 
an element of order 71 and negated by an element of order 94,
respectively. In :cite:`LPWW98` it is shown that only the neutral
element of the monster fixes both vectors, :math:`w_{71}` and 
:math:`w_{94}`. So we may check if an element is the identity in
the monster. Here the corresponding calculations can be done modulo 
any odd prime; and we may even use different primes for the 
operation on :math:`w_{71}` and on :math:`w_{94}`. 

The python script ``mmgroup.structures.find_order_vectors.py``
generates suitable vectors :math:`w_{71}` and :math:`w_{94}` at
random and stores data for the fast recomputation of these vectors
in file  ``order_vector_data.py``. Since we actually work in the
representation :math:`\rho` of dimension 1 + 196883, we have to 
make sure that the generated vectors are not fixed by the monster.

More precisely, we generate :math:`w_{71}` as a vector in the 
representation :math:`\rho_3` of the monster in characteristic
3, and  :math:`w_{94}` as a vector in the representation 
:math:`\rho_5` of the monster in characteristic 5. We combine
these two vectors to a vector in the representation 
:math:`\rho_{15}` (modulo 15) of the monster via Chinese 
remaindering. Note that a computation in :math:`\rho_{15}` is
faster than the combination of the corresponding computations 
in :math:`\rho_{3}`. We write :math:`w` for the vector
combined from :math:`w_{71}` and :math:`w_{94}`.

   
Testing membership in :math:`2^{1+24}.\mbox{Co}_1`
...................................................

For testing membership in :math:`2^{1+24}.\mbox{Co}_1` we impose
additional conditions on the vector  :math:`w_{71}` in the 
representation :math:`\rho_3` defined in the last section.

For an odd number :math:`p` let :math:`\rho_{A,p}` be the subspace of 
:math:`\rho_p` spanned by the basis vectors of :math:`\rho_p` with tag
``A``. Then :math:`\rho_{A,p}` is fixed by the subgroup :math:`G_{x0}` 
of structure :math:`2^{1+24}.\mbox{Co}_1` of the monster. A vector in
:math:`\rho_{A,p}` has a natural interpretation as a symmetric bilinear
form on the Leech lattice (modulo :math:`p`), see :cite:`Con85`. Let 
:math:`Q_{p}` be the symmetric bilinear form on the Leech lattice 
(modulo :math:`p`) given by the scalar product; then :math:`Q_{p}` 
is given by the the :math:`24 \times 24` unit matrix, up to a 
scalar factor.

Let :math:`w` be a vector in :math:`\rho_{15}` as in the last section. 
For  any :math:`w` in :math:`\rho_{p}`, let :math:`w_{(A,p)}` be the 
bilinear from on the Leech lattice (modulo :math:`p`) corresponding 
to the projection of :math:`w` onto :math:`\rho_{A,p}`.
We shall search for a vector :math:`w` satisfying the 
condition stated an the last section. In addition we require that 
for some integer :math:`d` (modulo 3) the bilinear form 
:math:`w_{(A,3)} - d \cdot Q_3` has a kernel of dimension one, and 
that this kernel contains a type-4 vector :math:`x_4` of the Leech 
lattice modulo 3.The chance of finding such a 
vector :math:`w` is discussed in :cite:`Seysen22`, Appendix B.
In average, finding a suitable vector :math:`w` take less
than a minute on the author's compute, when running with 12
processes in parallel.


Let :math:`\Omega` be the standard frame of
the Leech lattice (modulo 3). Then :math:`\Omega` is  of type 4
and consists of the vectors :math:`(0,\ldots,0,\pm 1,0,\ldots,0)`.
Once we have found a suitable vector :math:`w`, we 
apply a transformation in :math:`G_{x0}` to :math:`w` that
maps the frame in the Leech lattice given by  :math:`x_4` to 
:math:`\Omega`. Therefore we work in the Leech lattice modulo 2, as
described in section :ref:`computation-leech2`. Conversion from a
type-4 vector in the Leech lattice mod 3 to a type-4 vector in the 
Leech lattice mod 2 is easy; the C function ``gen_leech3to2_type4`` 
in file ``gen_leech3.c`` does that job. In the sequel we also write 
:math:`w` for the transformed vector. So we may assume that the 
kernel of :math:`w_{(A,3)} - d \cdot Q_3` is one dimensional and 
spanned by a vector in the standard frame :math:`\Omega` of the 
Leech lattice modulo 3. 

For doing linear algebra in :math:`\mathbb{F}_3` we store a vector 
of 16 integers modulo 3 in a 64-bit integer; and we use arithmetic, 
Boolean and shift operations. 
 
Assume that :math:`g` is in :math:`G_{x0}`. Then the kernel of 
:math:`w_{A,3} \cdot g - d \cdot Q_3` is spanned by a vector in the 
frame :math:`\Omega \cdot g`. So we may also compute 
:math:`\Omega \cdot g` in the Leech lattice modulo 2. Using the 
methods in section :ref:`computation-leech2` we can compute a 
:math:`h_1` in  :math:`G_{x0}` that maps :math:`\Omega \cdot g` to 
:math:`\Omega`. So :math:`gh_1` fixes :math:`\Omega`. Hence 
:math:`gh_1` is in the subgroup :math:`N_{x0}` of structure 
:math:`2^{1+24+11}.M_{24}` of :math:`G_{x0}`. 

The factor group  :math:`M_{24}` of  :math:`N_{x0}` acts on
the :math:`24 \times 24` matrix representing the bilinear form
:math:`w_{(A,p)}` by permuting the rows and the columns, up to
changes of signs that we will ignore at the moment. So we
store the multiset of the absolute values of the entries of each
row of the matrix :math:`w_{(A,p)}`. Then :math:`M_{24}` permutes
these 24 multisets by its natural action. So we may recover
the element of :math:`M_{24}` performing that action
and compute a preimage :math:`h_2` of the inverse of that
element in :math:`N_{x0}`. Thus :math:`gh_1 h_2` is in the
subgroup  :math:`2^{1+24+11}` of :math:`N_{x0}`. Note that, 
originally, the vector :math:`w` has been defined modulo 15, so 
that almost certainly we have enough information in the matrix
:math:`w_{(A,15)}` to compute the permutation in :math:`M_{24}`. In 
the implementation we compute a hash function on the 24 multisets
obtained from the rows of :math:`w_{(A,15)}`, and we generate a new 
vector :math:`w` if these 24 hash values are not mutually different.

Now that task of reducing an element of  :math:`2^{1+24+11}` is easy.
We first eliminate the factor :math:`2^{11}` be checking signs 
in :math:`w_{(A,15)}`, or signs of suitable entries of vector 
:math:`w` with tag ``A``, which is the same thing. Then we
eliminate the factor :math:`2^{24}` be checking signs of suitable 
entries of :math:`w` with tags ``B, C, T, X``, and, finally, the 
factor :math:`2^{1}` by checking a sign of an entry with 
tag ``Y`` or ``Z``.

So we obtain a decomposition :math:`g h_1 h_2 h_3 = 1`, with 
:math:`h_i` in a transversal of :math:`G_i / G_{i+1}` for
:math:`G_1 = G_{x0}`, :math:`G_2 = N_{x0}`,  
:math:`G_3 = 2^{1+24+11}`, :math:`G_4 = \{1\}`. If  
:math:`g \notin  G_{x0}` then this procedure fails in one of the 
steps mentioned above. We abort as early as possible in case of 
failure, in order to accelerate functions searching for an 
element of :math:`G_{x0}` that have a low probability of success.

Module ``mmgroup.dev.mm_reduce.find_order_vector.py`` contains
a function for finding a suitable vector  :math:`w`. Once a 
vector :math:`w` has been found, it creates a module
``mmgroup.dev.mm_reduce.order_vector_data.py`` containing
information for recomputing the same vector :math:`w`
more quickly. The code generation process creates a C file 
``mm_order_vector.c`` containing the vector :math:`w`
in its source code together with some data required for
using :math:`w` and for checking its correctness.
We believe that this method simplifies the integration of
our project into computer algebra systems like GAP or Magma,
since these systems may now simply call C functions for 
computing in the monster group, without bothering how to 
generate a suitable vector :math:`w`.




