

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
:math:`\mbox{Co}_1` on :math:`\Lambda / 2 \Lambda`,
where :math:`\Lambda` is the Leech lattice.

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
:math:`2^{11}.M_{24} \setminus \mbox{Co}_1`. We have 
:math:`N_{x0} \setminus G_{x0} \cong 2^{11}.M_{24} 
\setminus \mbox{Co}_1`.
So identifying the right coset :math:`N_{x0} h` for a
:math:`h \in G_{x0}` reduces to the computation of the frame
:math:`\Omega h`.

In our construction of the monster, the element :math:`h` may be
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
the *shapes* of the shortest vectors in the corrsponding orbit of
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
in :cite:`Iva99` assigns a name and a shape to each orbit of  
:math:`N_{x0}` on the vectors of type 4 in  
:math:`\Lambda / 2 \Lambda`. We repreduce this information
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
:math:`v \in  \Lambda / 2 \Lambda` has a unique decomposion 
:math:`v = \lambda_d + \lambda_\delta, d \in \mathcal{C}, 
\delta \in \mathcal{C}^*`. The subtype of a vector 
:math:`\lambda_d + \lambda_\delta \in \Lambda / 2\Lambda` can
be computed from :math:`d` and :math:`\delta`,
as indicated in the following table:


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
weight of a Golay code word :math:`d` and of (a shortest 
representative of) the Golay cocode element :math:`\delta`, 
respectively. Column :math:`\langle d, \delta \rangle` lists 
conditions on the scalar product of  :math:`|d|` and  
:math:`|\delta|`. All this information can easily be derived from  
:cite:`Iva99` and  :cite:`Con85` (or :cite:`Seysen20`).

The table provides enough information for effectively computing
the subtype of an element :math:`\lambda_d + \lambda_\delta` from
:math:`d` and :math:`\delta`. Function ``gen_leech2_type`` in file 
``gen_leech.c`` computes that subtype.

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
    x8 & \mbox{Vectors $\lambda_d+\lambda_\delta$ with $|d| = 24$; 
    $\delta = 0$} \\ 
    \hline 
    \end{array}



Operation of the group :math:`G_{x0}` on the Leech lattice
..........................................................

The generators :math:`x_d, x_\delta, y_d, x_\pi, \xi`, with
:math:`d \in \mathcal{C}, \delta \in \mathcal{C}^*, \pi \in
\mbox{Aut}_{St} \mathcal{P}` operate on the subgroup 
:math:`G_{x0}` of the monster and also on 
:math:`\Lambda / 2 \Lambda`. Here :math:`x_d, x_\delta`
operate trivially on :math:`\Lambda / 2 \Lambda`.
:math:`y_d` and :math:`x_\pi` act as permutations and sign 
changes on the coordinates of :math:`\Lambda`, respectively; 
they do not change the subtype. 

From the *Leech graph* in :cite:`Iva99`, section 4.7, we see
that generators :math:`\xi` and :math:`\xi^2` may change
the subtype of a type-4 vector in  :math:`\Lambda / 2 \Lambda` 
as follows:

.. only:: html

    
        * subtype 48 to 40
        
        * subtype 40 to 40, 48, 42, 44
        
        * subtype 42 to 40, 42, 44, 43
        
        * subtype 44 to 40, 42, 44, 43, 46
        
        * subtype 43 to 42, 44, 43, 46
        
        * subtype 46 to 44, 43, 46


 

.. raw:: latex

    \begin{figure}[h]
    \setlength{\unitlength}{1cm}
    \centering
    \begin{picture}(10,6)
	\put(0.5,3){\oval(1.0, 0.75)}
	\put(0.3,2.9){$48$}
	\put(3.5,3){\oval(1.0, 0.75)}
    \put(3.3,2.9){$40$}
    \put(1.0,3){\line(1,0){2}}
    \put(4,3.1){\line(1,1){1.5}}
	\put(6.0,4.7){\oval(1.0, 0.75)}
    \put(5.8,4.6){$42$}
    \put(4,2.9){\line(1,-1){1.5}}
    \put(6.0,1.3){\oval(1.0, 0.75)}
    \put(6.0,1.75){\line(0,1){2.5}}
    \put(5.8,1.2){$44$}
    \put(6.5,1.5){\line(1,1){2.8}}s
	\put(9.5,4.7){\oval(1.0, 0.75)}
    \put(9.3,4.6){$43$}
	\put(9.5,1.3){\oval(1.0, 0.75)}
    \put(9.3,1.2){$46$}
    \put(6.5,4.7){\line(1,0){2.5}}
    \put(6.5,1.2){\line(1,0){2.5}}
    \put(9.5,1.75){\line(0,1){2.5}}   
    \end{picture}
    \caption{A line between two subtypes means that $\xi$ and $\xi^2$ 
    may map one subtype to the other subtype. Every subtype apart from
    48 may also be mapped to itself.}    
    \end{figure}


In the next few subsections we show how to map a vector of subtype
4x to a vector of subtype 48 by a sequence of generators
:math:`x_\pi, \xi^k` of length at most 6. The proofs are a bit 
involved, but the C function ``gen_leech2_reduce_type4_vector``
in file ``gen_leech.c`` constructs that sequence very fast.

We will provide some more information about the operation of 
:math:`\xi^k` on :math:`\Lambda`. 

We say that a vector :math:`w \in \mathbb{Z}^n` is of shape
:math:`(m^\alpha 0^{n-\alpha} \bmod 2m)` if :math:`w` has
:math:`\alpha` coordinates equal to :math:`m` :math:`\pmod{2m}`
and :math:`m-\alpha` coordinates equal to 
:math:`0` :math:`\pmod{2m}`.
For a vector :math:`v = (v_0,\ldots, v_{23}) \in \mathbb{R}^{24}`
define 
:math:`(v_{4i}, v_{4i+1}, v_{4i+2}, v_{4i+3}) \in \mathbb{R}^4` 
to be the :math:`i`-th column of :math:`v`. This definition is
related to the Miracle Octad Generator (MOG) used for the
description of the Golay code and the Leech lattice, see
:cite:`CS99`, chapter 11.

Lemma 

Let :math:`v \in \sqrt{8} \Lambda`, and let 
:math:`v^{(k)} = v^{\xi^k}`.
Let :math:`w` be a column of :math:`v`, and let :math:`w^{(k)}`
be the corresponding column of :math:`v^{(k)}`. Then
:math:`|w^{(k)}| = |w|`. If :math:`w` has shape
:math:`(m^4 \bmod{2m})` then there is a :math:`k \in \{1,2\}`
such that :math:`w^{(k)}` has shape :math:`(0^4 \bmod{2m})`.
If :math:`w` has shape :math:`(0^4 \bmod{2m})` or 
:math:`(m^2 0^2 \bmod{2m})` then :math:`w^{(k)}` has 
the same shape (modulo :math:`2m`) for all
:math:`k \in \mathbb{Z}`.

Sketch proof

In :cite:`Seysen20`, section 9.5 the operation of :math:`\xi` 
on a column is given as a product of two orthogonal symmetric
:math:`4 \times 4` matrices :math:`\xi_{4a}` and
:math:`\xi_{4b}`. Matrix :math:`\xi_{4b}` negates the first
coordinate of a column, and the operation of matrix 
:math:`\xi_{4a}` corresponds to the subtraction of the halved
sum of all coordinates from each coordinate. With this
description of :math:`\xi_{4a}` and :math:`\xi_{4b}` it is
easy to show the Lemma.

q.e.d.




From subtype 46 to subtype 44
.............................

Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of 
subtype 46. Then :math:`d` is a dodecad. There is a permutation
:math:`x_\pi` that maps :math:`d` to the dodecad
:math:`e = \{0,1,2,3,8,9,14,15,18,19,22,23\}`.
Function ``mat24_perm_from_dodecads`` in file ``mat24_functions.c`` 
can be used for constructing a suitable permutation. 
:math:`x_\pi` maps :math:`\delta` to an even cocode element
:math:`\epsilon` with :math:`\langle e, \epsilon \rangle = 1``.

There is a vector :math:`v \in \sqrt{8}\Lambda` with 
:math:`v/\sqrt{8} = \lambda_e +\lambda_\epsilon \pmod{2\Lambda}` 
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
    \end{array}
    
up to sign. By the Lemma given above  there is a 
:math:`k \in \{1,2\}` such that column 0 of the the vector 
:math:`w = v \cdot {\xi^k}` has shape :math:`(4 \, 0^3)`. The other
columns of :math:`w` have the same shape as the corresponding 
columns of :math:`v`. Thus  :math:`w` has shape 
:math:`(4 \, 2^8 0^{15})`. So using the table given above we see
that  :math:`w` has type 44. 

With a little extra effort we can show 
:math:`k = 2 - \langle f, \epsilon \rangle` for
:math:`f = \{0,1,2,3,4,5,6,7\}`, and 
:math:`\langle ., . \rangle: \mathcal{C} \times \mathcal{C}^*
\rightarrow \{0,1\}` the scalar product.



From subtype 43 to subtype 42
.............................


Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of 
subtype 43. Then :math:`|\delta| = 3`.There is a permutation
:math:`x_\pi` that maps :math:`\delta` to 
:math:`\epsilon = \{0,1,2\}`, and  :math:`d` to a Golay code 
word :math:`e`.

There is a vector :math:`v \in \sqrt{8}\Lambda` with 
:math:`v/\sqrt{8} = \lambda_e +\lambda_\epsilon \pmod{2\Lambda}` 
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
    \end{array}

up to sign. By the Lemma given above  there is a 
:math:`k \in \{1,2\}` such that one column of
:math:`w = v \cdot {\xi^k}` (and hence all columns) have even 
entries. Column 0 of :math:`w` has squared norm 
:math:`44 = 6^2 + 2^2 + 2^2 + 0^2`. That decompostion of
:math:`44` into a sum of four even squares is unique, so 
column 0 has shape :math:`(6 \, 2^2 0)`. The other columns
have shape  :math:`(2 \, 0^3)`. So :math:`w` is of
shape  :math:`(6 \, 2^7 0^{15})` and hence of subtype 42.
 
With a little extra effort we can show 
:math:`k = 2 - \langle e, \varphi \rangle` for
:math:`\varphi = \{0,1,2,3\} \in \mathcal{C}^*`, and 
:math:`\langle ., . \rangle` as above. 


From subtype 44 to subtype 40
.............................

Let :math:`\lambda_r = \lambda_d + \lambda_\delta` be of 
subtype 44. Then :math:`d` is an octad or a complement of
an octad. We call that octad :math:`d'`. There is a permutation 
:math:`x_\pi` that maps :math:`d'` to the standard octad
:math:`e = \{0,1,2,3,4,5,6,7\}`. But a suitable :math:`x_\pi`
also has to deal with :math:`\delta`.

If the cocode word :math:`\delta` is a duad, let :math:`c_0, c_1` 
be the two points in that duad. Otherwise, let :math:`c` by any 
tetrad intersecting the octad :math:`d'` in two points, so that
:math:`c` is equal to the sextet :math:`\delta` modulo the cocode
:math:`\mathcal{C}^*`. Let :math:`c_0, c_1` be the two points in  
:math:`c` which are not in the octad  :math:`d'`. 

Select 5 points from octad :math:`d'`, and let :math:`h` be the 
heptad formed by these 5 points and :math:`c_0, c_1`. Then 
:math:`h` is an umbral heptad (i.e. an heptad not contained in an 
octad); and :math:`h` has a unique point :math:`o_0` such that 
the difference set :math:`h \setminus o_0`  is a standard hexad,
i.e. a hexad contained in an octad. Furthermore, :math:`o_0`  
is in the octad  :math:`d'`.

Select 4 points :math:`o_1, o_2, o_3, o_4` (differerent from 
:math:`o_0`) from the octad  :math:`d'` and map the umbral 
heptad :math:`(o_1, o_2, o_3, o_4, c_0, c_1, o_0)` to the umbral 
heptad :math:`(0,1,2,3,8,9,4)` via a permutation :math:`x_\pi`. 
Then :math:`x_\pi` maps :math:`d'` to :math:`e` and 
:math:`\delta` to :math:`\epsilon \in \mathcal{C}^*` such that 
a representative of :math:`\epsilon` is a subset of
:math:`\{i \in \mathbb{Z} \mid 0 \leq i < 10\}`.
Function ``mat24_perm_from_heptads``
in file ``math24_functions.c`` constructs mappings from
umbral heptads to umbral heptads.

There is a vector :math:`v \in \sqrt{8}\Lambda` with 
:math:`v/\sqrt{8} = \lambda_e +\lambda_\epsilon \pmod{2\Lambda}` 
that has MOG coordinates
  
.. math::
    \begin{array}{|c|c|c|c|c|c|}
    \hline 
    2 & 2 & 4 & 0 & 0 & 0 \\
    \hline 
    2 & 2 & 4 & 0 & 0 & 0  \\ 
    \hline 
    2 & 2 & 0 & 2 & 2 & 2 \\ 
    \hline 
    2 & 2 & 0 & 2 & 2 & 2  \\ 
    \hline 
    \end{array}

up to sign. Similar to the previous cases, we can show
there is a :math:`k \in \{1,2\}` such that 
:math:`w = v \cdot {\xi^k}` has shape 
:math:`(4^4 0^{20})`. So :math:`w` is of subtype 40.

One can show :math:`k = 2 - \sigma`, where
:math:`\sigma = |c \cap \{0,1,2,3\}| \pmod{2}`, 
:math:`\sigma \in \{0,1\}`. 


From subtype 42 to subtype 40
.............................

Here a simplified version of the construction for the 
reduction of subtype 44 can be used. More specifically, we 
may ignore the contibution of the cocode part 
:math:`\lambda_\delta` to the permutation :math:`x_\pi` .

From subtype 40 to subtype 48
.............................

Let :math:`\lambda_r` be of subtype 40. Then
:math:`\lambda_r = \lambda_\delta` for some 
:math:`\delta \in \mathcal{C}^*`, :math:`|\delta|=4`.
A suitable permutation :math:`x_\pi` maps :math:`\delta` 
to the standard tetrad  :math:`\omega` represented by
:math:`\{0,1,2,3\}`. 
Then :math:`\lambda_r \cdot x_\pi = \alpha \lambda_\Omega + 
\lambda_\omega`. A suitable power of  :math:`\xi` maps that
vector to the unque vector :math:`\lambda_\Omega` of type 48.
