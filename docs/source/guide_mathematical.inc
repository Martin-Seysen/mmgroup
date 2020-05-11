

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
  











