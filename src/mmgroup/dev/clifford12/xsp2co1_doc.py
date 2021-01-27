r"""



The functions in file ``xsp2co1.c`` implement the group operation 
of the subgroup 
:math:`G_{x0}` (of structure :math:`2^{1+24}.\mbox{Co}_1`)
of the monster. 

Represenation of :math:`G_{x0}` on the tensor product :math:`4096_x \otimes \Lambda`
....................................................................................

In :cite:`Seysen20`, section 7.4  and 9, the generators  
:math:`x_d, x_\delta, y_\delta, x_\pi, \xi` of :math:`G_{x0}`
are also defined as generators of a group  
:math:`G_{x1} \subset  G(4096_x) \times  G(24_x)`. Here 
:math:`G(4096_x)` is a subgroup of the real Clifford group
:math:`\mathcal{C}_{12}` operating on the 4096-dimensional 
real vector space :math:`4096_x`, and :math:`G(24_x)` is
the automorphism group  :math:`\mbox{Co}_0` of the Leech 
lattice :math:`\Lambda`. :math:`G_{x1}` is a preimage of 
:math:`G_{x0}` with :math:`|G_{x1}:G_{x0}| = 2`, and 
:math:`G_{x0}` operates faithfully on the tensor product 
:math:`4096_x \otimes \Lambda`. Let 
:math:`(x_g, L_g) \in G_{x1} \subset G(4096_x) \times  G(24x)`.
Component :math:`x_g` determines component :math:`L_g` of
the pair up to sign. So it suffices to store the image 
:math:`v_g = v_0 \cdot L_g` of a fixed shortest vector 
:math:`v_0 \in \Lambda / 3 \Lambda` instead 
of the whole automorphism :math:`L_g` of :math:`\Lambda`.
We put :math:`v_0 = (0,0,4,-4,0,...,0)`  in the standard basis 
of the Leech lattice.

We can reconstruct :math:`L_g` from   :math:`x_g` and 
:math:`v_g` as follows. :math:`G(4096_x)` has an extraspecial 
subgroup :math:`Q_{x1}` of structure :math:`2^{1+24}`. The 
quotient of :math:`Q_{x1}` by its center :math:`Z(Q_{x1})` is 
isomorphic to :math:`\Lambda / 2 \Lambda`. By definition of 
:math:`G_{x1}`, the element :math:`L_g` operates on 
:math:`\Lambda / 2 \Lambda` in the same way as  :math:`x_g` 
operates on :math:`Q_{x1}/Z(Q_{x1})` by conjugation.
So that operation of :math:`L_g` can be reconstructed from
:math:`x_g`. A short vector in  :math:`\Lambda / 2 \Lambda`
has precisely two short preimages  in  :math:`\Lambda` which
are opposite. Thus the image of one short vector is known
exactly, and the images of all short vectors are known up to 
sign. Since :math:`L_g` preserves the scalar product, the 
images of all short vectors in :math:`\Lambda` can be computed. 
Function ``xsp2co1_elem_to_leech_op`` in file ``xsp2co1.c``
constructs :math:`L_g`  from  :math:`x_g` and :math:`v_g`.
The functions in file ``gen_leech.c`` support operations on
:math:`\Lambda / 2 \Lambda` and :math:`\Lambda / 3 \Lambda`. 

The groups  :math:`G_{x0}` also has an extraspecial subgroup
:math:`Q_{x0}` of structure :math:`2^{1+24}`.

Embedding a group related to :math:`G_{x0}` into a Clifford group
.................................................................

The representation :math:`4096_x` of the subgroup  
:math:`G(4096_x)` of :math:`\mathcal{C}_{12}` extends to 
the standard representation of the Clifford group 
:math:`\mathcal{C}_{12}`. So we may use the theory in 
section :ref:`clifford-group-label` for computing in 
:math:`G(4096_x)`. In that section the basis vectors of the
standard representation  of :math:`\mathcal{C}_{12}` are 
numbered from :math:`0` to :math:`2^{12}-1`. We will 
assign compatible numbers to the basis vectors of the
representation :math:`4096_x` . 


Let  :math:`\mathcal{P}` be the Parker loop as defined in 
section :ref:`parker-loop-label`. The vector space 
:math:`4096_x` has basis vectors 
:math:`d_1^+, d_1^-, d \in \mathcal{P}`, with  relations 
:math:`(-d)_1^\pm = -d_1^\pm`, :math:`(\Omega  d)_1^+ = d_1^+`, 
:math:`(\Omega  d)_1^- = -d_1^-`. 
In the same section the elements of :math:`\mathcal{P}` are 
numbered from :math:`0` to :math:`2^{13}-1`. For 
:math:`0 \leq d < 2^{11}` we identify the 
basis vectors :math:`d_1^+` and :math:`d_1^-` of  :math:`4096_x`
with the basis vectors with numbers corresponding to :math:`d` 
and to :math:`d+2^{11}` of the representation of 
:math:`\mathcal{C}_{12}`, respectively. 



So we may represent an element :math:`x_g` of :math:`G(4096_x)`
as a **quadratic state matrix** in a structure of type 
``qstate12_type``, as defined in file ``clifford12.h``. Then
we may use the functions in file ``qmatrix12.c`` for computing
in :math:`G(4096_x)`.



The :math:`G_{x0}` representation 
.................................. 

We actually represent an element of :math:`G_{x0}` in **G_x0 
representation**. This is as an array ``elem`` of 26 integers 
of type ``uint64_t``. That array contains a pair
:math:`(x_g,v_g)` as described above with :math:`v_g` stored 
in the first entry ``elem[0]``, and  :math:`x_g` stored in the 
remaining 25 entries of array ``elem``.

Vector :math:`v_g \in \Lambda / 3 \Lambda` is stored in
``elem[0]`` as a 48-bit integer in **Leech lattice mod 3 
encoding** as described in the documentation of module 
``gen_leech.c``.

We store the inverse :math:`x_g^{-1}` of  :math:`x_g` in
the upper 25 entries of the array ``elem``. The rationale
for storing the inverse is that this simplifies the operation
of :math:`G(4096)_x` on :math:`Q_{x1}` by conjugation, using
function ``qstate12_pauli_conjugate`` in file ``qmatrix12.c``.
Note that :math:`x_g^{-1} = x_g^\top`, since  :math:`x_g`
is orthogonal.


In section :ref:`clifford-group-label` an element :math:`c`
of :math:`\mathcal{C}_{12}` is given as a real
:math:`2^{12} \times 2^{12}` matrix stored in a structure
of type ``qstate12_type``. In case :math:`c \in G(4096_x)`
that matrix has rational entries, where the denominators are 
powers of two. A structure of type ``qstate12_type`` 
representing a  :math:`c \in G(4096_x)` contains 
a triple :math:`(e,A,Q)`. There :math:`e` is a signed integral  
power of two, and :math:`(A,Q)` is a pair of bit 
matrices with up to 25 rows and up to 49 colums, where the 
first 24 columns belong to matrix :math:`A`, and the remaining 
columns belong to  matrix :math:`Q`.


Changing the sign of :math:`e` corresponds to negation of a 
matrix :math:`x_g^{-1} \in G(4096_x)` given by :math:`(e,A,Q)`
or to  multiplication by :math:`x_{-1}`. In the group 
:math:`G_{x0}`, changing the sign of :math:`v_g` has the same 
effect as changing the sign of :math:`x_g`; so  :math:`e` can 
always be made positive. The absolute value of :math:`e` is 
just there for scaling the operator norm of a matrix 
:math:`c` to 1, and can be omitted.

So we can represent one of the values :math:`\pm x_g^{-1}` as 
a pair :math:`(A,Q)`  with an implied positive scalar factor 
:math:`e` as above, and negate component :math:`v_g` of the 
pair :math:`(x_g,v_g)`, if necessary. 
 
Assuming that  :math:`x_g^{-1}` is stored in a structure 
``qs`` of type ``qstate_type``, we copy all valid entries 
``qs.data[i]`` to ``elem[i+1]``. This amounts to copying the
components  :math:`(A,Q)` of the triple  :math:`(e,A,Q)`
from ``qs`` to ``elem``.  We always reduce the quadratic
state matrix in ``qs`` before copying it, as indicated in 
section :ref:`clifford-group-label`. We fill unused entries 
and bits in the array ``elem`` with zeros. Thus  the 
representation of a :math:`g \in G_{x0}` in memory is unique.
``qs`` can easily be reconstructed from ``elem``, assuming
that the scalar factor :math:`e` in the triple
:math:`(e,A,Q)` is positive, and that trailing zero 
entries in the array ``elem`` are unused.

The normal subgroup :math:`Q_{x0}` of :math:`G_{x0}`
.....................................................


The group :math:`G_{x0}` has an extraspecial normal subgroup 
:math:`Q_{x0}` of structure :math:`2^{1+24}`.
We also represent an element of :math:`Q_{x0}` as a 25-bit
integer in **Leech lattice encoding** as described in 
section :ref:`mmgroup-generators-label`. 
We also represent an element of the subgroup :math:`Q_{x1}`
of :math:`G(4096_x)` in that way.


Changing the basis of the space :math:`4096_x`
...............................................


We also use the vectors :math:`(d'), d \in \mathcal{P}` as basis 
vectors of :math:`4096_x` , where:

.. math::
    (d)' = \frac{1}{\sqrt2} \left( d_1^+ + d_1^- \right) \, , \quad 
    d_1^\pm = \frac{1}{\sqrt2} \left( (d)' \pm (\Omega d)' \right) \, . 

    
The reason for introducing this basis is that the operation of 
:math:`G_{x0}` on  :math:`Q_{x0}` (by conjugation) is
easy when we use this basis.
Then we obviously have :math:`(-d)' = -(d)'`. We can transform 
the coordinates of a vector :math:`v` from the basis  given by
:math:`d^\pm`  to the basis given by  :math:`(d)'`  
by multiplying the coordinate vector with a matrix :math:`H`. 
Multiplication with matrix :math:`H` is equivalent to
applying the Hadamard gate to Qubit 11 (i.e. to the bit with 
valence :math:`2^{11}`) of the coordinates of  :math:`v` .  
Thus :math:`H` is an involution. 

The numbering of the 4096 basis vectors :math:`(d)'` corresponds 
to the numbering of the positive elemments of :math:`\mathcal{P}`. 
We may put :math:`(d \oplus 2^{12})' = -(d)'`, where 
:math:`\oplus` means bitwise addition modulo 2; then that 
correspondence holds for all values :math:`0 \leq d < 2^{13}`. 
The exact defnition of the operator :math:`\oplus` on the Parker 
loop :math:`\mathcal{P}` is given in section  :ref:`implement-gen-mm`.

In this basis the operation of the extraspecial group 
:math:`Q_{x0}` is very simple:

    * :math:`x_e x_{\theta(e)}` maps  :math:`(d)'` to   
      :math:`(d \oplus e)'` for :math:`e \in \mathcal{P}`.
      Here :math:`\theta` is the cocycle, see section 
      :ref:`parker-loop-label`. In the language of 
      quantum computing this corresponds to a sequence
      of commuting **not** gates.
      
    * :math:`x_\epsilon` maps    :math:`(d)'` to   
      :math:`(-1)^{\langle d, \epsilon\rangle}(d)'`
      for :math:`\epsilon \in \mathcal{C}^*`. Here 
      :math:`\mathcal{C}^*` is the Golay cocode. In 
      the language of quantum computing this corresponds 
      to a sequence of commuting **phase**-:math:`\pi`
      gates.


Since :math:`H` is in :math:`\mathcal{C}_{12}`, it operates on 
:math:`Q_{x0}` by conjugation. The following subtlety has to be
considered in function ``conv_pauli_vector_xspecial`` (in 
file ``xsp2cco1.c``) implementing that operation. 
:math:`H` exchanges the anticcommuting
elements of :math:`Q_{x0}` of with numbers ``0x800`` and  
``0x800000``  in Leech lattice encoding. Thus it negates 
the element with number ``0x800800``.
 



     
"""
