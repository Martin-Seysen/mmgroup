

Computation in the Subgroup G_x0 of the monster
------------------------------------------------

According to Conway's construction of the monster  :math:`\mathbb{M}`
in :cite:`Con85`, the subgroup :math:`G_{x0}` of structure 
:math:`2^{1+24}_+.\mbox{Co}_1` has a faithful rational representation
on the tensor product :math:`4096_x \otimes 24_x`. Here :math:`24_x`
is the representation of the group :math:`\mbox{Co}_1` as the 
automorphism group of the real Leech lattice, and :math:`4096_x` is 
the representation of a group  :math:`G(4096_x)` as a subgroup of 
the real Clifford group :math:`\mathcal{C}_{12}`, where the
representation of :math:`\mathcal{C}_{12}` is as in section
:ref:`qstate_matrix_label`. Note that :math:`G(4096_x)`
is also of structure :math:`2^{1+24}_+.\mbox{Co}_1`.

We use the ideas in section :ref:`qstate_matrix_label` for 
implementing the representation :math:`4096_x`. In principle,
the representation :math:`24_x` is straightforward. The part
:math:`4096_x` of a representation of an element :math:`g` of
:math:`G_{x0}` determines :math:`g` up to sign. So :math:`g` is
quniquely determined if we store the image under :math:`g` of a 
single vector in the representation :math:`24_x`. The space
:math:`24_x` is equivalent to the Leech lattice :math:`\Lambda`, 
and we store we store the image of a fixed shortest vector
:math:`\beta` in :math:`\Lambda/ 3 \Lambda`, which is the Leech 
lattice modulo 3.


These ideas lead to a very fast implementation of the subgroup
:math:`G_{x0}` of the monster  :math:`\mathbb{M}`. Details of
that implementationa are given in 
the **C interface of the mmgroup project** in section 
:ref:`c-functions-G-x0-label`.

Class ``mmgroup.Xsp2_Co1`` provides an implementation  of the
group :math:`G_{x0}` based on these ideas.

