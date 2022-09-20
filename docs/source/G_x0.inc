

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

We remark that the subgroup :math:`G_{x0}` of the monster has
five classes of involutions that map to 2B involutions in the
monster, see :cite:`Wilson13`. From Table 2 in :cite:`Nor98` we 
conclude that :math:`G_{x0}` has two classes of involutions
that map to 2A involutions in the monster. 
According to :cite:`Con85`, the restriction of the
196883-dimensional representation of the monster to
:math:`G_{x0}` leads to a representation of :math:`G_{x0}`
of shape 

.. math::
    299_x \oplus 98280_x \oplus 24_x \otimes 4096_x \, .

We can compute the characters of all these representations
of :math:`G_{x0}`, e.g. with the C function 
``xsp2co1_traces_fast`` in file ``xsp2co1_traces.c``, or with
function ``xsp2co1_traces_all`` in file ``xsp2co1_elem.c``.
The calculations in the python module
``mmgroup.tests.test_involutions.make_involution_samples``
show that the class of an involution in the group :math:`G_{x0}`
can be identified by computing the four characters given above. 


We use the ideas in section :ref:`qstate_matrix_label` for 
implementing the representation :math:`4096_x`. In principle,
the representation :math:`24_x` is straightforward. The part
:math:`4096_x` of a representation of an element :math:`g` of
:math:`G_{x0}` determines :math:`g` up to sign. So :math:`g` is
uniquely determined if we store the image under :math:`g` of a 
single vector in the representation :math:`24_x`. The space
:math:`24_x` is equivalent to the Leech lattice :math:`\Lambda`, 
and we store the image of a fixed shortest vector
:math:`\beta` in :math:`\Lambda/ 3 \Lambda`, which is the Leech 
lattice modulo 3, as follows.

Let :math:`\phi` be the natural homomorphism from the subgroup 
:math:`Q_{x0}` of :math:`G_{x0}` of structure :math:`2^{1+24}` onto  
:math:`\Lambda/ 2 \Lambda`, which is the Leech lattice modulo 2.

Let :math:`g \in G_{x0}` be represented as a pair
:math:`(g_1, g_2) \in G(4096_x) \times \mbox{Co}_1`.
Then :math:`g_2` can be reconstructed from the image
:math:`\beta \cdot g_2` of a single shortest vector :math:`\beta`
in the Leech lattice (modulo 3).

The operation of :math:`g_2` on a shortest vector 
:math:`\beta' \in \Lambda` is given by


.. math::
    \beta' \mapsto 
    \phi(g^{-1} \cdot \phi^{-1}(\beta' + 2 \Lambda) \cdot g)
    \in \Lambda / 2 \Lambda,
    

up to sign. Here a short vector in :math:`\Lambda / 2 \Lambda` 
has precisely two opposite shortest preimages in :math:`\Lambda`,
and we may take an arbitrary preimage in :math:`Q_{x0}` when
applying the inverse :math:`\phi^{-1}` of  :math:`\phi`. 

We compute the image :math:`\beta' \cdot g_2` of an arbitrary 
shortest vector :math:`\beta' \in \Lambda` as follows. Assume first
that :math:`\beta, \beta'` are shortest vectors in :math:`\Lambda` 
which are not perpendicular to each other, and assume that the
shortest Leech lattice vector :math:`\beta \cdot g_2` is known for
some  :math:`g_2 \in \mbox{Co}_1`. If  :math:`\beta' \cdot g_2` 
is known up to sign only, we can compute the sign of :math:`\beta'`
from the scalar products :math:`(\beta, \beta')`, and 
:math:`(\beta \cdot g_2, \beta' \cdot g_2)`, which must be equal.
Since  :math:`|(\beta, \beta')| \in \{1,2,4\}`, it suffices to 
compute the scalar products modulo 3. For perpendicular shortest
vectors  :math:`\beta, \beta'` we can easily find a shortest vector
:math:`\beta'' \in \Lambda` with :math:`(\beta, \beta'') \neq 0`
and :math:`(\beta'', \beta') \neq 0`. So we may compute first
:math:`\beta'' \cdot g_2` and then :math:`\beta' \cdot g_2` 
from these scalar products.

These ideas lead to a very fast implementation of the subgroup
:math:`G_{x0}` of the monster  :math:`\mathbb{M}`. Details of
that implementation are given in 
the **C interface of the mmgroup project** in section 
:ref:`c-functions-G-x0-label`.

Class ``mmgroup.Xsp2_Co1`` provides an implementation  of the
group :math:`G_{x0}` based on these ideas.

