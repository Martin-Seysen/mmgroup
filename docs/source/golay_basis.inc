
We chose a basis ``b_0``,..., ``b_11`` of the Golay code 
:math:`\mathcal{C}` such that the resulting Golay
code is compatible to :cite:`CS99`, Chapter 11. 

As usual in the Python or C language we number the indices 
of the basis vectors of the space :math:`\mathbb{F}_2^{24}`
from   ``0`` to ``23``. Similarly, we number the basis vectors 
of the Golay code from  ``0`` to ``11``. 

In  :cite:`CS99` the Golay code is described in terms
of the MOG (Miracle Octad Generator), which can be considered as
a :math:`4 \times 6` matrix for arranging the basis vectors 
of  :math:`\mathbb{F}_2^{24}`. Our numbering of these basis 
vectors corresponds to the MOG as follows:


    ====== ====== ====== ====== ====== ====== 
       0      4      8     12     16     20    
       1      5      9     13     17     21    
       2      6     10     14     18     22    
       3      7     11     15     19     23    
    ====== ====== ====== ====== ====== ====== 

For the cocode :math:`\mathcal{C}^*` we use the reciprocal
basis ``c_0``,..., ``c_11`` of the basis ``b_0``,..., ``b_11``
of :math:`\mathcal{C}`.  Then the scalar product 
``b_i & c_j`` of ``b_i`` and ``c_j`` is equal to ``1`` if
``i == j`` and to ``0`` otherwise. 
 


.. include:: golay_basis_data.inc

The cocycle 
:math:`\theta:  \mathcal{C} \times \mathcal{C} \rightarrow \mathbb{F}_2`
is used for the multiplication of two elements of the Parker loop,
see section :ref:`parker-loop-label`. 
For basis vectors  :math:`b_i`,
:math:`b_j` of :math:`\mathcal{C}` the column ``theta`` in the 
table above contains the value :math:`\theta(b_i, b_j)` in bit 
:math:`j` of row  :math:`i`. From these values other values of the 
cocycle :math:`\theta` can be computed by using


.. math::
     \theta(g_1 + g_2, g_3)   =   \theta(g_1, g_3) + \theta(g_2, g_3)
     + |g_1 \& g_2 \& g_3|,   \\
     \theta(g_1,  g_2 + g_3)  =  \theta(g_1, g_2) + \theta(g_1, g_3),
      \quad  g_1, g_2, g_3 \in  \mathcal{C} \; .


Here :math:`\&` means the bitwise ``and`` operation of bit 
vectors and :math:`|.|` means the bit weight of a bit vector.


The cocycle :math:`\theta` has been chosen in accordance with 
the requirements in  :cite:`Seysen20`. In that context the
*grey* subspace of :math:`\mathcal{C}` is spanned by 
:math:`b_0,\ldots,b_3,b_{10}, b_{11}` and the 
*coloured* subspace of :math:`\mathcal{C}` is spanned by 
:math:`b_4,\ldots,b_9`. Similarly, the
*grey* subspace of :math:`\mathcal{C}^*` is spanned by 
:math:`c_0,\ldots,c_3,c_{10}, c_{11}` and the 
*coloured* subspace of :math:`\mathcal{C}^*` is spanned by 
:math:`c_4,\ldots,c_9`.

These specific properties of our chosen basis and cocycle 
are vital for obtaining an effective implementation
of a complete generating set of the monster :math:`\mathbb{M}`.

