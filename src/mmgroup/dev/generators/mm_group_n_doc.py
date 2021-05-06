r"""We describe an implementation of the subgroup :math:`N_{0}` of the monster group.


The subgroup :math:`N_{0}` of the monster group of structure 
:math:`2^{2+11+2\cdot 11}.(\mbox{Sym}_3 \times\mbox{M}_{24})` has been 
described in :cite:`Con85`. Theorem 5.1 in :cite:`Seysen20` reduces the 
group operation in :math:`N_{0}` to easy calculations in the Parker loop 
:math:`\mathcal{P}`, the Colay cocode :math:`\mathcal{C}^*`, and the group
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` of standard automorphisms of
:math:`\mathcal{P}`. The loops :math:`\mathcal{P}`, :math:`\mathcal{C}^*`, 
and :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` are described in section
:ref:`basic_label`. Module ``mat24_functions.c`` provides the required
functions for computing in these loops.

Using the notation in section :ref:`mmgroup-label` we may describe an 
element of  :math:`N_{0}` as a product:

.. math::
    \tau^t y_f x_e x_\delta x_\pi \; ; \quad 0 \leq t < 3; \; 
    e, f \in   \mathcal{P}; \; \delta \in \mathcal{C}^*; \;
    \pi \in  {{\rm Aut}_{{\rm St}} \mathcal{P}} \; .


This representation is unique if we require :math:`f` to be in a
transversal of :math:`\mathcal{P} / Z(\mathcal{P})` and 
:math:`\pi` to be a standard representative in
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` as described in section 
:ref:`aut_ploop_label`.

We store an element of :math:`N_{0}` an a quintuple
:math:`t, f, e, \delta, \pi` of five integers of type ``uint32_t``.
For :math:`f, e` we use the numbering in class ``PLoop`` in module
``mmgroup``; for :math:`\delta` we use the numbering in class 
``Cocode`` in module ``mmgroup``. Here :math:`\pi` refers to a standard 
representative in :math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}`. These 
standard representatives correspond to the elements of the Mathieu 
group :math:`\mbox{M}_{24}`; and we use the numbering of the elements 
of :math:`\mbox{M}_{24}` described in class ``AutPL|`` in module 
``mmgroup``. 

Most functions in module ``mm_group_n.c`` take a pointer to a 5-tuple
:math:`(t, f, e, \delta, \pi)` representing an element :math:`g` 
of the :math:`N_{0}`  as their first argument. Then the tuple 
representing :math:`g` is modified  to a tuple representing an 
element :math:`g_2 = g \cdot g_1`, with the element :math:`g_1` of 
:math:`N` given by one or more subsequent arguments of the function.

The functions in module ``mm_group_n.c`` may cause some overhead due 
to the fact that element of the Mathieu group :math:`\mbox{M}_{24}`
is represented as an integer. But compared to an operation of the 
monster group on its 196884-dimensional that overhead is negligible. 


"""