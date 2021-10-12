
.. include:: definitions.rst

=========================
The mmgroup API reference
=========================


Introduction
============

The *mmgroup* package is a python implementation of Conway's 
construction :cite:`Con85` of the monster group :math:`\mathbb{M}`. 
In mathematics the *monster group* :math:`\mathbb{M}` is the largest 
sporadic finite simple group. The group :math:`\mathbb{M}` has order

   :math:`2^{46} \cdot 3^{20} \cdot 5^9 \cdot 7^6 \cdot 11^2 \cdot 13^3 \cdot 17 \cdot 19 \cdot 23 \cdot 29 \cdot 31 \cdot 41 \cdot 47 \cdot 59 \cdot 71` 
   = :math:`\small 808.017.424.794.512.875.886.459.904.961.710.757.005.754.368.000.000.000` . 

It has a rational representation :math:`\rho` of dimension 
:math:`196884` with coefficients in :math:`\mathbb{Z}[\frac{1}{2}]`,
see :cite:`Con85`. 
So we obtain a modular representation :math:`\rho_p` of 
:math:`\mathbb{M}` by reducing the rational representation 
:math:`\rho` modulo an arbitrary odd number :math:`p`. The current
version of the *mmgroup* package supports the representations 
:math:`\rho_p` of :math:`\mathbb{M}` for 
:math:`p = 3, 7, 15, 31, 127, 255`.

The *mmgroup* package uses highly optimized C programs for 
calculating in the representations of :math:`\rho_p` of
:math:`\mathbb{M}`, which are most efficient if :math:`p+1` 
is a small power of two. Most of these C programs have been
generated automatically by a program code generator which is
also part of this project, but not part of the public 
interface described in this document.

In the description of the *mmgroup* package we use the
notation in :cite:`Seysen20`, where an explicit description
of the representations :math:`\rho(g), g \in G` is given
for a generating subset  :math:`G` of 
:math:`\mathbb{M}`.


Installation
============

.. include:: installation.inc



.. _basic_label:

Basic structures
================

.. _golay-label:

The Golay code and its cocode
-----------------------------

.. automodule:: mmgroup.structures.gcode


.. autoclass:: mmgroup.GCode
   :members:

.. autoclass:: mmgroup.GcVector
   :members:

.. autoclass:: mmgroup.Cocode
   :members:

.. autoclass:: mmgroup.Parity
   :members:



.. _parker-loop-label:


The Parker loop
---------------


.. automodule:: mmgroup.structures.ploop

.. autoclass:: mmgroup.PLoop
   :members:

.. autofunction:: mmgroup.PLoopZ


.. _octads_label:

Octads and suboctads
--------------------

.. automodule:: mmgroup.structures.suboctad

.. autofunction:: mmgroup.Octad

.. autoclass:: mmgroup.SubOctad
   :members:

.. _aut_ploop_label:

Automorphisms of the Parker loop
--------------------------------

.. automodule:: mmgroup.structures.autpl


.. autoclass:: mmgroup.AutPL
   :members: 


.. _xleech2-label:


The Leech lattice modulo 2 and its preimage :math:`Q_{x0}`
-----------------------------------------------------------

.. automodule:: mmgroup.structures.xleech2

.. autoclass:: mmgroup.XLeech2
   :members: 


.. _basis-golay-label:

The basis of the Golay code and of its cocode
---------------------------------------------


.. include:: golay_basis.inc



.. _mmgroup-label:

The monster group
=================


.. automodule:: mmgroup.mm_group


Python classes implementing the monster group
--------------------------------------------- 


.. autoclass:: mmgroup.MM
   :members:  as_tuples, copy,  is_reduced, reduce, order, mmdata,
              half_order, in_G_x0, chi_G_x0, in_N_x0, in_Q_x0,
              as_Q_x0_atom, type_Q_x0, simplify


.. _mmrep-label:

The representation of the monster group
=======================================


.. automodule:: mmgroup.mm_space


Python classes implementing the representation of the monster group
------------------------------------------------------------------- 


.. autoclass:: mmgroup.MMVector
   :members:  as_sparse, as_tuples, projection,
              mul_exp

.. autoclass:: mmgroup.MMSpace
   :members:  tuple_to_index, index_to_tuple, index_to_short




Auxiliary functions for the representation of the monster group
--------------------------------------------------------------- 

.. autofunction:: mmgroup.characteristics

.. autofunction:: mmgroup.MMV


.. _clifford-group-label:

The subgroup :math:`2^{1+24}.Co_1` of the monster and the Clifford group
========================================================================

The section describes the fast computation in a certain subgroup 
:math:`G_{x0}` of the monster :math:`\mathbb{M}` in detail. A person
who simply wants do do calculations in the monster group need not 
read this section. 


Introduction
------------


In Conway's construction :cite:`Con85` the monster :math:`\mathbb{M}`
has a subgroup :math:`G_{x0}` of structure 
:math:`2^{1+24}_+.\mbox{Co}_1`.
There :math:`G_{x0}` is constructed as a diagonal product of the
two groups :math:`\mbox{Co}_0` of structure :math:`2.\mbox{Co}_1`
and :math:`G(4096_x)`. :math:`G(4096_x)` is also of structure
:math:`2^{1+24}_+.\mbox{Co}_1` but not isomorphic to  :math:`G_{x0}`.
Here :math:`2^{1+24}_+` means an extraspecial group of plus
type. :math:`\mbox{Co}_0`, and :math:`\mbox{Co}_1` are the automorphism
groups of the :math:`24`-dimensional Leech lattice and of the Leech 
lattice modulo :math:`2`, see e.g. :cite:`Atlas` or :cite:`CS99` for 
background. 

Computation in  :math:`\mbox{Co}_0` is easy since that group has a 
:math:`24`-dimensional rational representation. The smallest real
representation of the group :math:`G(4096_x)` has dimension
:math:`4096`, so naive computation in that representation is
rather inefficient.

The group :math:`G(4096_x)` is a subgroup of the real Clifford group
:math:`\mathcal{C}_{12}`. The real Clifford group :math:`\mathcal{C}_{n}`
of structure :math:`2^{1+2n}_+.\mbox{GO}^+_{2n}(2)` is defined e.g. in
:cite:`NRS01`. :math:`\mathcal{C}_{12}` is a subgroup
of the complex Clifford group :math:`\mathcal{X}_{n}`  of structure 
:math:`\frac{1}{2} ( 2^{1+2n}_+ \times Z_8 ). \mbox{Sp}_{2n}(2)`, 
which is also defined in :cite:`NRS01`.  Here a factor 
:math:`\frac{1}{2}` means identification of central subgroups of 
order :math:`2` of the factors of a direct product. 
:math:`\mbox{GO}^+_{2n}(2)` and :math:`\mbox{Sp}_{2n}(2)` are the
orthogonal subgroup (of :math:`+`-type) and the symplectic
subgroup of the group :math:`GL_{2n}(2)` of invertible 
:math:`2n \times 2n`-matrices with entries in :math:`\mathbb{F}_2`,
as defined in :cite:`Atlas`.

Effective computation in the group  :math:`\mathcal{X}_{n}` has received
a lot of attention from the theory of quantum computation, see e.g.
:cite:`AG04`. In the next section we present efficient algorithms for
computing in :math:`\mathcal{X}_{n}` and in its :math:`2^{n}` 
dimensional complex representation.




.. automodule:: mmgroup.dev.clifford12.qstate12


Class ``QStateMatrix`` modelling a quadratic state matrix
---------------------------------------------------------


.. autoclass:: mmgroup.structures.qs_matrix.QStateMatrix
   :members: conjugate, T, H, shape, reshape, copy,
             rot_bits, xch_bits,
             gate_not, gate_ctrl_not, gate_phi, gate_ctrl_phi, gate_h,
             extend_zero, extend, restrict_zero, restrict, sumup,
             lb_rank, lb_norm2, inv, power, order,
             pauli_vector, pauli_conjugate, show, trace

.. autofunction:: mmgroup.structures.qs_matrix.qs_unit_matrix






.. only:: html

   .. rubric:: **References**

.. bibliography:: references.bib



