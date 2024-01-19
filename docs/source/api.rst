
.. include:: definitions.rst

=========================
The mmgroup API reference
=========================


Introduction
============

In the area of mathematics known as group theory, the Monster group 
:math:`\mathbb{M}` is the largest finite sporadic simple group.
It has order 

   :math:`2^{46} \cdot 3^{20} \cdot 5^9 \cdot 7^6 \cdot 11^2 \cdot 13^3 \cdot 17 \cdot 19 \cdot 23 \cdot 29 \cdot 31 \cdot 41 \cdot 47 \cdot 59 \cdot 71` 
   = :math:`\small 808.017.424.794.512.875.886.459.904.961.710.757.005.754.368.000.000.000` . 

The Monster group has first been constructed by Griess :cite:`Gri82`.
That construction has been simplified by Conway :cite:`Con85`. 
For more information about the Monster group we refer to
:cite:`wiki:monster`.

The *mmgroup* package is a python implementation of Conway's 
construction :cite:`Con85` of the Monster group :math:`\mathbb{M}`. 
Its is the first implementation of the Monster group where arbitrary
operations in that group can effectively be performed. On the author's 
PC (Intel i7-8750H at 4 GHz running on 64-bit Windows) 
the group multiplication in :math:`\mathbb{M}` takes less than 30 ms.
This is more than five orders of magnitude faster than estimated 
in 2013 in :cite:`Wilson13`.

The Monster group :math:`\mathbb{M}` has a rational representation 
:math:`\rho` of dimension  :math:`196884`, see :cite:`Con85`. In that 
representation the denominators of the matrix coefficients are powers 
of two. So reducing these coefficients modulo a small odd prime 
:math:`p` leads to a representation :math:`\rho_p` of 
:math:`\mathbb{M}` over the finite field :math:`\mathbb{F}_p`. 
 
The *mmgroup* package uses highly optimized C programs for 
calculating in such representations :math:`\rho_p` of the Monster
:math:`\mathbb{M}`. The main ingredient for speeding up the
computation in :math:`\mathbb{M}` is the calculation and the
analysis of the images of certain vectors in :math:`\rho_p`
that are called 2A axes in :cite:`Con85`.

In the description of the *mmgroup* package we use the notation
in :cite:`Seysen20`, where an explicit generating set of the
Monster :math:`\mathbb{M}` is given. For a mathematical description
of the implementation we refer to :cite:`Seysen22`.


Installation and test
=====================

.. include:: installation.inc





.. _basic_label:

Basic structures
================

Conway's construction of the Monster group starts with the extended 
binary Golay code :math:`\mathcal{C}`, which is a 12-dimensional 
linear subspace of the 24-dimensional vector space 
:math:`\mathbb{F}_2^{24}`. The Golay code has Hamming distance 8. 
Its automorphism group is the Mathieu group :math:`M_{24}`, which is
a simple group operating as a permutation group on 24 points. 
These points are identified with the basis vectors of 
:math:`\mathbb{F}_2^{24}`.

Module ``mmgroup`` contains fast C routines for computing with the
Golay code  :math:`\mathcal{C}`, its cocode :math:`\mathcal{C^*}`,
and the Mathieu group :math:`M_{24}`. There are Python classes 
for a more convenient handling of these objects that wrap these
C functions; these classes are described in this section.

In this section we also describe Python classes modelling more
complicated mathematical structures. One of these structures is
the Parker loop :math:`\mathcal{P}`, which is a non-associative 
loop that can be constructed as a double cover of the Golay code
:math:`\mathcal{C}`. We also consider a group 
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` of automorphisms of
:math:`\mathcal{P}`, which has structure 
:math:`2^{12}.M_{24}`. Here the normal subgroup of 
:math:`{{\rm Aut}_{{\rm St}} \mathcal{P}}` or order :math:`2^{12}`
is isomorphic to the Golay cocode  :math:`\mathcal{C^*}`.

Another important ingredient of the construction of the Monster
is the Leech lattice :math:`\Lambda`, which is the densest lattice
in dimension 24. We also consider the Leech lattice modulo 2,
which we denote by :math:`\Lambda / 2 \Lambda`, and the automorphism
group  :math:`\mbox{Co}_1` of :math:`\Lambda / 2 \Lambda`, which is
simple.




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

.. autofunction:: mmgroup.SubOctad

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


.. autofunction:: mmgroup.leech2_orbits_raw



.. _basis-golay-label:

The basis of the Golay code and of its cocode
---------------------------------------------


.. include:: golay_basis.inc



.. _mmgroup-label:

The Monster group
=================



.. automodule:: mmgroup.mm_group


Python classes implementing the Monster group
--------------------------------------------- 

.. _monster-classes-label:

.. autoclass:: mmgroup.MM
   :members:  as_tuples, copy,  is_reduced, reduce, order, mmdata,
              half_order, in_G_x0, chi_G_x0, in_N_x0, in_Q_x0,
              half_order_chi, chi_powers, conjugate_involution,
              conjugate_involution_G_x0, as_int,
              as_Co1_bitmatrix


Functions dealing with elements of the Monster group
------------------------------------------------------

.. autofunction:: mmgroup.MM_from_int


.. _random_MM_label:

Generating random elements of certain subgroups of the Monster
--------------------------------------------------------------

.. include:: random_MM.inc


.. _mmrep-label:

The representation of the Monster group
=======================================


.. automodule:: mmgroup.mm_space


Python classes implementing the representation of the Monster group
------------------------------------------------------------------- 


.. autoclass:: mmgroup.MMVector
   :members:  as_sparse, as_tuples, projection,
              mul_exp

.. autoclass:: mmgroup.MMSpace
   :members:  tuple_to_index, index_to_tuple, index_to_short




Auxiliary functions for the representation of the Monster group
--------------------------------------------------------------- 

.. autofunction:: mmgroup.characteristics

.. autofunction:: mmgroup.MMV

.. autofunction:: mmgroup.mmv_scalprod

.. autofunction:: mmgroup.order_vector

.. _clifford-group-label:

The subgroup :math:`G_{x0}` of the Monster and the Clifford group
========================================================================

The section describes the fast computation in a certain subgroup 
:math:`G_{x0}` of structure :math:`2^{1+24}.Co_1` of the Monster 
:math:`\mathbb{M}` in detail. A person who simply wants do do 
calculations in the Monster group need not read this section. 


Introduction
------------


In Conway's construction :cite:`Con85` the Monster :math:`\mathbb{M}`
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


.. include:: qstate12.inc

.. _qstate_matrix_label:


Class ``QStateMatrix`` modelling a quadratic state matrix
---------------------------------------------------------


.. autoclass:: mmgroup.structures.qs_matrix.QStateMatrix
   :members: conjugate, T, H, shape, reshape, copy,
             rot_bits, xch_bits,
             gate_not, gate_ctrl_not, gate_phi, gate_ctrl_phi, gate_h,
             extend_zero, extend, restrict_zero, restrict, sumup,
             lb_rank, lb_norm2, inv, power, order, to_symplectic,
             pauli_vector, pauli_conjugate, show, trace

.. autofunction:: mmgroup.structures.qs_matrix.qs_unit_matrix


.. _group_g_x0_label:

.. include:: G_x0.inc


Class ``Xsp2_Co1`` modelling an element of the group :math:`G_{x0}`
-------------------------------------------------------------------

.. autoclass:: mmgroup.structures.xsp2_co1.Xsp2_Co1
   :members: order, type_Q_x0, conjugate_involution, subtype,
             conjugate_involution_G_x0, chi_G_x0,




.. _y555_label:


The Coxeter group :math:`Y_{555}` and the Bimonster
=============================================================

.. automodule:: mmgroup.bimm.readme


Implementing the Bimonster and the homomorphism from :math:`Y_{555}` to it
---------------------------------------------------------------------------

The projective plane over :math:`\mathbb{F}_3` and its automorphism group
..........................................................................

.. automodule:: mmgroup.bimm.inc_p3


.. autoclass:: mmgroup.bimm.P3_node
   :members:  y_name, name, ord


.. autoclass:: mmgroup.bimm.AutP3
   :members:  order, map, point_map, line_map

.. autofunction:: mmgroup.bimm.P3_incidences 

.. autofunction:: mmgroup.bimm.P3_incidence

.. autofunction:: mmgroup.bimm.P3_remaining_nodes

.. autofunction:: mmgroup.bimm.P3_is_collinear


Norton's presentation of the Monster group
...........................................


.. automodule:: mmgroup.bimm.p3_to_mm

.. autofunction::  mmgroup.bimm.Norton_generators


The Bimonster and its presentation  :math:`Y_{555}`
.....................................................

.. automodule::  mmgroup.bimm.bimm


.. autoclass::  mmgroup.bimm.BiMM
   :members:  order, decompose

.. autofunction::  mmgroup.bimm.P3_BiMM


.. autofunction::  mmgroup.bimm.AutP3_BiMM

Example: Checking the spider relation in the group :math:`Y_{555}` 
.......................................................................

.. automodule:: mmgroup.bimm.readme_example


Construction of the mapping from the Coxeter group into the Bimonster
---------------------------------------------------------------------

.. automodule:: mmgroup.bimm.readme_math





Version history
===============


.. list-table:: List of releases
   :widths: 11 17 72
   :header-rows: 1

   * - Version
     - Date
     - Description

   * - 0.0.1 
     - 2020-05-20 
     - First release

   * - 0.0.2
     - 2020-06-04
     - Order oracle added; bugfixes

   * - 0.0.3
     - 2020-06-10
     - Bugfixes in code generator

   * - 0.0.4     
     - 2020-06-15
     - MSVC compiler is now supported

   * - 0.0.5
     - 2021-08-02
     - Word shortening in monster implemented

   * - 0.0.6
     - 2021-12-01
     - Group operation accelerated

   * - 0.0.7
     - 2021-12-01
     - Bugfix in version generation

   * - 0.0.8
     - 2022-07-12
     - Performance improved

   * - 0.0.9
     - 2022-09-01
     - Performance improved

   * - 0.0.10
     - 2022-10-11
     - Support for cibuildwheel added

   * - 0.0.11
     - 2022-10-19
     - Bugfixes and macOS version added

   * - 0.0.12
     - 2023-01-09
     - Support for the group Y_555 and the Bimonster

   * - 0.0.13
     - 2023-10-27
     - Supports numbering elements of the Monster, and Python 3.12

   * - 0.0.14
     - 2024-01-19
     - Demonstration code for reduction in the Monster added




.. only:: html

   .. rubric:: **References**

.. bibliography:: references.bib



