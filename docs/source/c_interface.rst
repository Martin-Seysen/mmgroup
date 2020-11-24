

================================
The C interface
================================


Introduction
============

This *document* describes the functionality of the C modules in this
project. For most of these C modules there is also a python extension.




Description of the ``mmgroup.mat24`` extension
==============================================

.. automodule:: mmgroup.dev.mat24.mat24_doc


C interface
-----------

.. doxygenfile:: mat24_functions.c




Generating C code for the ``mmgroup.mat24`` extension
-----------------------------------------------------

In this section we give a brief overview over the modules in
``mmgroup.dev.mat24`` used for generating C code for the 
``mmgroup.mat24`` extension.



Module ``mat24_ref``
....................

.. automodule:: mmgroup.dev.mat24.mat24_ref

Module ``mat24tables``
......................

.. automodule:: mmgroup.dev.mat24.mat24tables

Module ``mat24aux``
...................

.. automodule:: mmgroup.dev.mat24.mat24aux


Module ``mat24heptad``
......................

.. automodule:: mmgroup.dev.mat24.mat24heptad

Module ``make_addition_table``
..............................

.. automodule:: mmgroup.dev.mat24.make_addition_table

Module ``make_mul_transp``
..........................

.. automodule:: mmgroup.dev.mat24.make_mul_transp


Description of the ``mmgroup.mat24_xi`` extension
=================================================

This extension is used for computing tables containing the monomial
part of the operation of the generators :math:`\xi^i, i=1,2`.
These elements operate monomially on the ``98280`` entries of the 
representation :math:`\rho_p` with tags ``B``, ``C``, ``T``, and ``X``. 
One table is generated for each exponent :math:`i`.

Module ``mmgroup.dev.mat24_xi.mat24_xi_ref`` is a pure python 
substitute for this extension; but calculating the tables for the
generators :math:`\xi^i` with that extension takes a rather long time.

More details will be documented in a future version of this project.


Description of the ``mmgroup.mm`` extension
===========================================


Module ``mmgroup.mm`` provides basic functions for the representation
of the monster. 

Module ``mmgroup.mm`` is implemented as an extension with ``Cython``.
The main source file for that extension is ``mm_basics.pyx`` in
directory ``src.mmgroup.dev.mm_basics``. Most functions exported by
this extension are in one-to-one correspondence with certain C
functions exported from several automatically generated .c files. 

Representation of a vector in :math:`\rho_p`
--------------------------------------------

.. automodule:: mmgroup.dev.mm_basics.mm_doc

Basic table-providing classes for module ``mmgroup.mm``
-------------------------------------------------------

.. automodule:: mmgroup.dev.mm_basics.mm_basics

.. autoclass:: mmgroup.dev.mm_basics.mm_basics.MM_Const
   :members: snippet



More details will be documented in a future version of this project.


Description of the ``mmgroup.mm<p>`` extensions
===============================================


Module ``mmgroup.mm<p>`` implements the representation :math:`\rho_p`
for a fixed modulus :math:`p = 2^k-1`. E.g  module ``mmgroup.mm3`` 
implements :math:`\rho_3`. 

Module ``mmgroup.mm<p>`` is implemented as an extension with 
``Cython``. The main source file for the ``mmgroup.mm<p>`` extension
is the automatically generated file ``mm_op<p>.pyx`` in directory 
``src.mmgroup.dev.pxd_files``. Each function exported by such an 
extension corresponds to a certain C function  exported from an 
automatically generated .c file.

For reasons of efficiency a dedicated set of C files is generated for 
each modulus ``p``. Corresponding C files for different moduli are 
generated from the same source file with extension ``.ske``. 


The basic table-providing class for ``mmgroup.mm<p>``
-----------------------------------------------------

.. automodule:: mmgroup.dev.mm_op.mm_op

.. autoclass:: mmgroup.dev.mm_op.mm_op.MM_Op


More details will be documented in a future version of this project.


Computation in the subgroup :math:`2^{1+24}.Co_1` of the monster
================================================================

In Conway's construction :cite:`Con85` the monster :math:`\mathbb{M}`
has a subgroup :math:`G_{x0}` of structure 
:math:`2^{1+24}_+.\mbox{Co}_1`.
Here :math:`G_{x0}` is constructed as a diagonal product of the
two groups :math:`\mbox{Co}_0` of structure :math:`2.\mbox{Co}_1`
and :math:`N(4096_x)`. :math:`N(4096_x)` is also of structure
:math:`2^{1+24}_+.\mbox{Co}_1` but not isomorphic to  :math:`G_{x0}`.
Computation in  :math:`\mbox{Co}_0` is easy since that group has a 
:math:`24`-dimensional rational representation. The smallest real
representation of the group :math:`N(4096_x)` has dimension
:math:`4096`, so naive computation in that representation is
rather inefficient.

The group :math:`N(4096_x)` is a subgroup of the real Clifford group
:math:`\mathcal{C}_{12}`. The real Clifford group :math:`\mathcal{C}_{n}`
of structure :math:`2^{1+2n}_+.\mbox{O}_{2n}(2)` is defined e.g. in
:cite:`NRS01`. :math:`\mathcal{C}_{12}` is a subgroup
of the complex Clifford group :math:`\mathcal{X}_{n}`, which is also
defined in :cite:`NRS01`.

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
             lb_norm2, inv,
             power, order

.. autofunction:: mmgroup.structures.qs_matrix.qs_unit_matrix


