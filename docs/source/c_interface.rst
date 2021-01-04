

======================================
The C interface of the mmgroup project
======================================


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


Description of the ``mmgroup.generators`` extension
===================================================

.. automodule:: mmgroup.dev.generators.generators_doc

Header file mmgroup_generators.h
--------------------------------

.. doxygenfile::  mmgroup_generators.h
 



.. automodule:: mmgroup.dev.generators.mm_group_n_doc

C interface for file mm_group_n.c
---------------------------------------

.. doxygenfile:: mm_group_n.c




.. automodule:: mmgroup.dev.generators.gen_leech_doc

C interface for file gen_leech.c
--------------------------------

.. doxygenfile:: gen_leech.c



.. automodule:: mmgroup.dev.generators.gen_xi_ref

C interface for file gen_xi_functions.c
---------------------------------------

.. doxygenfile:: gen_xi_functions.c


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


Description of the ``qstate12`` and ``clifford12`` extensions
=============================================================

C functions dealing with quadratic state vectors
------------------------------------------------


The C functions in modules ``qstate.c`` and ``qsmatrix.c`` perform 
operations on quadratic state matrices given by triples 
:math:`(e, A, Q)` as described in **API reference** the section 
:ref:`clifford-group-label`. 



A quadratic state vector :math:`v` of type
``qbstate12_type`` with component ``ncols = n`` models a complex 
vector in a vector space  :math:`V` of dimension :math:`2^n`, 
or an  :math:`2^{n-k} \times 2^k` matrix.

The basis of vector space :math:`V` is labelled by the elements of 
the  Boolean vector space :math:`\mathbb{F}_2^n`. In C and python 
programs we represent the element :math:`(x_{n-1}, \ldots, x_{0})` 
of :math:`\mathbb{F}_2^n` by the integer 
:math:`\sum_{0 \leq i < n} 2^i \cdot x_i`. This leads to a natural
representation of a vector :math:`v \in V` as a one-dimensional 
complex array of length :math:`2^n`, starting with index 
:math:`0`.

A quadratic state matrix is a quadratic shape vector augmented
by an information about its matrix shape. For more details we 
refer to the description of struct ``qbstate12_type`` in file
``clifford12.h``.

  
Header file ``clifford12.h``
----------------------------
   .. doxygenfile:: clifford12.h

Module ``qstate12.c``
---------------------
   
   
   .. doxygenfile:: qstate12.c

Module ``qmatrix12.c``
----------------------
   
   
   .. doxygenfile:: qmatrix12.c




