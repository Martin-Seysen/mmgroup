

================================
The mmgroup guide for developers
================================


Introduction
============

This *guide for developers* explains some selected topics of the
``mmgroup`` project which are relevant for developers who want to make 
changes or bugfixes. It also explains some mathematical aspects of the 
implementation which are not covered by :cite:`Seysen20`.

The *guide* is far from complete and we hope that we can provide 
more information in future versions of it.

Directory structure
===================

.. include:: guide_directory.inc


Installation from a source distribution
=======================================


.. include:: installation_from_source.inc


The build process
=================

.. include:: guide_build_process.inc


Some mathematical aspects of the implementation
===============================================

.. include:: guide_mathematical.inc


.. _code-generation-label:

Code generation
===============

In this section we describe the most important functions and classes
used for the automatic generation of C code.
  

.. automodule:: mmgroup.generate_c.make_c_tables_doc


Classes and functions provided by the code generator
----------------------------------------------------
 

.. autoclass:: mmgroup.generate_c.TableGenerator
   :members: generate, generate_pxd


.. autoclass:: mmgroup.generate_c.UserDirective

.. autoclass:: mmgroup.generate_c.UserFormat

          
.. autofunction::   mmgroup.generate_c.c_snippet  

.. autofunction::   mmgroup.generate_c.format_item

.. autofunction::   mmgroup.generate_c.make_doc

.. autofunction::   mmgroup.generate_c.make_table

.. autofunction::   mmgroup.generate_c.prepend_blanks
  
.. autofunction::   mmgroup.generate_c.pxd_to_function_list
  
.. autofunction::   mmgroup.generate_c.pxd_to_pyx

        


How the code generator is used
------------------------------

.. include:: guide_use_code_generation.inc




Description of the ``mmgroup.mat24`` extension
==============================================

.. automodule:: mmgroup.dev.mat24.mat24_doc




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

