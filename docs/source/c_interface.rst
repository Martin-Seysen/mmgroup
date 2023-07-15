

======================================
The C interface of the mmgroup project
======================================


Introduction
============

This *document* describes the functionality of the C modules in this
project. For most of these C modules there is also a python extension.


.. highlight:: none


Description of the ``mmgroup.mat24`` extension
==============================================

.. automodule:: mmgroup.dev.mat24.mat24_doc


C interface
-----------

Header files
............

.. doxygenfile:: mat24_functions.h


C functions for the Mathieu group :math:`M_{24}`
................................................

.. doxygenfile:: mat24_functions.c


C functions for generating random elements of :math:`M_{24}`
..............................................................


.. doxygenfile:: mat24_random.c


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


.. _mmgroup-generators-label:

Description of the ``mmgroup.generators`` extension
===================================================

.. automodule:: mmgroup.dev.generators.generators_doc


.. _header-mmgroup-generators-label:


Header file mmgroup_generators.h
--------------------------------

.. doxygenfile::  mmgroup_generators.h
 

C functions implementing the group :math:`N_{0}`
-------------------------------------------------


.. automodule:: mmgroup.dev.generators.mm_group_n_doc

C interface for file mm_group_n.c
..................................

.. doxygenfile:: mm_group_n.c


.. _op-leech-label:

C functions for the operation of :math:`G_{x0}` on the Leech lattice
----------------------------------------------------------------------


.. automodule:: mmgroup.dev.generators.gen_leech_doc


C interface for file gen_leech.c
.................................

.. doxygenfile:: gen_leech.c


C interface for file gen_leech_type.c
......................................

.. doxygenfile:: gen_leech_type.c


C interface for file gen_leech3.c
.................................

.. doxygenfile:: gen_leech3.c


C interface for file gen_leech_reduce.c
........................................

.. doxygenfile:: gen_leech_reduce.c


C interface for file gen_leech_reduce_n.c
.........................................

.. doxygenfile:: gen_leech_reduce_n.c


C functions for the generator  :math:`\xi` of the monster group
----------------------------------------------------------------   

.. automodule:: mmgroup.dev.generators.gen_xi_ref

C interface for file gen_xi_functions.c
........................................

.. doxygenfile:: gen_xi_functions.c


C functions implementing a random generator
--------------------------------------------

.. automodule:: mmgroup.dev.generators.gen_random_doc



C interface for file gen_random.c
..................................

.. doxygenfile:: gen_random.c





Description of the ``qstate12`` and ``clifford12`` extensions
=============================================================

Quadratic state vectors
-------------------------


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

C functions in  ``bitmatrix64.c``
---------------------------------------

   .. doxygenfile:: bitmatrix64.c

C functions in  ``uint_sort.c``
--------------------------------

   .. doxygenfile:: uint_sort.c


C functions in   ``qstate12.c``
--------------------------------
     
   .. doxygenfile:: qstate12.c

C functions in   ``qstate12io.c``
----------------------------------
     
   .. doxygenfile:: qstate12io.c

C functions in   ``qmatrix12.c``
---------------------------------
     
   .. doxygenfile:: qmatrix12.c


.. _c-functions-G-x0-label:

Computing in the subgroup :math:`G_{x0}` of the Monster
--------------------------------------------------------

.. automodule:: mmgroup.dev.clifford12.xsp2co1_doc



C functions in ``xsp2co1.c``
----------------------------

   .. doxygenfile:: xsp2co1.c


C functions in ``xsp2co1_elem.c``
----------------------------------

   .. doxygenfile:: xsp2co1_elem.c


C functions in ``leech3matrix.c``
----------------------------------

   .. doxygenfile:: leech3matrix.c


C functions in ``involutions.c``
---------------------------------

   .. doxygenfile:: involutions.c


C functions in ``xsp2co1_traces.c``
------------------------------------

   .. doxygenfile:: xsp2co1_traces.c


C functions in ``xsp2co1_map.c``
---------------------------------

   .. doxygenfile:: xsp2co1_map.c




Description of the ``mmgroup.mm`` extension
===========================================


Module ``mmgroup.mm`` provides basic functions for the representations
of the monster modulo small numbers ``p``.  For each supported
modulus ``p`` there is also a specific module ``mmgroup.mm_op`` 
containing highly optimized C programs dealing with the
representations modulo ``p``.

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



Header file mm_basics.h
--------------------------------

.. doxygenfile::  mm_basics.h
 

C interface for file mm_aux.c
---------------------------------------

.. doxygenfile:: mm_aux.c


C interface for file mm_tables.c
---------------------------------------

.. doxygenfile:: mm_tables.c

C interface for file mm_random.c
---------------------------------------

Module ``mm_random.c`` is deprecated an no longer in use.


Description of the ``mmgroup.mm_op`` extensions
===============================================


Module ``mmgroup.mm_op`` implements the 196884-dimensional rational
representation :math:`\rho_p` of the monster group modulo several
fixed odd moduli :math:`p = 2^k-1`.

Module ``mmgroup.mm_op`` is implemented as an extension with 
``Cython``. The main source file for the ``mmgroup.mm_op`` extension
is the automatically generated file ``mm_op.pyx`` in directory 
``src.mmgroup.dev.pxd_files``. Each function exported by such an 
extension corresponds to a certain C function  exported from an 
automatically generated .c file.

For reasons of efficiency a dedicated set of C files is generated for 
each modulus ``p``. Corresponding C files for different moduli are 
generated from the same source file with extension ``.ske``. 


Internal operation
------------------

The user of this module should call functions from files
``mm_op_p_vector.c`` and ``mm_op_p_axis.c`` only. Each function
in this module takes the modulus ``p`` as its first argument.

For reasons of efficiency a dedicated set of C files is generated for 
each modulus ``p``. Corresponding C files for different moduli are 
generated from the same source file with extension ``.ske``. 
The C functions in module ``mm_op_p_vector.c`` and ``mm_op_p_axis.c``
contain automatically-generated ``case`` statements that call the
working functions for the appropriate modulus ``p``.

E.g. the function ``mm_op_pi`` in module ``mm_op_p_vector.c` calls
working functions ``mm3_op_pi``, ``mm7_op_pi``, etc.  in 
case ``p = 3``, ``p = 7``, etc. for doing the actual work in the 
representation of the Monster modulo ``p``. Such working functions
are implemented in different C files, with one C file for each 
modulus. So functions ``mm3_op_pi`` and ``mm7_op_pi`` are implemented
in the C files ``mm3_op_pi.c`` and `mm7_op_pi.c``, respecively.
Here such a C file may implement several working functions for the
same modules ``p``.

All functions exported from file ``mm_op_p_vector.c`` support the
same set of moduli ``p``. Functions exported from file
``mm_op_p_axis.c`` deal with 2A axes; the usually support moduli
3 and 15 (or just one of these two values) only

The process of switching to the appropriate function for a given
modulus is completely invisible for Python and Cython. In older
versions of the *mmgroup* package this was not the case. In the
main directory of the package there are some legacy python
modules in the main emulating the old behaviour. 




The basic table-providing class for ``mmgroup.mm_op``
-----------------------------------------------------

.. automodule:: mmgroup.dev.mm_op.mm_op

.. autoclass:: mmgroup.dev.mm_op.mm_op.MM_Op



Header file mm_op_p.h
--------------------------------

.. doxygenfile::  mm_op_p.h


C interface for file mm_op_p_vector.c
---------------------------------------

.. doxygenfile:: mm_op_p_vector.c


C interface for file mm_op_p_axis.c
---------------------------------------

.. doxygenfile:: mm_op_p_axis.c




Description of the ``mm_reduce`` extension
==========================================

The functions in this module implement the fast reduction of an
element of the monster group described in :cite:`Seysen22`. There
we define a triple of vector :math:`(v_1, v^+, v^-)` in the
representation :math:`\rho_{15}` such that an element 
:math:`g` of the monster van be recognized from the triple  
:math:`(v_1 \cdot g, v^+ \cdot g, v^- \cdot g)`. A precomputed
vector :math:`v_1` is stored in file ``mm_order_vector.c``.
Module ``mm_order.c`` contains functions for computing the
order of an element of the monster. Module ``mm_reduce.c``
contains the function ``mm_reduce_M`` that implements the
fast reduction algorithm in the monster group.
 


Generating an **order vector**
------------------------------

Python module ``mmgroup.dev.mm_reduce.find_order_vector``
.........................................................

.. automodule:: mmgroup.dev.mm_reduce.find_order_vector


Header file mm_reduce.h
--------------------------------------

.. doxygenfile:: mm_reduce.h



C interface for file mm_order_vector.c
--------------------------------------

.. doxygenfile:: mm_order_vector.c


C interface for file mm_order.c
-------------------------------

.. doxygenfile:: mm_order.c



C interface for file mm_reduce.c
--------------------------------

.. doxygenfile:: mm_reduce.c


C interface for file mm_shorten.c
---------------------------------

.. doxygenfile:: mm_shorten.c


C interface for file mm_compress.c
-----------------------------------

.. doxygenfile:: mm_compress.c


C interface for file mm_vector_v1_mod3.c
-----------------------------------------

.. doxygenfile:: mm_vector_v1_mod3.c




Shared libraries and dependencies between Cython extensions
===========================================================

For each Cython extension the relevant C functions are 
implemented in a shared library (or a DLL in Windows). Names
of the corresponding DLLs in windows are given in the following 
table:


.. table:: Shared libraries implementing Cython extensions
    :widths: 35 65

    ====================== =======================================
    Cython extension       Implemented in shared library
    ====================== =======================================
    ``mm_reduce``          ``mmgroup_reduce``
    ``mm_op``              ``mmgroup_mm_op``
    ``mm``                 ``mmgroup_mm_basics``
    ``clifford12``         ``mmgroup_clifford12``
    ``generators, mat24``  ``mmgroup_mat24``
    ====================== =======================================

In Linux, the name of each shared library is prefixed by the
string ``lib``. So the Windows DLL  ``mmgroup_mm_op.dll``
corresponds to the Linux shared library ``libmmgroup_mm_op.so.``



The Cython extension ``mm_<p>`` depends on other Cython extensions
as shown in the following paragraph:


.. code-block:: none
    
    mm_reduce
       + mm_op
          + mm
             + clifford12
                + generators
                    + mat24
 




